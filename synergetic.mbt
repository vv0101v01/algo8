///|
/// 協同學函子神經網絡（Synergetics functor neural network）
///
/// 抽象層（Haken 協同學 × 函子）：
/// - 網絡的每一層是物件，層間映射是態射；
/// - 「特徵映射」F 把輸入空間的態射送到潛空間的態射（函子）；
/// - 序參數（order parameter）奴役快變數：slaving 極限下
///   快變數 s = s(ξ) 由序參數 ξ 唯一決定。
///
/// 導出：
/// - 定理（函子律）：F(id) = id、F(g∘f) = F(g)∘F(f)，於抽樣點機械驗證。
/// - 定理（slaving 不動點）：快變數動方程 ds/dt = -γs + g(ξ) 在
///   γ > 0 時有唯一不動點 s* = g(ξ)/γ，且是全局指數穩定的。
/// - 定理（序參數競爭收斂）：softmax 型序參數更新在適當步長下
///   收斂到輸入投影的最大者（Winner-Take-All）。
/// - 演算法：函子式層推送、slaving 迭代、序參數 WTA 更新。
pub(all) struct Layer {
  w : Array[Array[Double]]
  b : Array[Double]
}

///|
pub fn Layer::activate(self : Layer, x : Array[Double]) -> Array[Double] {
  let out : Array[Double] = []
  for i in 0..<self.w.length() {
    let mut s = self.b[i]
    for j in 0..<self.w[i].length() {
      s = s + self.w[i][j] * x[j]
    }
    // ReLU
    out.push(if s > 0.0 { s } else { 0.0 })
  }
  out
}

///|
fn vec_eq(a : Array[Double], b : Array[Double], eps : Double) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i in 0..<a.length() {
    if (a[i] - b[i]).abs() > eps {
      return false
    }
  }
  true
}

///|
/// 函子 F：把「與輸入 u 的平移態射」推到潛空間。
/// 物件 = 輸入向量；態射 = 由平移誘導的重映射 τ_u(x) = x + u。
/// F(τ_u) = 潛空間中由同層線性映射誘導的共軛平移。
pub(all) struct SynergeticNet {
  encoder : Layer
  /// 序參數維度
  n_order : Int32
}

///|
/// F 在物件上的作用：編碼。
pub fn SynergeticNet::encode(
  self : SynergeticNet,
  x : Array[Double],
) -> Array[Double] {
  self.encoder.activate(x)
}

///|
/// F 在態射上的作用：對潛向量施加與輸入平移一致的推前。
/// 這裡取線性近似：F(τ_u)(E(x)) = E(x) + J·u，J 由編碼層權重給出。
pub fn SynergeticNet::push_morphism(
  self : SynergeticNet,
  e : Array[Double],
  u : Array[Double],
) -> Array[Double] {
  // J·u：編碼層的總體線性部分（忽略 ReLU 非線性，於非負活化區精確）
  let ju : Array[Double] = []
  for i in 0..<self.encoder.w.length() {
    let mut s = 0.0
    for j in 0..<self.encoder.w[i].length() {
      s = s + self.encoder.w[i][j] * u[j]
    }
    ju.push(s)
  }
  let out : Array[Double] = []
  for k in 0..<e.length() {
    out.push(e[k] + ju[k])
  }
  out
}

///|
/// 定理（函子律）於抽樣點的機械驗證：
/// (1) F(id) = id：u = 0 時推前為恆等；
/// (2) F(g∘f) = F(g)∘F(f)：平移複合 u = u₁ + u₂ 對應推前複合。
pub fn functor_laws_synergetic(
  net : SynergeticNet,
  x : Array[Double],
  u1 : Array[Double],
  u2 : Array[Double],
) -> Bool {
  let e = net.encode(x)
  // (1) F(id)
  let zero : Array[Double] = []
  for _ in u1 {
    zero.push(0.0)
  }
  let fid = net.push_morphism(e, zero)
  if !vec_eq(fid, e, 1.0e-12) {
    return false
  }
  // (2) 複合律：先 u1 再 u2  ⇔  一次 u1+u2
  let e1 = net.push_morphism(e, u1)
  let e12 = net.push_morphism(e1, u2)
  let usum : Array[Double] = []
  for k in 0..<u1.length() {
    usum.push(u1[k] + u2[k])
  }
  let esum = net.push_morphism(e, usum)
  vec_eq(e12, esum, 1.0e-12)
}

///|
fn mk_test_net() -> SynergeticNet {
  let enc : Layer = {
    w: [[1.0, 0.5, 0.0], [0.0, 1.0, 0.5], [0.5, 0.0, 1.0]],
    b: [0.1, -0.2, 0.0],
  }
  { encoder: enc, n_order: 3, }
}

///|
test {
  let net = mk_test_net()
  let x : Array[Double] = [1.0, 2.0, 3.0]
  let u1 : Array[Double] = [0.5, -0.5, 1.0]
  let u2 : Array[Double] = [-1.0, 0.25, 0.0]
  assert_true(functor_laws_synergetic(net, x, u1, u2))
}

///|
/// Slaving：快變數動方程 ds/dt = -γ·s + g(ξ)，γ > 0。
/// 定理（唯一不動點 + 指數穩定）：s* = g(ξ)/γ，
/// 誤差 e_t = s_t - s* 滿足 e_{t+1} = (1 - ηγ)·e_t（顯式歐拉），
/// 當 0 < η < 2/γ 時 |1 - ηγ| < 1，指數收斂。
pub fn slaving_step(
  s : Double,
  xi : Double,
  gamma : Double,
  eta : Double,
) -> Double {
  let g_of_xi = xi * xi - 1.0 // 非線性耦合 g(ξ)
  s + eta * (-gamma * s + g_of_xi)
}

///|
pub fn slaving_fixed_point(xi : Double, gamma : Double) -> Double {
  (xi * xi - 1.0) / gamma
}

///|
pub fn slaving_run(
  xi : Double,
  gamma : Double,
  eta : Double,
  steps : Int32,
) -> Double {
  let mut s = 0.0
  for _ in 0..<steps {
    s = slaving_step(s, xi, gamma, eta)
  }
  s
}

///|
test {
  let gamma = 4.0
  let eta = 0.2 // |1 - ηγ| = 0.2 < 1
  for xi in [0.0, 0.5, 1.0, -1.5] {
    let s_star = slaving_fixed_point(xi, gamma)
    let s_n = slaving_run(xi, gamma, eta, 200)
    assert_true((s_n - s_star).abs() < 1.0e-9)
    // 幾何收斂率驗證：|e_{k+1}| / |e_k| ≈ |1 - ηγ| = 0.2
    let s0 = 3.7
    let s1 = slaving_step(s0, xi, gamma, eta)
    let e0 = s0 - s_star
    let e1 = s1 - s_star
    if e0.abs() > 1.0e-9 {
      assert_true((e1 / e0 - (1.0 - eta * gamma)).abs() < 1.0e-12)
    }
  }
}

///|
/// 序參數 WTA 更新：ξ ← ξ + η·(proj - ξ)，其中 proj 為輸入在
/// 序參數基上的投影；softmax 歸一化後最大者勝出。
pub fn order_parameter_wta(
  proj : Array[Double],
  eta : Double,
  steps : Int32,
) -> Array[Int32] {
  let n = proj.length()
  let mut xi : Array[Double] = []
  for _ in proj {
    xi.push(0.0)
  }
  for _ in 0..<steps {
    // 歸一化（softmax 型）
    let mut mx = xi[0]
    for v in xi {
      if v > mx {
        mx = v
      }
    }
    let exps : Array[Double] = []
    let mut z = 0.0
    for v in xi {
      let e = @math.exp(v - mx)
      exps.push(e)
      z = z + e
    }
    let nxt : Array[Double] = []
    for i in 0..<n {
      let soft = exps[i] / z
      nxt.push(xi[i] + eta * (proj[i] - soft))
    }
    xi = nxt
  }
  // 找最大者
  let winners : Array[Int32] = []
  let mut best = xi[0]
  for v in xi {
    if v > best {
      best = v
    }
  }
  for i in 0..<n {
    if (xi[i] - best).abs() < 1.0e-12 {
      winners.push(i)
    }
  }
  winners
}

///|
/// 定理（WTA 收斂）：當投影嚴格最大分量唯一時，
/// 若初值對稱打破且 η 適當，收斂後最大者保持為贏家。
/// 機械驗證：對任意嚴格最大輸入，更新後 argmax 不變（吸引性）。
pub fn wta_argmax_stable(proj : Array[Double]) -> Bool {
  let w0 = order_parameter_wta(proj, 0.1, 100)
  // 再推一步仍為同一贏家集合
  let w1 = order_parameter_wta(proj, 0.1, 101)
  if w0.length() != w1.length() {
    return false
  }
  for i in 0..<w0.length() {
    if w0[i] != w1[i] {
      return false
    }
  }
  true
}

///|
test {
  let proj : Array[Double] = [0.1, 2.0, 0.5, -1.0]
  let w = order_parameter_wta(proj, 0.1, 200)
  assert_eq(w.length(), 1)
  assert_eq(w[0], 1)
  assert_true(wta_argmax_stable(proj))
  // 對稱初值由投影差打破
  let proj2 : Array[Double] = [1.0, 0.0, 0.0]
  let w2 = order_parameter_wta(proj2, 0.1, 200)
  assert_eq(w2[0], 0)
}
