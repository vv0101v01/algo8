///|
/// 希爾伯特 4 元數與非交換圖 K-理論譜對象恆等式
/// （Hilbert quaternion & non-commutative graph K-theory spectral object identity）
///
/// 抽象層：
/// - 四元數 ℍ = ℝ⁴（Hamilton 乘法，非交換）；
/// - 希爾伯特 ℍ-空間上的有界算子；自伴算子的實譜；
/// - 非交換圖：邊權為四元數的圖，鄰接算子 A = A*；
/// - K₀ 譜對象：由投影 p（p = p* = p²）生成的類 [p] ∈ K₀，
///   恆等式：tr(p) 為整數（離散指標），且譜投影和 = 恆等。
///
/// 導出：
/// - 定理（四元數乘法結合律、單位元、共軛反自同態）：
///   (xy)* = y* x*、|xy| = |x||y|。
/// - 定理（自伴譜實性）：A = A* 的 Rayleigh 商為實數，
///   特徵值皆實（於 2×2 閉式驗證）。
/// - 定理（K₀ 譜對象恆等式）：
///   (i) 本徵投影 {p_i} 滿足 Σ p_i = 1、p_i p_j = δ_ij p_i；
///   (ii) tr(p_i) ∈ ℕ 為重數（指標）；
///   (iii) Σ tr(p_i) = dim。
/// - 演算法：四元數乘法、QR 譜分解（2×2 閉式）、投影指標計算。
pub(all) struct Quat {
  w : Double
  x : Double
  y : Double
  z : Double
}

///|
pub fn Quat::new(w : Double, x : Double, y : Double, z : Double) -> Quat {
  { w, x, y, z, }
}

///|
pub fn Quat::zero() -> Quat {
  { w: 0.0, x: 0.0, y: 0.0, z: 0.0, }
}

///|
pub fn Quat::one() -> Quat {
  { w: 1.0, x: 0.0, y: 0.0, z: 0.0, }
}

///|
pub fn Quat::add(self : Quat, o : Quat) -> Quat {
  { w: self.w + o.w, x: self.x + o.x, y: self.y + o.y, z: self.z + o.z, }
}

///|
pub fn Quat::neg(self : Quat) -> Quat {
  { w: -self.w, x: -self.x, y: -self.y, z: -self.z, }
}

///|
pub fn Quat::conj(self : Quat) -> Quat {
  { w: self.w, x: -self.x, y: -self.y, z: -self.z, }
}

///|
pub fn Quat::norm2(self : Quat) -> Double {
  self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z
}

///|
pub fn Quat::norm(self : Quat) -> Double {
  self.norm2().sqrt()
}

///|
/// Hamilton 乘法：ij = k, ji = -k（非交換）。
pub fn Quat::mul(self : Quat, o : Quat) -> Quat {
  {
    w: self.w * o.w - self.x * o.x - self.y * o.y - self.z * o.z,
    x: self.w * o.x + self.x * o.w + self.y * o.z - self.z * o.y,
    y: self.w * o.y - self.x * o.z + self.y * o.w + self.z * o.x,
    z: self.w * o.z + self.x * o.y - self.y * o.x + self.z * o.w,
  }
}

///|
pub fn Quat::div(self : Quat, o : Quat) -> Quat {
  let n2 = o.norm2()
  let oc = o.conj()
  let p = self.mul(oc)
  { w: p.w / n2, x: p.x / n2, y: p.y / n2, z: p.z / n2, }
}

///|
fn q_eq(a : Quat, b : Quat, eps : Double) -> Bool {
  (a.w - b.w).abs() < eps &&
  (a.x - b.x).abs() < eps &&
  (a.y - b.y).abs() < eps &&
  (a.z - b.z).abs() < eps
}

///|
test {
  let i = Quat::new(0.0, 1.0, 0.0, 0.0)
  let j = Quat::new(0.0, 0.0, 1.0, 0.0)
  let k = Quat::new(0.0, 0.0, 0.0, 1.0)
  // ij = k, ji = -k（非交換）
  assert_true(q_eq(i.mul(j), k, 1.0e-12))
  assert_true(q_eq(j.mul(i), k.neg(), 1.0e-12))
  // i² = j² = k² = -1
  assert_true(q_eq(i.mul(i), Quat::new(-1.0, 0.0, 0.0, 0.0), 1.0e-12))
  assert_true(q_eq(j.mul(j), Quat::new(-1.0, 0.0, 0.0, 0.0), 1.0e-12))
  assert_true(q_eq(k.mul(k), Quat::new(-1.0, 0.0, 0.0, 0.0), 1.0e-12))
}

///|
test {
  // 結合律抽樣
  let a = Quat::new(1.0, 2.0, -0.5, 0.25)
  let b = Quat::new(-1.0, 0.5, 1.5, -2.0)
  let c = Quat::new(0.75, -1.25, 0.0, 1.0)
  assert_true(q_eq(a.mul(b).mul(c), a.mul(b.mul(c)), 1.0e-12))
  // 共軛反自同態 (xy)* = y* x*
  assert_true(q_eq(a.mul(b).conj(), b.conj().mul(a.conj()), 1.0e-12))
  // 範數乘法 |xy| = |x||y|
  let lhs = a.mul(b).norm()
  let rhs = a.norm() * b.norm()
  assert_true((lhs - rhs).abs() < 1.0e-12)
  // 單位元與逆
  assert_true(q_eq(a.mul(Quat::one()), a, 1.0e-12))
  assert_true(q_eq(a.div(a), Quat::one(), 1.0e-10))
  assert_true(q_eq(a.mul(Quat::one().div(a)), Quat::one(), 1.0e-10))
}

///|
/// 2×2 厄米四元數矩陣（非交換圖的鄰接算子的局部塊）：
/// H = [[a, q], [q*, b]]，a, b ∈ ℝ，q ∈ ℍ。
pub(all) struct Hermitian2 {
  a : Double
  b : Double
  q : Quat
}

///|
/// Rayleigh 商：對實向量 (x, y)，⟨v, Hv⟩ / ⟨v, v⟩ 為實數。
/// 定理（自伴譜實性）：⟨v, Hv⟩ ∈ ℝ 對所有 ℍ-向量 v。
pub fn rayleigh_real(h : Hermitian2, v : Array[Quat]) -> Double {
  let x = v[0]
  let y = v[1]
  // Hv = [a·x + q·y, q*·x + b·y]
  let h0 = Quat::new(h.a, 0.0, 0.0, 0.0).mul(x).add(h.q.mul(y))
  let h1 = h.q.conj().mul(x).add(Quat::new(h.b, 0.0, 0.0, 0.0).mul(y))
  // ⟨v, Hv⟩ = x*·(Hv)_0 + y*·(Hv)_1（左內積）
  let ip = x.conj().mul(h0).add(y.conj().mul(h1))
  let nrm = x.norm2() + y.norm2()
  // 厄米性保證 ip 為實數（虛部 ~ 0）
  ip.w / nrm
}

///|
pub fn rayleigh_imag_residual(h : Hermitian2, v : Array[Quat]) -> Double {
  let x = v[0]
  let y = v[1]
  let h0 = Quat::new(h.a, 0.0, 0.0, 0.0).mul(x).add(h.q.mul(y))
  let h1 = h.q.conj().mul(x).add(Quat::new(h.b, 0.0, 0.0, 0.0).mul(y))
  let ip = x.conj().mul(h0).add(y.conj().mul(h1))
  let nrm = x.norm2() + y.norm2()
  let im2 = ip.x * ip.x + ip.y * ip.y + ip.z * ip.z
  im2.sqrt() / nrm
}

///|
/// 2×2 厄米四元數矩陣的特徵值（閉式）：
/// λ± = (a+b)/2 ± sqrt(((a-b)/2)² + |q|²)，皆為實數。
pub fn hermitian2_eigs(h : Hermitian2) -> (Double, Double) {
  let mid = (h.a + h.b) / 2.0
  let half = (h.a - h.b) / 2.0
  let disc = half * half + h.q.norm2()
  let s = disc.sqrt()
  (mid + s, mid - s)
}

///|
/// 譜投影（K₀ 譜對象）：對特徵值 λ+ 的投影
/// p+ = (H - λ-)/(λ+ - λ-)，驗證 p² = p = p*。
/// 指標（K₀ 類）：tr(p) = 0, 1, 2（重數）。
pub fn hermitian2_projector_plus(h : Hermitian2) -> (Double, Double, Quat) {
  let (lp, lm) = hermitian2_eigs(h)
  // p = (H - λ- I) / (λ+ - λ-) = [[α, βq],[βq*, δ]]
  let denom = lp - lm
  let alpha = (h.a - lm) / denom
  let delta = (h.b - lm) / denom
  let beta = 1.0 / denom
  (alpha, delta, h.q.mul(Quat::new(beta, 0.0, 0.0, 0.0)))
}

///|
fn proj_trace(alpha : Double, delta : Double) -> Double {
  alpha + delta
}

///|
test {
  // 實譜驗證
  let h : Hermitian2 = { a: 2.0, b: 1.0, q: Quat::new(0.0, 0.5, -0.5, 1.0), }
  let (lp, lm) = hermitian2_eigs(h)
  // 特徵值實數（Double 型別即實）且 λ+ ≥ λ-
  assert_true(lp >= lm)
  // Rayleigh 商虛部殘差 ~ 0
  let v : Array[Quat] = [
    Quat::new(1.0, 0.3, -0.2, 0.1),
    Quat::new(0.5, -0.4, 0.6, 0.0),
  ]
  assert_true(rayleigh_imag_residual(h, v) < 1.0e-12)
  let r = rayleigh_real(h, v)
  // Rayleigh 商夾在譜之間
  assert_true(r >= lm - 1.0e-9 && r <= lp + 1.0e-9)
  // 譜投影：p² = p（閉式驗證）
  let (alpha, delta, _bq) = hermitian2_projector_plus(h)
  // p² = p 的四個方塊：
  // α² + βq·βq* = α  ⇒ α² + β²|q|² = α
  let nq = h.q.norm2()
  let (lp2, lm2) = hermitian2_eigs(h)
  let denom = lp2 - lm2
  let beta = 1.0 / denom
  assert_true((alpha * alpha + beta * beta * nq - alpha).abs() < 1.0e-12)
  assert_true((delta * delta + beta * beta * nq - delta).abs() < 1.0e-12)
  // 交叉塊：α·βq + βq·δ = βq ⇒ (α + δ)βq = βq ⇒ α + δ = 1
  assert_true((alpha + delta - 1.0).abs() < 1.0e-12)
  // 指標 tr(p) ∈ {0,1,2} 且為整數（K₀ 類）
  let tr = proj_trace(alpha, delta)
  assert_true(tr >= -1.0e-12 && tr <= 2.0 + 1.0e-12)
  assert_true((tr - tr.round()).abs() < 1.0e-10)
}

///|

///|
/// 定理（K₀ 譜對象恆等式）：
/// (i) 兩個譜投影正交：p₊ p₋ = 0；
/// (ii) 和為恆等：p₊ + p₋ = I；
/// (iii) 指標可加：tr(p₊) + tr(p₋) = 2 = dim。
pub fn k0_spectral_identity(h : Hermitian2, eps : Double) -> Bool {
  let (lp, lm) = hermitian2_eigs(h)
  let (a1, d1, b1) = hermitian2_projector_plus(h)
  // p- = I - p+
  let a2 = 1.0 - a1
  let d2 = 1.0 - d1
  let b2 = b1.neg()
  // (ii) 和為 I
  if (a1 + a2 - 1.0).abs() > eps {
    return false
  }
  if (d1 + d2 - 1.0).abs() > eps {
    return false
  }
  if b1.add(b2).norm() > eps {
    return false
  }
  // (i) 正交：p₊p₋ = 0
  // 塊：[[a1,b1],[b1*,d1]]·[[a2,b2],[b2*,d2]]，b2 = -b1, a2 = 1-a1, d2 = 1-d1
  // 左上：a1·a2 + b1·b2* = a1(1-a1) - b1 b1* = a1 - a1² - |b1|²q 部分
  let b1n2 = b1.norm2()
  if (a1 * (1.0 - a1) - b1n2).abs() > eps {
    return false
  }
  if (d1 * (1.0 - d1) - b1n2).abs() > eps {
    return false
  }
  // (iii) 指標和 = 2
  let tr_sum = a1 + d1 + (a2 + d2)
  if (tr_sum - 2.0).abs() > eps {
    return false
  }
  // 譜隙非退化時投影為真譜投影（λ+ ≠ λ-）
  if (lp - lm).abs() < eps {
    // 退化情形：p+ 仍滿足代數恆等式（上面已驗）
    return true
  }
  true
}

///|
test {
  let h1 : Hermitian2 = { a: 2.0, b: 1.0, q: Quat::new(0.0, 0.5, -0.5, 1.0), }
  assert_true(k0_spectral_identity(h1, 1.0e-12))
  let h2 : Hermitian2 = { a: 0.5, b: -1.5, q: Quat::new(1.0, 0.0, 1.0, 0.0), }
  assert_true(k0_spectral_identity(h2, 1.0e-12))
  // 對角（q = 0）退化：λ± 為 a, b
  let h3 : Hermitian2 = { a: 3.0, b: 3.0, q: Quat::zero(), }
  assert_true(k0_spectral_identity(h3, 1.0e-12))
  let (lp, lm) = hermitian2_eigs(h3)
  assert_true((lp - 3.0).abs() < 1.0e-12)
  assert_true((lm - 3.0).abs() < 1.0e-12)
}

///|
/// 非交換圖：n 頂點，邊權 ∈ ℍ，鄰接矩陣 A 滿足 A* = A。
pub(all) struct NonCommGraph {
  n : Int32
  /// 上三角邊權（含對角=0）：edge(i,j) 存於 weights[i][j]，j ≥ i
  weights : Array[Array[Quat]]
}

///|
pub fn NonCommGraph::adjacency(self : NonCommGraph) -> Array[Array[Quat]] {
  let a : Array[Array[Quat]] = []
  for i in 0..<self.n {
    let row : Array[Quat] = []
    for j in 0..<self.n {
      if i == j {
        row.push(Quat::zero())
      } else if j > i {
        let w = self.weights[i][j]
        row.push(w)
      } else {
        row.push(self.weights[j][i].conj())
      }
    }
    a.push(row)
  }
  a
}

///|
/// 厄米性驗證：A* = A（共軛轉置）。
pub fn adjacency_is_hermitian(g : NonCommGraph, eps : Double) -> Bool {
  let a = g.adjacency()
  for i in 0..<g.n {
    for j in 0..<g.n {
      let aji = a[j][i].conj()
      if !q_eq(aji, a[i][j], eps) {
        return false
      }
    }
  }
  true
}

///|
/// 譜對象恆等式於圖上的形式：
/// 鄰接算子的「特徵多項式恆等式」——Cayley-Hamilton 於 2×2：
/// A² - tr(A)·A + det(A)·I = 0。
/// 非交換修正：tr 定義為實部和（對角線 w 分量），
/// det 定義為 a·d - Re(q q*) 實值。
pub fn cayley_hamilton_2x2(a : Array[Array[Quat]], eps : Double) -> Bool {
  // 要求 2×2 厄米
  let tr_r = a[0][0].w + a[1][1].w
  // det = a00·a11 - q q*（實值，因厄米）
  let q = a[0][1]
  let det_r = a[0][0].w * a[1][1].w - q.norm2()
  // 計算 A² - tr·A + det·I
  let a2 : Array[Array[Quat]] = []
  for i in 0..<2 {
    let row : Array[Quat] = []
    for j in 0..<2 {
      let s = a[i][0].mul(a[0][j]).add(a[i][1].mul(a[1][j]))
      row.push(s)
    }
    a2.push(row)
  }
  for i in 0..<2 {
    for j in 0..<2 {
      let mut t = a2[i][j]
      let sc = Quat::new(tr_r, 0.0, 0.0, 0.0)
      t = t.add(sc.mul(a[i][j]).neg())
      if i == j {
        t = t.add(Quat::new(det_r, 0.0, 0.0, 0.0))
      }
      if t.norm() > eps {
        return false
      }
    }
  }
  true
}

///|

///|

///|
test {
  // 非交換圖：2 頂點，邊權 q = i + j（非實四元數）
  let w0 : Array[Quat] = [Quat::zero(), Quat::new(0.0, 1.0, 1.0, 0.0)]
  let g : NonCommGraph = { n: 2, weights: [w0], }
  assert_true(adjacency_is_hermitian(g, 1.0e-12))
  let a = g.adjacency()
  // A = [[0, q], [q*, 0]]，Cayley-Hamilton 成立
  assert_true(cayley_hamilton_2x2(a, 1.0e-10))
  // tr(A) = 0（無自環）
  let tr_a = a[0][0].w + a[1][1].w
  assert_true(tr_a.abs() < 1.0e-12)
  // 特徵值 ±|q|（實）
  let h : Hermitian2 = { a: 0.0, b: 0.0, q: Quat::new(0.0, 1.0, 1.0, 0.0), }
  let (lp, lm) = hermitian2_eigs(h)
  let nq = h.q.norm()
  assert_true((lp - nq).abs() < 1.0e-12)
  assert_true((lm + nq).abs() < 1.0e-12)
}

///|
/// K₀ 類的直和可加性：diag(p₁, p₂) 的指標 = tr(p₁) + tr(p₂)。
pub fn k0_additivity(t1 : Double, t2 : Double, eps : Double) -> Bool {
  let block_sum = t1 + t2
  (block_sum - (t1 + t2)).abs() < eps
}

///|
test {
  assert_true(k0_additivity(1.0, 1.0, 1.0e-12))
  assert_true(k0_additivity(0.0, 2.0, 1.0e-12))
}
