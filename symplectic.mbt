///|
/// 神經偏微分辛相調恆等式方程（Neural-PDE symplectic phase-harmony identity）
///
/// 抽象層：
/// - 神經動力系統取 Hamilton 形式：dq/dt = ∂H/∂p, dp/dt = -∂H/∂q，
///   H 由小型神經網絡（線性 + tanh）表示；
/// - 「辛相調恆等式」：辛 Euler 映射 Φ_h 保持辛形式 ω = dq ∧ dp，
///   即 Jacobi 矩陣 J 的滿足 JᵀΩJ = Ω（Ω 為辛矩陣）；
/// - 對應幾何：相空間面積守恆（離散 Liouville 定理）。
///
/// 導出：
/// - 定理（辛恆等式）：對可分離 H = T(p) + V(q)，辛 Euler
///   p' = p - h∇V(q), q' = q + h∇T(p') 的變分矩陣滿足
///   MᵀΩM = Ω，精確成立（非近似）。
/// - 定理（面積守恆的離散版本）：det(∂(q',p')/∂(q,p)) = 1。
/// - 演算法：辛 Euler 積分器、變分方程追蹤、長期能量漂移量測。
pub(all) struct NeuralH {
  /// 動量項權重（T(p) = Σ wᵢ pᵢ²/2）
  tw : Array[Double]
  /// 位置勢能網絡：V(q) = Σ aⱼ·tanh(Σ Wⱼₖ qₖ + bⱼ)
  va : Array[Double]
  vw : Array[Array[Double]]
  vb : Array[Double]
}

///|
pub fn NeuralH::kinetic(self : NeuralH, p : Array[Double]) -> Double {
  let mut s = 0.0
  for i in 0..<p.length() {
    s = s + self.tw[i] * p[i] * p[i] / 2.0
  }
  s
}

///|
pub fn NeuralH::potential(self : NeuralH, q : Array[Double]) -> Double {
  let mut s = 0.0
  for j in 0..<self.va.length() {
    let mut z = self.vb[j]
    for k in 0..<q.length() {
      z = z + self.vw[j][k] * q[k]
    }
    s = s + self.va[j] * @math.tanh(z)
  }
  s
}

///|
pub fn NeuralH::energy(
  self : NeuralH,
  q : Array[Double],
  p : Array[Double],
) -> Double {
  self.kinetic(p) + self.potential(q)
}

///|
/// ∇V(q)
pub fn NeuralH::grad_potential(
  self : NeuralH,
  q : Array[Double],
) -> Array[Double] {
  let g : Array[Double] = []
  for k in 0..<q.length() {
    let mut s = 0.0
    for j in 0..<self.va.length() {
      let mut z = self.vb[j]
      for m in 0..<q.length() {
        z = z + self.vw[j][m] * q[m]
      }
      s = s + self.va[j] * (1.0 - @math.tanh(z) * @math.tanh(z)) * self.vw[j][k]
    }
    g.push(s)
  }
  g
}

///|
/// ∂T/∂p = diag(tw)·p
pub fn NeuralH::grad_kinetic(
  self : NeuralH,
  p : Array[Double],
) -> Array[Double] {
  let g : Array[Double] = []
  for i in 0..<p.length() {
    g.push(self.tw[i] * p[i])
  }
  g
}

///|
/// 辛 Euler 一步（可分離 Hamilton）：
///   p' = p - h·∇V(q)
///   q' = q + h·∇T(p')
pub fn symplectic_euler_step(
  h : NeuralH,
  q : Array[Double],
  p : Array[Double],
  dt : Double,
) -> (Array[Double], Array[Double]) {
  let gv = h.grad_potential(q)
  let p2 : Array[Double] = []
  for i in 0..<p.length() {
    p2.push(p[i] - dt * gv[i])
  }
  let gt = h.grad_kinetic(p2)
  let q2 : Array[Double] = []
  for k in 0..<q.length() {
    q2.push(q[k] + dt * gt[k])
  }
  (q2, p2)
}

///|
/// 一般（非辛）顯式 Euler，作為對照：
///   p' = p - h·∇V(q), q' = q + h·∇T(p)
pub fn explicit_euler_step(
  h : NeuralH,
  q : Array[Double],
  p : Array[Double],
  dt : Double,
) -> (Array[Double], Array[Double]) {
  let gv = h.grad_potential(q)
  let gt = h.grad_kinetic(p)
  let p2 : Array[Double] = []
  for i in 0..<p.length() {
    p2.push(p[i] - dt * gv[i])
  }
  let q2 : Array[Double] = []
  for k in 0..<q.length() {
    q2.push(q[k] + dt * gt[k])
  }
  (q2, p2)
}

///|
/// 變分矩陣：對 (q,p) 的 2n×2n Jacobi，以有限差分計算。
fn variational_jacobian(
  step : (Array[Double], Array[Double], Double) -> (
    Array[Double],
    Array[Double],
  ),
  q : Array[Double],
  p : Array[Double],
  eps : Double,
) -> Array[Array[Double]] {
  let n = q.length()
  let dim = 2 * n
  let cols : Array[Array[Double]] = []
  for c in 0..<dim {
    let dq : Array[Double] = []
    let dp : Array[Double] = []
    for i in 0..<n {
      dq.push(if i == c && c < n { eps } else { 0.0 })
      dp.push(if i == c - n && c >= n { eps } else { 0.0 })
    }
    // 中心差分對每個基方向
    let qp : Array[Double] = []
    let pp : Array[Double] = []
    for i in 0..<n {
      qp.push(q[i] + dq[i])
      pp.push(p[i] + dp[i])
    }
    let (q1, p1) = step(qp, pp, eps)
    let qm : Array[Double] = []
    let pm : Array[Double] = []
    for i in 0..<n {
      qm.push(q[i] - dq[i])
      pm.push(p[i] - dp[i])
    }
    let (q0, p0) = step(qm, pm, eps)
    let col : Array[Double] = []
    for i in 0..<n {
      col.push((q1[i] - q0[i]) / (2.0 * eps))
    }
    for i in 0..<n {
      col.push((p1[i] - p0[i]) / (2.0 * eps))
    }
    cols.push(col)
  }
  // 轉置：jac[r][c] = cols[c][r]
  let jac : Array[Array[Double]] = []
  for r in 0..<dim {
    let row : Array[Double] = []
    for c in 0..<dim {
      row.push(cols[c][r])
    }
    jac.push(row)
  }
  jac
}

///|
/// 定理（辛相調恆等式）：MᵀΩM = Ω，其中
/// Ω = [[0, I], [-I, 0]]。以 Frobenius 範數驗證殘差。

///|
/// 共用：辛形式殘差 ‖MᵀΩM − Ω‖_F，M 為 step 的變分 Jacobi。
fn symplectic_form_residual(
  step : (Array[Double], Array[Double], Double) -> (
    Array[Double],
    Array[Double],
  ),
  q : Array[Double],
  p : Array[Double],
) -> Double {
  let n = q.length()
  let dim = 2 * n
  let m = variational_jacobian(step, q, p, 1.0e-5)
  // Ω = [[0, I], [-I, 0]]
  let omega : Array[Array[Double]] = []
  for r in 0..<dim {
    let row : Array[Double] = []
    for c in 0..<dim {
      let v = if r < n && c >= n && c - r == n {
        1.0
      } else if r >= n && c < n && r - c == n {
        -1.0
      } else {
        0.0
      }
      row.push(v)
    }
    omega.push(row)
  }
  // Mᵀ
  let mt : Array[Array[Double]] = []
  for r in 0..<dim {
    let row : Array[Double] = []
    for c in 0..<dim {
      row.push(m[c][r])
    }
    mt.push(row)
  }
  // MᵀΩ
  let om : Array[Array[Double]] = []
  for r in 0..<dim {
    let row : Array[Double] = []
    for c in 0..<dim {
      let mut acc = 0.0
      for k in 0..<dim {
        acc = acc + mt[r][k] * omega[k][c]
      }
      row.push(acc)
    }
    om.push(row)
  }
  // MᵀΩM
  let mtom : Array[Array[Double]] = []
  for r in 0..<dim {
    let row : Array[Double] = []
    for c in 0..<dim {
      let mut acc = 0.0
      for k in 0..<dim {
        acc = acc + om[r][k] * m[k][c]
      }
      row.push(acc)
    }
    mtom.push(row)
  }
  // ‖MᵀΩM − Ω‖_F
  let mut resid = 0.0
  for r in 0..<dim {
    for c in 0..<dim {
      let d = mtom[r][c] - omega[r][c]
      resid = resid + d * d
    }
  }
  resid.sqrt()
}

///|
/// 定理（辛相調恆等式）：MᵀΩM = Ω。以 Frobenius 範數驗證殘差。
pub fn symplectic_identity_residual(
  h : NeuralH,
  q : Array[Double],
  p : Array[Double],
  dt : Double,
) -> Double {
  let step = fn(
    qq : Array[Double],
    pp : Array[Double],
    _dt : Double,
  ) -> (Array[Double], Array[Double]) {
    symplectic_euler_step(h, qq, pp, dt)
  }
  symplectic_form_residual(step, q, p)
}

///|
/// 離散面積守恆：det Jacobi = 1（n = 1 的 2×2 情形直接公式）。
pub fn symplectic_det_area(
  q : Array[Double],
  p : Array[Double],
  dt : Double,
) -> Double {
  let h = default_neural_h()
  let step = fn(
    qq : Array[Double],
    pp : Array[Double],
    _dt : Double,
  ) -> (Array[Double], Array[Double]) {
    symplectic_euler_step(h, qq, pp, dt)
  }
  let m = variational_jacobian(step, q, p, 1.0e-5)
  det_general(m)
}

///|
/// 一般 n×n 行列式（高斯消去，部分樞軸）。
fn det_general(a_in : Array[Array[Double]]) -> Double {
  let n = a_in.length()
  let a : Array[Array[Double]] = []
  for row in a_in {
    let r : Array[Double] = []
    for x in row {
      r.push(x)
    }
    a.push(r)
  }
  let mut det = 1.0
  for col in 0..<n {
    let mut piv = col
    for i in col..<n {
      if a[i][col].abs() > a[piv][col].abs() {
        piv = i
      }
    }
    if a[piv][col].abs() < 1.0e-300 {
      return 0.0
    }
    if piv != col {
      let tmp = a[col]
      a[col] = a[piv]
      a[piv] = tmp
      det = -det
    }
    det = det * a[col][col]
    for i in (col + 1)..<n {
      let f = a[i][col] / a[col][col]
      if f.abs() > 0.0 {
        for j in col..<n {
          a[i][j] = a[i][j] - f * a[col][j]
        }
      }
    }
  }
  det
}

///|
pub fn default_neural_h() -> NeuralH {
  {
    tw: [1.0, 1.0],
    va: [0.6, -0.4],
    vw: [[1.2, -0.7], [0.5, 1.1]],
    vb: [0.1, -0.3],
  }
}

///|
test {
  let h = default_neural_h()
  let q : Array[Double] = [0.7, -0.4]
  let p : Array[Double] = [0.3, 0.9]
  // 辛恆等式殘差應 ~ 0（有限差分精度內）
  let r = symplectic_identity_residual(h, q, p, 0.1)
  assert_true(r < 1.0e-6)
  // 面積 det = 1
  let d = symplectic_det_area(q, p, 0.1)
  assert_true((d - 1.0).abs() < 1.0e-6)
}

///|
/// 定理（非辛對照）：顯式 Euler 的 Jacobi 不保持 Ω（殘差顯著非零）。

///|
/// 定理（非辛對照）：顯式 Euler 的 Jacobi 不保持 Ω（殘差顯著非零）。
pub fn explicit_identity_residual(
  h : NeuralH,
  q : Array[Double],
  p : Array[Double],
  dt : Double,
) -> Double {
  let step = fn(
    qq : Array[Double],
    pp : Array[Double],
    _dt : Double,
  ) -> (Array[Double], Array[Double]) {
    explicit_euler_step(h, qq, pp, dt)
  }
  symplectic_form_residual(step, q, p)
}

///|
test {
  let h = default_neural_h()
  let q : Array[Double] = [0.7, -0.4]
  let p : Array[Double] = [0.3, 0.9]
  let r = explicit_identity_residual(h, q, p, 0.1)
  assert_true(r > 1.0e-3)
}

///|
/// 定理（長期能量穩定）：辛積分器長時間平均能量漂移
/// 遠小於顯式 Euler（幾何積分器的結構保持優勢）。
pub fn energy_drift(symp : Bool, steps : Int32, dt : Double) -> Double {
  let h = default_neural_h()
  let q0 : Array[Double] = [0.7, -0.4]
  let p0 : Array[Double] = [0.3, 0.9]
  let e0 = h.energy(q0, p0)
  let mut q = q0
  let mut p = p0
  for _ in 0..<steps {
    let (nq, np) = if symp {
      symplectic_euler_step(h, q, p, dt)
    } else {
      explicit_euler_step(h, q, p, dt)
    }
    q = nq
    p = np
  }
  let e1 = h.energy(q, p)
  (e1 - e0).abs()
}

///|
test {
  let dt = 0.05
  let steps = 2000
  let d_s = energy_drift(true, steps, dt)
  let d_e = energy_drift(false, steps, dt)
  // 辛方法漂移顯著更小
  assert_true(d_s < d_e)
  assert_true(d_s < 0.05)
}
