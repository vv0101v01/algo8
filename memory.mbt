///|
/// 長期記憶細胞狀態轉移（Long-term memory cell state transition）
///
/// 抽象層：
/// 1. 閘控記憶細胞 `MemoryCell`：狀態轉移映射 step : Cell × Gates → Cell，
///    常數閘下構成可解的幾何級數（轉移半群）。
/// 2. 零散相位馬可夫鏈 `Phase`：記憶相位之間的機率轉移，
///    涵蓋遍歷鏈（平穩分佈）與吸收鏈（吸收時間）兩種結構。
pub(all) struct Gates {
  forget : Double
  input : Double
  output : Double
  candidate : Double
}

///|
pub(all) struct MemoryCell {
  c : Double
  h : Double
}

///|
/// 一步狀態轉移（LSTM 型閘控更新）。
pub fn MemoryCell::step(self : MemoryCell, g : Gates) -> MemoryCell {
  let c2 = g.forget * self.c + g.input * g.candidate
  { c: c2, h: g.output * @math.tanh(c2), }
}

///|
fn ipow(f : Double, t : Int32) -> Double {
  let mut acc = 1.0
  for _ in 0..<t {
    acc = acc * f
  }
  acc
}

///|
/// 定理（閉式解）：常數閘 f、i·g 下，第 t 步細胞電位為
///   c_t = f^t · c_0 + i·g · (1 - f^t)/(1 - f)   (f ≠ 1)
/// 即轉移半群的幾何級數求和公式。
pub fn memory_closed_form(
  c0 : Double,
  f : Double,
  ig : Double,
  t : Int32,
) -> Double {
  if (f - 1.0).abs() < 1.0e-14 {
    c0 + ig * t.to_double()
  } else {
    let ft = ipow(f, t)
    ft * c0 + ig * (1.0 - ft) / (1.0 - f)
  }
}

///|
/// 幾何閉式解與逐步迭代的一致性（轉移半群定理的機械驗證）。
test {
  let g : Gates = { forget: 0.9, input: 0.8, output: 1.0, candidate: 0.5, }
  let ig = g.input * g.candidate
  let mut cell : MemoryCell = { c: 1.0, h: 0.0, }
  for t in 1..<=200 {
    cell = cell.step(g)
    if t <= 20 {
      let closed = memory_closed_form(1.0, g.forget, ig, t)
      assert_true(approx(cell.c, closed, 1.0e-9))
    }
  }
  // 不動點定理：|f| < 1 時 c_t → i·g/(1-f)
  let fixed = ig / (1.0 - g.forget)
  assert_true(approx(cell.c, fixed, 1.0e-6))
}

///|
pub(all) enum Phase {
  ShortTerm
  Encoding
  LongTerm
  Retrieval
  Lost
}

///|
pub fn phase_index(p : Phase) -> Int32 {
  match p {
    ShortTerm => 0
    Encoding => 1
    LongTerm => 2
    Retrieval => 3
    Lost => 4
  }
}

///|
/// 遍歷轉移矩陣（ShortTerm/Encoding/LongTerm/Retrieval 四相通訊）。
pub fn memory_ergodic_matrix() -> Array[Array[Double]] {
  [
    [0.4, 0.4, 0.1, 0.1],
    [0.2, 0.4, 0.3, 0.1],
    [0.1, 0.2, 0.5, 0.2],
    [0.3, 0.2, 0.2, 0.3],
  ]
}

///|
/// 吸收轉移矩陣（Lost 為吸收態）。
pub fn memory_absorbing_matrix() -> Array[Array[Double]] {
  [
    [0.5, 0.4, 0.05, 0.0, 0.05],
    [0.1, 0.4, 0.4, 0.05, 0.05],
    [0.05, 0.1, 0.7, 0.1, 0.05],
    [0.2, 0.3, 0.3, 0.15, 0.05],
    [0.0, 0.0, 0.0, 0.0, 1.0],
  ]
}

///|
fn row_times_matrix(
  row : Array[Double],
  m : Array[Array[Double]],
) -> Array[Double] {
  let n = m.length()
  let out : Array[Double] = []
  for j in 0..<n {
    let mut s = 0.0
    for i in 0..<n {
      s = s + row[i] * m[i][j]
    }
    out.push(s)
  }
  out
}

///|
fn l1_rows(a : Array[Double], b : Array[Double]) -> Double {
  let mut s = 0.0
  for i in 0..<a.length() {
    s = s + (a[i] - b[i]).abs()
  }
  s
}

///|
/// 演算法：平穩分佈的冪迭代 π_{k+1} = π_k P。
pub fn stationary_power_iteration(
  p : Array[Array[Double]],
  steps : Int32,
) -> Array[Double] {
  let n = p.length()
  let mut pi : Array[Double] = []
  for _ in 0..<n {
    pi.push(1.0 / n.to_double())
  }
  for _ in 0..<steps {
    pi = row_times_matrix(pi, p)
  }
  pi
}

///|
/// 定理（遍歷鏈）：本鏈為素性（primitive）時存在唯一平穩分佈，
/// 滿足 πP = π，且與初始分佈無關；冪迭代收斂到它。
test {
  let p = memory_ergodic_matrix()
  // 素性驗證：P^3 全為正。
  let mut m = p
  for _ in 0..<3 {
    let out : Array[Array[Double]] = []
    for row in m {
      out.push(row_times_matrix(row, p))
    }
    m = out
  }
  for row in m {
    for x in row {
      assert_true(x > 0.0)
    }
  }
  let pi = stationary_power_iteration(p, 200)
  let pi2 = row_times_matrix(pi, p)
  assert_true(l1_rows(pi, pi2) < 1.0e-12)
  // 初始分佈無關性：不同初值收斂到同一 π。
  let n = p.length()
  let mut other : Array[Double] = []
  for i in 0..<n {
    other.push(if i == 0 { 1.0 } else { 0.0 })
  }
  for _ in 0..<200 {
    other = row_times_matrix(other, p)
  }
  assert_true(l1_rows(pi, other) < 1.0e-9)
  let mut total = 0.0
  for x in pi {
    total = total + x
  }
  assert_true(approx(total, 1.0, 1.0e-12))
  // 相位建構子全被使用。
  let all_phases = [ShortTerm, Encoding, LongTerm, Retrieval, Lost]
  assert_eq(phase_index(all_phases[4]), 4)
  assert_eq(phase_index(all_phases[0]), 0)
}

///|
/// 演算法：吸收時間的兩種算法——
/// (1) 基本矩陣線性方程組 (I-Q)t = 1；
/// (2) 生存機率級數 E[T] = Σ_n P(T > n)。
/// 定理：兩者相等。
pub fn absorption_time_linear_solve(
  q : Array[Array[Double]],
  r : Array[Double],
) -> Array[Double] {
  let n = q.length()
  // 建立 (I-Q) 與增廣欄（就地高斯-約當消去）。
  let a : Array[Array[Double]] = []
  for i in 0..<n {
    let row : Array[Double] = []
    for j in 0..<n {
      row.push((if i == j { 1.0 } else { 0.0 }) - q[i][j])
    }
    a.push(row)
  }
  let rhs = r
  for col in 0..<n {
    let mut piv = col
    for i in col..<n {
      if a[i][col].abs() > a[piv][col].abs() {
        piv = i
      }
    }
    let tmp = a[col]
    a[col] = a[piv]
    a[piv] = tmp
    let tv = rhs[col]
    rhs[col] = rhs[piv]
    rhs[piv] = tv
    let d = a[col][col]
    for j in col..<n {
      a[col][j] = a[col][j] / d
    }
    rhs[col] = rhs[col] / d
    for i in 0..<n {
      if i != col {
        let f = a[i][col]
        if f.abs() > 1.0e-15 {
          for j in col..<n {
            a[i][j] = a[i][j] - f * a[col][j]
          }
          rhs[i] = rhs[i] - f * rhs[col]
        }
      }
    }
  }
  rhs
}

///|
pub fn absorption_time_series(
  start : Int32,
  p : Array[Array[Double]],
  n_transient : Int32,
  tol : Double,
) -> Double {
  let row : Array[Double] = []
  for i in 0..<p.length() {
    row.push(if i == start { 1.0 } else { 0.0 })
  }
  let mut cur = row
  let mut acc = 0.0
  for _ in 0..<100000 {
    // P(T > step) = 尚在暫態的總機率。
    let mut transient_mass = 0.0
    for i in 0..<n_transient {
      transient_mass = transient_mass + cur[i]
    }
    acc = acc + transient_mass
    if transient_mass < tol {
      break
    }
    cur = row_times_matrix(cur, p)
  }
  acc
}

///|
test {
  let p = memory_absorbing_matrix()
  // 列和為 1（隨機矩陣）。
  for row in p {
    let mut s = 0.0
    for x in row {
      s = s + x
    }
    assert_true(approx(s, 1.0, 1.0e-12))
  }
  // 暫態子鏈 Q 與逃逸向量 r。
  let q : Array[Array[Double]] = []
  let r : Array[Double] = []
  for i in 0..<4 {
    let qrow : Array[Double] = []
    for j in 0..<4 {
      qrow.push(p[i][j])
    }
    q.push(qrow)
    r.push(p[i][4])
  }
  // 吸收機率為 1：Q^n → 0。
  let mut m = q
  for _ in 0..<500 {
    let out : Array[Array[Double]] = []
    for row in m {
      out.push(row_times_matrix(row, q))
    }
    m = out
  }
  for row in m {
    for x in row {
      assert_true(x.abs() < 1.0e-9)
    }
  }
  // 算法 (1)：解 (I-Q)t = 1。
  let ones : Array[Double] = [1.0, 1.0, 1.0, 1.0]
  let t_lin = absorption_time_linear_solve(q, ones)
  // 算法 (2)：級數求和。
  let t_ser = absorption_time_series(0, p, 4, 1.0e-13)
  assert_true(approx(t_lin[0], t_ser, 1.0e-6))
  assert_true(t_lin[0] > 1.0)
}
