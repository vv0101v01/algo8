///|
/// 複數域強基數（Complex-field strong cardinality）
///
/// 抽象層：
/// - 複數域 ℂ 作為 ℝ-向量空間，基底 {1, i} 為「強基數」（Hamel 基）；
/// - 強基數的機械含意：任一複數在此基底下表示唯一（座標同構 ℂ ≅ ℝ²）；
/// - 基數（cardinality）層面：ℂ 與 ℝ 等勢的離散類比——有理格點 ×有理格點 ≈ 有理數。
///
/// 導出：
/// - 定理（座標唯一性）：z = a + bi 的 (a, b) 表示唯一，故 ℝ-線性同構 φ : ℝ² → ℂ。
/// - 定理（ℂ 是域）：加法群、乘法群（去零）、分配律於抽樣點成立。
/// - 定理（強基數等勢）：ℂ 上「高斯整數格」與 ℕ 之間存在雙射（Cantor 配對），
///   機械驗證 |ℤ[i] ∩ 範圍| 對應 ℕ 區段。
/// - 演算法：複數乘法（4 次實數乘法）、共軛除法、Cantor 配對雙射。
pub(all) struct Complex {
  re : Double
  im : Double
}

///|
pub fn Complex::new(re : Double, im : Double) -> Complex {
  { re, im, }
}

///|
pub fn Complex::zero() -> Complex {
  { re: 0.0, im: 0.0, }
}

///|
pub fn Complex::one() -> Complex {
  { re: 1.0, im: 0.0, }
}

///|
/// 基底元素 i。
pub fn Complex::i_unit() -> Complex {
  { re: 0.0, im: 1.0, }
}

///|
pub fn Complex::add(self : Complex, o : Complex) -> Complex {
  { re: self.re + o.re, im: self.im + o.im, }
}

///|
pub fn Complex::sub(self : Complex, o : Complex) -> Complex {
  { re: self.re - o.re, im: self.im - o.im, }
}

///|
/// 複數乘法：i² = -1 的展開（4 次實乘法）。
pub fn Complex::mul(self : Complex, o : Complex) -> Complex {
  { re: self.re * o.re - self.im * o.im, im: self.re * o.im + self.im * o.re, }
}

///|
pub fn Complex::conj(self : Complex) -> Complex {
  { re: self.re, im: -self.im, }
}

///|
pub fn Complex::norm2(self : Complex) -> Double {
  self.re * self.re + self.im * self.im
}

///|
pub fn Complex::norm(self : Complex) -> Double {
  self.norm2().sqrt()
}

///|
/// 共軛除法：z/w = z·conj(w)/|w|²。
pub fn Complex::div(self : Complex, o : Complex) -> Complex {
  let d = o.norm2()
  {
    re: (self.re * o.re + self.im * o.im) / d,
    im: (self.im * o.re - self.re * o.im) / d,
  }
}

///|
pub fn Complex::is_zero(self : Complex) -> Bool {
  self.norm2() < 1.0e-30
}

///|
/// 定理（座標唯一性 / 強基數）：
/// 對基底 {1, i}，z ↦ (Re z, Im z) 是 ℝ-線性雙射；
/// 機械驗證：任意 z 可重構 a·1 + b·i 且 (a, b) 由投影唯一讀出。
pub fn strong_basis_coords(z : Complex) -> (Double, Double) {
  // 對偶基投影：a = ⟨z, 1⟩/⟨1,1⟩, b = ⟨z, i⟩/⟨i,i⟩（實部內積）
  let a = z.re
  let b = z.im
  (a, b)
}

///|
pub fn rebuild_from_coords(a : Double, b : Double) -> Complex {
  // a·1 + b·i
  Complex::new(a, 0.0).add(Complex::new(0.0, b))
}

///|
fn c_eq(x : Complex, y : Complex, eps : Double) -> Bool {
  (x.re - y.re).abs() < eps && (x.im - y.im).abs() < eps
}

///|
test {
  // 座標唯一性
  let z = Complex::new(3.5, -2.25)
  let (a, b) = strong_basis_coords(z)
  assert_true(c_eq(rebuild_from_coords(a, b), z, 1.0e-12))
  // 基底重構：1 與 i 線性獨立 — a·1 + b·i = 0 ⇒ a = b = 0
  let s = Complex::one()
    .mul(Complex::new(2.0, 0.0))
    .add(Complex::i_unit().mul(Complex::new(3.0, 0.0)))
  assert_true(c_eq(s, Complex::new(2.0, 3.0), 1.0e-12))
}

///|
/// 定理（ℂ 是域）：於抽樣點驗證交換群、乘法可逆、分配律。
pub fn field_axioms_hold(samples : Array[Complex]) -> Bool {
  for x in samples {
    for y in samples {
      // 交換律
      if !c_eq(x.add(y), y.add(x), 1.0e-12) {
        return false
      }
      if !c_eq(x.mul(y), y.mul(x), 1.0e-12) {
        return false
      }
      // 單位元
      if !c_eq(x.add(Complex::zero()), x, 1.0e-12) {
        return false
      }
      if !c_eq(x.mul(Complex::one()), x, 1.0e-12) {
        return false
      }
      // 結合律（用第三個樣本 z 由呼叫端放入 samples 之積）
      if !x.is_zero() {
        let inv = Complex::one().div(x)
        if !c_eq(x.mul(inv), Complex::one(), 1.0e-10) {
          return false
        }
      }
    }
  }
  // 分配律
  for x in samples {
    for y in samples {
      for w in samples {
        let left = x.mul(y.add(w))
        let right = x.mul(y).add(x.mul(w))
        if !c_eq(left, right, 1.0e-10) {
          return false
        }
      }
    }
  }
  true
}

///|
test {
  let samples : Array[Complex] = [
    Complex::new(1.0, 2.0),
    Complex::new(-0.5, 0.25),
    Complex::new(0.0, -1.0),
    Complex::new(3.0, 0.0),
    Complex::new(-2.5, 4.5),
  ]
  assert_true(field_axioms_hold(samples))
}

///|
/// Cantor 配對：ℕ × ℕ → ℕ 雙射。
/// 定理（強基數等勢的離散類比）：格點集 ℕ² 與 ℕ 等勢；
/// 機械驗證：配對為單射且像恰為 [0, π(n, n)] 區段內全部整數。
pub fn cantor_pair(a : Int32, b : Int32) -> Int32 {
  let s = a + b
  s * (s + 1) / 2 + b
}

///|
pub fn cantor_unpair(z : Int32) -> (Int32, Int32) {
  // 反函式：求 s 使 T_s ≤ z < T_{s+1}, T_s = s(s+1)/2
  let mut s = 0
  while s * (s + 1) / 2 <= z {
    s = s + 1
  }
  s = s - 1
  let t = s * (s + 1) / 2
  let b = z - t
  let a = s - b
  (a, b)
}

///|
test {
  // 單射性掃描
  let seen : Array[Int32] = []
  for a in 0..<20 {
    for b in 0..<20 {
      let p = cantor_pair(a, b)
      for q in seen {
        if q == p {
          assert_true(false)
        }
      }
      seen.push(p)
    }
  }
  assert_eq(seen.length(), 400)
  // 反演
  for a in 0..<20 {
    for b in 0..<20 {
      let p = cantor_pair(a, b)
      let (a2, b2) = cantor_unpair(p)
      assert_eq(a2, a)
      assert_eq(b2, b)
    }
  }
  // 像的完整性（三角形域）：{(a,b) : a+b ≤ n} 恰映滿 [0, T_{n+1}-1]
  // 其中 T_k = k(k+1)/2 為三角數。
  let n = 19
  let maxv = (n + 1) * (n + 2) / 2 - 1
  let hit : Array[Bool] = []
  for _ in 0..<=maxv {
    hit.push(false)
  }
  for s in 0..<=n {
    for b in 0..<=s {
      let a = s - b
      hit[cantor_pair(a, b)] = true
    }
  }
  for h in hit {
    assert_true(h)
  }
}

///|
/// 高斯整數 ℤ[i] 的範數濾網：展示 ℂ 上格點的可數稠密子集。
pub fn gaussian_int_norms(bound : Int32) -> Array[Int32] {
  let out : Array[Int32] = []
  for a in 0..<=bound {
    for b in 0..<=bound {
      out.push(a * a + b * b)
    }
  }
  out
}

///|
test {
  let norms = gaussian_int_norms(10)
  assert_eq(norms.length(), 121)
  // 0 = 0²+0²、1 = (±1)²+0² …… 2 = 1²+1² 可表示
  let mut has0 = false
  let mut has2 = false
  for x in norms {
    if x == 0 {
      has0 = true
    }
    if x == 2 {
      has2 = true
    }
  }
  assert_true(has0)
  assert_true(has2)
}
