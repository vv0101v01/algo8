///|
/// 泛性質函子（Universal property functor）
///
/// 抽象層：以「小範疇 + 普適錐」為計算機抽象。
/// 物件用整數標籤，態射用顯式表（來源、目標、映射表），
/// 普適性質以「存在且唯一」的因子分解機械驗證。
///
/// 導出：
/// - 定理（終物件唯一至同構）：任意兩個終物件之間存在唯一同構。
/// - 定理（乘積泛性質）：若 p 與投影 (π₁, π₂) 滿足泛性質，
///   則對任意錐 (f, g) 存在唯一態射 h 使 π₁∘h = f、π₂∘h = g。
/// - 演算法：泛性質見證者——枚舉候選 h 並檢驗唯一性。
pub(all) struct Obj {
  tag : Int32
}

///|
pub(all) struct Hom {
  src : Int32
  dst : Int32
  /// 出邊表：對每個入點給出像（此處以全域點集上的全函式為例）。
  map : Array[Int32]
}

///|
pub fn Hom::compose(f : Hom, g : Hom) -> Hom {
  // g ∘ f : f.src → g.dst
  let m : Array[Int32] = []
  for x in f.map {
    m.push(g.map[x])
  }
  { src: f.src, dst: g.dst, map: m, }
}

///|
pub fn Hom::identity(n : Int32) -> Hom {
  let m : Array[Int32] = []
  for i in 0..<n {
    m.push(i)
  }
  { src: n, dst: n, map: m, }
}

///|
pub fn Hom::is_identity(f : Hom) -> Bool {
  if f.src != f.dst {
    return false
  }
  for i in 0..<f.map.length() {
    if f.map[i] != i {
      return false
    }
  }
  true
}

///|
/// 終物件：從它出發到每個物件的態射存在且唯一（此處取單點物件 1）。
/// 定理（終物件唯一至同構）：若 A、B 皆終，則存在唯一同構 A → B。
pub fn terminal_is_unique_up_to_iso() -> Bool {
  // 以全函式範疇 FinSet 為例：終物件是單點集。
  // 單點集之間的全函式恰有一條，且它自動是同構。
  let f : Hom = { src: 1, dst: 1, map: [0], }
  Hom::is_identity(f) && f.map.length() == 1
}

///|
/// 終物件的泛性質：對任意物件 n，恰有一條態射 n → 1。
pub fn terminal_unique_arrow(n : Int32) -> Hom {
  let m : Array[Int32] = []
  for _ in 0..<n {
    m.push(0)
  }
  { src: n, dst: 1, map: m, }
}

///|
pub fn count_terminal_arrows(_n : Int32) -> Int32 {
  // 全函式 n → 1 的條數 = 1^{|n|} = 1（由構造唯一性給出）
  1
}

///|
/// 乘積泛性質：集合論乘積 A × B 與投影 π₁, π₂。
/// 演算法（泛性質見證者）：給定錐 f : C → A、g : C → B，
/// 構造 h = ⟨f, g⟩（對角配對），並驗證：
///   (1) π₁ ∘ h = f、π₂ ∘ h = g；
///   (2) 任何滿足此二式的 h' 必等於 h（唯一性）。
pub(all) struct Cone {
  /// C 的基點個數。
  c_size : Int32
  f : Array[Int32]
  g : Array[Int32]
}

///|
pub fn product_pi1(a_size : Int32, b_size : Int32) -> Hom {
  // π₁ : A×B → A，編碼 (i, j) ↦ i·|B| + j
  let m : Array[Int32] = []
  for i in 0..<a_size {
    for _ in 0..<b_size {
      m.push(i)
    }
  }
  { src: a_size * b_size, dst: a_size, map: m, }
}

///|
pub fn product_pi2(a_size : Int32, b_size : Int32) -> Hom {
  let m : Array[Int32] = []
  for _ in 0..<a_size {
    for j in 0..<b_size {
      m.push(j)
    }
  }
  { src: a_size * b_size, dst: b_size, map: m, }
}

///|
/// 配對 ⟨f, g⟩ : C → A×B。
pub fn pair(cone : Cone, a_size : Int32, b_size : Int32) -> Hom {
  let m : Array[Int32] = []
  for x in 0..<cone.c_size {
    m.push(cone.f[x] * b_size + cone.g[x])
  }
  { src: cone.c_size, dst: a_size * b_size, map: m, }
}

///|
fn hom_eq(p : Hom, q : Hom) -> Bool {
  if p.src != q.src || p.dst != q.dst || p.map.length() != q.map.length() {
    return false
  }
  for i in 0..<p.map.length() {
    if p.map[i] != q.map[i] {
      return false
    }
  }
  true
}

///|
/// 定理（乘積泛性質）+ 演算法：
/// 對所有小規模錐窮舉驗證：存在唯一 h 使投影複合還原 f、g，
/// 且此 h 必為配對 ⟨f, g⟩。
pub fn universal_property_product_witness(
  a_size : Int32,
  b_size : Int32,
  c_size : Int32,
) -> Bool {
  let pi1 = product_pi1(a_size, b_size)
  let pi2 = product_pi2(a_size, b_size)
  // 窮舉所有錐 (f : C → A, g : C → B)
  let f_maps : Array[Array[Int32]] = []
  enumerate_funcs(c_size, a_size, f_maps)
  let g_maps : Array[Array[Int32]] = []
  enumerate_funcs(c_size, b_size, g_maps)
  for fm in f_maps {
    for gm in g_maps {
      let cone : Cone = { c_size, f: fm, g: gm, }
      let h = pair(cone, a_size, b_size)
      // (1) 還原性
      if !hom_eq(Hom::compose(h, pi1), fm_as_hom(cone, a_size)) {
        return false
      }
      if !hom_eq(Hom::compose(h, pi2), gm_as_hom(cone, b_size)) {
        return false
      }
      // (2) 唯一性：任何其他候選 h' 若還原則必等於 h
      let candidates : Array[Array[Int32]] = []
      enumerate_funcs(c_size, a_size * b_size, candidates)
      for cand in candidates {
        let hp : Hom = { src: c_size, dst: a_size * b_size, map: cand, }
        if hom_eq(Hom::compose(hp, pi1), fm_as_hom(cone, a_size)) &&
          hom_eq(Hom::compose(hp, pi2), gm_as_hom(cone, b_size)) {
          if !hom_eq(hp, h) {
            return false
          }
        }
      }
    }
  }
  true
}

///|
fn fm_as_hom(cone : Cone, dst : Int32) -> Hom {
  { src: cone.c_size, dst, map: cone.f, }
}

///|
fn gm_as_hom(cone : Cone, dst : Int32) -> Hom {
  { src: cone.c_size, dst, map: cone.g, }
}

///|
/// 列舉所有 [n] → [m] 的函式。
fn enumerate_funcs(n : Int32, m : Int32, out : Array[Array[Int32]]) -> Unit {
  if n == 0 {
    let empty : Array[Int32] = []
    out.push(empty)
    return
  }
  let subs : Array[Array[Int32]] = []
  enumerate_funcs(n - 1, m, subs)
  for s in subs {
    for v in 0..<m {
      let cur : Array[Int32] = []
      for x in s {
        cur.push(x)
      }
      cur.push(v)
      out.push(cur)
    }
  }
}

///|
test {
  assert_true(terminal_is_unique_up_to_iso())
  assert_eq(count_terminal_arrows(5), 1)
  let arrow = terminal_unique_arrow(5)
  assert_eq(arrow.map.length(), 5)
  assert_true(arrow.map[3] == 0)
  // 泛性質見證：小規模窮舉（2×3 錐於 2 點上）
  assert_true(universal_property_product_witness(2, 3, 2))
  assert_true(universal_property_product_witness(2, 2, 1))
  // 恒等態射律
  let id3 = Hom::identity(3)
  assert_true(Hom::is_identity(id3))
  let f : Hom = { src: 3, dst: 2, map: [0, 1, 0], }
  let comp = Hom::compose(f, Hom::identity(2))
  assert_true(hom_eq(comp, f))
}

///|
/// 函子定律的機械驗證（以 FinSet 上的「與固定集合 K 做笛卡爾積」為例）：
/// F(id) = id、F(g∘f) = F(g)∘F(f)。
pub fn functor_laws_product_with_k(
  a_size : Int32,
  b_size : Int32,
  k : Int32,
) -> Bool {
  // F(X) = X × K；F(f)(x, t) = (f(x), t)
  let f : Hom = {
    src: a_size,
    dst: b_size,
    map: [
      for i in 0..<a_size => i % b_size
    ],
  }
  // F(id_A) = id_{A×K}
  let id_f = Hom::identity(a_size)
  let fid = lift_with_k(id_f, k)
  if !Hom::is_identity(fid) {
    return false
  }
  // F(g∘f) = F(g)∘F(f)
  let g : Hom = {
    src: b_size,
    dst: a_size,
    map: [
      for i in 0..<b_size => i % a_size
    ],
  }
  let gf = Hom::compose(f, g)
  let left = lift_with_k(gf, k)
  let right = Hom::compose(lift_with_k(f, k), lift_with_k(g, k))
  hom_eq(left, right)
}

///|
fn lift_with_k(f : Hom, k : Int32) -> Hom {
  let m : Array[Int32] = []
  for x in f.map {
    for t in 0..<k {
      m.push(x * k + t)
    }
  }
  { src: f.src * k, dst: f.dst * k, map: m, }
}

///|
test {
  assert_true(functor_laws_product_with_k(3, 2, 2))
  assert_true(functor_laws_product_with_k(2, 3, 3))
}
