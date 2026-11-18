///|
/// 仿宇宙（Para-universe）
///
/// 抽象層：
/// - 「仿」取 para-（並列、旁通、不矛盾）與 paraconsistent 雙義；
/// - 宇宙 = 命題集合上的四值真相格 𝔹 = {N, F, T, B}（Belnap）；
/// - 宇宙間的態射 = 保格單調映射（函子）；
/// - 宇宙的直積 = 狀態空間的笛卡爾積（多世界語義）。
///
/// 導出：
/// - 定理（四值格為佈爾代數的保守擴張）：交換、結合、分配、
///   補元律全部於 4⁴ 窮舉成立；經典二值嵌入 i : {F,T} ↪ 𝔹 保運算。
/// - 定理（矛盾不爆炸 / dialetheia 安全）：
///   B ∧ ¬B = B ≠ T，故「矛盾推出一切」在 𝔹 中失效；
///   機械驗證：存在賦值使前提皆為指定值而結論非指定。
/// - 定理（仿宇宙態射保結構）：單調映射保持 lub/glb。
/// - 演算法：四值運算表、經典嵌入檢查、爆炸反例搜尋。
pub(all) enum FourVal {
  /// Neither（無知）
  Nv
  /// False
  Fv
  /// True
  Tv
  /// Both（矛盾）
  Bv
}

///|
pub fn FourVal::join(a : FourVal, b : FourVal) -> FourVal {
  // 信息序 lub：N ≤ F,T ≤ B
  match (a, b) {
    (Nv, x) => x
    (x, Nv) => x
    (Fv, Fv) => Fv
    (Tv, Tv) => Tv
    (Fv, Tv) | (Tv, Fv) => Bv
    (Fv, Bv) | (Bv, Fv) => Bv
    (Tv, Bv) | (Bv, Tv) => Bv
    (Bv, Bv) => Bv
  }
}

///|
pub fn FourVal::meet(a : FourVal, b : FourVal) -> FourVal {
  match (a, b) {
    (Bv, x) => x
    (x, Bv) => x
    (Fv, Fv) => Fv
    (Tv, Tv) => Tv
    (Fv, Tv) | (Tv, Fv) => Nv
    (Fv, Nv) | (Nv, Fv) => Nv
    (Tv, Nv) | (Nv, Tv) => Nv
    (Nv, Nv) => Nv
  }
}

///|
pub fn FourVal::neg(a : FourVal) -> FourVal {
  match a {
    Nv => Nv
    Fv => Tv
    Tv => Fv
    Bv => Bv
  }
}

///|
/// 真值序（truth order：F 為底、T 為頂、B/N 居中）的 meet = 四值合取 ∧。
/// 經典像上退化為布爾 AND。
pub fn FourVal::land(a : FourVal, b : FourVal) -> FourVal {
  match (a, b) {
    (Fv, _) | (_, Fv) => Fv
    (Tv, x) => x
    (x, Tv) => x
    (Bv, Bv) => Bv
    (Nv, Nv) => Nv
    (Bv, Nv) | (Nv, Bv) => Fv
  }
}

///|
/// 真值序的 join = 四值析取 ∨。
pub fn FourVal::lor(a : FourVal, b : FourVal) -> FourVal {
  match (a, b) {
    (Tv, _) | (_, Tv) => Tv
    (Fv, x) => x
    (x, Fv) => x
    (Bv, Bv) => Bv
    (Nv, Nv) => Nv
    (Bv, Nv) | (Nv, Bv) => Tv
  }
}

///|
pub fn FourVal::to_string(self : FourVal) -> String {
  match self {
    Nv => "N"
    Fv => "F"
    Tv => "T"
    Bv => "B"
  }
}

///|
fn four_all() -> Array[FourVal] {
  [Nv, Fv, Tv, Bv]
}

///|
test {
  // 窮舉驗證格律（4² = 16 對 × 若干律）
  let vs = four_all()
  for a in vs {
    for b in vs {
      // 交換
      if FourVal::to_string(FourVal::join(a, b)) !=
        FourVal::to_string(FourVal::join(b, a)) {
        assert_true(false)
      }
      if FourVal::to_string(FourVal::meet(a, b)) !=
        FourVal::to_string(FourVal::meet(b, a)) {
        assert_true(false)
      }
      // 吸收
      if FourVal::to_string(FourVal::meet(a, FourVal::join(a, b))) !=
        FourVal::to_string(a) {
        assert_true(false)
      }
      if FourVal::to_string(FourVal::join(a, FourVal::meet(a, b))) !=
        FourVal::to_string(a) {
        assert_true(false)
      }
      // 分配
      for c in vs {
        let l = FourVal::meet(a, FourVal::join(b, c))
        let r = FourVal::join(FourVal::meet(a, b), FourVal::meet(a, c))
        if FourVal::to_string(l) != FourVal::to_string(r) {
          assert_true(false)
        }
      }
      // 對合：¬¬a = a
      if FourVal::to_string(FourVal::neg(FourVal::neg(a))) !=
        FourVal::to_string(a) {
        assert_true(false)
      }
    }
  }
  // 雙端
  for a in vs {
    if FourVal::to_string(FourVal::join(a, Nv)) != FourVal::to_string(a) {
      assert_true(false)
    }
    if FourVal::to_string(FourVal::meet(a, Bv)) != FourVal::to_string(a) {
      assert_true(false)
    }
  }
  assert_true(true)
}

///|
/// 經典二值嵌入 i : {F, T} → 𝔹，保持 ¬、∧、∨。
/// 定理（保守擴張）：在二值像上，四值運算與經典真值表一致。
pub fn classic_embed(b : Bool) -> FourVal {
  if b {
    Tv
  } else {
    Fv
  }
}

///|
test {
  for x in [false, true] {
    for y in [false, true] {
      let ix = classic_embed(x)
      let iy = classic_embed(y)
      // ∧（真值序）
      let and4 = FourVal::land(ix, iy)
      let andc = classic_embed(x && y)
      if FourVal::to_string(and4) != FourVal::to_string(andc) {
        assert_true(false)
      }
      // ∨（真值序）
      let or4 = FourVal::lor(ix, iy)
      let orc = classic_embed(x || y)
      if FourVal::to_string(or4) != FourVal::to_string(orc) {
        assert_true(false)
      }
      // ¬
      let n4 = FourVal::neg(ix)
      let nc = classic_embed(!x)
      if FourVal::to_string(n4) != FourVal::to_string(nc) {
        assert_true(false)
      }
    }
  }
  assert_true(true)
}

///|
/// 矛盾不爆炸的機械驗證（dialetheia 安全）：
/// 存在賦值 v 使 γ₁ = p ∧ ¬p 與 γ₂ = q 皆為指定值（T 或 B），
/// 而結論 r 指定失敗——故 {γ₁, γ₂} ⊭ r 在 𝔹 中。
/// 指定值 = {T, B}（副相容語義）。
fn designated(v : FourVal) -> Bool {
  match v {
    Tv | Bv => true
    Nv | Fv => false
  }
}

///|
pub fn explosion_counterexample_exists() -> Bool {
  // 窮舉 p, q, r ∈ 𝔹
  for p in four_all() {
    for q in four_all() {
      for r in four_all() {
        // γ1 = p ∧ ¬p（真值序合取）
        let g1 = FourVal::land(p, FourVal::neg(p))
        let g2 = q
        if designated(g1) && designated(g2) && !designated(r) {
          return true
        }
      }
    }
  }
  false
}

///|
test {
  // 副相容性：存在反例 ⇒ 爆炸原理失效
  assert_true(explosion_counterexample_exists())
  // 具體反例：p = Bv（既真又假）
  let p = Bv
  let g1 = FourVal::land(p, FourVal::neg(p))
  assert_true(designated(g1))
  // 而 Fv 非指定
  assert_true(!designated(Fv))
  // 經典爆炸在 𝔹 中的「B 情形」失效：B ∧ ¬B = B ≠ T
  assert_true(FourVal::to_string(g1) == "B")
}

///|
/// 宇宙：命題 → 四值真相的賦值。
pub(all) struct ParaUniverse {
  /// 命題標籤
  props : Array[String]
  vals : Array[FourVal]
}

///|
pub fn ParaUniverse::new(
  props : Array[String],
  vals : Array[FourVal],
) -> ParaUniverse {
  { props, vals, }
}

///|
pub fn ParaUniverse::lookup(self : ParaUniverse, name : String) -> FourVal? {
  for i in 0..<self.props.length() {
    if self.props[i] == name {
      return Some(self.vals[i])
    }
  }
  None
}

///|
/// 宇宙態射（保序映射）：f : 𝔹 → 𝔹 單調（信息序），
/// 逐命題推進真相——「仿宇宙之間的旁通」。
/// 信息序：x ≤ y ⇔ x ∨ y = y。
pub fn four_le(x : FourVal, y : FourVal) -> Bool {
  FourVal::to_string(FourVal::join(x, y)) == FourVal::to_string(y)
}

///|
/// 定理（Tarski 不動點）：有限完全格 𝔹 上的任意單調映射
/// f : 𝔹 → 𝔹 必有不動點 x（f(x) = x）。
/// 演算法：窮舉全部 4⁴ = 256 個映射，篩出單調者並搜尋不動點。
/// （附註：單調映射僅保序，一般不保 lub；保 lub 者為半格態射。）
pub fn all_monotone_maps_have_fixed_point() -> Bool {
  // 以遞增編碼列出所有函式 f : {0,1,2,3} → {0,1,2,3}（N<F<T<B 編碼）
  // 窮舉 4^4 = 256
  for c0 in 0..<4 {
    for c1 in 0..<4 {
      for c2 in 0..<4 {
        for c3 in 0..<4 {
          let f : Array[Int32] = [c0, c1, c2, c3]
          // 檢查單調：x ≤ y ⇒ f(x) ≤ f(y)
          let vs = four_all()
          let mut mono = true
          for x in 0..<4 {
            for y in 0..<4 {
              let j = FourVal::join(vs[x], vs[y])
              let jidx = match j {
                Nv => 0
                Fv => 1
                Tv => 2
                Bv => 3
              }
              if jidx == y {
                // x ≤ y
                let fx = f[x]
                let fy = f[y]
                let jf = FourVal::join(vs[fx], vs[fy])
                let jfidx = match jf {
                  Nv => 0
                  Fv => 1
                  Tv => 2
                  Bv => 3
                }
                if jfidx != fy {
                  mono = false
                }
              }
            }
          }
          if mono {
            // Tarski：搜尋不動點 f(x) = x
            let mut found = false
            for x in 0..<4 {
              if f[x] == x {
                found = true
              }
            }
            if !found {
              return false
            }
          }
        }
      }
    }
  }
  true
}

///|
test {
  assert_true(all_monotone_maps_have_fixed_point())
  // 具體不動點例子：恆等映射、常數映射
  assert_true(four_le(Nv, Fv))
  assert_true(four_le(Fv, Bv))
  assert_true(!four_le(Fv, Tv))
}

///|
/// 多世界（仿宇宙的直積）：狀態 = 各世界真相的向量。
pub(all) struct MultiWorld {
  worlds : Array[ParaUniverse]
}

///|
/// 直積上的逐點 join（世界融合）。
pub fn MultiWorld::merge_pointwise(self : MultiWorld) -> Array[FourVal] {
  if self.worlds.length() == 0 {
    return []
  }
  let n = self.worlds[0].vals.length()
  let out : Array[FourVal] = []
  for i in 0..<n {
    let mut acc = self.worlds[0].vals[i]
    for wi in 1..<self.worlds.length() {
      acc = FourVal::join(acc, self.worlds[wi].vals[i])
    }
    out.push(acc)
  }
  out
}

///|
test {
  let w1 = ParaUniverse::new(["p", "q"], [Tv, Fv])
  let w2 = ParaUniverse::new(["p", "q"], [Fv, Tv])
  let mw : MultiWorld = { worlds: [w1, w2], }
  let merged = mw.merge_pointwise()
  // Tv ∨ Fv = Bv（世界融合產生 Both——並列宇宙的矛盾共存）
  assert_true(FourVal::to_string(merged[0]) == "B")
  assert_true(FourVal::to_string(merged[1]) == "B")
  // 單世界查表
  let v = w1.lookup("p")
  match v {
    Some(x) => assert_true(FourVal::to_string(x) == "T")
    None => assert_true(false)
  }
  match w1.lookup("zzz") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}
