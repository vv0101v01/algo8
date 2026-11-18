///|
/// 展示程式：以數學結構為計算機抽象，導出定理與演算法。
/// Run with `moon run cmd/main` from the project root.
fn main {
  println("=== 數學結構計算機抽象：定理與演算法展示 ===")
  println("")

  // 1. 長期記憶細胞狀態轉移
  println("[1] 長期記憶細胞狀態轉移")
  let g : @lib.Gates = { forget: 0.9, input: 0.8, output: 1.0, candidate: 0.5, }
  let ig = g.input * g.candidate
  let mut cell : @lib.MemoryCell = { c: 1.0, h: 0.0, }
  for _ in 0..<10 {
    cell = cell.step(g)
  }
  let closed = @lib.memory_closed_form(1.0, g.forget, ig, 10)
  println("  迭代 c₁₀ = \{cell.c}")
  println(
    "  閉式 c₁₀ = \{closed}   ← 定理：幾何級數轉移半群",
  )
  let pi = @lib.stationary_power_iteration(@lib.memory_ergodic_matrix(), 200)
  println(
    "  平穩分佈 π = [\{pi[0]}, \{pi[1]}, \{pi[2]}, \{pi[3]}]  ← 冪迭代演算法",
  )
  let t_lin = @lib.absorption_time_linear_solve(
    [
      [0.5, 0.4, 0.05, 0.0],
      [0.1, 0.4, 0.4, 0.05],
      [0.05, 0.1, 0.7, 0.1],
      [0.2, 0.3, 0.3, 0.15],
    ],
    [1.0, 1.0, 1.0, 1.0],
  )
  println("  吸收時間 E[T] = \{t_lin[0]}   ← (I-Q)t=1 與級數法同值")
  println("")

  // 2. 泛性質函子
  println("[2] 泛性質函子")
  println(
    "  終物件唯一箭數 = \{@lib.count_terminal_arrows(7)}  ← 定理：1^{|n|} = 1",
  )
  let ok_prod = @lib.universal_property_product_witness(2, 2, 2)
  println(
    "  乘積泛性質見證（2×2 錐） = \{ok_prod}  ← 存在唯一因子分解",
  )
  let ok_fun = @lib.functor_laws_product_with_k(3, 2, 2)
  println("  函子律 F(id)=id, F(g∘f)=F(g)∘F(f) = \{ok_fun}")
  println("")

  // 3. 複數域強基數
  println("[3] 複數域強基數")
  let z = @lib.Complex::new(3.5, -2.25)
  let (a, b) = @lib.strong_basis_coords(z)
  let back = @lib.rebuild_from_coords(a, b)
  println("  z = \{z.re} + \{z.im}i，基底座標 (a, b) = (\{a}, \{b})")
  println(
    "  重構 a·1 + b·i = (\{back.re}, \{back.im})  ← 定理：座標唯一（Hamel 強基數）",
  )
  println(
    "  Cantor 配對 π(3, 5) = \{@lib.cantor_pair(3, 5)}  ← ℕ² ≈ ℕ 等勢雙射",
  )
  let (pa, pb) = @lib.cantor_unpair(@lib.cantor_pair(3, 5))
  println("  反演 π⁻¹(π(3,5)) = (\{pa}, \{pb})")
  println("")

  // 4. 協同學函子神經網絡
  println("[4] 協同學函子神經網絡")
  let enc : @lib.Layer = {
    w: [[1.0, 0.5, 0.0], [0.0, 1.0, 0.5], [0.5, 0.0, 1.0]],
    b: [0.1, -0.2, 0.0],
  }
  let net : @lib.SynergeticNet = { encoder: enc, n_order: 3, }
  let x : Array[Double] = [1.0, 2.0, 3.0]
  let u1 : Array[Double] = [0.5, -0.5, 1.0]
  let u2 : Array[Double] = [-1.0, 0.25, 0.0]
  let fl = @lib.functor_laws_synergetic(net, x, u1, u2)
  println("  函子律（編碼層推前） = \{fl}")
  let s_star = @lib.slaving_fixed_point(1.5, 4.0)
  let s_n = @lib.slaving_run(1.5, 4.0, 0.2, 200)
  println(
    "  slaving 不動點 s* = \{s_star}，迭代 s₂₀₀ = \{s_n}  ← 指數收斂",
  )
  let w = @lib.order_parameter_wta([0.1, 2.0, 0.5, -1.0], 0.1, 200)
  println("  序參數 WTA 贏家索引 = \{w[0]}  ← Winner-Take-All 收斂")
  println("")

  // 5. 神經偏微分辛相調恆等式方程
  println("[5] 神經偏微分辛相調恆等式方程")
  let h = @lib.default_neural_h()
  let q : Array[Double] = [0.7, -0.4]
  let p : Array[Double] = [0.3, 0.9]
  let r = @lib.symplectic_identity_residual(h, q, p, 0.1)
  let d = @lib.symplectic_det_area(q, p, 0.1)
  println(
    "  MᵀΩM − Ω 殘差 = \{r}  ← 定理：辛相調恆等式（精確成立）",
  )
  println(
    "  Jacobi 行列式 det = \{d}   ← 面積守恆（離散 Liouville）",
  )
  let ds = @lib.energy_drift(true, 2000, 0.05)
  let de = @lib.energy_drift(false, 2000, 0.05)
  println(
    "  2000 步能量漂移：辛 = \{ds}，顯式 = \{de}  ← 幾何積分優勢",
  )
  println("")

  // 6. 湧現算子叢
  println("[6] 湧現算子叢")
  let gab : Array[Double] = [2.0, 0.0, 0.0, 0.5]
  let gbc : Array[Double] = [3.0, 0.0, 0.0, 1.0 / 3.0]
  let gac : Array[Double] = [6.0, 0.0, 0.0, 1.0 / 6.0]
  let coc = @lib.cocycle_condition(gab, gbc, gac, 1.0e-12)
  println("  轉移函數上邊條件 g_αβ g_βγ = g_αγ = \{coc}")
  let aa : Array[Double] = [1.0, 2.0, 0.0, -1.0]
  let bb : Array[Double] = [0.0, 1.0, 1.0, 0.0]
  let gg : Array[Double] = [1.0, 1.0, 0.0, 1.0]
  let inv = @lib.emergent_commutator_trace(aa, bb, gg, 1.0e-10)
  println("  交換子跡共軛不變（湧現不變量） = \{inv}")
  let gap = @lib.emergent_spectral_gap([4.0, 0.0, 0.0, 1.0])
  println("  拼接全域算子譜隙 = \{gap}  ← 局部拼湧現整體譜")
  println("")

  // 7. 希爾伯特 4 元數與非交換圖 K-理論譜對象
  println("[7] 希爾伯特4元數與非交換圖K-理論譜對象恆等式")
  let i = @lib.Quat::new(0.0, 1.0, 0.0, 0.0)
  let j = @lib.Quat::new(0.0, 0.0, 1.0, 0.0)
  let k = i.mul(j)
  println(
    "  ij = \{k.w} + \{k.x}i + \{k.y}j + \{k.z}k   ← 非交換：ji = −k",
  )
  let hh : @lib.Hermitian2 = {
    a: 2.0,
    b: 1.0,
    q: @lib.Quat::new(0.0, 0.5, -0.5, 1.0),
  }
  let (lp, lm) = @lib.hermitian2_eigs(hh)
  println(
    "  厄米 2×2 譜：λ₊ = \{lp}, λ₋ = \{lm}  ← 自伴算子實譜",
  )
  let k0 = @lib.k0_spectral_identity(hh, 1.0e-12)
  println(
    "  K₀ 譜對象恆等式（正交投影和=1、指標可加） = \{k0}",
  )
  println("")

  // 8. 仿宇宙
  println("[8] 仿宇宙")
  let bth = @lib.FourVal::Bv
  let g1 = @lib.FourVal::land(bth, @lib.FourVal::neg(bth))
  println(
    "  B ∧ ¬B = \{@lib.FourVal::to_string(g1)}  ← 矛盾不爆炸（dialetheia 安全）",
  )
  let expl = @lib.explosion_counterexample_exists()
  println("  爆炸原理反例存在 = \{expl}  ← 副相容語義")
  let w1 = @lib.ParaUniverse::new(["p", "q"], [
    @lib.FourVal::Tv,
    @lib.FourVal::Fv,
  ])
  let w2 = @lib.ParaUniverse::new(["p", "q"], [
    @lib.FourVal::Fv,
    @lib.FourVal::Tv,
  ])
  let mw : @lib.MultiWorld = { worlds: [w1, w2], }
  let merged = mw.merge_pointwise()
  println(
    "  多世界逐點融合 p: \{@lib.FourVal::to_string(merged[0])}, q: \{@lib.FourVal::to_string(merged[1])}",
  )
  let tarski = @lib.all_monotone_maps_have_fixed_point()
  println("  Tarski 不動點（全體單調映射） = \{tarski}")
  println("")
  println("=== 全部定理於測試中機械驗證：moon test ===")
}
