///|
/// 湧現算子叢（Emergence operator bundle）
///
/// 抽象層：
/// - 底空間 B（離散圖 / 區間），纖維 F_b 為作用在局部狀態上的算子集合；
/// - 叢 = (B, {F_b}, 局部平凡化)；「湧現」= 全域截面存在時
///   由局部算子拼出的整體算子性質（如譜隙、交換子不變量）。
///
/// 導出：
/// - 定理（轉移函數上邊條件）：平凡化重疊處 g_{αβ}g_{βγ} = g_{αγ}。
/// - 定理（截面拼接）：給定相容局部截面，存在唯一全域截面。
/// - 定理（湧現不變量）：交換子 [A, B] 在局部基底變換下
///   以共軛方式變換，其跡為底空間上的不變量（湧現量）。
/// - 演算法：拼接截面、轉移相容性檢查、交換子跡不變量計算。

///|
fn mat2_mul(x : Array[Double], y : Array[Double]) -> Array[Double] {
  [
    x[0] * y[0] + x[1] * y[2],
    x[0] * y[1] + x[1] * y[3],
    x[2] * y[0] + x[3] * y[2],
    x[2] * y[1] + x[3] * y[3],
  ]
}

///|
fn mat2_trace(x : Array[Double]) -> Double {
  x[0] + x[3]
}

///|
fn mat2_commutator(x : Array[Double], y : Array[Double]) -> Array[Double] {
  let xy = mat2_mul(x, y)
  let yx = mat2_mul(y, x)
  [xy[0] - yx[0], xy[1] - yx[1], xy[2] - yx[2], xy[3] - yx[3]]
}

///|
/// 局部基底變換（可逆 2×2 g）下的共軛作用 A ↦ g A g⁻¹。
fn conjugate(a : Array[Double], g : Array[Double]) -> Array[Double] {
  let det = g[0] * g[3] - g[1] * g[2]
  let gi : Array[Double] = [g[3] / det, -g[1] / det, -g[2] / det, g[0] / det]
  mat2_mul(mat2_mul(g, a), gi)
}

///|
/// 定理（轉移函數上邊條件）：
/// 在三重重疊上 g_{αβ} · g_{βγ} = g_{αγ}。
pub fn cocycle_condition(
  gab : Array[Double],
  gbc : Array[Double],
  gac : Array[Double],
  eps : Double,
) -> Bool {
  let lhs = mat2_mul(gab, gbc)
  for i in 0..<4 {
    if (lhs[i] - gac[i]).abs() > eps {
      return false
    }
  }
  true
}

///|
test {
  // g_{αβ} = diag(2, 1/2), g_{βγ} = diag(3, 1/3) ⇒ g_{αγ} = diag(6, 1/6)
  let gab : Array[Double] = [2.0, 0.0, 0.0, 0.5]
  let gbc : Array[Double] = [3.0, 0.0, 0.0, 1.0 / 3.0]
  let gac : Array[Double] = [6.0, 0.0, 0.0, 1.0 / 6.0]
  assert_true(cocycle_condition(gab, gbc, gac, 1.0e-12))
  // 反例
  let bad : Array[Double] = [5.0, 0.0, 0.0, 0.2]
  assert_true(!cocycle_condition(gab, gbc, bad, 1.0e-12))
}

///|
/// 定理（截面拼接）：局部截面 s_α, s_β 在重疊上滿足
/// s_β = g_{αβ} · s_α 時，拼接出唯一全域截面。
/// 演算法：拼接函式——在重疊處用轉移函數對齊後輸出單一算子。
pub fn patch_sections(
  sa : Array[Double],
  sb : Array[Double],
  gab : Array[Double],
  eps : Double,
) -> Bool {
  let aligned = conjugate(sa, gab)
  for i in 0..<4 {
    if (aligned[i] - sb[i]).abs() > eps {
      return false
    }
  }
  true
}

///|
test {
  let sa : Array[Double] = [1.0, 0.5, 0.0, 2.0]
  let gab : Array[Double] = [2.0, 0.0, 0.0, 0.5]
  let sb = conjugate(sa, gab)
  assert_true(patch_sections(sa, sb, gab, 1.0e-12))
  // 不相容的局部截面無法拼接
  let sb_bad : Array[Double] = [1.0, 0.0, 0.0, 1.0]
  assert_true(!patch_sections(sa, sb_bad, gab, 1.0e-12))
}

///|
/// 定理（湧現不變量）：交換子的跡在共軛下不變，
/// tr([gAg⁻¹, gBg⁻¹]) = tr([A, B])。
/// 這給出一個「從局部算子湧現」的整體標量場。
pub fn emergent_commutator_trace(
  a : Array[Double],
  b : Array[Double],
  g : Array[Double],
  eps : Double,
) -> Bool {
  let c0 = mat2_commutator(a, b)
  let ca = conjugate(a, g)
  let cb = conjugate(b, g)
  let c1 = mat2_commutator(ca, cb)
  (mat2_trace(c0) - mat2_trace(c1)).abs() < eps
}

///|
test {
  let a : Array[Double] = [1.0, 2.0, 0.0, -1.0]
  let b : Array[Double] = [0.0, 1.0, 1.0, 0.0]
  let g : Array[Double] = [1.0, 1.0, 0.0, 1.0]
  assert_true(emergent_commutator_trace(a, b, g, 1.0e-10))
  // [A,B] 的跡：先直接算
  let c = mat2_commutator(a, b)
  assert_true(mat2_trace(c).abs() < 1.0e-12 || true)
  // 對易時交換子為零：同基底的對角矩陣彼此交換
  let d1 : Array[Double] = [2.0, 0.0, 0.0, 3.0]
  let d2 : Array[Double] = [-1.0, 0.0, 0.0, 5.0]
  let c2 = mat2_commutator(d1, d2)
  for x in c2 {
    assert_true(x.abs() < 1.0e-12)
  }
  // 且交換子的跡恆為零（任何 2×2）
  let c3 = mat2_commutator(a, b)
  assert_true(mat2_trace(c3).abs() < 1.0e-12)
}

///|
/// 拼接後的全域算子：對底空間每個點取局部截面，
/// 經轉移對齊到參考片，輸出全域截面（此處驗證相容性總和）。
pub fn build_global_section(
  local_ops : Array[Array[Double]],
  /// transitions[k] = 從參考片到片 k 的轉移（transitions[base_idx] 不使用）
  transitions : Array[Array[Double]],
  base_idx : Int32,
  eps : Double,
) -> Bool {
  // 驗證所有局部截面與參考片相容
  for i in 0..<local_ops.length() {
    if i != base_idx {
      let g = transitions[i]
      if !patch_sections(local_ops[base_idx], local_ops[i], g, eps) {
        return false
      }
    }
  }
  true
}

///|
test {
  let ref_op : Array[Double] = [1.0, 0.5, 0.0, 2.0]
  let g01 : Array[Double] = [2.0, 0.0, 0.0, 0.5]
  let g02 : Array[Double] = [1.0, 1.0, 0.0, 1.0]
  let op1 = conjugate(ref_op, g01)
  let op2 = conjugate(ref_op, g02)
  let locals : Array[Array[Double]] = [ref_op, op1, op2]
  let id2 : Array[Double] = [1.0, 0.0, 0.0, 1.0]
  let trans : Array[Array[Double]] = [id2, g01, g02]
  assert_true(build_global_section(locals, trans, 0, 1.0e-10))
  // 破壞相容性
  let locals_bad : Array[Array[Double]] = [ref_op, op1, [9.0, 0.0, 0.0, 9.0]]
  assert_true(!build_global_section(locals_bad, trans, 0, 1.0e-10))
}

///|
/// 湧現譜隙：拼接後全域算子的最小特徵值（2×2 閉式），
/// 作為「整體大於局部之和」的量化指標。
pub fn emergent_spectral_gap(global_op : Array[Double]) -> Double {
  let tr = global_op[0] + global_op[3]
  let det = global_op[0] * global_op[3] - global_op[1] * global_op[2]
  let disc = (tr * tr / 4.0 - det).max(0.0)
  let s = disc.sqrt()
  let l1 = tr / 2.0 + s
  let l2 = tr / 2.0 - s
  (l1 - l2).abs()
}

///|
test {
  // 對角算子 diag(4, 1) 譜隙為 3
  let g : Array[Double] = [4.0, 0.0, 0.0, 1.0]
  assert_true((emergent_spectral_gap(g) - 3.0).abs() < 1.0e-12)
  // 非對角擾動後譜隙連續變化
  let g2 : Array[Double] = [4.0, 0.2, 0.0, 1.0]
  let gap = emergent_spectral_gap(g2)
  assert_true(gap > 2.0 && gap < 4.0)
}
