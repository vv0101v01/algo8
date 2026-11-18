///|
/// 共用數值工具：近似判斷與邏輯斯迴歸閘函數。
pub fn approx(a : Double, b : Double, eps : Double) -> Bool {
  (a - b).abs() < eps
}

///|
