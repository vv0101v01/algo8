# math-structures

以**數學結構為計算機抽象**的 MoonBit 專案：把 8 個數學名詞落實為型別與運算，
並由抽象過程導出 31 條定理與對應演算法（全部於 `moon test` 機械驗證）。

## 環境

- MoonBit 穩定版：`moon 0.1.20261104` / `moonc v0.10.20`（2026-11-05）
- 安裝：`curl -fsSL https://cli.moonbitlang.com/install/unix.sh | bash`

## 模組與導出

| 檔案 | 數學名詞 | 導出定理 |
|---|---|---|
| `memory.mbt` | 長期記憶細胞狀態轉移 | 轉移半群閉式解、平穩分佈唯一、吸收時間同值 |
| `universal.mbt` | 泛性質函子 | 終物件唯一至同構、乘積泛性質、函子律 |
| `complex_field.mbt` | 複數域強基數 | 座標唯一（Hamel 基）、域公理、Cantor 等勢 |
| `synergetic.mbt` | 協同學函子神經網絡 | 函子律、slaving 唯一不動點、WTA 穩定 |
| `symplectic.mbt` | 神經偏微分辛相調恆等式方程 | MᵀΩM=Ω 辛恆等式、面積守恆、能量穩定 |
| `bundle.mbt` | 湧現算子叢 | 轉移上邊條件、截面拼接、交換子跡不變 |
| `quaternion_ktheory.mbt` | 希爾伯特4元數與非交換圖K-理論譜對象恆等式 | 實譜、K₀ 譜投影恆等式、Cayley-Hamilton |
| `para_universe.mbt` | 仿宇宙 | 四值分配格、矛盾不爆炸、Tarski 不動點 |

完整定理陳述、證明思路與演算法複雜度見 [THEOREMS.md](THEOREMS.md)。

## 使用

```bash
moon test            # 31 個定理機械驗證
moon run cmd/main    # 展示全部定理的計算結果
moon check --deny-warn
moon fmt && moon info
```
