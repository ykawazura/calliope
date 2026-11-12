# Smooth_wave: 円偏光 Alfvén 波と斜めの平面問題

論文の本文に載せる 2 つの多次元の試験、有限振幅の円偏光 Alfvén 波（CP 波）と斜めの平面問題の入力を置く。
どちらも、解が分かっている問題で 2D と 3D の演算子を確かめる試験である。
CP 波は厳密解と、斜めの平面問題は同じ問題の 1D の実行と比べる。
平面問題は衝撃波などの不連続を含み、滑らかな波ではないが、CP 波と同じ Stage 4 の試験なのでここに置いた。

結果と判定は `../stages/REPORT_stage4.md`（Stage 4）にある。
入力は ps1 の `~/work/calliope_isoth/runs/stage4/{alfven,cpdiag,planar}/` から写し、`runs/` の下のパスのまま置いた。

## 最終版で走らせ直していないこと

この試験は Stage 4 の実行ファイル `bin/calliope_isoth_s4` で走らせ、最終版（`bin/calliope_isoth_s13d`）では走らせ直していない。
Stage 14 で走らせ直したのは、1D MHD と 2D OT だけである。
`calliope.in` に書いていない変数は、走らせる版の既定値をとる。

## 多次元の実行の組み立て

REPORT_stage4 の 3 章による。

- 初期値は `../stages/gen_ic.py` が Fourier 係数として restart の形式で `ic/` に書き、calliope は `init_type = 'restart'` で読む。
  書く前に m と B に 2/3 マスクを掛け、B を div B = 0 に射影する。ln ρ は、滑らかな初期値だけにマスクを掛ける。
- 2D の実行はすべて x-z 面（nly = 2、`dim_decomp = 1`）で行う。
- dt は（問題、N）ごとに固定で、`reset_method = 'stop'` である。
- 平面問題の 1D の参照は、`profile1d.dat` から読む（`init_type = 'profile1d'`）。
- ケースには `ic/` と `profile1d.dat` を置かず、その md5 を `generated.md5` に書いた。

## 円偏光 Alfvén 波（`alfven`、`cpdiag`）

設定は REPORT_stage4 の 5.1 節による。

- 一様な ρ = 1、B0 = 1、c_s = 1 の上に、δB⊥/B0 = 0.1 と 1 の CP 波を置く。
- 波数は 2D（x-z 面）の k = (1,1)、(1,2) と、3D の k = (1,1,2) とその成分の入れ替え (1,2,1)、(2,1,1) である。
- N = 16–64（2D）、16–32（3D）で、1 周期と 10 周期を走らせる。dt = T/n_p（T = 2π/|k|）で、n_p は 200–3200 である。
- 係数をすべて 0 にした理想の実行（`ideal`）と、既定の設定（`default`。C_shock = 2、ES-SVV、C_eta = 2）を比べる。

ケース名は `<2d|3d>_k<波数>_a<01|1>_N<N>_<設定>_p<n_p>[_np<ランク数>]` で、`a01` が振幅 0.1、`a1` が 1 である。

| batch | 数 | 内容 | ランク | REPORT_stage4 |
|---|---|---|---|---|
| `alfven` | 63 | 2D の k = (1,1)、(1,2) と 3D の k = (1,1,2) を各振幅で。3D の (1,2,1)、(2,1,1) は振幅 1、N = 16 だけ | 1、4、8、16 | 5.1–5.3 節 |
| `cpdiag` | 16 | 振幅 1 の 4 つの波（2D の 2 つ、3D の N = 16 と 32）で、散逸の項を 1 つずつ入れる（`etaonly`、`noeta`、`nosvv`、`svvonly`） | 1、8 | 5.4 節 |

振幅 1 の 10 周期では丸め誤差が成長する。
REPORT_stage4 では、これを CP 波の線形不安定で説明した（5.3 節）。
線形安定性の計算は手元で `python3 tests/MHD_COMP_ISOTH/stages/cp_stability.py check` を走らせて求める。

## 斜めの平面問題（`planar`）

設定は REPORT_stage4 の 6.1 節と 6.2 節による。

- Stage 3 の 1D の分布 F（N 点、周期 L）を、x-z 格子の点 (i, l) に p(i, l) = mod(a i + b l, N) で写す。法線は n = (a, 0, b)/|k| で、実行の長さは t_end/|k| にする。
- 問題は Stage 3 の F1、C2、TDS_w1、RDh_w4、RDpc で（`../MHD1D/README.md`）、F1 と C2 は C_eta = 0 と 2（`_eta2`）の両方を走らせる。
- k = (1,2) では 2/3 マスクが残す帯域が半分になる。参照には、同じ分布を N/2 点に取った 1D の実行（`_1dh`）を加えた。

| 種類 | 数 | 内容 | 基準 |
|---|---|---|---|
| `ctl_*` | 8 | k = (1,0) の 2D の実行を、同じ dt の 1D の実行と比べる | 差 ≤ 1e-11 |
| `eqv_*` | 10 | k = (1,1)、SVV なし、200 step の 2D の実行を、縮めた 1D の実行（shock_nmax = shock_nsmooth = 2）と比べる | 差 ≤ 1e-10 |
| `phy_*` | 70 | 既定の設定で k = (1,1)（`_1d`、`_2d`）と k = (1,2)（`_1d`、`_1dh`、`_2d`）を N = 256、512 で走らせ、L1 誤差を 1D と比べる | 並進対称性の破れ ≤ 1e-12、恒等式の残差 ≤ 1e-10 など（6.1 節、6.3 節） |

ランクは 1、16、64 である。
88 本すべてが基準を満たした（REPORT_stage4 の 6.3 節）。

## 走らせ方

手順は `../README.md` の「走らせ方」による。この試験に固有の点は次のとおりである。

- ケースを ps1 の `runs/` の下の同じパスに写し、`python3 tests/MHD_COMP_ISOTH/stages/regen_inputs.py <case> ...` で `ic/` または `profile1d.dat` を書く。
- ケースに `exe` がないので、`qsub -v` の `EXE=` で実行ファイルを渡す。Stage 4 は `EXE=bin/calliope_isoth_s4` で走らせた。
- `alfven`、`cpdiag` の全ケースと、`planar` の `ctl_*`、`eqv_*` には `fftw_estimate` があり、FFT の計画を ESTIMATE に固定する。

## 解析

`python3 tests/MHD_COMP_ISOTH/stages/analyze_stage4.py <batch> [base]`（`<batch>` は `alfven`、`cpdiag`、`planar`。base の既定値は `runs/stage4`）で、判定と数値を標準出力と `<base>/<batch>/summary.json` に書く。
`planar` は、base の隣にある Stage 3 の 1D の実行（`runs/stage3_f`）も読むので、それがない環境では物理の比較の一部を使えない。

## 図

`plot_smooth_wave.py` が `fig/<名前>.pdf` を書く。コマンドと fetch は `../README.md` の「図」の節による。
入力は、上の解析が書いた `stage4/{alfven,cpdiag,planar}/summary.json` と `stage4/alfven/cp_rates.json`、平面問題の k = (1,2)、N = 256 の最後の分布（`runs/figdata/Smooth_wave/`）である。
`fetch` は、要約を ps1 から写し、分布を ps1 の出力から取り出す（合わせて 0.8 MB）。
C2 の FV の参照は、手元の `runs/stage3/fvref` から読む。
要約と FV の参照は、このフォルダの同じパス（`stage4/alfven/summary.json`、`stage3/fvref` など）に置いた symlink からもたどれる（`../README.md` の「図」の節）。
図の実行は、上に書いたとおり最終版で走らせ直していない。

成長率 γ は、ω0 = |k0| v_A を単位とする。

| 図 | 中身 | REPORT_stage4 |
|---|---|---|
| `cp_error` | (a) 振幅 0.1、ideal、N = 16 の max\|B − B_ex\|/A と周期の数 m。固有モードを RK3 で解いた解（破線）と、それとの差（白抜き、丸め誤差）を添える。(b) 振幅 1 の 2D k = (1,2)（N = 64）と 3D k = (1,1,2)（N = 32）の、ideal と default の RK3 の解との差と、ideal に当てはめた exp(2πγm) | 5.2、5.3 節 |
| `cp_dt` | default の 1 周期後の誤差と 1 周期のステップ数 n_p（振幅 0.1、1）。RK3 が 1 周期で失う振幅 1 − \|R\|^{n_p} ≈ (2π)⁴/(24 n_p³) を添える | 5.2 節 |
| `cp_growth` | 振幅 1、10 周期の γ と N（3 つの k0）。ideal（丸）、Galerkin の最大の成長率（破線）、default（四角。成長が測れないものは下端の三角）と、格子の上の物理の最大の成長率とその 1.25 倍（灰色。N ≥ 32 の default の上限） | 5.3、5.4 節 |
| `cp_diag` | 振幅 1、10 周期の γ を、同じ散逸を入れた Galerkin の支配的な族の成長率で割った比。散逸の項を 1 つずつ入れた `cpdiag` の実行 | 5.4 節 |
| `planar_profiles` | F1、C2、TDS_w1、RDpc の t_end の分布を法線に沿って描く。2D の k = (1,2)、N = 256 と、N/2 点（`_1dh`）と N 点（`_1d`）の 1D の実行に、厳密解（F1）か FV の参照（C2）を添える | 6.2、6.3 節 |

- `cp_dt` の誤差は時間積分の誤差だけで決まり、k0 にも振幅にもよらない。3D ではその 0.91 倍で、これは (1 − k_x²/|k0|²)^{1/2} = 0.913 にあたる（誤差は成分の最大で測るため）。
- N ≥ 32 の default で、γ が上限（格子の上の物理の最大の 1.25 倍）を超えたのは、3D の k = (1,1,2)、N = 32 だけである（γ = 0.132、物理の最大 0.098）。REPORT_stage4 の 5.4 節は、その原因を η_shock の項と記録し、判定を緩めていない。
- `planar_profiles` は、6 変数の L1 誤差を計算し直し、`planar/summary.json` の値と並べて標準出力に書く。
