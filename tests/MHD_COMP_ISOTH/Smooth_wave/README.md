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
