# MHD1D: 1D MHD の Riemann 問題と収束の試験

論文の本文に載せる 1D MHD の Riemann 問題と、その収束の試験の入力を置く。
問題の組と判定の量は Stage 3（`../stages/REPORT_stage3.md`）で決めた。
論文の数値は、最終版の実行ファイル `bin/calliope_isoth_s13d` で走らせ直した Stage 14（`../stages/REPORT_stage14.md`）と Stage 16（`../stages/REPORT_stage16.md`）の実行から取る。
入力は ps1 の `~/work/calliope_isoth/runs/stage14/mhd1d/` と `runs/stage16/conv/` から写し、`runs/` の下のパスのまま置いた。

C2（B_t が反転する複合問題）は、論文では限界の例として扱う（REPORT_stage3 の 1 章）。
suite と eta_suite の `C2_N*` がそれにあたる。

## 問題の設定

設定は REPORT_stage3 の 4.1 節と 4.2 節による。

- 単位は c_s = 1（β0 = 2）、箱の長さ L = 1（lx = 1/(2π)）で、境界は周期である。格子点は x_j = jL/N に置く。
- 不連続は、とくに断らなければ格子点の上に H = 1/2 で置く。
- 問題は次の 32 通りで、W = (ρ, u_x, u_y, u_z, B_y, B_z) と書く。

| 族 | 問題 | 初期値 | B_x | t_end |
|---|---|---|---|---|
| 垂直 MHD | P1a、P1b、P1c | (ρ, u_x, B_y) = (1, ±U, 1) の衝突。U = 1、2、4 | 0 | 0.11、0.09、0.06 |
| | P2 | (10, 0, 1) : (1, 0, 1) | 0 | 0.09 |
| | P3 | (1, 0, 1) : (0.125, 0, −1) | 0 | 0.055 |
| 共面の単一衝撃波 | F1、F2（fast）、S1、S2（slow）、I1（2→3 の中間） | 上流 (1, 0, 0, 1) を (x_a, x_b) に、下流を外側に置く。x_a = 0.125 | 1 | 0.16、0.105、0.25、0.25、0.24 |
| 接線不連続 | TDS（静止）、TDA（u_x = 1 で一周） | B_y = √3 : 1、c²ρ + B_y²/2 = 2。tanh の幅 w = 0、1、4 格子。`pre` は B_y に 2/3 則のマスクを掛けてから ρ を決める | 0 | 1 |
| 回転不連続 | RDh（Δψ = π/2）、RDp（Δψ = π） | ρ = \|B_t\| = 1、u_x = 0、u_t = −B_t。ψ を幅 w = 0、4、8 格子で混ぜて一周させる | 1 | 1 |
| | RDpc | Δψ = π。跳びを格子点の間に置き、B_z = u_z = 0 を厳密に保つ（共面） | 1 | 1 |
| 複合 | C1 | (1, 0, 0, 1) : (0.125, 0, 0, 1)（B_t の反転なし） | 0.75 | 0.1 |
| | C2 | (1, 0, 0, 1) : (0.125, 0, 0, −1)（B_t の反転あり） | 0.75 | 0.1 |
| | C3 | (1, ±1, 0, 1) の衝突 | 1 | 0.1 |

- 垂直 MHD と複合問題は、右の状態を (L/4, 3L/4) に置き、2 つの界面が出す波を重ね合わせる。
- dt は固定で、dt = t_end/nstep とする。nstep は 20 の倍数で、dt は同じ設定の適応 dt の予備実行での min(dt_adv, dt_visc) の最小値の 0.7 倍以下である。値は `../stages/stage3_nstep.json` にあり、各ケースの `case.json` にも書いてある。dt が変わる事態では `reset_method = 'stop'` で止める。
- 初期値は `profile1d.dat` として読む（`init_type = 'profile1d'`）。ケースには `profile1d.dat` を置かず、その md5 を `generated.md5` に書いた。

入力は、Stage 14 と Stage 16 で生成スクリプトから作り直さず、Stage 3 のケース（`runs/stage3_f`）から写した。
作り直すと、Stage 4 以降に足した既定値で `calliope.in` の字面が変わりうるからである（REPORT_stage14 の 2 章）。

## batch とケース

| batch | 数 | 内容 | ランク | REPORT |
|---|---|---|---|---|
| `stage14/mhd1d/suite` | 131 | 32 問題 × N = 256、512、1024、2048 と、C1、C2、C3 の N = 4096。C_η = 0 | 1、2、4 | Stage 14 の 5.1 節（設定は Stage 3） |
| `stage14/mhd1d/eta_suite` | 51 | 13 問題の C_η = 2 | 1、2、4 | 同上 |
| `stage14/mhd1d/decomp` | 15 | C3、F1、RDh_w4 の N = 256 を、x 方向の 1、2、4 ランク、y 方向の 2 ランク、z 方向の 4 ランクで | 1、2、4 | Stage 14 の 5.1 節（D4） |
| `stage14/mhd1d/eta_decomp` | 10 | C3、I1 の N = 256、C_η = 2 を同じ 5 通りで | 1、2、4 | 同上 |
| `stage16/conv` | 8 | I1 と RDpc の N = 4096、8192 を C_η = 0、2 で | 4（N = 4096）、8（8192） | Stage 16 の 2.2 節と 3 章（設定は Stage 3 の 10.6 節） |

eta_suite の 13 問題は、P1b、P2、F1、F2、S1、S2、I1（各 N = 256–2048）、C1、C2、C3（各 256–4096）、RDh_w4、RDpc（各 256–1024）、TDS_w1（256、512）である。
Stage 3 は C_η = 0.5、8 でも走らせたが、Stage 14 で走らせ直したのは既定値の C_η = 2 だけである。

Stage 14 の判定（REPORT_stage14 の 3 章で事前に固定）は、D1（完走）、D2（試作との差 ≤ 1e-11）、D3（Stage 3 の Fortran の実行との差 ≤ 2e-11）、D4（分割による差 ≤ 1e-13）が合格した。
D5（Stage 3 の記録との比較）は不合格である。
原因は、実行の継ぎ目で書く行の違い、再開ごとにやり直す λ_svv の冪乗法、FFT の計画の違いで、コードの違いによる差は見つからなかった（5.2 節）。
Stage 16 の conv は、完走（CV1）し、Stage 3 の実行と全 21 行でビット単位で一致した（CV2）。

## 走らせ方

手順は `../README.md` の「走らせ方」による。この試験に固有の点は次のとおりである。

- ケースを ps1 の `runs/` の下の同じパスに写し、`python3 tests/MHD_COMP_ISOTH/stages/regen_inputs.py <case> ...` で `profile1d.dat` を書く。
- 各ケースの `exe` が `bin/calliope_isoth_s13d` を指すので、`EXE=` は要らない。
- conv の 8 本には `fftw_estimate` があり、FFT の計画を ESTIMATE に固定する。Stage 3 の実行と同じ条件にするためである。
- Stage 14 は、207 本の一覧を `PAR=128` で投げ、予算で止まったケースを `PAR=28` の 2 本目で続けた（REPORT_stage14 の 10 章）。
  Stage 16 の conv は、`PAR=48` で投げた後、残ったケースのランクの和に合わせて `PAR` を下げ、同じ一覧で投げ直した（REPORT_stage16 の 13 章）。

## 解析

- 表と図: `python3 tests/MHD_COMP_ISOTH/stages/analyze_stage3.py table runs/stage14/mhd1d/suite runs/stage14/mhd1d/suite_rec` で、各ケースの記録と `tables.md`、図を書く。eta_suite も同じである。
- 最後の状態: `analyze_stage3.py last runs/stage14/mhd1d/suite runs/stage14/mhd1d/suite_last`（eta_suite も同じ）。
- 分割: `analyze_stage3.py decomp runs/stage14/mhd1d/decomp`。
- conv: `MPLBACKEND=Agg python3 tests/MHD_COMP_ISOTH/stages/analyze_stage16.py conv` で完走と旧との一致を調べ、`analyze_stage16.py table` で N = 256–8192 の表を作る。
  後者は、N ≤ 2048 を `runs/stage14/mhd1d/*_last` から、4096 と 8192 を `runs/stage16/conv_last` から読む。
- 旧との比較（`analyze_stage14.py mhd1d`、`proto`、`records`）は、ps1 にある Stage 3 の実行（`runs/stage3_f`）と試作の実行（`runs/stage3`）を読むので、それらがない環境では使えない。
