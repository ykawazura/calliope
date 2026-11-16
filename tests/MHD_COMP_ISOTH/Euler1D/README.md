# Euler1D: 1D 等温 Euler の衝撃波と滑らかな解

論文の本文に載せる 1D 等温 Euler（B = 0）の試験の入力を置く。
結果と判定は `../stages/REPORT_stage1.md`（Stage 1）にある。
入力は ps1 の `runs/stage1_es/` から写した（`runs/` は `tests/MHD_COMP_ISOTH/runs/`。`../README.md` の「`runs/` の場所」の節）。
各ケースは `runs/` の下のパスのまま置いてある（例: `stage1_es/shock/U2_N1024/`）。

## 最終版との照合

この試験は Stage 1 の実行ファイル `bin/calliope_isoth_es`（commit 2126274 のソース）で走らせた。
最終版 `bin/calliope_isoth_s13d`（07716cc）で同じ結果になることを、2026-11-17 に `check_s13d.py` で確かめた（job 497715）。
判定の基準は実行の前に `../../../tasks/MHD_COMP_ISOTH/todo.md` に書き、`check_s13d.py` の docstring に写した。

- 対象は 1 本の実行で終わる 99 ケースである。restart でつないだ 3 ケース（`shock/U4_N4096`、`dtconv/U4_N2048_cflsvv1`、
  `smoke/collide_U2_N512_est_chain`）は、つなぎ目の位置が実時間で決まり、行がそろわないので外した。
- 各ケースを es と s13d で、同じ `calliope.in`（md5 で原本と一致）、同じ rank 数、FFTW_ESTIMATE（`fftw_shim.so`）で走らせた。
  `calliope.in` に書いていない変数は s13d の既定値をとり、その既定値で s13d は es と同じ経路を通る。
- 99 対すべてで、終了状態（exit 0 が 92 対、exit 1 が 4 対、dt で止まった `dt_stop` が 3 対）がそろい、
  解のファイル（out1d、restart、shock_selftest.dat）がバイト単位で一致した。
- conservation.dat、cfl.dat、budget.dat の es が書く列（18、16、26 列）は、96 対で全行の文字列が一致した。
  budget.dat の列 10（eps_ideal。s13d は eps_eta も引く）も、この 96 対では一致した（収支のない `fast/` の 3 対を除く）。
- 残る 4 rank の 3 対（`decomp/U2_N512_x_np4`、`decomp/U2_N512_z_np4`、`smoke/selftest_np4`）では、budget.dat の大域和の列が
  丸め誤差の大きさで違った。差は、和の列（dEdt_sum、eps_shock、eps_svv、-D_svv）で列の最大値の 3.66e-16 以下、
  和の差の列（err_t、eps_ideal、R_zeta、R_svv）で絶対値 2.27e-13 以下である（和の列の最大値は 806）。
  残差の比の列（14、25、26）は、分子の残差が丸め誤差の大きさなので、相対では 0.167 まで変わる。
- この差の原因は MPI_Allreduce の足し算の順序である。s13d は収支の和を 1 回の Allreduce で集める配列を 16 個から 19 個に増やした。
  Allreduce を recursive doubling に固定して（`I_MPI_ADJUST_ALLREDUCE=1`、job 497729）走らせ直すと、3 対はすべての列で一致した。
  長さの変わらない conservation.dat の和（7 個）は、既定の Allreduce のままでも 99 対で一致した。
  同じ実行ファイルでも、既定と固定とでは conservation.dat と budget.dat の和の列が違い、解のファイルは一致した。
  2 rank の和は 1 回の足し算なので、順序によらない。
- 負の対照（入力の違う `scan/U4_N1024_C1` の es と `scan/U4_N1024_C4` の s13d）では差を検出した。
- 終了状態の判定は、実行の後で 1 か所直した。Stage 1 のジョブは dt で止まった 3 本にも `exit 0` を書いたが、今のジョブスクリプトは
  `dt_stop` を書く（3b784b1）。Stage 1 の calliope.out.std にも止まったことの表示があるので、期待する状態を今のスクリプトの規則で
  読み直した（`check_s13d.expected`）。
- 出力は ps1 の `runs/euler1d_s13d/`（`compare.json`、`judge.out`、再実行は `adjust/`）にある。
  費用は 19.33 core-h（48 コア × 23 分 59 秒と 24 コア × 22 秒）。

## 問題の設定

設定は REPORT_stage1 の 7 章の冒頭と 7.1 節、7.3 節による。

- 格子は N × 2 × 2 で、横方向を一様にして 1D の問題を解く（`pruned = F`。2/3 フィルタで横方向は k = 0 だけが残る）。
- c_s = 1（beta0 = 2）、lx = 1（箱の長さ L = 2π）で、境界は周期である。
- Riemann 問題は、右の状態を (L/4, 3L/4) に置いた周期問題で、界面は 2 枚ある。初期の界面は幅 1 格子の tanh でつなぐ。
  - Sod: ρ_L/ρ_R = 10、u = 0。t_end = 0.6。
  - U1、U2、U4: ρ = 1、u = ±U の衝突流。t_end は U1 と U2 が 0.8、U4 が 0.4。
  - どの t_end も、2 つの界面から出た波が出会う前である。
- 滑らかな解は、音波（振幅 1e-3、モード 1、1 周期、dt = 2π/4000）と、砕ける前の単純波である。
- dt は固定である。試作 `../stages/proto1d.py` を既定設定で走らせて min(dt_adv, dt_visc) の全ステップの最小値を問題と N ごとに求め、それに 0.7 を掛けた。

初期条件は calliope に組み込まれている（`init_type` が `riemann`、`sound_wave`、`simple_wave`、`single_mode`）。
生成するファイルがないので、ケースに `case.json` と `generated.md5` はなく、`regen_inputs.py` を当てる必要もない。

## batch とケース

| batch | 数 | 内容 | REPORT_stage1 |
|---|---|---|---|
| `shock` | 20 | Sod、U1、U2、U4 を N = 256、512、1024、2048、4096 で（既定設定） | 7.1 節 |
| `scan` | 31 | Sod、U2、U4 の N = 1024 で、係数と定式化を変える（C_shock、svv_eps0、SVV なし、`linear`、`filter_lnrho = T` など） | 7.2 節 |
| `smooth` | 19 | 音波（N = 32、64）と単純波（N = 32–512、C_shock = 0 と 2、SVV の有無） | 7.3 節 |
| `dtconv` | 5 | U2（N = 512）で dt を 4 通りに変えたものと、U4（N = 2048）の 1 本 | 7.4 節 |
| `cflmax` | 6 | U2（N = 512）で dt の倍率を 1.0–3.0 に変える | 7.4 節 |
| `decomp` | 9 | U2（N = 512）を x、y、z の各方向の分割と 1、2、4 ランクで | 7.4 節 |
| `smoke` | 9 | 単体試験（`selftest_np1/2/4`）、非有限値で止まること（`nan`）、音波（`sound64`）、U2 の衝突流と試作との照合、restart の試験 | 7.5 節 |
| `fast` | 3 | 16³ の 3D で B ≠ 0 の速進波の固有モード（`ideal`、`default`、`default_linear`） | 7.5 節 |

`fast` は 1D でも B = 0 でもないが、Stage 1 の同じ batch で走らせ、REPORT_stage1 の 7.5 節で扱っているので、ここに置いた。
REPORT_stage1 の 11.2 節にある `smoke/collide_U2_N512_chain` は、restart の試験を `FFTW_ESTIMATE` で作り直す前の版で、判定に使っていない。ps1 にも残っていないので、ここにもない。

## 走らせ方

手順は `../README.md` の「走らせ方」による。この試験に固有の点は次のとおりである。

- ケースに `exe` がないので、`qsub -v` の `EXE=` で実行ファイルを渡す。Stage 1 は `EXE=bin/calliope_isoth_es` で走らせた。
- restart の試験は、`smoke/collide_U2_N512_est` と、それを 1 ステップ目で止めて restart から続ける `smoke/collide_U2_N512_est_chain` を比べる。
  Stage 1 では後者に `terminate` ファイルを置いて投入した。
  calliope はこのファイルを見つけると消して止まるので、ps1 の元には残っておらず、ここにもない。
  走らせ直すときは、投入の前に `touch tests/MHD_COMP_ISOTH/runs/stage1_es/smoke/collide_U2_N512_est_chain/terminate` とする。
- ランクは 1、2、4 である。Stage 1 では、restart の試験の 2 本を除く全ケースを 1 つの一覧にまとめ、`PAR=64` で投入した。止まった U4 のケースは、残りの一覧で続けた（REPORT_stage1 の 11.2 節）。
- 判定は `python3 tests/MHD_COMP_ISOTH/stages/judge_stage1.py tests/MHD_COMP_ISOTH/runs/stage1_es` で、各ケースに `summary.json` を書き、判定表を出す。
- fast 波の判定、試作との照合（最初の 1000 ステップ）、SVV の冪乗法と dt の統計は、後処理 `post_stage1_es.sh` で求めた。
  これと fast 波の解析 `fast_wave.py` は ca9a205 で消したので、`git show cf9bb60:tests/MHD_COMP_ISOTH/<file> > <file>` で取り出して使う。
  `post_stage1_es.sh` は `T=tests/MHD_COMP_ISOTH` の下の `proto1d.py` と `fast_wave.py` を呼ぶので、取り出した 2 本は `stages/` に置き、`T` を `stages/` に直す。

## 図

`plot_euler1d.py` が `fig/<名前>.pdf` を書く。コマンドと fetch は `../README.md` の「図」の節による。
入力は、Stage 1 の判定表 `runs/stage1_es/judge.json`（`judge_stage1.py` が書く）と、`shock` の N = 512、2048 の t_end の out1d の行である。
判定表は、このフォルダの同じパス（`stage1_es/judge.json`）に置いた symlink からもたどれる（`../README.md` の「図」の節）。
`fetch` は、手元にないものを ps1 から写す（合わせて 2.7 MB）。
図の実行は es のものである。s13d との照合は「最終版との照合」の節による。

| 図 | 中身 | REPORT_stage1 |
|---|---|---|
| `euler_profiles` | Sod、U1、U2、U4 の t_end の ρ（対数軸）と u。N = 512、2048 と厳密解 | 7.1 節 |
| `euler_l1` | 4 問題の ρ の L1 誤差と N（256–4096）。N⁻¹ の線を添える | 7.1 節 |
| `euler_simple` | 砕ける前の単純波（t = t_b/2）の max\|u − u_exact\| と N（32–512）。C_shock = 0 の SVV なしとあり、C_shock = 2 の 3 通り。N⁻² の線を添える | 7.3 節 |
| `euler_dt` | U2（N = 512、t = 0.2）の質量の誤差 max\|M/M_0 − 1\| と、最小の dt の実行との ρ の差を、dt に対して描く。dt³ の線を添える | 7.4 節 |

- `euler_profiles` は、ρ の L1 誤差と中間状態の密度を `analyze_stage1.riemann_analysis` で計算し直し、`judge.json` の値と並べて標準出力に書く。
- 単純波の C_shock = 0 の誤差は、N = 256 で 2.1e-11 に達し、N = 512 でも下がらない。この値は時間積分の誤差で決まる。C_shock = 2 の誤差は N⁻² で下がる。
