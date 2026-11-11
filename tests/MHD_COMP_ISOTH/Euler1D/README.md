# Euler1D: 1D 等温 Euler の衝撃波と滑らかな解

論文の本文に載せる 1D 等温 Euler（B = 0）の試験の入力を置く。
結果と判定は `../stages/REPORT_stage1.md`（Stage 1）にある。
入力は ps1 の `~/work/calliope_isoth/runs/stage1_es/` から写した。
各ケースは `runs/` の下のパスのまま置いてある（例: `stage1_es/shock/U2_N1024/`）。

## 最終版で走らせ直していないこと

この試験は Stage 1 の実行ファイル `bin/calliope_isoth_es`（commit 2126274 のソース）で走らせ、最終版（`bin/calliope_isoth_s13d`）では走らせ直していない。
Stage 14 で走らせ直したのは、1D MHD と 2D OT だけである。
`calliope.in` に書いていない変数は、走らせる版の既定値をとる。
現在の src はこれらの `calliope.in` の変数をすべて宣言している（`../README.md`）が、現在の版で Stage 1 と同じ結果になるかは確かめていない。

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
  走らせ直すときは、投入の前に `touch runs/stage1_es/smoke/collide_U2_N512_est_chain/terminate` とする。
- ランクは 1、2、4 である。Stage 1 では、restart の試験の 2 本を除く全ケースを 1 つの一覧にまとめ、`PAR=64` で投入した。止まった U4 のケースは、残りの一覧で続けた（REPORT_stage1 の 11.2 節）。
- 判定は `python3 tests/MHD_COMP_ISOTH/stages/judge_stage1.py runs/stage1_es` で、各ケースに `summary.json` を書き、判定表を出す。
- fast 波の判定、試作との照合（最初の 1000 ステップ）、SVV の冪乗法と dt の統計は、後処理 `post_stage1_es.sh` で求めた。
  これと fast 波の解析 `fast_wave.py` は ca9a205 で消したので、`git show cf9bb60:tests/MHD_COMP_ISOTH/<file> > <file>` で取り出して使う。
  `post_stage1_es.sh` は `T=tests/MHD_COMP_ISOTH` の下の `proto1d.py` と `fast_wave.py` を呼ぶので、取り出した 2 本は `stages/` に置き、`T` を `stages/` に直す。
