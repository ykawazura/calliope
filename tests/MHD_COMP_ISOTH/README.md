# tests/MHD_COMP_ISOTH

等温圧縮性 MHD のモデル `src/model/MHD_COMP_ISOTH` の検証に使うスクリプト、入力、報告を置く。
2026-11-11 に、論文で使う Stage のものだけを残すよう整理した。整理の前の状態は cf9bb60 にある。
計算の出力は git の管理外で、ps1 の `~/work/calliope_isoth/runs/` と手元の `runs/` に置く。
各スクリプトの使い方は、その docstring に書いてある。

## 論文での位置

論文に載せる試験と Stage の対応を次の表に示す。
分類は、検討の記録（worktree の `chat.tex`、git の管理外）の「論文に乗せるテストはどれ？」への答えに従い、その後に行った Stage 16 を加えた。

| 論文での位置 | 試験 | Stage |
|---|---|---|
| 本文 | 1D 等温 Euler の衝撃波と単純波 | 1 |
| 本文 | 1D 等温 MHD の Riemann 問題 | 3、14、16（conv） |
| 本文 | 多次元の滑らかな解 | 4 |
| 本文 | 2D Orszag–Tang 渦 | 4、14 |
| 本文 | 3D 減衰乱流（頑健さ、収支、スペクトル、密度 PDF、ESf の停止） | 6、7、15、16（K） |
| 付録 | MPI の分割によらないこと | 6、15 |
| 付録 | 係数の感度 | 8、15（ES25） |
| 限界 | M_s ≥ 2 | 4、5 |
| 限界 | B_t が反転する複合問題 C2 | 3 |
| 限界 | 質量の項を外した ES-SVV（svv_mass = F） | 13 |

`chat.tex` の分類は、1D 等温 Euler に Stage 14 も挙げている。
しかし Stage 14 で走らせ直したのは Stage 3 の 1D MHD と Stage 4 の OT だけで、Stage 1 の Euler の試験は最終版（s13d）で走らせ直していない。
Stage 7 のうち論文に載せるのは密度 PDF の裾で、密度の穴の機構は載せない。
Stage 6 と Stage 13 の batch `regress`（Fortran の変更の前後のビット一致を見る回帰試験）も、論文には載せない。

## Stage ごとのファイル

| Stage | 内容 | 生成 | 解析 | 報告 |
|---|---|---|---|---|
| 1 | 1D 等温 Euler（B = 0） | `gen_stage1.py` | `analyze_stage1.py`、`judge_stage1.py` | `REPORT_stage1.md` |
| 3 | 1D 等温 MHD | `stage3.py`（試作の実行）、`gen_stage3.py` | `analyze_stage3.py`、`stage3_table.py`、`stage3_aux.py` | `REPORT_stage3.md` |
| 4 | 2D と 3D の等温 MHD | `gen_stage4.py` | `analyze_stage4.py`、`cp_stability.py` | `REPORT_stage4.md` |
| 5 | M_s ≈ 4 の密度の穴、128³ の解像度試験、64³ での M_s の上限 | `gen_stage5.py` | `analyze_stage5.py`、`svv3d.py` | `REPORT_stage5.md` |
| 6 | dt を dt_cfl に追従させる 'track' と、turb1 の 256³ | `gen_stage6.py` | `analyze_stage6.py` | `REPORT_stage6.md` |
| 7 | turb1 の密度の穴の診断 | なし（既存の実行を読む） | `analyze_stage7.py` | `REPORT_stage7.md` |
| 8 | 密度の穴の散逸係数への依存（128³） | `gen_stage8.py` | `analyze_stage8.py` | `REPORT_stage8.md` |
| 13 | svv_mass = F（128³） | `gen_stage13.py` | `analyze_stage13.py` | `REPORT_stage13.md` |
| 14 | Stage 3 の 1D と Stage 4 の OT の、最終版での再実行 | `gen_stage14.py` | `analyze_stage14.py` | `REPORT_stage14.md` |
| 15 | 3D 減衰乱流での ES-SVV と標準の散逸の比較 | `gen_stage15.py` | `analyze_stage15.py` | `REPORT_stage15.md` |
| 16 | conv の再実行と、時刻をそろえた K の判定 | `gen_stage16.py` | `analyze_stage16.py`、`check_report_stage16.py` | `REPORT_stage16.md` |

## 共通の部品

- `isoth_case.py`: 試験のケースの `calliope.in` を作る。
- `gen_ic.py`: restart の形式で初期条件を書き、出力を読む。
- `riemann_exact.py`: 1D 等温 Euler の Riemann 問題の厳密解。
- `mhd_exact.py`: 1D 等温 MHD の厳密解（Stage 3）。
- `fv_mhd1d.py`: 1D 等温 MHD の 2 次精度の有限体積法。Stage 3 の参照解を作る。
- `proto1d.py`: 離散化の 1D の numpy 試作（B = 0）。Stage 2 で比べた定式化の変種も、このファイルに残っている。
- `proto_mhd1d.py`: B を含む 1D の試作（Stage 3）。
- `svv3d.py`: 3D の状態での ES-SVV。Stage 5 の静的な解析のために作り、Stage 7 と Stage 13 も使う。
- `check_seams.py`: 実行の継ぎ目（再開）をまたぐエネルギー収支を調べる（Stage 4–6）。
- `job.pbs-ps-fes-dev`: ps1 のジョブスクリプト。
- `fftw_shim.c`: p3dfft が FFTW に渡す計画の flag を stderr に書き出し、`FFTW_SHIM_ESTIMATE=1` のときは FFTW_ESTIMATE に置き換える LD_PRELOAD の shim。FFTW_MEASURE では同じ入力の 2 回の実行が 1e-13 程度ずれるので、ビット一致の試験に使う。
- `stage3_nstep.json`: Stage 3 の試作の実行が t_end に届いたときの step 数と dt。`gen_stage3.py nstep` が書き、Stage 3 の Fortran のケースはこの dt で走る。
- `incomp_regress.in`: MHD_INCOMP の回帰試験の入力。`gen_stage6.py` と `gen_stage13.py` の batch `regress` が使う。

## 整理で消したファイル

次の 19 本を消した。どれも cf9bb60 から取り出せる。

```
git show cf9bb60:tests/MHD_COMP_ISOTH/<file> > <file>
```

| 群 | ファイル | 内容 |
|---|---|---|
| どこからも使われない診断 | `diag_u4.py`、`fast_wave.py`、`void_probe.py`、`post_stage1_es.sh`、`stage0/fast_wave.in`、`stage0/fast_wave_new_default.in`、`stage0/fast_wave_new_ideal.in` | U4 の衝突流（Stage 1）と 3D の実行の密度の最小点の診断、Stage 1 の ps1 での後処理、Stage 0 の速進波の固有モードの入力と解析 |
| Stage 2 | `stage2.py`、`stage2_table.py`、`REPORT_stage2.md` | 数値定式化の比較 |
| Stage 9–12 | `analyze_stage9.py`、`analyze_stage10.py`、`analyze_stage11.py`、`analyze_stage12.py`、`gen_stage10.py`、`REPORT_stage9.md`–`REPORT_stage12.md` | 密度の穴の質量収支と、ES-SVV の項の svv_eps0 への依存 |

残したスクリプトが使っていた定義は、消す前に次のように移した。
移した定義の AST は元と同じである。
ただし `check_control` からは、Stage 10 の対照に固有の既定値 `allowed=ALLOWED` を外し、`allowed` を必ず渡す引数にした。

| 定義 | 元（消した） | 移した先 |
|---|---|---|
| `shells`、`spectrum` | `analyze_stage12.py` | `analyze_stage13.py` |
| `md5s`、`check_control` | `gen_stage10.py` | `gen_stage13.py` |
| `COLOR`、`MARK`、`GRAY`、`e`、`stats`、`status`、`gmax`、`order`、`table` | `stage2_table.py` | `stage3_table.py` |

`incomp_regress.in` は `tests/MHD_INCOMP_regression/calliope.in` から内容を変えずに移した。

残した REPORT と一部の docstring は、消したファイルや `tests/MHD_INCOMP_regression` に言及している。
これらは書いた時点の記録なので変えていない。
言及されたファイルは、上の `git show` で取り出せる。
ps1 の `~/work/calliope_isoth/tests` は、整理の前の配置のままである。
