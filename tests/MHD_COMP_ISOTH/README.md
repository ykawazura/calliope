# tests/MHD_COMP_ISOTH

等温圧縮性 MHD のモデル `src/model/MHD_COMP_ISOTH` の検証の入力、スクリプト、報告を置く。
論文の試験ごとのフォルダ（`Euler1D`、`MHD1D`、`Smooth_wave`、`OT2D`、`Turb3D`）に、その試験の全ケースの入力、説明（`README.md`）、図のスクリプト（`plot_<試験>.py`）、図が読む要約への symlink を置く。
Stage ごとの生成と解析のスクリプト、報告（`REPORT_stage*.md`）、共通の部品は `stages/` に置く。
計算の出力は git の管理外である。
ps1 の `~/work/calliope_isoth/tests/MHD_COMP_ISOTH` は、このフォルダと同じ配置で、論文の 589 ケースの実行（入力と出力）を持つ（「ps1 の配置」の節）。
解析の要約、図の入力、論文に使わない実行は、このフォルダの `runs/` に置く（「`runs/` の場所」の節）。
各スクリプトの使い方は、その docstring に書いてある。

## `runs/` の場所

このフォルダの `runs/` のうち、手元のものは解析の要約と図の入力を持ち、ps1 のものは、論文に使わない実行、batch の要約、試験のフォルダのケースを指す symlink を持つ（「ps1 の配置」の節）。
どちらも git の管理外である（`.gitignore`）。

2026-11-13 より前は、これをリポジトリの根の `runs/` に置いていた。
根の `runs/` は、試験とは別の production run（`runs/<モデル>/<問題>/...`）のために空けた。
REPORT、docstring、検討の記録（`chat.tex`）、ケースの一覧、JSON の記録にある `runs/...` は、移す前の根の `runs/` を指し、いまの `tests/MHD_COMP_ISOTH/runs/...` にあたる。
これらは書いた時点の記録なので、書き換えていない。
この README と試験のフォルダの README でも、`runs/` はこのフォルダの `runs/` を指す。

スクリプトが引数を省いたときに読み書きする場所は、`stages/paths.py` の `RUNS` が決める。
`RUNS` は、走らせた場所からこのフォルダの `runs/` への相対パスで、リポジトリの根で走らせると `tests/MHD_COMP_ISOTH/runs` になる。
生成のスクリプトが表示するケースのパス（`CASES`）もこの形になり、根で投入するジョブにそのまま渡せる。

## 論文の試験とフォルダ

論文に載せる試験と Stage、フォルダの対応を次の表に示す。
分類は、検討の記録（worktree の `chat.tex`、git の管理外）の「論文に乗せるテストはどれ？」への答えに従い、その後に行った Stage 16 を加えた。

| 論文での位置 | 試験 | Stage | フォルダ |
|---|---|---|---|
| 本文 | 1D 等温 Euler の衝撃波と単純波 | 1 | `Euler1D` |
| 本文 | 1D 等温 MHD の Riemann 問題 | 3、14、16（conv） | `MHD1D` |
| 本文 | 多次元の滑らかな解（円偏光 Alfvén 波、斜めの平面問題） | 4 | `Smooth_wave` |
| 本文 | 2D Orszag–Tang 渦 | 4、14 | `OT2D` |
| 本文 | 3D 減衰乱流（頑健さ、収支、スペクトル、密度 PDF、ESf の停止） | 6、7、15、16（K） | `Turb3D` |
| 付録 | MPI の分割によらないこと | 6、15 | `Turb3D` |
| 付録 | 係数の感度 | 8、15（ES25） | `Turb3D` |
| 付録 | 理想部の 3 つの定式化 A、B、C の比較（1D の試作） | 2 | なし（`runs/stage2`） |
| 限界 | M_s ≥ 2 | 4、5 | `Turb3D` |
| 限界 | B_t が反転する複合問題 C2 | 3 | `MHD1D` |
| 限界 | 質量の項を外した ES-SVV（svv_mass = F） | 13 | `Turb3D` |

Stage 2 の行は、論文の付録に載せると決めた 2026-11-17 に足した（「Stage 2 の復元」の節）。
`chat.tex` の分類は、1D 等温 Euler に Stage 14 も挙げている。
しかし Stage 14 で走らせ直したのは Stage 3 の 1D MHD と Stage 4 の OT だけで、Stage 1 の Euler の試験は最終版で走らせ直していない。
Stage 7 のうち論文に載せるのは密度 PDF の裾で、密度の穴の機構は載せない。
Stage 6 と Stage 13 の batch `regress`（Fortran の変更の前後のビット一致を見る回帰試験）も、論文には載せない。

フォルダごとのケースの数と、そのケースを走らせた実行ファイルを次に示す。
最終版は Stage 13 の 4 つ目のビルド `bin/calliope_isoth_s13d`（md5 026c510bac1f2e8811cdaf734281d15e）で、ソースは 07716cc である。
それ以降、Fortran のソースは変えていない（`src/time.F90` のコメント 1 行を除く）。

| フォルダ | ケース | 実行ファイル（`bin/calliope_isoth_*`） |
|---|---|---|
| `Euler1D` | 102 | `es`（Stage 1） |
| `MHD1D` | 215 | `s13d`（Stage 14、16 で走らせ直した） |
| `Smooth_wave` | 167 | `s4`（Stage 4） |
| `OT2D` | 15 | `s13d`（Stage 14 で走らせ直した） |
| `Turb3D` | 90 | `s4`、`s4b`、`s6`（Stage 4–8）、`s13d`（Stage 13、15、16） |

最終版で走らせ直していないケース（Euler1D、Smooth_wave、Turb3D の Stage 4–8）の扱いは、各フォルダの README に書いた。

## ケースのファイル

各ケースは、ps1 の `runs/` の下のパスのまま、試験のフォルダの下に置いた（例: `runs/stage14/ot/ot_N256_eta0` は `OT2D/stage14/ot/ot_N256_eta0/`）。
ps1 の元のディレクトリをパスでたどれるようにするためである。
ps1 では、元のディレクトリそのものをこの位置に移し、`runs/` の元のパスには移した先を指す symlink を置いた（「ps1 の配置」の節）。
ケースのディレクトリには次のファイルがある。

| ファイル | 中身 | あるケース |
|---|---|---|
| `calliope.in` | calliope の入力 | すべて |
| `nproc` | MPI のランク数 | すべて |
| `case.json` | 問題の設定（問題名、N、dt、ステップ数など）。`regen_inputs.py` が読む | Euler1D 以外 |
| `generated.md5` | ps1 の元の `profile1d.dat` または `ic/*.dat` の md5（`md5sum` の形式） | Euler1D 以外 |
| `exe` | 実行ファイルのパス（作業ディレクトリからの相対）。ジョブスクリプトの `EXE` より優先する | MHD1D、OT2D、Turb3D の Stage 13–16 |
| `fftw_estimate` | FFTW の計画を ESTIMATE に固定する印（`stages/fftw_shim.c`） | Euler1D の restart の試験、MHD1D の conv、Smooth_wave の `alfven`、`cpdiag` と `planar` の一部、Turb3D の Stage 15 の分割の試験 |
| `t_end` | 終了時刻。ジョブスクリプトがこの時刻まで restart で続ける | Turb3D の 'decrement' と 'track' のケース |

写したのは ps1 の元にある入力のファイルだけで、589 ケースの入力 2151 本の md5 が ps1 の元と一致した。
生成した初期値（`profile1d.dat` と `ic/`）と実行の出力は写していない。
Stage 1 の restart の試験で投入の前に置いた `terminate` は、calliope が読んで消したので元に残っておらず、ここにもない（`Euler1D/README.md`）。

## ps1 の配置

ps1 の `~/work/calliope_isoth/tests/MHD_COMP_ISOTH` は、2026-11-12 から手元と同じ配置である。
git の管理下のファイルは、手元から `tar` で送って置いた（ps1 では git を使えない）。手元で変えたら、同じように送り直す。
589 ケースのディレクトリは元の実行そのもので、入力、生成した初期値（`profile1d.dat`、`ic/`）、出力を持つ。
元は `runs/<stage>/<batch>/[<sub>/]<case>` にあり、同じファイルシステムの中で mv で移した（複写していない）。
ケースに足したのは git の `generated.md5` だけで、git のほかの入力は ps1 の元のものと md5 が一致した。

- `runs/` の元のパスは、移した先を指す相対の symlink である。REPORT と解析のスクリプトの `runs/...` のパスは、そのまま使える。
- batch の要約（`summary.json` など）と、論文に使わないケースは、`runs/` に残した。図が読む要約は、試験のフォルダの symlink からもたどれる（「図」の節）。
- 手元の Stage 16 の解析で書いた `runs/stage16/conv_table.json`、`judge_tm.json`、`judge1_hf.json` は、ps1 になかったので、2026-11-13 に手元から写した。`check_report_stage16.py` は 4 本の要約を読み、ps1 でも手元と同じ結果を出す。
- `du runs/...` は symlink の先を数えない。ケースの量は `du -L` か、試験のフォルダの側で数える。
- ケースのディレクトリで投入すると、元の出力を上書きする（「走らせ方」の節）。

移す前の平らな `tests/MHD_COMP_ISOTH` と `tests/MHD_INCOMP_regression` は、`~/work/calliope_isoth/attic/tests_flat_20261112/` に残した（ps1 では何も消していない）。
平らな配置の `fig/`（Stage 2–5 の解析が書いた PNG）だけは、`stages/fig/` に移した。
同じ attic に、移したケースの一覧（`moves.txt`）、移した記録（`restructure.log`）、元に戻すスクリプト（`rollback.sh`。`DRY=0` を付けたときだけ動かす）がある。
移した後、ケースの中の 61506 項目の inode、大きさ、mtime は移す前と一致した。
変わったのは、`generated.md5` を足した 487 ケースの最上位のディレクトリの mtime だけである。

2026-11-13 に、根の `runs/` を `tests/MHD_COMP_ISOTH/runs/` に移した（同じファイルシステムの中の rename）。
`runs/` の中の symlink 1084 本のうち、移すと切れる 719 本の指し先を書き換えた。

- ケースを相対パスで指す 589 本は、指し先から `tests/MHD_COMP_ISOTH/` を抜いた。
- `runs/` の中を絶対パスで指す 130 本（`runs/stage14/{seams,dry,old}` の下）は、指し先の `/calliope_isoth/runs/` を `/calliope_isoth/tests/MHD_COMP_ISOTH/runs/` に置き換えた。
- `runs/` の中を相対パスで指す 365 本は、書き換えずにそのまま解決する。

書き換えた後、1084 本が指す実体は、移す前と一致した。
書き換えた一覧（`links.tsv`）と元に戻すスクリプト（`rollback.py`。`DRY=0` を付けたときだけ動かす）は、`~/work/calliope_isoth/attic/runs_move_20261113/` にある。
2026-11-12 の `rollback.sh` は根の `runs/` を前提にするので、`tests/MHD_COMP_ISOTH/runs` があれば止まる行を足した（元は同じ場所の `rollback.sh.orig`）。
2026-11-12 の再配置まで戻すときは、先に `rollback.py` で `runs/` を根に戻す。

`runs/` にある過去のスクリプト（`runs/stage14/ana_*.sh`、`runs/stage4/*.py` など）は、平らな配置のパス `tests/MHD_COMP_ISOTH/<file>` を参照するので、そのままでは走らない。
`runs/` にあるケースの一覧（`cases_*.txt` など）とジョブのスクリプトは、根からのパス `runs/...` を書いているので、投げ直すときは `tests/MHD_COMP_ISOTH/runs/...` に直す。
これらは記録なので書き換えていない。
attic の平らな配置のファイルのうち `analyze_stage16.py` と `REPORT_stage5.md` は、git のどの版とも違う。どちらも git の版より古い版である。

## 走らせ方

ps1 で 1 つのケースの組を走らせる手順を示す。
試験に固有の点（`EXE=`、ランク数、ジョブの数）は、各フォルダの README と REPORT の再現手順の章にある。

1. 送る。ps1 は外部のネットワークにつながらず、git も使えないので、手元から `tar` と `ssh` で送る。
   - 作業ディレクトリには、ps1 の `~/work/calliope_isoth` でなく、別の作業ディレクトリ（`bin/` と `tests/` を持つ calliope のコピー）を使う。`~/work/calliope_isoth` のケースのディレクトリには元の実行の出力があり、そこで投入すると上書きするからである。
   - `tests/MHD_COMP_ISOTH/` を、作業ディレクトリの同じパスに置く。
   - ケースを、作業ディレクトリの `tests/MHD_COMP_ISOTH/runs/` の下の同じパスに置く。生成と解析のスクリプトが読み書きするパスの既定値が、そこだからである（「`runs/` の場所」の節）。
2. 初期値を書く。`python3 tests/MHD_COMP_ISOTH/stages/regen_inputs.py <case> ...` で、`case.json` から `profile1d.dat` または `ic/` を書き、`generated.md5` と照らす（Euler1D は組み込みの初期値なので不要）。
3. 実行ファイルを用意する。最終版は次のとおりビルドした（REPORT_stage13 の 12 章）。
   ```sh
   module load intel/2025.1 fftw/3.3.10/oneapi2025.1 hdf5/1.14.6/oneapi2025.1-parallel \
     netcdf-fortran/4.6.2/oneapi2025.1-parallel
   make clean MODEL=MHD_COMP_ISOTH && make MODEL=MHD_COMP_ISOTH && cp calliope bin/calliope_isoth_s13d
   ```
4. 投入する。
   ```sh
   qsub -l select=1:ncpus=<cores>:mem=<mem>:mpiprocs=<cores>:ompthreads=1 \
     -v CASEFILE=<list>,PAR=<cores>[,EXE=bin/calliope_isoth_<ver>][,BUDGET=<min>] \
     tests/MHD_COMP_ISOTH/stages/job.pbs-ps-fes-dev
   ```
   `CASEFILE` はケースのパス（作業ディレクトリからの相対）を 1 行に 1 つ並べたファイルで、`CASES="tests/MHD_COMP_ISOTH/runs/a tests/MHD_COMP_ISOTH/runs/b"` でも渡せる。
   `PAR` は複数のケースを同時に走らせるときのコア数、`BUDGET` はジョブの開始から実行を止めるまでの分数（既定 27）である。
   予算で止まったケースは restart を保存して `exit_status` が `incomplete` になり、同じ一覧で投げ直すと続きから走る。
   ジョブスクリプトの使い方の詳細は、その冒頭のコメントにある。
5. 解析する。各フォルダの README の「解析」の節による。

## 生成する初期値の書き直し

`stages/regen_inputs.py` は、`case.json` の内容から生成の関数（`gen_stage3.write_profile`、`gen_stage4` の場の関数と `gen_ic.write_restart`）を選んで初期値を書き、`generated.md5` と照らす。
既にあるファイルは書き直さず、照らすだけにする。

md5 がビット単位で一致するのは、ps1 の実行と同じ numpy（2.0.2、x86_64）を使ったときである。
ps1 では、`profile1d.dat` の 266 ケースと、`ic/` の md5 の組 57 から 1 ケースずつ（256³ を含む）の、合わせて 323 ケースがすべて一致した（書いた内容はメモリで照らし、ps1 には書いていない）。
手元（numpy 1.26.4、macOS の arm64）では、`profile1d.dat` の 53 ケースと `ic/` の全 57 組が一致しなかった。
差を調べた `profile1d.dat` の 3 ケースでは絶対差が 7.5e-16 以下、`ic/` の 5 ケースでは場の最大値に対して 5.2e-16 以下で、FFT と初等関数の丸めの違いと読める。

`profile1d.dat` と `ic/*.dat` は、リポジトリの根の `.gitignore` の `*.dat*` で git から外れる。
試験のフォルダの中で書いても、誤って commit されることはない。
256³ の `ic/` は 1 ケース約 0.95 GB である。

589 本の `calliope.in` が使う namelist の変数は、いまの `src` がすべて宣言している（静的な照合）。
いまの版で走らせたときに元と同じ結果になることは、最終版で走らせたケースのほかは確かめていない。

## 図

各試験のフォルダの `plot_<試験>.py`（`Euler1D/plot_euler1d.py` など）が、論文に載せる図をそのフォルダの `fig/<名前>.pdf` に書く。
図の一覧と中身は、各フォルダの README の「図」の節にある。
スタイルは gallope の strong scaling の図（`~/Dropbox/Codes/scaling/gallope/flow/strong scaling/202610/plot_flow.py`）に合わせ、`stages/figstyle.py` にまとめた。
LaTeX の Times（newtx）、30 pt の字、内向きの目盛り、Okabe–Ito の色で、PDF に書く。

コマンドは、リポジトリの根で次のとおりである（例は Euler1D で、ほかのフォルダも同じ形である）。

```sh
python3 tests/MHD_COMP_ISOTH/Euler1D/plot_euler1d.py list              # 図と入力の一覧
python3 tests/MHD_COMP_ISOTH/Euler1D/plot_euler1d.py fetch             # ps1 にしかない入力を runs/ に持ち帰る
python3 tests/MHD_COMP_ISOTH/Euler1D/plot_euler1d.py [plot] [名前 ...]  # 図を書く（名前を省くとすべて）
```

`--runs` で入力の根（既定は、走らせた場所によらずこのフォルダの `runs/`）を、`--host` と `--remote` で fetch の相手（既定 `ps1` の `~/work/calliope_isoth`）を変えられる。
入力のない図は、足りないファイルを表示して飛ばし、最後に飛ばした図を並べる。

- fetch は ps1 に何も書かず、手元にあるファイルを上書きしない（あれば飛ばす）。
- 小さいファイル（要約の JSON、out1d の行）は、ssh と tar で、ps1 と同じパスで手元の `runs/` に写す。
- 大きい出力（2D と 3D の面、`cfl.dat`、平面問題の分布）からは、要る配列だけを取り出す。`stages/figdata.py` を ssh の標準入力で ps1 に送って走らせ、npz を標準出力で `runs/figdata/<試験>/` に受け取る。
- 写すものも取り出すものも、ps1 の `tests/MHD_COMP_ISOTH/runs/` からのパスで指す。
- 取り出しは、ps1 の作業ディレクトリにある解析のモジュール（`analyze_stage1.py`、`analyze_stage4.py`、`gen_ic.py`、`gen_stage4.py`、`riemann_exact.py`、`cp_stability.py`）を import する。モジュールは ps1 の `tests/MHD_COMP_ISOTH/stages/` から読み、md5 が手元の `stages/` のものと違えば止まる。
- fetch が持ち帰る量は、Euler1D が 2.7 MB、Smooth_wave が 0.8 MB、OT2D が 20 MB、Turb3D が 14 MB である。解析の要約（Stage 3 の FV の参照、Stage 7–8、13–16 の JSON など）は、各 Stage の解析で手元の `runs/` に置いたものを読み、fetch の対象にしていない。MHD1D の図は、これだけで描ける。

図が読む要約（Stage の判定表、batch の `summary.json`、FV の参照など。5 つのフォルダで合わせて 24 項目）は、試験のフォルダからもたどれるように symlink を置いた。
symlink は、試験のフォルダの下で、要約の `runs/` の下のパスと同じパスに置き、このフォルダの `runs/` の同じパスを相対で指す（例: `Euler1D/stage1_es/judge.json` は `../../runs/stage1_es/judge.json` を指す）。
図のスクリプトは、symlink でなく `--runs` の下を読む。
symlink は git で管理する。手元と ps1 は試験のフォルダと `runs/` の相対の位置が同じなので、同じ symlink がどちらでも解決する。
`runs/` は git の外にあるので、`runs/` のない場所（新しい clone など）では、どの symlink も切れている。
fetch が取り出す npz（`runs/figdata/<試験>/`）と、ケースの出力（Euler1D の out1d の行など）には、symlink を置いていない。npz は手元にしかなく、ケースの出力は ps1 では試験のフォルダのケースの中に実物があるからである。

図は手元でしか描けない。newtx が要り、ps1 にはないからである。
各スクリプトは、描いた量の主な数値（収束の次数、最大値、判定の値）を標準出力に書く。REPORT の表と照らすためである。
PDF は、リポジトリの根の `.gitignore` の `*.pdf` で git から外れる。`fig/` はスクリプトが作る。

## stages/

### Stage ごとのファイル

| Stage | 内容 | 生成 | 解析 | 報告 |
|---|---|---|---|---|
| 1 | 1D 等温 Euler（B = 0） | `gen_stage1.py` | `analyze_stage1.py`、`judge_stage1.py` | `REPORT_stage1.md` |
| 2 | 理想部の定式化（A、B、C）、衝撃波粘性の形、フィルタの位置、SVV の形の比較（1D の試作、B = 0） | `stage2.py`（試作の実行） | `stage2.py check`、`stage2_table.py` | `REPORT_stage2.md` |
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

Stage 2 と Stage 3 の表と図のスクリプト（`stage2_table.py`、`stage3_table.py`、`stage3_aux.py`）は、図を `stages/fig/` に書く。このディレクトリは git の管理外である。
`stage2_table.py` は、表の書式と色の定義を `stage3_table.py` から読む。

### 共通の部品

- `regen_inputs.py`: 試験のフォルダのケースの初期値を書き直し、`generated.md5` と照らす（上の節）。
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
- `fftw_shim.c`: p3dfft が FFTW に渡す計画の flag を stderr に書き出し、`FFTW_SHIM_ESTIMATE=1` のときは FFTW_ESTIMATE に置き換える LD_PRELOAD の shim。FFTW_MEASURE では同じ入力の 2 回の実行が 1e-13 程度ずれるので、ビット一致の試験に使う。ケースに `fftw_estimate` があると、ジョブスクリプトがこれをビルドして使う。
- `stage3_nstep.json`: Stage 3 の試作の実行が t_end に届いたときの step 数と dt。`gen_stage3.py nstep` が書き、Stage 3 の Fortran のケースはこの dt で走る。
- `incomp_regress.in`: MHD_INCOMP の回帰試験の入力。`gen_stage6.py` と `gen_stage13.py` の batch `regress` が使う。
- `figstyle.py`: 試験のフォルダの図のスタイル（上の「図」の節）。rc の設定、色、線の太さ、パネルと地図の配置、`fig/` への書き出しを持つ。
- `figdata.py`: 図のスクリプトが共有するコマンドライン（`plot`、`fetch`、`list`）と fetch。ps1 で走らせて出力から配列を取り出す部分も、このファイルにある。

## 整理の記録

### 試験ごとのフォルダへの移動（2026-11-11）

`tests/MHD_COMP_ISOTH/` の直下にあったスクリプト、REPORT、共通の部品を、中身を変えずに `stages/` へ移した。
変えたのは、移したことで壊れるパスだけである。

- `fv_mhd1d.py` の `REF_DIR` と `stage3.py` の `OUT`（ファイルの位置からリポジトリの根をたどる相対パス）に `'..'` を 1 つ足した。
- `job.pbs-ps-fes-dev` が `fftw_shim.c` をビルドするパスと、使い方のコメントのパス。
- 使い方の docstring の `python3 tests/MHD_COMP_ISOTH/<file>`（`analyze_stage5.py`、`check_report_stage16.py`、`check_seams.py`、`gen_stage14.py`、`gen_stage16.py`）。
- `src/time.F90` のコメントにある `REPORT_stage6.md` のパス（コメントだけで、コードは変わらない）。

REPORT の本文は記録なので変えていない。
REPORT の再現手順などにあるパス `tests/MHD_COMP_ISOTH/<file>` は、いまは `tests/MHD_COMP_ISOTH/stages/<file>` にある。

### 論文で使う Stage への絞り込み（2026-11-11、ca9a205）

論文で使う Stage のものだけを残し、次の 19 本を消した。
どれも cf9bb60 から取り出せる（移動の前のパスで指定する）。

```
git show cf9bb60:tests/MHD_COMP_ISOTH/<file> > <file>
```

| 群 | ファイル | 内容 |
|---|---|---|
| どこからも使われない診断 | `diag_u4.py`、`fast_wave.py`、`void_probe.py`、`post_stage1_es.sh`、`stage0/fast_wave.in`、`stage0/fast_wave_new_default.in`、`stage0/fast_wave_new_ideal.in` | U4 の衝突流（Stage 1）と 3D の実行の密度の最小点の診断、Stage 1 の ps1 での後処理、Stage 0 の速進波の固有モードの入力と解析 |
| Stage 2 | `stage2.py`、`stage2_table.py`、`REPORT_stage2.md` | 数値定式化の比較（2026-11-17 に `stages/` に戻した。「Stage 2 の復元」の節） |
| Stage 9–12 | `analyze_stage9.py`、`analyze_stage10.py`、`analyze_stage11.py`、`analyze_stage12.py`、`gen_stage10.py`、`REPORT_stage9.md`–`REPORT_stage12.md` | 密度の穴の質量収支と、ES-SVV の項の svv_eps0 への依存 |

残したスクリプトが使っていた定義は、消す前に次のように移した。
移した定義の AST は元と同じである。
ただし `check_control` からは、Stage 10 の対照に固有の既定値 `allowed=ALLOWED` を外し、`allowed` を必ず渡す引数にした。

| 定義 | 元（消した） | 移した先 |
|---|---|---|
| `shells`、`spectrum` | `analyze_stage12.py` | `analyze_stage13.py` |
| `md5s`、`check_control` | `gen_stage10.py` | `gen_stage13.py` |
| `COLOR`、`MARK`、`GRAY`、`e`、`stats`、`status`、`gmax`、`order`、`table` | `stage2_table.py` | `stage3_table.py`（戻した `stage2_table.py` も、ここから読む） |

`incomp_regress.in` は `tests/MHD_INCOMP_regression/calliope.in` から内容を変えずに移した。

残した REPORT と一部の docstring は、消したファイルや `tests/MHD_INCOMP_regression` に言及している。
これらは書いた時点の記録なので変えていない。
言及されたファイルは、上の `git show` で取り出せる。

### 試験の `runs/` の移動（2026-11-13）

リポジトリの根の `runs/` を `tests/MHD_COMP_ISOTH/runs/` に移し、根の `runs/` を production run のために空けた（「`runs/` の場所」の節）。
スクリプトで変えたのは、`runs/` の場所を決めていた箇所だけである。

- `stages/paths.py` を足し、場所を `RUNS` の 1 箇所で決めるようにした。
- 生成と解析のスクリプトで、引数を省いたときの既定値と module の定数（`os.path.join('runs', ...)`）を `paths.RUNS` から作るようにした。`stage3.py` の `OUT` と `fv_mhd1d.py` の `REF_DIR` も同じである。
- `figdata.py` の `--runs` の既定値を `paths.RUNS` にした。fetch は ps1 の `tests/MHD_COMP_ISOTH/runs/` で写し、取り出しもそこで走らせる。そのため、図のスクリプトが取り出しに渡すケースのパスから、先頭の `runs/` を外した。
- `job.pbs-ps-fes-dev` の使い方のコメントにある `CASES` の例を直した。
- 試験のフォルダの 24 本の symlink の指し先を、2 段浅くした。
- `.gitignore` に `/runs/` と `tests/MHD_COMP_ISOTH/runs/` を足した。

変えた後、手元の 30 枚の図は変える前と画素が一致し、`check_report_stage16.py` の出力も一致した。

### Stage 2 の復元（2026-11-17）

論文の付録に理想部の定式化の比較（Stage 2）を載せることにしたので、ca9a205 で消した次の 3 本を cf9bb60 から `stages/` に戻した。

- `REPORT_stage2.md`: 内容を変えていない。書いた時点の記録であり、論文の出典の行番号はこの版の行を指す。
- `stage2.py`: 出力の場所 `OUT` を `paths.RUNS` から作るようにした。ほかは変えていない。
- `stage2_table.py`: ca9a205 で `stage3_table.py` に移した 9 つの定義を消し、`stage3_table.py` から読むようにした。移した定義の AST は、`stage2_table.py` の元の定義と同じである。

戻した後、`stage2.py table` と `stage2.py check` を手元で走らせた（Python 3.12.2、numpy 1.26.4、matplotlib 3.8.4）。
`runs/stage2/tables.md`、`stages/fig/stage2_*.png` の 2 枚、check の出力は、戻す前の記録（`check` の出力は `runs/stage2/check.txt`）とバイト単位で一致した。
174 本の実行（`stage2.py run`）は走らせ直していない。
