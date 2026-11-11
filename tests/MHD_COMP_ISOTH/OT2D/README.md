# OT2D: 2D の等温 Orszag–Tang 渦

論文の本文に載せる 2D の等温 Orszag–Tang 渦（OT）の入力を置く。
問題の設定と判定の量は Stage 4（`../stages/REPORT_stage4.md` の 7 章）で決めた。
論文の数値は、最終版の実行ファイル `bin/calliope_isoth_s13d` で走らせ直した Stage 14（`../stages/REPORT_stage14.md` の 6 章）の実行から取る。
入力は ps1 の `~/work/calliope_isoth/runs/stage14/ot/` から写し、`runs/` の下のパスのまま置いた。
Stage 14 の入力は、Stage 4 の実行（`runs/stage4/ot`）の `calliope.in`、`nproc`、`case.json`、`ic/` を複製したものである（REPORT_stage14 の 2 章）。
生成スクリプトで作り直さなかったのは、Stage 4 より後に足した既定値で `calliope.in` の字面が変わりうるからである。

## 問題の設定

設定は REPORT_stage4 の 7.1 節による。

- Mignone（2007）の 5.5 節の等温 OT 渦（音速 a = 2）を、u と B を a で割り t に a を掛けて c_s = 1 に直した。
- x-z 面（nly = 2、`dim_decomp = 1`）で、ρ = 1、u = (−sin z, 0, sin x)、B = B0 (−sin z, 0, sin 2x)、B0 = √(3/5)（β = 10/3）である。箱は (2π)² で、境界は周期である。
- t_end = 2π で、面を Δt = 0.25 ごとに書く（`write_intvl_planes`）。Mignone の t = 1.5、3 はここでの t = 3、6 にあたるので、t_end を計画の π から 2π に延ばした。
- 係数は C_shock = 2、ES-SVV（svv_eps0 = 0.5）、C_eta = 0 と 2 である。
- dt は N ごとに固定で、`reset_method = 'stop'` である。N = 128、256 は dt = dt_adv(t = 0)/4 とし、N = 512 と 1024 はその下の N の実行の min(dt_adv, dt_visc) から決めた（7.1 節の式）。値は各ケースの `case.json` にある。
- 初期値は `../stages/gen_ic.py` が restart の形式で `ic/` に書き、calliope は `init_type = 'restart'` で読む。ケースには `ic/` を置かず、その md5 を `generated.md5` に書いた。

## ケース

`stage14/ot/` の下に 15 本がある。

| ケース | 数 | 内容 | ランク |
|---|---|---|---|
| `ot_N{128,256,512}_eta{0,2}` | 6 | t = 2π まで | 4（N = 128）、16（256）、128（512） |
| `ot_N1024_eta0_p256` | 1 | N = 1024、C_eta = 0 を t = π まで（25146 ステップ） | 256 |
| `dec_N256_np{4,16,64,256}` | 4 | N = 256、C_eta = 2 を t = 1 まで（1511 ステップ）。MPI の分割による差を見る | 4、16、64、256 |
| `dec_N256_np{4,16,64,256}_n10` | 4 | 同じものを 10 ステップだけ | 4、16、64、256 |

N = 1024 は、計画の t_end（π）までにとどめた。
ε の最大（t ≈ 1.68）を含み、t ≥ 4 では N = 512 までで各点の量が収束しないからである（REPORT_stage4 の 7.6 節、7.7 節）。
`case.json` の `nstep`（50292）は t = 2π に達するステップ数で、`nrun` と `calliope.in` の `nstep`（25146）が実際に走らせたステップ数である。
C_eta = 2 の N = 1024 は、費用のため走らせていない。

Stage 14 の判定（REPORT_stage14 の 3 章で事前に固定）は、O1（7 本の完走と収支）、O2（分割による 10 ステップ後の差 ≤ 1e-12）、O3（収束の判定の合否が Stage 4 と同じ）、O4（t ≤ π の Stage 4 との相対差 ≤ 1e-6）のすべてが合格した。
収束の判定（O3）そのものは、C_eta = 0、2 の両方で不合格の項目を含む。
t ≥ 4 ではプラズモイドができ、各点の量が N を上げても収束しないためで、理由は REPORT_stage4 の 7.5–7.7 節にある。

## 走らせ方

手順は `../README.md` の「走らせ方」による。この試験に固有の点は次のとおりである。

- ケースを ps1 の `runs/` の下の同じパスに写し、`python3 tests/MHD_COMP_ISOTH/stages/regen_inputs.py <case> ...` で `ic/` を書く。
- 各ケースの `exe` が `bin/calliope_isoth_s13d` を指すので、`EXE=` は要らない。
- Stage 14 は 3 つの一覧（N = 128、256 と分割の 8 本、N = 512、N = 1024）を、それぞれ 256 コアのジョブ（`-l select=1:ncpus=256:mem=745gb:mpiprocs=256:ompthreads=1`、`PAR=256`）で投げ、予算で止まったら同じ一覧で投げ直した（REPORT_stage14 の 10 章）。N = 512 は 2 本、N = 1024 は 3 本のジョブで終わった。

## 解析

- `python3 tests/MHD_COMP_ISOTH/stages/analyze_stage4.py ot runs/stage14` で、判定（停止、保存則、分割、収束）と数値を標準出力と `runs/stage14/ot/summary.json` に書き、図を同じディレクトリに置く。
- Stage 4 との比較（O1–O4）は `python3 tests/MHD_COMP_ISOTH/stages/analyze_stage14.py ot runs/stage14` で行う。Stage 4 の実行を `analyze_stage4.py ot` で解析した `runs/stage14/old/ot/summary.json` と、その実行の出力を読むので、ps1 の `runs/stage4/ot` がない環境では使えない。

## 図

`plot_ot2d.py` が `fig/<名前>.pdf` を書く。コマンドと fetch は `../README.md` の「図」の節による。
`ot_eps` と `ot_spec` は、手元の `runs/stage14/ot/summary.json`（上の `analyze_stage4.py ot` が書く）を読む。
`ot_rho` と `ot_J` は、C_eta = 0 の実行の t = 3、6 の面の ρ、B_x、B_z（float32）を読む。`fetch` が ps1 の出力から `runs/figdata/OT2D/` に取り出す（合わせて 20 MB）。
どの図も最終版（s13d）の実行である。

| 図 | 中身 | REPORT |
|---|---|---|
| `ot_eps` | ε(t) = −dE/dt。C_eta = 0（N = 128–1024。1024 は t = π まで）と C_eta = 2（N = 128–512） | Stage 4 の 7.4 節、Stage 14 の 6.2 節 |
| `ot_spec` | t = 1、3 の E(k) = E_{√ρu}(k) + E_B(k)（殻の和）。C_eta = 0、N = 128–1024。k = N/3 の点線と k⁻² の線を添える | Stage 4 の 7.5 節、7.7 節 |
| `ot_rho` | N = 512、C_eta = 0 の t = 3、6 の ρ の地図。共通の対数の色の目盛り | Stage 4 の 7.6 節 |
| `ot_J` | C_eta = 0 の J_y = ∂_z B_x − ∂_x B_z の地図。上の行は t = 3 の N = 256、512、1024、下の行は t = 6 の N = 128、256、512。色は \|J_y\| = 50 で飽和させる | Stage 4 の 7.6 節、7.7 節 |

- `ot_J` は、取り出した面から求めた max\|J_y\| と ρ の範囲を、`summary.json` の `late`（ps1 で同じ記録から float64 で求めた値）と照らして標準出力に書く。
- N = 512、t = 6 の ρ の最大 6.3 は、箱の中央の 1 か所の値である。ρ > 3 の点は 99 点で、13 × 9 格子の範囲に収まる。
