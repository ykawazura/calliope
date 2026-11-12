# Turb3D: 3D 等温 MHD の減衰乱流

論文に載せる 3D 等温 MHD の減衰乱流の入力を置く。
本文では、頑健さ、収支、スペクトル、密度 PDF と、ES-SVV と標準の散逸の比較を示す（Stage 6、7、15、16）。
付録では MPI の分割によらないこと（Stage 6、15）と係数の感度（Stage 8、15）を、限界の例では M_s ≥ 2（Stage 4、5）と質量の項を外した ES-SVV（svv_mass = F、Stage 13）を示す。
入力は ps1 の `~/work/calliope_isoth/runs/` の下から写し、`runs/` の下のパスのまま置いた。

## 最終版で走らせ直していないこと

Stage 13、15、16 の実行は、最終版の実行ファイル `bin/calliope_isoth_s13d` で走らせた。
Stage 4、5、6、8 の実行は、それぞれの Stage の実行ファイル（`bin/calliope_isoth_s4`、`_s4b`、`_s6`）で走らせ、最終版では走らせ直していない。
s4b は、restart の継ぎ目で収支の行が欠けないように s4 の終了処理を直した版である（REPORT_stage4 の 3.3 節）。
版の間で既定の設定の結果が変わらないことは、小さい回帰のケースで確かめた。
ISOTH の回帰の 7 組で、s6 は s4b と、s13d（svv_mass = T）は s6 と、ビット単位で一致した（REPORT_stage6 の 5 章、REPORT_stage13 の 4 章と 5.2 節の V1）。
ここに置いた乱流の実行そのものを、版の間で比べてはいない。
`calliope.in` に書いていない変数は、走らせる版の既定値をとる。

## 問題の設定

設定は REPORT_stage4 の 8.1 節による。

- 箱は (2π)³ で、境界は周期である。ρ = 1、B0 = 1 ẑ、c_s = 1（β0 = 2）とする。
- 初期の u と δB は、スペクトル E(k) ∝ k⁴ exp(−2k²/k_p²)（k_p = 2）を 1 ≤ |k| ≤ 4 に持つ乱数の場である。乱数の種 20261103 で 16³ の格子に作り、どの N にも同じ Fourier 係数を埋め込む。
- 問題は M_s = rms|u|（c_s = 1）と rms|δB|/B0 の組で、次の 6 つである。どの問題も同じ乱数の場を使い、u と δB の振幅だけが違う。

| 問題 | M_s | rms\|δB\|/B0 | 論文での扱い |
|---|---|---|---|
| `turb1` | 1 | 1 | 本文、付録、svv_mass = F |
| `turb05` | 0.5 | 0.5 | Stage 15 |
| `turb4` | 4 | 4 | M_s ≥ 2（Stage 4、5） |
| `turb4b` | 4 | 1 | M_s ≥ 2（Stage 4、5） |
| `turb2b`、`turb3b` | 2、3 | 1 | M_s の上限（Stage 5 の 12 章） |

- 時間の単位は τ0 = π/M_s で、t_end は 3τ0（Stage 4–6）、2τ0（Stage 15、16）、1.2τ0（Stage 8、13）である。
- 既定の係数は C_shock = 2、ES-SVV（svv_eps0 = 0.5）、C_eta = 0 である。
- 初期値は `../stages/gen_ic.py` が restart の形式で `ic/` に書き、calliope は `init_type = 'restart'` で読む。ケースには `ic/` を置かず、その md5 を `generated.md5` に書いた。
- dt の決め方は Stage で違う。
  - Stage 4 は dt を固定し、dt_cfl を超えたら `reset_method = 'stop'` で止める。止まった `turb4_N64` と `turb1_N128` は `'decrement'`（dt を下げて続ける）で走らせ直し、名前に `_a` を付けた。`turb4b_N64` と Stage 5 の実行は、はじめから `'decrement'` で走らせた。
  - Stage 6 以降は、dt を dt_cfl に追従させる `'track'`（有効数字 3 桁への切り下げで MPI の分割によらない）で走らせた（REPORT_stage6 の 3 章）。
  - `'decrement'` の実行では、ケースのファイル `t_end` が終了時刻を決め、ジョブスクリプトがその時刻まで restart で続ける。
  - `'track'` の実行は `calliope.in` の `t_end` で止まり、`nstep` は上限である。ケースのファイル `t_end` にも同じ値を置いた。`stage6q/track/dec_N64_np*` は 400 ステップで終わるので、`t_end` を持たない。
  - Stage 15 の段階 0 の dt と分割の試験は、dt を固定して決めたステップ数だけ走らせる（`'stop'`）。

## 方式の記号

Stage 15 と 16 のケース名の `<方式>` は、次の散逸の設定を指す（REPORT_stage15 の 2.2 節、REPORT_stage16 の 2.4 節）。

| 記号 | 設定 | 役割 |
|---|---|---|
| ES | 既定（svv_eps0 = 0.5、filter_lnrho = F） | 本番の方式 |
| ESf | ES に filter_lnrho = T を足す | ln ρ のフィルタの効果を切り分ける対照 |
| ES25 | svv_eps0 = 0.25 | 係数の感度（64³ だけ） |
| R0 | SVV なし、ν = η = 0、filter_lnrho = T | SVV が要るかの対照 |
| L(r) | SVV なし、filter_lnrho = T、m と B のラプラシアン（nu = eta、lmd = 0） | 標準の散逸 |
| H(r) | SVV なし、filter_lnrho = T、m、B、ln ρ の k⁸ の超粘性 | 標準の散逸 |
| HF(r) | H で filter_lnrho = F にしたもの | Stage 16 の 2 本目の参照の候補 |

r は L と H の係数を 1 つの量でそろえる。
k_c = N/3 での減衰率を γ_c = r ε_ref^{1/3} k_c^{2/3} とおき、L は nu = eta = (27/4)γ_c、H と HF は nu_h = eta_h = lmd_h = (27/4)⁴γ_c とする。
ε_ref は 64³ の ES のエネルギー散逸率の最大で、M_s = 1 で 0.1424333、M_s = 0.5 で 0.01471188 である。
L と H は、ES と比べるための参照であり、ES を安定にするための粘性ではない。

## batch とケース

ケースは 90 本である。
下の表の「実行ファイル」が s4、s4b、s6 の batch にはケースに `exe` がなく、s13d の batch では各ケースの `exe` がその実行ファイルを指す。

### 本文

| batch | 数 | 内容 | ランク | 実行ファイル | REPORT |
|---|---|---|---|---|---|
| `stage6q/turb` | 1 | `turb1_N256`。turb1 の 256³ を 3τ0 まで | 256 | s6 | Stage 6 の 8、9 章 |
| `stage6q/track` | 2 | `turb1_N64_t`、`turb1_N128_t`。'track' の turb1 を Stage 4 の実行と比べる | 64、128 | s6 | Stage 6 の 7 章 |
| `stage4/turb` | 2 | `turb1_N64`、`turb1_N128_a`。turb1 の 64³ と 128³ | 64、128 | s4 | Stage 4 の 8 章、Stage 6 の 8 章 |
| `stage15/pre` | 2 | `turb1_N64_ES`、`turb05_N64_ES`。ES の 64³（ε_ref を決める） | 64 | s13d | Stage 15 の 2.3 節 |
| `stage15/a` | 30 | 64³ で、M_s ごとに ESf、ES25、R0、L(1/4–4)、H(1/4–16) | 64 | s13d | Stage 15 |
| `stage15/b` | 14 | 128³ で、turb1 の ES、ESf、H(2、4、8)、L(1、2) と turb05 の ES、H(1/4、1/2、1)、L(1、2)、R0 | 128 | s13d | Stage 15 |
| `stage15/c` | 2 | 256³ の turb1 の H(8)、H(32) | 256 | s13d | Stage 15 の 4.1 節と 5 章 |
| `stage16/hf64`、`hf128`、`hf256` | 7、3、1 | turb1 の HF を 64³（r = 1/4–16）、128³（r = 2、4、8）、256³（r = 8）で | 64、128、256 | s13d | Stage 16 の 5 章 |

- 密度 PDF（Stage 7）は、新しい実行を作らず、`stage4/turb` の 2 本と `stage6q/turb/turb1_N256` の out3d を読む。
- Stage 15 の K の判定は、256³ の参照に `stage6q/turb/turb1_N256`（ES-256）を使う。128³ の ES として Stage 8 の `stage8/coef/e0.5_C2` も読む。
- `stage15/c` の 2 本と `stage16/hf256` は t_end の前に止まった。止まったこと自体が結果で（REPORT_stage15 の 4.1 節と 5 章、REPORT_stage16 の 5.2 節）、256³ の参照は ES-256 の 1 本だけになった。
- Stage 15 の K の事前の判定は「差が決まらない」、Stage 16 で時刻をそろえた事後の判定は「ES が優る」で、どちらも参照が ES-256 だけであることの印を付けた。どちらの判定でも、事前の規則では ES に費用対効果はない。

### 付録

| batch | 数 | 内容 | ランク | 実行ファイル | REPORT |
|---|---|---|---|---|---|
| `stage6q/track` | 3 | `dec_N64_np{1,8,64}`。turb1 の 64³ を 'track' で 400 ステップ。MPI の分割による dt と場の差を見る | 1、8、64 | s6 | Stage 6 の 6.2 節 |
| `stage15/pre` | 4 | `dec_turb1_N64_{H_r1,L_r1}_np{8,64}`。L(1)、H(1) の 10 ステップの分割の試験 | 8、64 | s13d | Stage 15 の 2.3 節 |
| `stage15/pre` | 6 | `dt{1,2}_turb1_N64_{ES,H_r1,L_r1}`。dt と dt/2 での収支の残差 | 16 | s13d | Stage 15 の 2.3 節 |
| `stage8/coef` | 5 | `e<svv_eps0>_C<C_shock>`。turb1 の 128³ の対照 `e0.5_C2` と、svv_eps0 と C_shock を 2 倍、1/2 倍にした 4 本 | 128 | s6 | Stage 8 |
| `stage15/a` | 2 | `turb{1,05}_N64_ES25`（上の 30 本に含む） | 64 | s13d | Stage 15 |

`stage8/coef` の判定は、svv_eps0 の軸も C_shock の軸も「鈍い」（応答の絶対値は最大 0.083、閾値は 0.1）である（REPORT_stage8 の 1 章）。

### 限界の例

| batch | 数 | 内容 | ランク | 実行ファイル | REPORT |
|---|---|---|---|---|---|
| `stage4/turb` | 3 | `turb4_N64`（'stop' で止まった最初の実行）、`turb4_N64_a`、`turb4b_N64` | 64 | s4 | Stage 4 の 8.3 節 |
| `stage5/res` | 2 | `turb4_N128_a`、`turb4b_N128`。M_s = 4 の 128³ | 64 | s4b | Stage 5 の 6 章 |
| `stage5/ms` | 2 | `turb2b_N64`、`turb3b_N64`。M_s = 2、3 の 64³ | 64 | s4b | Stage 5 の 12 章 |
| `stage13d/nomass` | 1 | `e0.5_C2`。turb1 の 128³ を svv_mass = F で | 128 | s13d | Stage 13 |

- M_s = 4 では、ρ が格子 1 点の幅の穴で 0 に向かい、dt が 2e-8 以下に縮んで t_end に届かない。128³ でも、穴は格子 1 点の厚さで深まった（REPORT_stage5 の 6 章）。
- M_s = 2、3 は dt が 3e-6 まで縮み、ユーザーの判断で 1 本目のジョブで打ち切った。この初期場を 64³ で既定の設定で 3τ0 まで走らせきれる M_s の上限は、1 と 2 の間にある（REPORT_stage5 の 12 章）。
- svv_mass = F の実行は 0.466τ0 で不安定になった。1 つのセルで ρ が 0 に向かい、dt が ρ_min^1.5 に比例して縮み、ジョブの予算で止まった。判定のラベルは出ていない。対照は `stage8/coef/e0.5_C2` である（REPORT_stage13 の 1 章と 7 章）。

## 走らせ方

手順は `../README.md` の「走らせ方」による。この試験に固有の点は次のとおりである。

- ケースを ps1 の `runs/` の下の同じパスに写し、`python3 tests/MHD_COMP_ISOTH/stages/regen_inputs.py <case> ...` で `ic/` を書く。256³ の `ic/` は 1 ケース約 0.95 GB である。
- `exe` のない batch は、`qsub -v` の `EXE=` で実行ファイルを渡す。`stage4` は `bin/calliope_isoth_s4`、`stage5` は `_s4b`、`stage6q` と `stage8` は `_s6` で走らせた。
- `stage15/pre` の `dec_*` には `fftw_estimate` があり、FFT の計画を ESTIMATE に固定する。分割の試験でビット単位の一致を見るためである。
- 256³ は 256 ランクで、1 ノードを丸ごと要求する（`-l select=1:ncpus=256:mem=745gb:mpiprocs=256:ompthreads=1`）。`turb1_N256` は 30 分のジョブ 16 本で 3τ0 に達した。同じ一覧で投げ直すと restart から続く。
- 各 Stage の qsub の行（コア数、`PAR`、`BUDGET`）は、REPORT の再現手順の章にある（Stage 5 の 10 章、Stage 6 の 12 章、Stage 8 の 10 章、Stage 13 の 12 章、Stage 15 と 16 の 13 章）。そこにあるパス `tests/MHD_COMP_ISOTH/<file>` は、いまは `tests/MHD_COMP_ISOTH/stages/<file>` にある。

## 解析

ps1 の `~/work/calliope_isoth` で走らせ、`runs/` の下の出力を読む。
Stage 16 の `rejudge` と `judge1` だけは、ps1 から持ち帰った JSON を手元で読む。
ROOT は `runs` である。

| Stage | コマンド（`tests/MHD_COMP_ISOTH/stages/` の下） | 読む実行 |
|---|---|---|
| 4 | `analyze_stage4.py turb` | `stage4/turb` |
| 5 | `analyze_stage5.py res runs`、`analyze_stage5.py ms runs` | `stage4/turb`、`stage5/res`、`stage5/ms` |
| 6 | `analyze_stage6.py track runs/stage6q`、`analyze_stage6.py turb runs/stage6q` | `stage6q`、`stage4/turb` |
| 7 | `analyze_stage7.py tail ROOT N`、`hole ROOT N`、`judge ROOT` | `stage4/turb`、`stage6q/turb` |
| 8 | `analyze_stage8.py rhomin ROOT`、`tail`、`hole`、`judge` | `stage8/coef` |
| 13 | `analyze_stage13.py rhomin ROOT`、`fields`、`judge`、`instab` | `stage13d/nomass`、`stage8/coef/e0.5_C2` |
| 15 | `analyze_stage15.py pre`、`runs <batch>`、`spec`、`judge` | `stage15`、`stage6q/turb`、`stage8/coef/e0.5_C2` |
| 16 | `analyze_stage16.py hf <batch>`、`spec_hf`、`rejudge`、`judge1` | `stage16/hf*`、Stage 15 の `spec.json` |

各コマンドの引数と出力は、スクリプトの docstring と REPORT の再現手順の章にある。

## 図

`plot_turb3d.py` が `fig/<名前>.pdf` を書く。コマンドと fetch は `../README.md` の「図」の節による。
とくに断らなければ turb1（M_s = 1）の図で、時刻は τ0 を単位とする。

入力の多くは、上の解析が書いて手元の `runs/` に置いた要約である（`stage6q/turb/summary.json`、`stage7/tail_N*.json`、`stage8/{rhomin,summary}.json`、`stage13d/rhomin.json`、`stage15/{spec,judge}.json`、`stage16/judge_tm.json`）。
`fetch` は、残りの 3 つを ps1 から持ち帰る（合わせて 14 MB）。
`turb_slice` の面（256³、t = τ0、z = 0 の ρ、ζ、D）と `turb_track` の `cfl.dat` の表（4 本）は ps1 の出力から取り出し、`turb_ms` の `stage5/ana/{res,ms}.json` は写す。
上に挙げた要約と `stage5/ana/{res,ms}.json` は、このフォルダの同じパスに置いた symlink からもたどれる（`../README.md` の「図」の節）。

| 図 | 論文 | 中身 | REPORT |
|---|---|---|---|
| `turb_eps` | 本文 | ε(t) = −dE/dt（実線）と、そのうちの衝撃波粘性の分（破線）。N = 64、128、256 | Stage 6 の 8.3 節 |
| `turb_spec` | 本文 | t = 1、3 の E(k) = E_{√ρu}(k) + E_B(k)（殻の和）。N = 64、128、256。k = N/3 の点線を添える | Stage 6 の 8.4 節 |
| `turb_slice` | 本文 | 256³ の t = 1 の z = 0 の面の ln ρ と、衝撃波粘性の係数 ζ = C_shock Δ² S[max(0, −div u)] | Stage 6 の 8.6 節 |
| `turb_track` | 本文 | dt/dt_cfl（dt_cfl = min(dt_adv, dt_visc)、毎ステップ）。N = 64、128 の Stage 4 の実行（固定の dt と 'decrement'）と 'track' の実行。点線は track_safety（0.8）と track_raise（0.7） | Stage 6 の 7 章 |
| `turb_rho_pdf` | 本文 | t = 0.5、1 の s = ln ρ の体積 PDF。N = 64、128、256。破線は N = 256 の s の平均と標準偏差をもつ Gauss 分布 | Stage 7 の 6.2 節 |
| `turb_K_spec` | 本文 | 窓 W1 = [0.75, 1.25]、W2 = [1.5, 2] で平均した E_B(n) と E_m(n)（m = ρu）。128³ の ES、H(8)（B*）、L(1)、L(2) と、ES-256（灰色） | Stage 15 の 6.1 節 |
| `turb_K` | 本文 | W1、W2 の ln T(n)、T = E_128(n)/E_ES-256(n)。灰色の帯は \|ln T\| ≤ ln 1.25、丸は ES と B* の k_× | Stage 15 の 6.1 節 |
| `turb_kappa` | 本文 | ES と B* の記録ごとの κ。同じ時刻の ES-256 の記録に対して求めた（時刻をそろえた事後の判定）。横棒は窓で求めた κ | Stage 16 の 4.2 節 |
| `turb_coef` | 付録 | Stage 8 の 128³ の 5 本の ρ_min(t) と ln(ρ_min/ρ_min,対照)。灰色は窓 [0.6, 1.2]、横棒は R = L − L_対照（L は窓での ln ρ_min の時間平均）、点線は R = ±0.1 | Stage 8 の 5.1 節、6.3 節 |
| `turb_nomass` | 限界 | svv_mass = F と対照（128³）の ρ_min(t) と、dt と ρ_min の関係。ρ_min^1.5 の線を添える | Stage 13 の 7 章 |
| `turb_ms` | 限界 | 64³ の M_s = 1–4（turb1、turb2b、turb3b、turb4b）と 128³ の M_s = 4（turb4b、破線）の ρ_min と dt/dt0（dt0 は各実行の最初の dt）。点線は ρ_min = 1e-3 と dt0/10（判定 H と S の閾値） | Stage 5 の 6 章、12 章 |

- K の量は REPORT_stage15 の 3 章による。T(n) の |ln T| が初めて ln 1.25 を超える殻を k_× とし、κ = k_×/k_c（k_c = 128/3）とする。判定と、参照が ES-256 だけであることの印は、上の「本文」の節のとおりである。
- 最終版（s13d）の実行は、`turb_K_spec`、`turb_K`、`turb_kappa` の 128³ の実行と、`turb_nomass` の svv_mass = F の実行である。ES-256 と `turb_nomass` の対照は s6、ほかの図は s4、s4b、s6 の実行である。
