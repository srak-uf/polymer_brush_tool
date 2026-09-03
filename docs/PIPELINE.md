# polymer-brush-tool 処理内容の解説

表面グラフト型 polymer brush の全原子 MD 初期構造を、量子化学計算の出力から GROMACS で走らせられる状態まで自動で組み上げるツールです。本書は「何を、どのツールで、どういう理屈で」やっているかを、共同研究者が中身を追えるようにまとめたものです。

対象読者: GROMACS / AMBER をある程度使ったことがある人。コードを読まなくても各ステップの意味と出力ファイルが分かることを目標にしています。

---

## 1. 全体像

```
Gaussian (構造最適化 + RESP 用 ESP)          ← 手作業
      │  .log
      ▼
antechamber  → mpc.ac (GAFF2 原子タイプ + RESP 電荷)   ← 手作業
      │
      ▼  ここから pbuild が自動化
[1] prepgen      HEAD / MID / TAIL 残基の .prepi を作る
[2] tleap        残基を重合して 1 本の鎖 (chain.prmtop / chain.pdb)
[3] sander       末端間距離を拘束して最小化 → 伸びた鎖 (chain_min_pull.pdb)
[4] PACKMOL      鎖を z 軸方向に立てる (aligned_chain.pdb)
[5] PACKMOL      nx × ny 本を格子状に配置 (grafted_chain.pdb)
[6] tleap+ParmEd 多本系に力場を割り当て GROMACS 形式へ (grafted_chain.top / .gro)
[7] Python       リンカー原子を z=1.5 Å（linker_height）に揃え、position restraints を追加
[8] gmx          真空中で最小化 + 短い NVT（鎖の重なりをほどく）
[9] gmx solvate  TIP3P 水で充填、ボックス高さを調整
      │
      ▼
grafted_chain_water_box.gro + *_restraint.top + tip3p.itp  → 02_MD_template で本計算
```

化学的な設定: 鎖は **基板—HEAD—(MID)ₙ—TAIL** という 3 種の残基で表現します。同じモノマー（例: MPC）でも、結合する側の水素を落として prepgen に渡すため、末端用 HEAD / TAIL と内部用 MID を別残基名（`hmp` / `mmp` / `tmp`）で定義します。

基板は明示的な原子としては置かず、GROMACS の **wall ポテンシャル**（`pbc = xy`, `nwall = 2`, `wall-type = 10-4`）と、リンカー原子の **position restraint** で表現します。

トポロジーは 2 種類:

| topology | 形 | グラフト点 | 特徴 |
|---|---|---|---|
| `linear` | 基板—HEAD—(MID)ₙ—TAIL | HEAD 末端の 1 点 | 鎖を伸長鎖長まで引っ張って立てる |
| `loop` | 基板—HEAD—(MID)ₙ—TAIL—基板 | HEAD と TAIL の 2 点 | 末端間距離を短く拘束し、中点を持ち上げてアーチにする |

---

## 2. 前準備（手作業）

### 2.1 Gaussian

1. 構造最適化: `wb97xd/6-311+g(2d,p) opt freq`。水系なら PCM / SMD で水の誘電率を使う。
2. RESP 用 ESP 計算: GAFF の規約に合わせ `hf/6-31g(d) pop=mk iop(6/33=2,6/42=6)`。

### 2.2 antechamber

```bash
antechamber -fi gout -i mpc_pcm_resp.log -fo ac  -o mpc.ac  -c resp -pf y -at gaff2
antechamber -fi gout -i mpc_pcm_resp.log -fo pdb -o mpc.pdb -c resp -pf y -at gaff2
```

`mpc.ac` が以降の唯一の化学的入力です。`mpc.pdb` は VESTA で開いて、**どの炭素で隣と結合するか、その際に落とす水素はどれか、末端原子はどれか** を原子名で確認するために使います。

---

## 3. 各ステップの詳細

以下、`work_dir` は `pbuild ... --work-dir` で指定した作業ディレクトリです。すべての外部ツールはこのディレクトリをカレントにして実行されます。

### [1] prepgen — 残基定義 (`ff/fragments.py: write_fragment`)

各残基について `.chain` ファイルを書き、`prepgen` で `.prepi` を生成します。

```
HEAD_NAME       C11     ← 前の残基と結合する原子（HEAD にはない）
PRE_HEAD_TYPE   c3      ← 前の残基側の GAFF 原子タイプ
TAIL_NAME       C2      ← 次の残基と結合する原子（TAIL にはない）
POST_TAIL_TYPE  c3
OMIT_NAME       H23     ← 結合形成のために除去する水素（結合点 1 つにつき 1 個）
OMIT_NAME       H24
CHARGE 0.000
```

設定ファイルとの対応:

| 残基 | `headname`/`pre_headtype` | `tailname`/`post_tailtype` | `omitnames` の個数 | `termname` |
|---|---|---|---|---|
| HEAD (`hmp`) | なし | 必須 | 1 | グラフト点の原子（基板側末端） |
| MID (`mmp`) | 必須 | 必須 | 2 | 不要 |
| TAIL (`tmp`) | 必須 | なし | 1 | 自由末端の原子 |

除去した水素の電荷は prepgen 側で残基内に再配分されます（`CHARGE 0.000` は残基全体の目標電荷）。

### [2] tleap — 鎖の重合 (`write_tleap`)

```
source leaprc.gaff2
loadamberprep hmp.prepi / mmp.prepi / tmp.prepi
chain = sequence {HMP MMP MMP ... MMP TMP}
saveAmberParm chain chain.prmtop chain.inpcrd
savepdb chain chain.pdb
```

`n_mid_repeat_units` 個の MID を挟みます。出力の `chain.pdb` は ASE で読み直して書き戻し、書式を正規化しています。

### [3] sander — 拘束付き最小化で鎖を伸ばす (`ff/minimize.py`)

tleap が作る鎖は丸まっているので、**末端間距離に平底ポテンシャルを掛けて気相で最小化**し、伸びた形にします。

AMBER の NMR 拘束（`nmropt=1`, `DISANG`）を使います:

```
&rst
  iat=<head>,<tail>, r1=0., r2=0.8*r3, r3=<d>, r4=1.2*r3, rk2=1000.7, rk3=1000.7,
/
```

- `r2 ≤ r ≤ r3` で力ゼロ、外側で調和ポテンシャル。実質「末端間距離を d の ±20 % に収めろ」という拘束です。
- `head` / `tail` は `chain.pdb` の中で `head.termname` / `tail.termname` に一致する原子の 1-based index。

**d の決め方（`d_polymer`）**

- linear: `d = d_cc × n_cc_all × 0.8`（Å）。`n_cc_all = (head_n_cc+1) + n×(mid_n_cc+1) + (tail_n_cc+1) + 2` は主鎖 C–C 結合数の見積もり、0.8 は結合角による短縮の係数です。ほぼ全伸長に近い長さです。
- loop: **`d_polymer` は 2 つのグラフト点の間隔**（旧スクリプトでは 14.9 Å）で、必ず明示指定します。伸長鎖長を入れると環になりません。加えて、中点原子（index が head と tail の中間の原子）と両末端の距離を `loop_height = d_cc × n_cc_all × 0.8 / 2` に拘束し、アーチを立ち上げます。

最小化は `maxcyc=50000, ncyc=1000, cut=999`（カットオフなし、周期境界なし）。結果は `ambpdb` で `chain_min_pull.pdb` に変換します。

### [4] PACKMOL — 鎖を z 軸に立てる (`structure/graft.py: align_chain_z`)

`chain_min_pull.pdb` を **底面 x × (x·xyratio)、高さ 1000 Å の細長い箱** に押し込み、かつ「基板側の原子」を z = 1.2 Å の平面に固定するよう PACKMOL に要求します。x を 2.1 Å から 0.1 Å 刻みで増やし、**初めて PACKMOL が成功した x を採用**します。細い箱で成功する = 鎖がほぼ真っすぐ z 方向に立っている、という理屈です。

「基板側の原子」の index（`bottom_atom_index`）は、`chain_min_pull.pdb` を VESTA で見て人間が決めます。設定ファイルに書いておけば対話なしで進みます。HEAD の `termname` と同じ原子を指すことが多いですが、最小化後の PDB での番号なので必ず確認してください。

成功判定は PACKMOL のログに `ERROR` が含まれないことで行っています（元スクリプトの挙動をそのまま踏襲）。

### [5] PACKMOL — 格子状にグラフト (`graft_brush`)

グラフト密度 `rho` (chains/nm²) と本数 `nx`, `ny` から底面サイズを決めます:

```
box_x = sqrt(nx·ny / rho) × 10   [Å]
box_y = box_x                     （常に正方形。xyratio は関与しない）
```

各鎖は格子セルの中心 `((ix+0.5)·box_x/nx, (iy+0.5)·box_y/ny)` に、回転も z 移動もなし（`fixed x y 0 0 0 0`）で置かれ、PDB の chain ID を A, B, C, … と振ります。箱の高さは鎖の z 方向の広がり + 30 Å。

`xyratio` は [4] の PACKMOL 整列箱の縦横比にのみ使われ、シミュレーション箱の寸法やグラフト密度は変えません（旧スクリプトと同じ挙動）。loop 例で nx=1, ny=2 とすると、正方形の箱の中で 1 本あたりの占有面積が y 方向に半分の長方形になります。

**注意:** `gmx grompp` は箱の最短辺が `2 × rlist`（同梱 mdp では 2.8 nm）より短いとエラーになります。`sqrt(nx·ny/rho)` が 2.8 nm 以上になるよう `nx`, `ny`, `rho` を選んでください。

### [6] tleap + ParmEd — 多本系の力場と GROMACS 変換 (`write_tleap_grafted`)

PACKMOL の出力 `grafted_chain.pdb`（CRYST1 に箱情報あり）を tleap で読み直し、GAFF2 を割り当てて `grafted_chain.prmtop / .inpcrd` を作ります。ParmEd で `.top` / `.gro` に変換し、`gmx editconf -c` で箱の中心に寄せます。

ParmEd は同一分子をまとめるため、`[ moleculetype ]` は 1 つで `[ molecules ]` に `mol  nx·ny` と書かれます。これは次のステップの前提になります。

### [7] リンカー原子の固定と position restraints (`step_position_restraints`)

1. **リンカー原子**（基板に結合している原子）を決めます。既定は HEAD の `termname`。loop では HEAD と TAIL の両方を `linker_atoms` で指定します。
2. 系全体の最小 z を求め、リンカー原子の z をその値に揃えてから、全体を平行移動してリンカー原子が z = `linker_height`（既定 1.5 Å）に来るようにします（`grafted_chain_shifted.gro`）。
   GROMACS の wall は z = 0 にあるので、`linker_height` は「グラフト原子と基板面の距離」です。基板と共有結合しているとみなして結合長相当の 1.5 Å を既定にしています。
   z = 0 に置いても `wall-r-linpot = 0.3` のおかげで数値的には壊れませんが、リンカーに結合した C 原子（1〜2 Å 上）が 10-4 wall ポテンシャルの斥力芯に入り、
   常に約 300 kJ/mol/nm の反発を受け続けるため、soft restraint の本計算では鎖の根元が持ち上がる方向に働きます。1.5 Å ならこれを避けられます。
3. `[ position_restraints ]` を `[ system ]` の直前に挿入し、2 種類の topology を書きます:

| ファイル | 力の定数 (kJ/mol/nm²) | 用途 |
|---|---|---|
| `hardrest_grafted_chain_restraint.top` | 1,000,000 | 真空中緩和・溶媒充填時にグラフト点を動かさない |
| `grafted_chain_restraint.top` | 10,000 | 本計算用の緩い拘束 |

**注意（ParmEd の分子まとめ由来）**: `[ position_restraints ]` の原子番号は moleculetype 内の番号（1 本目の鎖の番号）でなければなりません。コードは全系の gro から得たリンカー index のうち、`[ atoms ]` に存在する番号だけを書くことでこれを満たしています。鎖が全て同一なら 1 本目のリンカー番号が全鎖に適用されるので正しく動きますが、鎖ごとに異なる分子を混ぜる場合はこの前提が崩れます。

### [8] GROMACS 真空中緩和 (`step_vacuum_relax`)

PACKMOL が並べただけの鎖同士の接触や歪みをほどくため、hard restraint 付きで

1. `min_vac.mdp`: steepest descent、`emtol = 100`, 最大 40,960 ステップ
2. `nvt_vac.mdp`: `dt = 0.5 fs`, 100,000 ステップ (50 ps), v-rescale 300 K

を実行します。両方とも `pbc = xy`, `nwall = 2` の slab 設定です。

### [9] 溶媒充填 (`step_solvate`)

1. `gmx solvate -cp nvt_vac.gro -p hardrest_...top` で TIP3P を充填（`grafted_chain_water_raw.gro`）。`-p` を渡した topology には `SOL  N` が追記されます。
   `gmx solvate` は wall を知らないので、wall (z = 0) とグラフト原子の間のスラブにも水を置きます。O が `solvent_min_z`（既定 3.0 Å ≈ 水と c3 wall の LJ 接触距離）より低い水分子を削除して `grafted_chain_water.gro` とし、`SOL  N` を書き直します。
2. 同じ `[ molecules ]` ブロックを soft 側 topology にもコピーします（**旧スクリプトではここが抜けており、soft 側 topology は原子数不一致で grompp が通らない状態でした**）。
3. 箱の高さを +4 Å して `grafted_chain_water_box.gro` を作ります。`gmx editconf -box` は既定で系を箱の中心に寄せ直してしまう（リンカーが z = `linker_height` から浮く）ため、`-noc` を付けて座標は動かしません。
4. 両 topology の `[ moleculetype ]` の直前に `#include "tip3p.itp"` を挿入し、`tip3p.itp` を作業ディレクトリにコピーします。

---

## 4. 出力ファイルと本計算へのつなぎ

作業ディレクトリに最終的に必要なのは次の 4 点です。

```
grafted_chain_water_box.gro                  座標（溶媒込み、wall z=0、リンカー z=linker_height）
grafted_chain_water_restraint.top            soft restraint (10,000)
hardrest_grafted_chain_water_restraint.top   hard restraint (1,000,000)
tip3p.itp
```

これを `02_MD_template/` の mdp と同じ場所に置き、

| mdp | 内容 |
|---|---|
| `min.mdp` | 溶媒込み最小化、all-bonds 拘束 |
| `nvt.mdp` | NVT 平衡化 (dt 1 fs, 200 ps) |
| `anneal.mdp` | 300→600→300 K のアニール、残基種ごとに温度制御 |
| `npt.mdp` | 本計算 NPT (dt 2 fs, 200 ns), Nosé–Hoover + Parrinello–Rahman, z 方向のみ圧縮 (semi-isotropic) |

の順で流します。`relax_local.sh`（ローカル）と `mdrun.sh`（SGE, Intel MPI 36 rank）が例です。まず `grafted_chain_water_box.gro` を VMD で見て、鎖が立っているか・水が抜けていないかを確認してください。

---

## 5. 設定ファイル（YAML）の項目

| キー | 意味 | linear 既定 | loop 既定 |
|---|---|---|---|
| `topology` | `linear` / `loop` | linear | loop |
| `n_mid_repeat_units` | MID の繰り返し数 | 12 | 26 |
| `rho` | グラフト密度 chains/nm² | 0.45 | 0.225 |
| `nx`, `ny` | x, y 方向の鎖数 | 2, 2 | 1, 2 |
| `xyratio` | PACKMOL 整列箱の y/x 比（シミュレーション箱には影響しない） | 1.0 | 0.5 |
| `d_cc` | 主鎖 C–C 結合長 Å | 1.54 | 1.54 |
| `d_polymer` | sander の HEAD–TAIL 拘束距離 Å | null（自動: 伸長鎖長） | **14.9（必須）** |
| `t_mpi`, `t_omp` | `gmx mdrun -ntmpi/-ntomp` | 8, 1 | 8, 1 |
| `head`/`mid`/`tail` | 残基定義（§3 [1] 参照） | | |
| `bottom_atom_index` | [4] で基板側に置く原子（1-based） | 対話 | 対話 |
| `linker_atoms` | [7] で基板面に固定・拘束する原子 | HEAD の termname | 対話（HEAD と TAIL を指定） |
| `linker_height` | [7] でリンカー原子を置く wall からの高さ Å | 1.5 | 1.5 |
| `solvent_min_z` | [9] でこの高さ (Å) 未満の水を削除 | 3.0 | 3.0 |

`examples/linear_config.yaml`, `examples/loop_config.yaml` にコメント付きの完全な例があります。`pbuild init --topology loop` で同じものが生成されます。

---

## 6. 実行方法

```bash
pip install -e .                       # ase, numpy, parmed, pyyaml
pbuild init --topology linear -o my_brush.yaml
# my_brush.yaml を編集、mpc.ac を作業ディレクトリへ
pbuild linear --config my_brush.yaml --work-dir ./out_linear
```

途中で失敗したら、原因を直して `--start <step名>` で再開できます（例: PACKMOL のスキャンが失敗した後に `bottom_atom_index` を直して `--start step_align_z`）。ステップ名は `step_prepgen, step_build_chain, step_amber_minimize, step_align_z, step_graft, step_assign_ff_grafted, step_position_restraints, step_vacuum_relax, step_solvate`。

Python から:

```python
from polymer_brush_tool.config import BrushConfig
from polymer_brush_tool.workflows import LinearBrushWorkflow

cfg = BrushConfig.from_yaml("my_brush.yaml")
wf = LinearBrushWorkflow(cfg, work_dir="out_linear")
wf.run()                      # 全部
wf.run(start="step_graft")    # 途中から
wf.step_align_z()             # 1 ステップだけ
```

旧来の `01_FF_template/prep_chain_linear.py` / `prep_chain_loop.py` も残してあり、冒頭の変数を編集して実行できます（中身はライブラリ呼び出し）。

---

## 7. 既知の注意点・限界

1. **`bottom_atom_index` は人が決める必要があります。** 最小化後の PDB での番号なので、モノマーや鎖長を変えたら毎回確認してください。
2. **loop の `d_polymer` は末端間距離です。** 自動計算に任せると伸長鎖長になり、ループになりません（設定読み込み時にエラーで止まります）。
3. **PACKMOL の成功判定はログ中の `ERROR` 文字列**で行っています。バージョンによって "ENDED WITHOUT PERFECT PACKING" のような出力で終わる場合、成功扱いになる可能性があります。`aligned_chain.pdb` を目視確認してください。
4. **鎖は全て同一分子である前提**（§3 [7] 参照）。組成の異なる鎖を混ぜる場合は position restraints の付け方を変える必要があります。
5. **`(head_idx + tail_idx) // 2` を「中点原子」としています。** 原子番号は残基順に振られるので概ね鎖の中央ですが、厳密に中央のモノマーではありません。
6. **電荷**: prepgen に `CHARGE 0.000` を渡しているため、各残基は中性に再配分されます。荷電モノマーを扱う場合はここを変える必要があります。
7. 外部ツール（AmberTools, PACKMOL, GROMACS）が PATH にあること、`parmed` が import できること（AmberTools 同梱版か conda 推奨）が前提です。

---

## 8. 元スクリプトからの変更点（レビューで修正したもの）

ライブラリ化に伴い、挙動を変えた・直した箇所です。

| 箇所 | 元の挙動 | 現在 |
|---|---|---|
| soft restraint topology | `gmx solvate -p` を hard 側にしか掛けておらず、soft 側に `SOL N` が無い | hard 側の `[ molecules ]` を soft 側へコピー |
| 外部コマンド失敗 | `os.system` で無視して続行 | 非ゼロ終了で例外、`--start` で再開 |
| 対話入力 | 必須 | YAML に書けば無人実行可 |
| loop の `d_polymer` | スクリプト内で 14.9 を直書き | YAML で必須項目として明示（未指定はエラー） |
| PACKMOL 入力のパス | カレント前提 | 作業ディレクトリ相対に統一（長い絶対パスは PACKMOL の文字数制限に掛かるため） |
| コード重複 | linear / loop で 6 関数を丸ごとコピー | 共通化。差分は sander 拘束の有無のみ |
| リンカーの高さ | z = 0（wall 上）に置き、溶媒化後の `editconf -box` の再センタリングで約 3 Å 浮いていた | `linker_height`（既定 1.5 Å）に明示配置し、`editconf -noc` で溶媒化後も保持 |
| wall 直上の水 | `gmx solvate` が wall とグラフト点の間に置いた水がそのまま残っていた | `solvent_min_z`（既定 3.0 Å）未満の水を削除し `SOL N` を修正 |
