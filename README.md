# polymer-brush-tool

表面グラフト型 polymer brush の初期構造を作成し、GROMACS で MD 計算できる状態まで準備するツールです。
Gaussian の RESP 電荷計算結果から、AMBER (antechamber / prepgen / tleap / sander)、PACKMOL、ParmEd、GROMACS を順に呼び出して、溶媒化済みの `.gro` / `.top` を出力します。

対応トポロジー:

- **linear** — 基板—HEAD—(MID)ₙ—TAIL（片末端グラフト）
- **loop** — 基板—HEAD—(MID)ₙ—TAIL—基板（両末端グラフト）

## インストール

```bash
git clone <this repo>
cd polymer_brush_tool
pip install -e .          # 開発用途なら pip install -e ".[dev]"
```

Python 依存: numpy, ase, parmed, pyyaml
外部ツール（PATH に必要）: AmberTools (`antechamber`, `prepgen`, `tleap`, `sander`, `ambpdb`), `packmol`, GROMACS (`gmx`)

## ディレクトリ構成

```
src/polymer_brush_tool/
├── config.py          BrushConfig / MonomerSpec（YAML 読み込み、box 寸法・鎖長の計算）
├── runner.py          外部コマンド実行ラッパー（失敗時に ExternalToolError）
├── cli.py             pbuild コマンド
├── ff/
│   ├── fragments.py   prepgen 用 .chain ファイル、tleap スクリプト生成
│   ├── minimize.py    sander による pull 拘束付き最小化
│   └── topology.py    GROMACS .top へ tip3p include / position_restraints 挿入
├── structure/
│   ├── atoms.py       末端原子・リンカー原子のインデックス検索
│   └── graft.py       PACKMOL による z 軸整列とグリッド配置
└── workflows/
    ├── base.py        共通パイプライン（step_* メソッド）
    ├── linear.py      LinearBrushWorkflow
    └── loop.py        LoopBrushWorkflow
md_template/           本計算用 mdp (min, nvt, anneal, npt) とジョブスクリプト
examples/              mpc.ac, Gaussian log, 設定ファイルのサンプル
tests/                 pytest（外部ツール不要）
```

## 手順

### 1. Gaussian 構造最適化

水系の brush なら PCM / SMD で水の誘電率を使うのが妥当です。

```
%chk=hoge.chk
#P wb97xd/6-311+g(2d,p)
opt freq

Gaussian input

1 1
C   1.1592610000   0.0057950000   0.0020800000
...
```

### 2. 電荷計算と antechamber

GAFF では HF/6-31G(d) の RESP 電荷を使います。

```
%chk=hoge_resp.chk
#p hf/6-31g(d) pop=mk iop(6/33=2,6/42=6)
```

```bash
antechamber -fi gout -i mpc_pcm_resp.log -fo ac  -o mpc.ac  -c resp -pf y -at gaff2
antechamber -fi gout -i mpc_pcm_resp.log -fo pdb -o mpc.pdb -c resp -pf y -at gaff2
```

`mpc.pdb` を VESTA で開き、結合点や除去する水素の原子名を確認しておきます。

### 3. 設定ファイルの作成

```bash
pbuild init --topology linear --output my_brush.yaml   # loop なら --topology loop
```

`examples/linear_config.yaml` / `examples/loop_config.yaml` にコメント付きの例があります。主な項目:

| キー | 意味 |
|---|---|
| `n_mid_repeat_units` | MID モノマーの繰り返し数 |
| `rho`, `nx`, `ny` | グラフト密度 (chains/nm²) と x, y 方向の鎖数。box 寸法はここから自動計算 |
| `xyratio` | PACKMOL で鎖を z 軸に整列させる際の細長い箱の y/x 比。シミュレーション箱は常に正方形（box_y = box_x）で密度は `rho` のまま |
| `t_mpi`, `t_omp` | `gmx mdrun` の並列数（積が CPU コア数） |
| `d_polymer`, `d_cc` | sander 最小化で HEAD–TAIL 末端間に掛ける拘束距離 (Å)。linear は `null` で伸長鎖長を自動計算。**loop では 2 つのグラフト点の間隔を明示指定（旧スクリプトは 14.9）** |
| `linker_height` | リンカー原子を置く wall (z = 0) からの高さ (Å)。基板と共有結合しているとみなし既定は結合長相当の 1.5 |
| `solvent_min_z` | `gmx solvate` 後にこの高さ (Å) より低い水を削除。既定 3.0（水と c3 wall の LJ 接触距離）。0 で削除しない |
| `head` / `mid` / `tail` | 各モノマーの `resname`, `ac_file`, 結合原子名 (`headname`, `tailname`), 除去水素 (`omitnames`), 隣接 GAFF 型 (`pre_headtype`, `post_tailtype`), 末端原子 (`termname`), 主鎖 C-C 数 (`n_cc`) |
| `bottom_atom_index` | 最小化後の鎖で基板側に置く原子の 1-based index。省略時は実行中に対話で入力 |
| `linker_atoms` | 基板に固定する原子のリスト `[{resname: HMP, atomname: H1}]`。省略時は HEAD の `termname` を既定値として対話で確認 |

### 4. 実行

`mpc.ac` を作業ディレクトリに置いて実行します。

```bash
pbuild linear --config my_brush.yaml --work-dir ./out_linear
pbuild loop   --config my_loop.yaml  --work-dir ./out_loop
```

`--mdp-dir` を省略した場合はパッケージ同梱の `min_vac.mdp`, `nvt_vac.mdp`, `tip3p.itp`（`src/polymer_brush_tool/templates/`）を使います。

途中で外部ツールが失敗した場合は、原因を直してから `--start <step>` で途中から再開できます（例: `--start step_graft`）。ステップ名は `pbuild linear --help` を参照してください。

処理内容の詳細な解説は [docs/PIPELINE.md](docs/PIPELINE.md) にあります。

実行中の対話（設定ファイルで指定していない場合のみ）:

1. `chain_min_pull.pdb` を VESTA で確認し、基板側の原子 index を入力
2. リンカー原子の指定。linear は既定値 (Enter または `y`) で通常問題なし。loop は HEAD と TAIL の両末端を `HMP,H1,TMP,H24` のように指定

完了すると次のファイルが出力されます（`tip3p.itp` は作業ディレクトリにコピーされます）。`md_template/` の mdp と同じディレクトリに置いて本計算に進みます。

```
grafted_chain_water_box.gro
grafted_chain_water_restraint.top            (soft restraint 10,000 kJ/mol/nm²)
hardrest_grafted_chain_water_restraint.top   (hard restraint 1,000,000 kJ/mol/nm²)
tip3p.itp
```

まず `grafted_chain_water_box.gro` を VMD で可視化して構造を確認してください。

### Python API

```python
from polymer_brush_tool.config import BrushConfig
from polymer_brush_tool.workflows import LinearBrushWorkflow

config = BrushConfig.from_yaml("my_brush.yaml")
wf = LinearBrushWorkflow(config, work_dir="./out_linear")
wf.run()

# 個別ステップも呼べます
# wf.step_prepgen(); wf.step_build_chain(); wf.step_amber_minimize(); ...
```

## テスト

外部ツールを必要としない純 Python 部分（設定読み込み、ファイル生成、トポロジー編集）を pytest で検証します。

```bash
pip install -e ".[dev]"
pytest
```

## パイプラインの概要

1. `write_fragment` で HEAD / MID / TAIL の `.chain` を書き、`prepgen` で `.prepi` を生成
2. `write_tleap` + `tleap` で鎖を重合し `chain.prmtop` / `chain.inpcrd` を出力
3. `amber_min_with_pull` で末端間距離拘束を掛けて `sander` 最小化（loop は中点拘束を追加）
4. `align_chain_z` で PACKMOL を使い鎖を z 軸に整列
5. `graft_brush` で nx × ny のグリッドに鎖を配置
6. `write_tleap_grafted` + `tleap` で力場を再割り当て、ParmEd で GROMACS 形式へ変換
7. リンカー原子を z = `linker_height`（既定 1.5 Å）に揃え、`insert_restraint_top` で position restraints を追加
8. GROMACS 真空中で最小化 → NVT
9. `gmx solvate` で TIP3P を充填、wall とグラフト点の間（z < `solvent_min_z`）の水を削除、`insert_tip3p_top` でトポロジーに include を追加
