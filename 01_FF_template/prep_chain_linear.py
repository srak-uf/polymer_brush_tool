"""
prep_chain_linear.py  —  後方互換ラッパー (backward-compatible wrapper)
=======================================================================

このスクリプトは旧来の直接実行方式との互換性のために残されています。
新規利用には pbuild CLI または Python API を推奨します。

新しい使い方 (recommended):
    # 設定ファイルを生成してから編集
    pbuild init --topology linear --output my_config.yaml
    pbuild linear --config my_config.yaml

Python API:
    from polymer_brush_tool.config import BrushConfig
    from polymer_brush_tool.workflows import LinearBrushWorkflow

    config = BrushConfig.from_yaml("my_config.yaml")
    LinearBrushWorkflow(config, work_dir="./output").run()

以下のパラメータをここで直接編集して `python3 prep_chain_linear.py` でも実行できます。
"""

import sys
from pathlib import Path

# ============================================================
# ユーザー設定 — 必要に応じて編集してください
# ============================================================

# GROMACS 並列数 (t_mpi × t_omp = CPU コア数)
t_mpi = 8
t_omp = 1

# x-y 比（正方形なら 1.0）
xyratio = 1.0

# ループ topology を使う場合は True
flag_loop = False

# HEAD モノマー定義
head_tailname   = "C2"
head_omitnames  = ["H24"]
head_post_tailtype = "c3"
head_termname   = "H1"
head_resname    = "hmp"
head_acfile     = "mpc.ac"
head_n_cc       = 1

# MID モノマー定義
mid_headname    = "C11"
mid_tailname    = "C2"
mid_omitnames   = ["H23", "H24"]
mid_pre_headtype = "c3"
mid_post_tailtype = "c3"
n_mid_repeat_units = 12
mid_resname     = "mmp"
mid_acfile      = "mpc.ac"
mid_n_cc        = 1

# TAIL モノマー定義
tail_headname   = "C11"
tail_omitnames  = ["H23"]
tail_pre_headtype = "c3"
tail_termname   = "H24"
tail_resname    = "tmp"
tail_acfile     = "mpc.ac"
tail_n_cc       = 1

# グラフト密度
rho = 0.45      # chains/nm²
nx = 2
ny = 2

# チェーン長（None で自動計算）
d_polymer = None
d_cc      = 1.54

# 対話プロンプトを省略する場合は値を設定（None で対話モード）
bottom_atom_index = None  # 例: 37
linker_atoms      = None  # 例: [{"resname": "HMP", "atomname": "H1"}]

# ============================================================
# 実行（編集不要）
# ============================================================

if __name__ == "__main__":
    try:
        from polymer_brush_tool.config import BrushConfig, MonomerSpec
        from polymer_brush_tool.workflows.linear import LinearBrushWorkflow
    except ImportError:
        print(
            "ERROR: polymer_brush_tool パッケージがインストールされていません。\n"
            "  pip install -e <repo_root> を実行してください。",
            file=sys.stderr,
        )
        sys.exit(1)

    config = BrushConfig(
        topology="linear",
        n_mid_repeat_units=n_mid_repeat_units,
        rho=rho,
        nx=nx,
        ny=ny,
        xyratio=xyratio,
        t_mpi=t_mpi,
        t_omp=t_omp,
        d_polymer=d_polymer,
        d_cc=d_cc,
        bottom_atom_index=bottom_atom_index,
        linker_atoms=linker_atoms,
        head=MonomerSpec(
            resname=head_resname,
            ac_file=head_acfile,
            termname=head_termname,
            tailname=head_tailname,
            omitnames=head_omitnames,
            post_tailtype=head_post_tailtype,
            n_cc=head_n_cc,
        ),
        mid=MonomerSpec(
            resname=mid_resname,
            ac_file=mid_acfile,
            headname=mid_headname,
            tailname=mid_tailname,
            omitnames=mid_omitnames,
            pre_headtype=mid_pre_headtype,
            post_tailtype=mid_post_tailtype,
            n_cc=mid_n_cc,
        ),
        tail=MonomerSpec(
            resname=tail_resname,
            ac_file=tail_acfile,
            termname=tail_termname,
            headname=tail_headname,
            omitnames=tail_omitnames,
            pre_headtype=tail_pre_headtype,
            n_cc=tail_n_cc,
        ),
    )

    mdp_dir = Path(__file__).parent  # MDP templates are in this directory
    LinearBrushWorkflow(config, work_dir=".", mdp_dir=mdp_dir).run()
