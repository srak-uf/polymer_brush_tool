"""Command-line interface for polymer_brush_tool.

Entry point: ``pbuild`` (defined in pyproject.toml [project.scripts]).

Sub-commands
------------
linear  Run the linear brush preparation pipeline.
loop    Run the loop brush preparation pipeline.
init    Write a template YAML configuration file.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


# --------------------------------------------------------------------------
# Sub-command handlers
# --------------------------------------------------------------------------

def _run_workflow(args: argparse.Namespace) -> int:
    from polymer_brush_tool.config import BrushConfig

    config_path = Path(args.config)
    if not config_path.exists():
        print(f"ERROR: config file not found: {config_path}", file=sys.stderr)
        return 1

    try:
        config = BrushConfig.from_yaml(config_path)
    except Exception as exc:
        print(f"ERROR: failed to load config: {exc}", file=sys.stderr)
        return 1

    # Override topology from sub-command name so the YAML topology field
    # is optional when the sub-command already specifies it.
    config.topology = args.topology

    work_dir = Path(args.work_dir)
    mdp_dir = Path(args.mdp_dir) if args.mdp_dir else None

    try:
        if args.topology == "linear":
            from polymer_brush_tool.workflows.linear import LinearBrushWorkflow
            workflow = LinearBrushWorkflow(config, work_dir=work_dir, mdp_dir=mdp_dir)
        else:
            from polymer_brush_tool.workflows.loop import LoopBrushWorkflow
            workflow = LoopBrushWorkflow(config, work_dir=work_dir, mdp_dir=mdp_dir)

        workflow.run()
    except Exception as exc:
        print(f"ERROR: workflow failed: {exc}", file=sys.stderr)
        return 1

    return 0


_LINEAR_CONFIG_TEMPLATE = """\
# polymer-brush-tool configuration — linear topology
# -------------------------------------------------------
# Edit this file and run: pbuild linear --config <this_file>

topology: linear
n_mid_repeat_units: 12

# Graft density and box dimensions
rho: 0.45            # chains/nm²
nx: 2                # chains along x
ny: 2                # chains along y
xyratio: 1.0         # box_y / box_x

# GROMACS parallelism (t_mpi × t_omp = CPU cores)
t_mpi: 8
t_omp: 1

# HEAD monomer (substrate-facing end)
head:
  resname: hmp
  ac_file: mpc.ac
  termname: H1           # terminal atom (grafting point)
  tailname: C2           # atom bonded to MID
  omitnames: [H24]       # hydrogens removed at the bonding point
  post_tailtype: c3      # GAFF atom type on the MID side

# MID monomer (repeat unit)
mid:
  resname: mmp
  ac_file: mpc.ac
  headname: C11
  tailname: C2
  omitnames: [H23, H24]
  pre_headtype: c3
  post_tailtype: c3
  n_cc: 1

# TAIL monomer (free end)
tail:
  resname: tmp
  ac_file: mpc.ac
  termname: H24          # terminal atom (free end)
  headname: C11
  omitnames: [H23]
  pre_headtype: c3
  n_cc: 1

# Optionally override the polymer length (Å).
# Leave as null to auto-calculate from bond counts.
d_polymer: null
d_cc: 1.54

# Set these to skip interactive prompts during the run.
# bottom_atom_index: 37    # 1-based index in chain_min_pull.pdb
# linker_atoms:
#   - {resname: HMP, atomname: H1}
"""

_LOOP_CONFIG_TEMPLATE = """\
# polymer-brush-tool configuration — loop topology
# -------------------------------------------------------
# Edit this file and run: pbuild loop --config <this_file>

topology: loop
n_mid_repeat_units: 26    # longer chain for loop geometry

# Graft density and box dimensions
rho: 0.225               # chains/nm² (halved for loop = 0.45/2)
nx: 2
ny: 2
xyratio: 0.5             # rectangular box for loop topology

# GROMACS parallelism
t_mpi: 8
t_omp: 1

# HEAD monomer
head:
  resname: hmp
  ac_file: mpc.ac
  termname: H1
  tailname: C2
  omitnames: [H24]
  post_tailtype: c3

# MID monomer
mid:
  resname: mmp
  ac_file: mpc.ac
  headname: C11
  tailname: C2
  omitnames: [H23, H24]
  pre_headtype: c3
  post_tailtype: c3
  n_cc: 1

# TAIL monomer
tail:
  resname: tmp
  ac_file: mpc.ac
  termname: H24
  headname: C11
  omitnames: [H23]
  pre_headtype: c3
  n_cc: 1

d_polymer: null
d_cc: 1.54

# Set these to skip interactive prompts.
# bottom_atom_index: 37
# linker_atoms:
#   - {resname: HMP, atomname: H1}
#   - {resname: TMP, atomname: H24}
"""


def _cmd_init(args: argparse.Namespace) -> int:
    topology = args.topology
    output = Path(args.output)

    if topology == "linear":
        content = _LINEAR_CONFIG_TEMPLATE
    else:
        content = _LOOP_CONFIG_TEMPLATE

    if output.exists() and not args.force:
        print(
            f"ERROR: {output} already exists. Use --force to overwrite.",
            file=sys.stderr,
        )
        return 1

    output.write_text(content)
    print(f"Template configuration written to: {output}")
    print("Edit the file and run:")
    print(f"  pbuild {topology} --config {output}")
    return 0


# --------------------------------------------------------------------------
# Argument parser
# --------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pbuild",
        description=(
            "polymer-brush-tool: build GROMACS-ready polymer brush systems.\n\n"
            "Workflow: antechamber → prepgen/tleap → AMBER minimize → "
            "PACKMOL graft → ParmEd convert → GROMACS solvate"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {_get_version()}",
    )

    sub = parser.add_subparsers(dest="command", metavar="COMMAND")
    sub.required = True

    # -- linear -----------------------------------------------------------
    p_linear = sub.add_parser(
        "linear",
        help="Run the linear brush preparation pipeline.",
    )
    p_linear.add_argument(
        "--config", "-c",
        required=True,
        metavar="YAML",
        help="Path to the YAML configuration file.",
    )
    p_linear.add_argument(
        "--work-dir", "-d",
        default=".",
        metavar="DIR",
        help="Working directory for all output files (default: current dir).",
    )
    p_linear.add_argument(
        "--mdp-dir",
        default=None,
        metavar="DIR",
        help="Directory containing min_vac.mdp and nvt_vac.mdp templates.",
    )
    p_linear.set_defaults(func=_run_workflow, topology="linear")

    # -- loop -------------------------------------------------------------
    p_loop = sub.add_parser(
        "loop",
        help="Run the loop brush preparation pipeline.",
    )
    p_loop.add_argument("--config", "-c", required=True, metavar="YAML")
    p_loop.add_argument("--work-dir", "-d", default=".", metavar="DIR")
    p_loop.add_argument("--mdp-dir", default=None, metavar="DIR")
    p_loop.set_defaults(func=_run_workflow, topology="loop")

    # -- init -------------------------------------------------------------
    p_init = sub.add_parser(
        "init",
        help="Write a template YAML configuration file.",
    )
    p_init.add_argument(
        "--topology", "-t",
        choices=["linear", "loop"],
        default="linear",
        help="Topology for the template (default: linear).",
    )
    p_init.add_argument(
        "--output", "-o",
        default="brush_config.yaml",
        metavar="FILE",
        help="Output YAML file name (default: brush_config.yaml).",
    )
    p_init.add_argument(
        "--force",
        action="store_true",
        help="Overwrite output file if it already exists.",
    )
    p_init.set_defaults(func=_cmd_init)

    return parser


def _get_version() -> str:
    try:
        from polymer_brush_tool import __version__
        return __version__
    except ImportError:
        return "unknown"


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
