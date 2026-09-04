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
from importlib import resources
from pathlib import Path

from polymer_brush_tool import __version__
from polymer_brush_tool.config import BrushConfig, ConfigError
from polymer_brush_tool.runner import ExternalToolError
from polymer_brush_tool.workflows.base import BrushWorkflowBase


def _template_text(name: str) -> str:
    return (resources.files("polymer_brush_tool") / "templates" / name).read_text()


# --------------------------------------------------------------------------
# linear / loop
# --------------------------------------------------------------------------

def _run_workflow(args: argparse.Namespace) -> int:
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"ERROR: config file not found: {config_path}", file=sys.stderr)
        return 1

    try:
        config = BrushConfig.from_yaml(config_path)
    except ConfigError as exc:
        print(f"ERROR: invalid config: {exc}", file=sys.stderr)
        return 1

    if config.topology != args.topology:
        print(
            f"ERROR: config declares topology='{config.topology}' but the "
            f"'{args.topology}' sub-command was used. Use 'pbuild {config.topology}'.",
            file=sys.stderr,
        )
        return 1

    if args.topology == "linear":
        from polymer_brush_tool.workflows.linear import LinearBrushWorkflow as WF
    else:
        from polymer_brush_tool.workflows.loop import LoopBrushWorkflow as WF

    try:
        WF(config, work_dir=args.work_dir, mdp_dir=args.mdp_dir).run(start=args.start)
    except ExternalToolError as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        print(
            "Fix the problem, then resume with --start <step> "
            f"(steps: {', '.join(BrushWorkflowBase.STEPS)})",
            file=sys.stderr,
        )
        return 1
    return 0


# --------------------------------------------------------------------------
# init
# --------------------------------------------------------------------------

def _cmd_init(args: argparse.Namespace) -> int:
    output = Path(args.output)
    if output.exists() and not args.force:
        print(f"ERROR: {output} already exists. Use --force to overwrite.", file=sys.stderr)
        return 1
    output.write_text(_template_text(f"{args.topology}_config.yaml"))
    print(f"Template configuration written to: {output}")
    print(f"Edit it, place your .ac file next to it, then run:\n  pbuild {args.topology} --config {output}")
    return 0


# --------------------------------------------------------------------------
# parser
# --------------------------------------------------------------------------

def _add_run_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--config", "-c", required=True, metavar="YAML", help="YAML configuration file.")
    p.add_argument("--work-dir", "-d", default=".", metavar="DIR",
                   help="Working directory; must contain the .ac file(s) (default: cwd).")
    p.add_argument("--mdp-dir", default=None, metavar="DIR",
                   help="Directory with min_vac.mdp / nvt_vac.mdp / tip3p.itp (default: bundled).")
    p.add_argument("--start", default=None, metavar="STEP", choices=BrushWorkflowBase.STEPS,
                   help="Resume from this step, reusing earlier outputs in the work dir.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pbuild",
        description=(
            "polymer-brush-tool: build GROMACS-ready polymer brush systems.\n\n"
            "antechamber(.ac) -> prepgen/tleap -> sander stretch -> PACKMOL align+graft\n"
            "-> tleap/ParmEd -> GROMACS vacuum relax -> gmx solvate"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="COMMAND", required=True)

    p_linear = sub.add_parser("linear", help="Run the linear brush pipeline.")
    _add_run_args(p_linear)
    p_linear.set_defaults(func=_run_workflow, topology="linear")

    p_loop = sub.add_parser("loop", help="Run the loop brush pipeline.")
    _add_run_args(p_loop)
    p_loop.set_defaults(func=_run_workflow, topology="loop")

    p_init = sub.add_parser("init", help="Write a template YAML configuration file.")
    p_init.add_argument("--topology", "-t", choices=["linear", "loop"], default="linear")
    p_init.add_argument("--output", "-o", default="brush_config.yaml", metavar="FILE")
    p_init.add_argument("--force", action="store_true", help="Overwrite an existing file.")
    p_init.set_defaults(func=_cmd_init)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
