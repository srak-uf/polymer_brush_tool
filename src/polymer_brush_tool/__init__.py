"""
polymer_brush_tool
==================

Polymer brush structure modeling and GROMACS MD preparation toolkit.

Workflow:
  1. Gaussian geometry optimization + RESP charge calculation (manual)
  2. Force-field preparation: fragment files → tleap → AMBER minimization
  3. Structure generation: z-axis alignment → PACKMOL graft placement
  4. Topology conversion: AMBER → GROMACS, solvation, position restraints
  5. GROMACS vacuum relaxation → solvated MD-ready system

Usage (CLI)::

    pbuild linear --config my_mpc_linear.yaml
    pbuild loop   --config my_mpc_loop.yaml
    pbuild init   --topology linear --output my_config.yaml

Usage (Python API)::

    from polymer_brush_tool.config import BrushConfig
    from polymer_brush_tool.workflows import LinearBrushWorkflow

    config = BrushConfig.from_yaml("my_config.yaml")
    workflow = LinearBrushWorkflow(config, work_dir="./output")
    workflow.run()
"""

__version__ = "0.1.0"
