"""Workflow for linear surface-grafted polymer brushes."""

from __future__ import annotations

from pathlib import Path

from polymer_brush_tool.config import BrushConfig
from polymer_brush_tool.ff import minimize
from polymer_brush_tool.structure.atoms import find_atom_index
from polymer_brush_tool.workflows.base import BrushWorkflowBase


class LinearBrushWorkflow(BrushWorkflowBase):
    """Complete preparation pipeline for a linear polymer brush system.

    The chain topology is: substrate──HEAD──(MID)ₙ──TAIL (one grafting end).

    Parameters
    ----------
    config:
        Brush configuration.  ``config.topology`` must be ``"linear"``.
    work_dir:
        Output directory.
    mdp_dir:
        Directory containing the GROMACS MDP template files.
    """

    def __init__(
        self,
        config: BrushConfig,
        work_dir: str | Path = ".",
        mdp_dir: str | Path | None = None,
    ) -> None:
        if config.topology != "linear":
            raise ValueError(
                f"LinearBrushWorkflow requires topology='linear', "
                f"got {config.topology!r}."
            )
        super().__init__(config, work_dir, mdp_dir)

    # ------------------------------------------------------------------
    # Step 3 – AMBER minimisation (linear variant)
    # ------------------------------------------------------------------

    def step_amber_minimize(self) -> None:
        """Run AMBER sander with a single end-to-end pull constraint."""
        from ase.io import read

        cfg = self.config
        print("\n=== Step 3: AMBER minimisation with pull constraint (linear) ===")

        chain_pdb = self.work_dir / "chain.pdb"
        atoms = read(str(chain_pdb))

        head_idx_0 = find_atom_index(
            atoms, cfg.head.termname, cfg.head.resname.upper()
        )
        tail_idx_0 = find_atom_index(
            atoms, cfg.tail.termname, cfg.tail.resname.upper()
        )

        minimize.amber_min_with_pull(
            self.work_dir / "chain.prmtop",
            self.work_dir / "chain.inpcrd",
            file_prefix="chain",
            head_idx=head_idx_0 + 1,   # 1-based for AMBER
            tail_idx=tail_idx_0 + 1,
            polymer_length=cfg.polymer_length(),
            work_dir=self.work_dir,
        )

    # ------------------------------------------------------------------
    # Full pipeline
    # ------------------------------------------------------------------

    def run(self) -> None:
        """Execute the complete linear brush preparation pipeline."""
        self.step_prepgen()
        self.step_build_chain()
        self.step_amber_minimize()
        self.step_align_z()
        self.step_graft()
        self.step_assign_ff_grafted()
        self.step_position_restraints()
        self.step_vacuum_relax()
        self.step_solvate()
