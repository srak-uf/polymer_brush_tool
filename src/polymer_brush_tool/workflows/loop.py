"""Workflow for loop polymer brushes (both ends grafted)."""

from __future__ import annotations

from pathlib import Path

from polymer_brush_tool.config import BrushConfig
from polymer_brush_tool.workflows.base import BrushWorkflowBase


class LoopBrushWorkflow(BrushWorkflowBase):
    """Preparation pipeline for a loop brush: substrate—HEAD—(MID)ₙ—TAIL—substrate.

    The sander stretch uses three restraints: HEAD–TAIL at
    ``config.d_polymer`` (the separation of the two grafting points, which
    must be given explicitly) and HEAD–MID / MID–TAIL at
    ``config.loop_height()`` to raise the arch.  Both terminal atoms should
    be listed in ``config.linker_atoms``.
    """

    def __init__(
        self,
        config: BrushConfig,
        work_dir: str | Path = ".",
        mdp_dir: str | Path | None = None,
    ) -> None:
        if config.topology != "loop":
            raise ValueError(
                f"LoopBrushWorkflow requires topology='loop', got {config.topology!r}."
            )
        super().__init__(config, work_dir, mdp_dir)
