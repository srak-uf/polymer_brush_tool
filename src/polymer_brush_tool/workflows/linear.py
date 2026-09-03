"""Workflow for linear surface-grafted polymer brushes."""

from __future__ import annotations

from pathlib import Path

from polymer_brush_tool.config import BrushConfig
from polymer_brush_tool.workflows.base import BrushWorkflowBase


class LinearBrushWorkflow(BrushWorkflowBase):
    """Preparation pipeline for a linear brush: substrate—HEAD—(MID)ₙ—TAIL.

    The sander stretch uses a single HEAD–TAIL restraint at the extended
    contour length (``config.polymer_length()``).  Everything else is shared
    with :class:`~polymer_brush_tool.workflows.base.BrushWorkflowBase`.
    """

    def __init__(
        self,
        config: BrushConfig,
        work_dir: str | Path = ".",
        mdp_dir: str | Path | None = None,
    ) -> None:
        if config.topology != "linear":
            raise ValueError(
                f"LinearBrushWorkflow requires topology='linear', got {config.topology!r}."
            )
        super().__init__(config, work_dir, mdp_dir)
