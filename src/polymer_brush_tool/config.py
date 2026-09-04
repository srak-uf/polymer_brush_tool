"""Configuration dataclasses and YAML loading for polymer brush workflows."""

from __future__ import annotations

from dataclasses import dataclass, fields
from pathlib import Path
from typing import Optional

import numpy as np
import yaml

TOPOLOGIES = ("linear", "loop")


class ConfigError(ValueError):
    """Raised when a configuration is missing required fields or is inconsistent."""


@dataclass
class MonomerSpec:
    """Specification for one monomer residue (HEAD, MID, or TAIL)."""

    resname: str
    """3-letter residue name used by tleap (e.g. 'hmp', 'mmp', 'tmp')."""

    ac_file: str
    """Path (relative to the working directory) to the antechamber AC file."""

    omitnames: list[str]
    """Atom names of hydrogens removed at bonding points (one per bond endpoint)."""

    n_cc: int = 1
    """Number of backbone C-C bonds in this monomer unit."""

    headname: Optional[str] = None
    """Atom bonded to the previous residue (None for HEAD)."""

    tailname: Optional[str] = None
    """Atom bonded to the next residue (None for TAIL)."""

    pre_headtype: Optional[str] = None
    """GAFF atom type of the previous residue's tail atom (None for HEAD)."""

    post_tailtype: Optional[str] = None
    """GAFF atom type of the next residue's head atom (None for TAIL)."""

    termname: Optional[str] = None
    """Terminal atom used to locate the chain end (HEAD and TAIL only)."""


@dataclass
class BrushConfig:
    """
    Complete configuration for a polymer brush simulation system.

    Parameters are equivalent to the hard-coded variables at the top of the
    legacy ``prep_chain_linear.py`` / ``prep_chain_loop.py`` scripts.
    """

    # ------------------------------------------------------------------
    # Chain topology
    # ------------------------------------------------------------------
    topology: str = "linear"
    """Chain topology: 'linear' or 'loop'."""

    n_mid_repeat_units: int = 12
    """Number of MID monomer repeat units in the chain."""

    # ------------------------------------------------------------------
    # Graft density and box geometry
    # ------------------------------------------------------------------
    rho: float = 0.45
    """Graft density in chains/nm²."""

    nx: int = 2
    """Number of chains along the x direction."""

    ny: int = 2
    """Number of chains along the y direction."""

    xyratio: float = 1.0
    """Aspect ratio of the PACKMOL *alignment* box used in ``align_chain_z``
    (y extent = x extent × xyratio).  It does **not** change the simulation
    box: ``box_y == box_x`` always, so the graft density stays exactly ``rho``
    (legacy prep_chain_*.py behaviour)."""

    # ------------------------------------------------------------------
    # Monomer definitions
    # ------------------------------------------------------------------
    head: Optional[MonomerSpec] = None
    mid: Optional[MonomerSpec] = None
    tail: Optional[MonomerSpec] = None

    # ------------------------------------------------------------------
    # Chain length
    # ------------------------------------------------------------------
    d_polymer: Optional[float] = None
    """
    Target HEAD–TAIL terminal distance in Ångström used as the sander
    restraint during the pull minimisation.

    * linear: leave ``None`` to auto-calculate the extended contour length
      ``d_cc * n_cc_all * 0.8``.
    * loop:   set explicitly to the desired separation of the two grafting
      points (the legacy script used 14.9 Å).  The arch height is derived
      separately by :meth:`loop_height`.
    """

    d_cc: float = 1.54
    """C-C backbone bond length in Ångström."""

    # ------------------------------------------------------------------
    # Interactive prompt overrides
    # ------------------------------------------------------------------
    bottom_atom_index: Optional[int] = None
    """
    1-based index of the grafting (bottom) atom in ``chain_min_pull.pdb``.
    When None the workflow prompts interactively (legacy behaviour).
    """

    linker_atoms: Optional[list[dict]] = None
    """
    Atoms pinned to the substrate plane and position-restrained, e.g.
    ``[{"resname": "HMP", "atomname": "H1"}]``.  When None the HEAD
    terminal atom is offered as the default in an interactive prompt.
    """

    linker_height: float = 1.5
    """
    Height (Å) of the linker atoms above the GROMACS wall at z = 0.

    The linker is treated as covalently bonded to the substrate, so the
    default is a bond length (1.5 Å) rather than 0.  Sitting exactly on the
    wall is numerically harmless thanks to ``wall-r-linpot``, but it puts
    the linker's bonded neighbours inside the repulsive core of the 10-4
    wall potential; ~1.5 Å keeps them in the attractive region.  This value
    is preserved through solvation (``gmx editconf -noc``), so the position
    restraint reference and the wall geometry stay consistent.
    """

    solvent_min_z: float = 3.0
    """
    Water molecules whose oxygen lies below this height (Å) after
    ``gmx solvate`` are deleted.

    ``gmx solvate`` fills the whole box, including the slab between the wall
    (z = 0) and the grafting atoms.  3.0 Å is roughly the LJ contact distance
    of a TIP3P oxygen on the c3 wall (σ ≈ 3.27 Å); nothing physical fits
    between the wall and that height, so such water would only be pushed out
    through the brush during equilibration.  Set to 0 to keep all water.
    """

    # ------------------------------------------------------------------
    # GROMACS thread parallelism
    # ------------------------------------------------------------------
    t_mpi: int = 8
    """Thread-MPI ranks for ``gmx mdrun -ntmpi``."""

    t_omp: int = 1
    """OpenMP threads per rank for ``gmx mdrun -ntomp``."""

    # ------------------------------------------------------------------
    # Loading / validation
    # ------------------------------------------------------------------

    @classmethod
    def from_yaml(cls, path: str | Path) -> "BrushConfig":
        """Load a BrushConfig from a YAML file and validate it."""
        with open(path) as fh:
            raw = yaml.safe_load(fh) or {}

        known = {f.name for f in fields(cls)}
        unknown = sorted(set(raw) - known)
        if unknown:
            raise ConfigError(
                f"Unknown key(s) in {path}: {', '.join(unknown)}. "
                f"Valid keys: {', '.join(sorted(known))}"
            )

        for key in ("head", "mid", "tail"):
            if isinstance(raw.get(key), dict):
                try:
                    raw[key] = MonomerSpec(**raw[key])
                except TypeError as exc:
                    raise ConfigError(f"Invalid '{key}' block in {path}: {exc}") from exc

        cfg = cls(**raw)
        cfg.validate()
        return cfg

    def validate(self) -> None:
        """Raise :class:`ConfigError` if the configuration cannot drive a workflow."""
        if self.topology not in TOPOLOGIES:
            raise ConfigError(
                f"topology must be one of {TOPOLOGIES}, got {self.topology!r}"
            )
        for name in ("head", "mid", "tail"):
            spec = getattr(self, name)
            if spec is None:
                raise ConfigError(f"'{name}' monomer block is required")
        if self.head.termname is None:
            raise ConfigError("head.termname is required (grafting-end atom)")
        if self.tail.termname is None:
            raise ConfigError("tail.termname is required (free-end atom)")
        if self.head.tailname is None or self.head.post_tailtype is None:
            raise ConfigError("head.tailname and head.post_tailtype are required")
        if self.mid.headname is None or self.mid.tailname is None:
            raise ConfigError("mid.headname and mid.tailname are required")
        if self.mid.pre_headtype is None or self.mid.post_tailtype is None:
            raise ConfigError("mid.pre_headtype and mid.post_tailtype are required")
        if self.tail.headname is None or self.tail.pre_headtype is None:
            raise ConfigError("tail.headname and tail.pre_headtype are required")
        if self.rho <= 0 or self.nx < 1 or self.ny < 1:
            raise ConfigError("rho must be > 0 and nx, ny must be >= 1")
        if self.n_mid_repeat_units < 0:
            raise ConfigError("n_mid_repeat_units must be >= 0")
        if self.topology == "loop" and self.d_polymer is None:
            raise ConfigError(
                "loop topology requires an explicit d_polymer (HEAD–TAIL grafting "
                "point separation in Å, e.g. 14.9); the auto-calculated contour "
                "length is not meaningful for a loop"
            )
        if self.linker_atoms is not None:
            for i, spec in enumerate(self.linker_atoms):
                if not {"resname", "atomname"} <= set(spec):
                    raise ConfigError(
                        f"linker_atoms[{i}] must have 'resname' and 'atomname' keys"
                    )
        if self.linker_height < 0:
            raise ConfigError("linker_height must be >= 0 (Å above the wall at z = 0)")
        if self.solvent_min_z < 0:
            raise ConfigError("solvent_min_z must be >= 0 (Å above the wall at z = 0)")

    # ------------------------------------------------------------------
    # Derived geometry
    # ------------------------------------------------------------------

    def box_x(self) -> float:
        """Simulation box length along x in Ångström."""
        return float(np.sqrt(self.nx * self.ny / self.rho) * 10)

    def box_y(self) -> float:
        """Simulation box length along y in Ångström.

        Equal to :meth:`box_x`; the box cross-section is always square so
        that ``nx * ny / (box_x * box_y)`` equals ``rho``.  ``xyratio`` only
        shapes the PACKMOL alignment box, exactly as in the legacy scripts.
        """
        return self.box_x()

    def n_cc_all(self) -> int:
        """Total number of backbone C-C bonds in the full chain."""
        h = self.head.n_cc if self.head else 1
        m = self.mid.n_cc if self.mid else 1
        t = self.tail.n_cc if self.tail else 1
        return (h + 1) + self.n_mid_repeat_units * (m + 1) + (t + 1) + 2

    def polymer_length(self) -> float:
        """HEAD–TAIL restraint distance in Ångström (see :attr:`d_polymer`)."""
        if self.d_polymer is not None:
            return self.d_polymer
        return self.d_cc * self.n_cc_all() * 0.8

    def loop_height(self) -> float:
        """Arch height for loop topology in Ångström: half the extended contour length."""
        return self.d_cc * self.n_cc_all() * 0.8 / 2
