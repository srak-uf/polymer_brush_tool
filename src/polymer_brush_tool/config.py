"""Configuration dataclasses and YAML loading for polymer brush workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml


@dataclass
class MonomerSpec:
    """Specification for one monomer residue (HEAD, MID, or TAIL)."""

    resname: str
    """3-letter residue name used by tleap (e.g. 'hmp', 'mmp', 'tmp')."""

    ac_file: str
    """Path to the AMBER AC file produced by antechamber."""

    omitnames: list[str]
    """Atom names of hydrogens to remove at bonding points (one per bond endpoint)."""

    n_cc: int = 1
    """Number of backbone C-C bonds in this monomer unit."""

    headname: Optional[str] = None
    """Atom name of the HEAD-side bonding carbon (None for HEAD monomer)."""

    tailname: Optional[str] = None
    """Atom name of the TAIL-side bonding carbon (None for TAIL monomer)."""

    pre_headtype: Optional[str] = None
    """GAFF atom type of the residue that bonds to this monomer's head (None for HEAD)."""

    post_tailtype: Optional[str] = None
    """GAFF atom type of the residue that bonds to this monomer's tail (None for TAIL)."""

    termname: Optional[str] = None
    """Atom name of the terminal atom (used to locate the end of the chain; HEAD and TAIL)."""


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
    """Aspect ratio box_y / box_x (default 1.0 = square)."""

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
    Extended chain length in Ångström.  When None it is auto-calculated
    from the C-C bond count: ``d_cc * n_cc_all * 0.8``.
    """

    d_cc: float = 1.54
    """C-C backbone bond length in Ångström (default 1.54 Å)."""

    # ------------------------------------------------------------------
    # Interactive prompt overrides
    # ------------------------------------------------------------------
    bottom_atom_index: Optional[int] = None
    """
    1-based index of the grafting (bottom) atom in the minimised chain PDB.
    When None, the script will prompt interactively (legacy behaviour).
    """

    linker_atoms: Optional[list[dict]] = None
    """
    Explicit list of linker atoms at the substrate attachment points.
    Each entry is ``{"resname": "HMP", "atomname": "H1"}``.
    When None, the script will use the HEAD terminal atom as the default
    or prompt interactively.
    """

    # ------------------------------------------------------------------
    # GROMACS MPI/OMP parallelism
    # ------------------------------------------------------------------
    t_mpi: int = 8
    """Number of thread-MPI ranks for gmx mdrun."""

    t_omp: int = 1
    """Number of OpenMP threads per rank for gmx mdrun."""

    # ------------------------------------------------------------------
    # Class methods
    # ------------------------------------------------------------------

    @classmethod
    def from_yaml(cls, path: str | Path) -> "BrushConfig":
        """Load a BrushConfig from a YAML file.

        Parameters
        ----------
        path:
            Path to the YAML configuration file.

        Returns
        -------
        BrushConfig
        """
        with open(path) as fh:
            raw = yaml.safe_load(fh)

        # Convert nested dicts for monomer specs
        for key in ("head", "mid", "tail"):
            if isinstance(raw.get(key), dict):
                raw[key] = MonomerSpec(**raw[key])

        return cls(**raw)

    def box_x(self) -> float:
        """Simulation box length along x in Ångström."""
        import numpy as np
        return float(np.sqrt(self.nx * self.ny / self.rho) * 10)

    def box_y(self) -> float:
        """Simulation box length along y in Ångström."""
        return self.box_x() * self.xyratio

    def n_cc_all(self) -> int:
        """Total number of backbone C-C bonds in the full chain."""
        h = self.head.n_cc if self.head else 1
        m = self.mid.n_cc if self.mid else 1
        t = self.tail.n_cc if self.tail else 1
        return (h + 1) + self.n_mid_repeat_units * (m + 1) + (t + 1) + 2

    def polymer_length(self) -> float:
        """Extended chain contour length in Ångström."""
        if self.d_polymer is not None:
            return self.d_polymer
        return self.d_cc * self.n_cc_all() * 0.8

    def loop_height(self) -> float:
        """Half-chain height for loop topology in Ångström."""
        return self.d_cc * self.n_cc_all() * 0.8 / 2
