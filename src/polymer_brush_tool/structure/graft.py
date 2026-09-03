"""PACKMOL-based chain alignment and graft placement on a 2-D grid.

Functions
---------
align_chain_z
    Iteratively call PACKMOL to align the chain bottom atom along z = 0.
graft_brush
    Write and execute a PACKMOL input that places chains on a regular grid.

All PACKMOL inputs reference files by paths *relative to work_dir* and
PACKMOL is executed with ``cwd=work_dir``; PACKMOL has a fixed-width
string buffer, so long absolute paths are avoided on purpose.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from ase.io import read

from polymer_brush_tool import runner

_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


def _rel(path: str | Path, work_dir: Path) -> str:
    return os.path.relpath(Path(path).resolve(), work_dir.resolve())


def align_chain_z(
    chain_pdb: str | Path,
    bottom_atom_index: int,
    xyratio: float = 1.0,
    *,
    scan_start: float = 2.1,
    scan_stop: float = 100.1,
    scan_step: float = 0.1,
    work_dir: Path = Path("."),
) -> Path:
    """Align the polymer chain along the z-axis using PACKMOL.

    The chain is constrained to a box ``[0,x] × [0,x·xyratio] × [0,1000]``
    with *bottom_atom_index* pinned to the plane z = 1.2 Å.  Starting from a
    very narrow box, x is increased until PACKMOL succeeds; the first
    success therefore yields the most upright orientation of the chain.

    Parameters
    ----------
    chain_pdb:
        Minimised chain PDB (e.g. ``chain_min_pull.pdb``).
    bottom_atom_index:
        1-based index of the atom that should sit at the substrate.
    xyratio:
        Aspect ratio box_y / box_x.
    scan_start, scan_stop, scan_step:
        Range of x-box sizes (Å) to try.
    work_dir:
        Working directory for PACKMOL input/output files.

    Returns
    -------
    Path
        ``<work_dir>/aligned_chain.pdb``.

    Raises
    ------
    RuntimeError
        If PACKMOL fails for every box size in the scan range.
    """
    work_dir = Path(work_dir).resolve()
    chain_rel = _rel(chain_pdb, work_dir)
    in_path = work_dir / "align_z.in"
    out_path = work_dir / "align_z.out"

    for x in np.arange(scan_start, scan_stop, scan_step):
        with open(in_path, mode="w") as fh:
            fh.write("tolerance 2.0\n")
            fh.write("output aligned_chain.pdb\n")
            fh.write("filetype pdb\n")
            fh.write("add_amber_ter\n")
            fh.write(f"structure {chain_rel}\n")
            fh.write("  number 1\n")
            fh.write(f"  inside box 0. 0. 0. {x:.2f} {x * xyratio:.2f} 1000.\n")
            fh.write(f"  atoms {bottom_atom_index}\n")
            fh.write(f"  inside box 0. 0. 1.2 {x:.2f} {x * xyratio:.2f} 1.2\n")
            fh.write("  end atoms\n")
            fh.write("end structure\n")

        print(f"Trying z-alignment with box_x = {x:.2f} Å", flush=True)
        # PACKMOL's exit status is not reliable across versions; the
        # original script checked the log for "ERROR" and so do we.
        runner.run("packmol < align_z.in > align_z.out", cwd=work_dir, check=False, verbose=False)

        if out_path.exists() and "ERROR" not in out_path.read_text():
            print(f"  -> success at box_x = {x:.2f} Å")
            return work_dir / "aligned_chain.pdb"

    raise RuntimeError(
        f"PACKMOL z-alignment failed for all box sizes in "
        f"[{scan_start}, {scan_stop}) Å. See {out_path}."
    )


def graft_brush(
    chain_pdb: str | Path,
    box_x: float,
    box_y: float,
    nx: int,
    ny: int,
    *,
    work_dir: Path = Path("."),
) -> Path:
    """Place ``nx * ny`` copies of a chain on a regular 2-D grid using PACKMOL.

    Each copy is placed at the centre of its grid cell with no rotation and
    no translation in z (``fixed x y 0 0 0 0``), and given a unique chain ID.
    The box height is the chain extent plus 30 Å of headroom.

    Parameters
    ----------
    chain_pdb:
        z-aligned single-chain PDB.
    box_x, box_y:
        Lateral box dimensions in Å.
    nx, ny:
        Number of chains along x and y.
    work_dir:
        Working directory for PACKMOL input/output files.

    Returns
    -------
    Path
        ``<work_dir>/grafted_chain.pdb``.
    """
    if nx * ny > len(_CHAIN_IDS):
        raise ValueError(
            f"nx*ny = {nx*ny} exceeds the {len(_CHAIN_IDS)} available PDB chain IDs"
        )

    work_dir = Path(work_dir).resolve()
    chain_rel = _rel(chain_pdb, work_dir)

    z = read(str(chain_pdb)).get_positions()[:, 2]
    box_z = z.max() - z.min() + 30.0

    with open(work_dir / "graft_brush.in", mode="w") as fh:
        fh.write("tolerance 2.0\n")
        fh.write("output grafted_chain.pdb\n")
        fh.write("filetype pdb\n")
        fh.write(f"pbc {box_x:.4f} {box_y:.4f} {box_z:.4f}\n")
        chain_id = 0
        for ix in range(nx):
            for iy in range(ny):
                x = (ix + 0.5) * (box_x / nx)
                y = (iy + 0.5) * (box_y / ny)
                fh.write(f"structure {chain_rel}\n")
                fh.write("  resnumbers 1\n")
                fh.write("  number 1\n")
                fh.write(f"  chain {_CHAIN_IDS[chain_id]}\n")
                fh.write(f"  fixed {x:.4f} {y:.4f} 0.0 0.0 0.0 0.0\n")
                fh.write("end structure\n\n")
                chain_id += 1

    runner.run("packmol < graft_brush.in > graft_brush.out", cwd=work_dir)
    return work_dir / "grafted_chain.pdb"
