"""PACKMOL-based chain alignment and graft placement on a 2-D grid.

Functions
---------
align_chain_z
    Iteratively call PACKMOL to align the chain bottom atom along z = 0.
graft_brush
    Write and execute a PACKMOL input that places chains on a regular grid.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from ase.io import read

from polymer_brush_tool import runner


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

    Iterates over increasing box x-sizes until PACKMOL succeeds in placing
    the chain with its bottom atom constrained to z ≈ 0.  The smallest
    successful box is used.

    Parameters
    ----------
    chain_pdb:
        Path to the minimised chain PDB (e.g. ``chain_min_pull.pdb``).
    bottom_atom_index:
        1-based index of the atom that should sit at the substrate (z = 0).
    xyratio:
        Aspect ratio box_y / box_x (default 1.0 = square cross-section).
    scan_start, scan_stop, scan_step:
        Range of x-box sizes (Å) to try.
    work_dir:
        Working directory for PACKMOL input/output files.

    Returns
    -------
    Path
        Path to the aligned PDB (``aligned_chain.pdb``).

    Raises
    ------
    RuntimeError
        If PACKMOL fails for every box size in the scan range.
    """
    work_dir = Path(work_dir)
    aligned_pdb = work_dir / "aligned_chain.pdb"

    for x in np.arange(scan_start, scan_stop, scan_step):
        in_path = work_dir / "align_z.in"
        out_path = work_dir / "align_z.out"

        with open(in_path, mode="w") as fh:
            fh.write("tolerance 2.0\n")
            fh.write(f"output {aligned_pdb}\n")
            fh.write("filetype pdb\n")
            fh.write("add_amber_ter\n")
            fh.write(f"structure {chain_pdb}\n")
            fh.write("  number 1\n")
            fh.write(f"  inside box 0. 0. 0. {x:.2f} {x * xyratio:.2f} 1000.\n")
            fh.write(f"  atoms {bottom_atom_index}\n")
            fh.write(f"  inside box 0. 0. 1.2 {x:.2f} {x * xyratio:.2f} 1.2\n")
            fh.write("  end atoms\n")
            fh.write("end structure\n")

        print(f"Trying z-alignment with box_x = {x:.2f} Å")
        runner.run(
            f"packmol < {in_path} > {out_path}",
            cwd=work_dir,
            check=False,  # We check the output file for errors instead
        )

        if out_path.exists() and "ERROR" not in out_path.read_text():
            print(f"  → Success at box_x = {x:.2f} Å")
            return aligned_pdb

    raise RuntimeError(
        f"PACKMOL z-alignment failed for all box sizes "
        f"in [{scan_start}, {scan_stop}) Å. "
        f"Inspect {work_dir / 'align_z.out'} for details."
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
    """Place ``nx * ny`` polymer chains on a regular 2-D grid using PACKMOL.

    Each chain is placed at the centre of its grid cell with its bottom
    fixed at z = 0 (``fixed x y 0 0 0 0``).  Each chain receives a unique
    chain ID character.

    Parameters
    ----------
    chain_pdb:
        Path to the z-aligned single-chain PDB.
    box_x, box_y:
        Simulation box dimensions in Ångström.
    nx, ny:
        Number of chains in x and y directions.
    work_dir:
        Working directory for PACKMOL input/output files.

    Returns
    -------
    Path
        Path to the multi-chain grafted PDB (``grafted_chain.pdb``).
    """
    work_dir = Path(work_dir)
    chain_pdb = Path(chain_pdb)

    atoms = read(str(chain_pdb))
    min_z = atoms.get_positions()[:, 2].min()
    max_z = atoms.get_positions()[:, 2].max()
    box_z = max_z - min_z + 30.0  # leave 30 Å headroom above chain

    output_pdb = work_dir / "grafted_chain.pdb"
    in_path = work_dir / "graft_brush.in"
    out_path = work_dir / "graft_brush.out"

    alphabet = (
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    )

    with open(in_path, mode="w") as fh:
        fh.write("tolerance 2.0\n")
        fh.write(f"output {output_pdb}\n")
        fh.write("filetype pdb\n")
        fh.write(f"pbc {box_x:.4f} {box_y:.4f} {box_z:.4f}\n")

        chain_id = 0
        for ix in range(nx):
            for iy in range(ny):
                x = (ix + 0.5) * (box_x / nx)
                y = (iy + 0.5) * (box_y / ny)
                fh.write(f"structure {chain_pdb}\n")
                fh.write("  resnumbers 1\n")
                fh.write("  number 1\n")
                fh.write(f"  chain {alphabet[chain_id]}\n")
                fh.write(f"  fixed {x:.4f} {y:.4f} 0.0 0.0 0.0 0.0\n")
                fh.write("end structure\n")
                fh.write("\n")
                chain_id += 1

    runner.run(
        f"packmol < {in_path} > {out_path}",
        cwd=work_dir,
    )

    return output_pdb
