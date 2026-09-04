"""AMBER sander energy minimisation with pull / loop constraints.

Functions
---------
amber_min_with_pull
    Write restraint + minimisation input files and run ``sander``.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from polymer_brush_tool import runner


def _rel(path: str | Path, work_dir: Path) -> str:
    """Path of *path* relative to *work_dir* (all tools run with cwd=work_dir)."""
    return os.path.relpath(Path(path).resolve(), work_dir.resolve())


def amber_min_with_pull(
    prmtop: str | Path,
    inpcrd: str | Path,
    file_prefix: str,
    head_idx: int,
    tail_idx: int,
    polymer_length: float,
    *,
    loop_height: Optional[float] = None,
    work_dir: Path = Path("."),
) -> Path:
    """Run AMBER sander gas-phase minimisation with distance restraints.

    Writes an NMR-style distance restraint file and a sander minimisation
    input, runs sander, then converts the result to PDB with ``ambpdb``.

    Parameters
    ----------
    prmtop, inpcrd:
        AMBER topology / coordinates produced by tleap.
    file_prefix:
        Prefix for intermediate files (``<prefix>_min_pull.{in,out,rst7,pdb}``).
    head_idx, tail_idx:
        1-based atom indices of the HEAD and TAIL terminal atoms.
    polymer_length:
        Target HEAD–TAIL distance r3 in Å.  The flat-bottom well is
        r2 = 0.8·r3 … r4 = 1.2·r3 with rk2 = rk3 = 1000.7 kcal/mol/Å².
    loop_height:
        When given, two extra restraints (HEAD–MID and MID–TAIL, where MID
        is the atom with index ``(head_idx + tail_idx) // 2``) pull the
        chain centre away from the grafting points to form an arch.
    work_dir:
        Working directory; all file names inside the inputs are written
        relative to it because sander is executed with ``cwd=work_dir``.

    Returns
    -------
    Path
        Path to ``<prefix>_min_pull.pdb``.
    """
    work_dir = Path(work_dir).resolve()

    def _rst_block(fh, i, j, r3):
        r2, r4 = 0.8 * r3, 1.2 * r3
        fh.write("&rst\n")
        fh.write(
            f"  iat={i},{j}, r1=0., r2={r2:.4f}, r3={r3:.4f}, r4={r4:.4f},"
            f" rk2=1000.7, rk3=1000.7,\n"
        )
        fh.write("/\n")

    # -- Distance restraint file ------------------------------------------
    restraint_name = "pull_termination.restraint"
    with open(work_dir / restraint_name, mode="w") as fh:
        _rst_block(fh, head_idx, tail_idx, polymer_length)
        if loop_height is not None:
            mid_idx = (head_idx + tail_idx) // 2
            _rst_block(fh, mid_idx, tail_idx, loop_height)
            _rst_block(fh, head_idx, mid_idx, loop_height)

    # -- sander minimisation input ----------------------------------------
    min_in = f"{file_prefix}_min_pull.in"
    with open(work_dir / min_in, mode="w") as fh:
        fh.write("Minimize\n")
        fh.write("&cntrl\n")
        fh.write("  imin=1, ntb=0, ntx=1, irest=0,\n")
        fh.write("  maxcyc=50000, ncyc=1000,\n")
        fh.write("  ntpr=100, ntwx=0, cut=999.0,\n")
        fh.write("  nmropt=1,\n")
        fh.write("/\n")
        fh.write("&wt type='END' /\n")
        fh.write(f"DISANG={restraint_name}\n")

    prmtop_rel = _rel(prmtop, work_dir)
    inpcrd_rel = _rel(inpcrd, work_dir)
    rst7 = f"{file_prefix}_min_pull.rst7"
    pdb = f"{file_prefix}_min_pull.pdb"

    runner.run(
        f"sander -O -i {min_in} -o {file_prefix}_min_pull.out"
        f" -p {prmtop_rel} -c {inpcrd_rel} -r {rst7}",
        cwd=work_dir,
    )
    runner.run(f"ambpdb -p {prmtop_rel} -c {rst7} > {pdb}", cwd=work_dir)

    return work_dir / pdb
