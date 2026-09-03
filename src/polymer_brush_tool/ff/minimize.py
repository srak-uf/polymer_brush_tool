"""AMBER sander energy minimisation with pull / loop constraints.

Functions
---------
amber_min_with_pull
    Write restraint + minimisation input files and run ``sander``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from polymer_brush_tool import runner


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
) -> None:
    """Run AMBER sander gas-phase minimisation with end-to-end distance restraints.

    Writes a NMR-style distance restraint file and a sander minimisation
    input file, then executes sander and converts the output to PDB with
    ``ambpdb``.

    Parameters
    ----------
    prmtop:
        AMBER topology file produced by tleap.
    inpcrd:
        AMBER coordinate file produced by tleap.
    file_prefix:
        Prefix used for all intermediate file names (e.g. ``"chain"``
        produces ``chain_min_pull.in``, ``chain_min_pull.rst7``, etc.).
    head_idx:
        1-based atom index of the HEAD terminal atom.
    tail_idx:
        1-based atom index of the TAIL terminal atom.
    polymer_length:
        Target end-to-end distance in Ångström for the linear stretch
        restraint (r3).  r2 = 0.8 * r3, r4 = 1.2 * r3.
    loop_height:
        When provided, additional mid-point restraints are added to model
        a loop topology.  The value sets the target height in Ångström.
    work_dir:
        Working directory for all files.
    """
    work_dir = Path(work_dir)

    # -- Distance restraint file ------------------------------------------
    restraint_path = work_dir / "pull_termination.restraint"
    with open(restraint_path, mode="w") as fh:
        r3 = polymer_length
        r2 = 0.8 * r3
        r4 = 1.2 * r3
        fh.write("&rst\n")
        fh.write(
            f"iat={head_idx},{tail_idx},"
            f"r1=0., r2={r2:.4f}, r3={r3:.4f}, r4={r4:.4f},"
            f"rk2=1000.7, rk3=1000.7,\n"
        )
        fh.write("/\n")

        if loop_height is not None:
            r3_loop = loop_height
            r2_loop = 0.8 * r3_loop
            r4_loop = 1.2 * r3_loop
            mid_idx = (head_idx + tail_idx) // 2

            fh.write("&rst\n")
            fh.write(
                f"iat={mid_idx},{tail_idx},"
                f"r1=0., r2={r2_loop:.4f}, r3={r3_loop:.4f}, r4={r4_loop:.4f},"
                f"rk2=1000.7, rk3=1000.7,\n"
            )
            fh.write("/\n")

            fh.write("&rst\n")
            fh.write(
                f"iat={head_idx},{mid_idx},"
                f"r1=0., r2={r2_loop:.4f}, r3={r3_loop:.4f}, r4={r4_loop:.4f},"
                f"rk2=1000.7, rk3=1000.7,\n"
            )
            fh.write("/\n")

    # -- sander minimisation input file ------------------------------------
    min_in_path = work_dir / f"{file_prefix}_min_pull.in"
    with open(min_in_path, mode="w") as fh:
        fh.write("Minimize\n")
        fh.write("&cntrl\n")
        fh.write("imin=1,\n")
        fh.write("ntb=0,\n")
        fh.write("ntx=1,\n")
        fh.write("irest=0,\n")
        fh.write("maxcyc=50000,\n")
        fh.write("ncyc=1000,\n")
        fh.write("ntpr=100,\n")
        fh.write("ntwx=0,\n")
        fh.write("cut=999.0,\n")
        fh.write("nmropt=1,\n")
        fh.write("/\n")
        fh.write("&wt type='END' /\n")
        fh.write(f"DISANG=pull_termination.restraint\n")

    # -- Run sander -------------------------------------------------------
    rst7_path = work_dir / f"{file_prefix}_min_pull.rst7"
    runner.run(
        f"sander -O"
        f" -i {min_in_path}"
        f" -o {work_dir / (file_prefix + '_min_pull.out')}"
        f" -p {prmtop}"
        f" -c {inpcrd}"
        f" -r {rst7_path}",
        cwd=work_dir,
    )

    # -- Convert rst7 → PDB -----------------------------------------------
    pdb_path = work_dir / f"{file_prefix}_min_pull.pdb"
    runner.run(
        f"ambpdb -p {prmtop} -c {rst7_path} > {pdb_path}",
        cwd=work_dir,
    )
