"""Write AMBER fragment chain files and tleap scripts.

Functions
---------
write_fragment
    Write a ``prepgen`` chain connectivity file for one monomer.
write_tleap
    Write a tleap script to polymerise HEAD + MID*n + TAIL.
write_tleap_grafted
    Write a tleap script to assign force-fields to a grafted multi-chain PDB.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional


def write_fragment(
    resname: str,
    *,
    headname: Optional[str] = None,
    tailname: Optional[str] = None,
    omitnames: Optional[list[str]] = None,
    pre_headtype: Optional[str] = None,
    post_tailtype: Optional[str] = None,
    work_dir: Path = Path("."),
) -> Path:
    """Write a ``prepgen`` chain connectivity file for one monomer.

    Parameters
    ----------
    resname:
        3-letter residue name (lower-case). The output file is
        ``<work_dir>/<resname>.chain``.
    headname:
        Atom name of the HEAD-side bonding atom.  Requires *pre_headtype*.
        Omit for the HEAD monomer (which has no upstream residue).
    tailname:
        Atom name of the TAIL-side bonding atom.  Requires *post_tailtype*.
        Omit for the TAIL monomer.
    omitnames:
        Atom names of hydrogens to remove at bonding points.  Must contain
        exactly one entry per specified head/tail name.
    pre_headtype:
        GAFF atom type of the residue that bonds at the head side.
    post_tailtype:
        GAFF atom type of the residue that bonds at the tail side.
    work_dir:
        Directory where the ``.chain`` file is written (default: cwd).

    Returns
    -------
    Path
        Path to the written ``.chain`` file.

    Raises
    ------
    ValueError
        If *pre_headtype* / *post_tailtype* is missing when the corresponding
        head/tail name is provided, or if the length of *omitnames* does not
        match the number of bonding endpoints.
    TypeError
        If *omitnames* is not a string or list.
    """
    n_omits = 0
    out_path = work_dir / f"{resname}.chain"

    with open(out_path, mode="w") as fh:
        if headname is not None:
            n_omits += 1
            if pre_headtype is None:
                raise ValueError(
                    "pre_headtype must be provided when headname is given."
                )
            fh.write(f"HEAD_NAME  {headname}\n")
            fh.write(f"PRE_HEAD_TYPE  {pre_headtype}\n")

        if tailname is not None:
            n_omits += 1
            if post_tailtype is None:
                raise ValueError(
                    "post_tailtype must be provided when tailname is given."
                )
            fh.write(f"TAIL_NAME  {tailname}\n")
            fh.write(f"POST_TAIL_TYPE  {post_tailtype}\n")

        if omitnames is not None:
            if isinstance(omitnames, str):
                omitnames = [omitnames]
            elif not isinstance(omitnames, list):
                raise TypeError("omitnames must be a string or a list of strings.")
            if len(omitnames) != n_omits:
                raise ValueError(
                    f"Length of omitnames ({len(omitnames)}) does not match "
                    f"the number of head/tail names provided ({n_omits})."
                )
            for name in omitnames:
                fh.write(f"OMIT_NAME  {name}\n")

        fh.write("CHARGE 0.000\n")

    return out_path


def write_tleap(
    head_resname: str,
    mid_resname: str,
    tail_resname: str,
    n_mid_repeat_units: int,
    *,
    work_dir: Path = Path("."),
) -> Path:
    """Write a tleap script to polymerise HEAD + MID*n + TAIL.

    The generated script sources GAFF2, loads the three ``.prepi`` files,
    and assembles the polymer sequence.  Output topology and coordinates are
    saved as ``chain.prmtop`` / ``chain.inpcrd``.

    Parameters
    ----------
    head_resname, mid_resname, tail_resname:
        3-letter residue names (lower-case).
    n_mid_repeat_units:
        Number of MID monomer repeat units.
    work_dir:
        Directory where ``build_chain.tleap`` is written.

    Returns
    -------
    Path
        Path to the written tleap script.
    """
    out_path = work_dir / "build_chain.tleap"
    with open(out_path, mode="w") as fh:
        fh.write("source leaprc.gaff2\n")
        fh.write(f"loadamberprep {head_resname}.prepi\n")
        fh.write(f"loadamberprep {mid_resname}.prepi\n")
        fh.write(f"loadamberprep {tail_resname}.prepi\n")
        fh.write(f"chain = sequence {{{head_resname.upper()} ")
        for _ in range(n_mid_repeat_units):
            fh.write(f"{mid_resname.upper()} ")
        fh.write(f" {tail_resname.upper()}}}\n")
        fh.write("saveAmberParm chain chain.prmtop chain.inpcrd\n")
        fh.write("savepdb chain chain.pdb\n")
        fh.write("quit\n")
    return out_path


def write_tleap_grafted(
    head_resname: str,
    mid_resname: str,
    tail_resname: str,
    box_x: float,
    box_y: float,
    box_z: float,
    *,
    pdb_file: str = "grafted_chain.pdb",
    work_dir: Path = Path("."),
) -> Path:
    """Write a tleap script to assign force-fields to a grafted multi-chain PDB.

    After PACKMOL places all chains, this script re-assigns GAFF2 parameters
    to the combined system and exports AMBER topology + GROMACS files.

    Parameters
    ----------
    head_resname, mid_resname, tail_resname:
        3-letter residue names (lower-case).
    box_x, box_y, box_z:
        Periodic box dimensions in Ångström.
    pdb_file:
        Filename of the multi-chain PDB produced by PACKMOL.
    work_dir:
        Directory where ``grafted_chain.tleap`` is written.

    Returns
    -------
    Path
        Path to the written tleap script.
    """
    out_path = work_dir / "grafted_chain.tleap"
    with open(out_path, mode="w") as fh:
        fh.write("source leaprc.gaff2\n")
        fh.write(f"loadamberprep {head_resname}.prepi\n")
        fh.write(f"loadamberprep {mid_resname}.prepi\n")
        fh.write(f"loadamberprep {tail_resname}.prepi\n")
        fh.write(f"mol = loadpdb {pdb_file}\n")
        fh.write(f"set mol box {{{box_x} {box_y} {box_z}}}\n")
        fh.write("saveAmberParm mol grafted_chain.prmtop grafted_chain.inpcrd\n")
        fh.write("quit\n")
    return out_path
