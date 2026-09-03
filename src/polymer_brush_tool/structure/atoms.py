"""Atom-lookup helpers for ASE Atoms objects loaded from polymer PDB/GRO files."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms


def find_atom_index(
    atoms: "Atoms",
    atomtype: str,
    resname: str,
) -> int:
    """Return the 0-based index of the first atom matching *atomtype* and *resname*.

    Parameters
    ----------
    atoms:
        ASE Atoms object loaded from a PDB or GRO file that carries
        ``"atomtypes"`` and ``"residuenames"`` in ``atoms.arrays``.
    atomtype:
        Atom name (e.g. ``"H1"``, ``"C2"``).
    resname:
        Residue name in UPPER CASE as stored in the PDB (e.g. ``"HMP"``).

    Returns
    -------
    int
        0-based atom index.

    Raises
    ------
    ValueError
        If no matching atom is found.
    """
    atypes = atoms.arrays["atomtypes"]
    rnames = atoms.arrays["residuenames"]
    for i, (at, rn) in enumerate(zip(atypes, rnames)):
        if at == atomtype and rn == resname:
            return i
    raise ValueError(
        f"No atom with atomtype={atomtype!r} and resname={resname!r} found."
    )


def find_linker_indices(
    atoms: "Atoms",
    linker_atoms: list[dict],
) -> list[int]:
    """Return 0-based indices for a list of ``{"resname": ..., "atomname": ...}`` specs.

    Parameters
    ----------
    atoms:
        ASE Atoms object.
    linker_atoms:
        List of dicts with keys ``"resname"`` (UPPER CASE) and
        ``"atomname"``.

    Returns
    -------
    list[int]
        0-based atom indices, one per linker spec (may include duplicates if
        the same atom matches multiple chains).
    """
    atypes = atoms.arrays["atomtypes"]
    rnames = atoms.arrays["residuenames"]
    indices: list[int] = []
    for spec in linker_atoms:
        aname = spec["atomname"]
        rname = spec["resname"]
        matched = [
            i
            for i, (at, rn) in enumerate(zip(atypes, rnames))
            if at == aname and rn == rname
        ]
        indices.extend(matched)
    return indices


def pin_linkers_to_substrate(
    atoms: "Atoms",
    linker_indices: list[int],
    height: float = 1.5,
) -> "Atoms":
    """Translate the system so the linker atoms sit at ``z = height`` (Å).

    Every linker atom is first moved onto the lowest z of the whole system,
    then the system is shifted rigidly so that plane lands at ``z = height``.
    Nothing else about the geometry changes.  The GROMACS wall is at z = 0,
    so *height* is the distance between the grafting atoms and the wall.

    Parameters
    ----------
    atoms:
        ASE Atoms object; modified in place and returned.
    linker_indices:
        0-based indices of the grafting atoms (one or more per chain).
    height:
        Target z (Å) of the linker atoms.  Must be >= 0.

    Returns
    -------
    Atoms
        The same object, for chaining.
    """
    if not linker_indices:
        raise ValueError("linker_indices must not be empty")
    if height < 0:
        raise ValueError("height must be >= 0")
    pos = atoms.get_positions()
    min_z = pos[:, 2].min()
    pos[linker_indices, 2] = min_z
    pos[:, 2] += height - min_z
    atoms.set_positions(pos)
    return atoms
