"""Patch GROMACS topology files.

Functions
---------
insert_tip3p_top
    Inject ``#include "tip3p.itp"`` before each ``[ moleculetype ]`` section.
insert_restraint_top
    Append ``[ position_restraints ]`` for linker atoms at two force constants.
"""

from __future__ import annotations

from pathlib import Path


def insert_tip3p_top(
    gmx_top: str | Path,
    output_top: str | Path,
) -> None:
    """Inject a TIP3P water model include before each ``[ moleculetype ]`` section.

    GROMACS topologies produced by ParmEd do not include the TIP3P water
    parameters.  This function patches the topology by inserting
    ``#include "tip3p.itp"`` just before every ``[ moleculetype ]`` line.

    Parameters
    ----------
    gmx_top:
        Path to the input GROMACS ``.top`` file.
    output_top:
        Path for the patched output file.
    """
    lines = Path(gmx_top).read_text().splitlines(keepends=True)

    with open(output_top, "w") as fh:
        for i, line in enumerate(lines):
            fh.write(line)
            if i < len(lines) - 1 and lines[i + 1].strip() == "[ moleculetype ]":
                fh.write('#include "tip3p.itp"\n')
                fh.write("\n")


def insert_restraint_top(
    gmx_top: str | Path,
    output_top: str | Path,
    linker_indices: list[int],
) -> None:
    """Append position restraints for linker atoms to a GROMACS topology.

    Two output files are written:

    * ``<output_top>`` — soft restraints at 10,000 kJ/mol/nm² (used for
      initial solvated equilibration).
    * ``hardrest_<output_top>`` — hard restraints at 1,000,000 kJ/mol/nm²
      (used for vacuum relaxation to prevent chain collapse).

    The ``[ position_restraints ]`` section is inserted just before
    ``[ system ]``.

    Parameters
    ----------
    gmx_top:
        Path to the input GROMACS ``.top`` file.
    output_top:
        Path for the soft-restraint output file.  The hard-restraint file
        is written in the same directory with the ``hardrest_`` prefix.
    linker_indices:
        1-based atom indices of the linker atoms to restrain.
    """
    output_top = Path(output_top)
    hardrest_top = output_top.parent / ("hardrest_" + output_top.name)

    lines = Path(gmx_top).read_text().splitlines(keepends=True)

    # Collect valid atom indices from the [ atoms ] section
    valid_indices: set[int] = set()
    in_atoms = False
    for line in lines:
        if line.strip() == "[ atoms ]":
            in_atoms = True
            continue
        if in_atoms:
            stripped = line.strip()
            if stripped == "" or stripped.startswith("["):
                break
            if stripped.startswith(";"):
                continue
            parts = stripped.split()
            if parts:
                valid_indices.add(int(parts[0]))

    def _write(path: Path, force_constant: int) -> None:
        with open(path, "w") as fh:
            for i, line in enumerate(lines):
                fh.write(line)
                if i < len(lines) - 1 and lines[i + 1].strip() == "[ system ]":
                    fh.write("\n[ position_restraints ]\n")
                    fh.write("; atom  type      fx      fy      fz\n")
                    for idx in linker_indices:
                        if idx in valid_indices:
                            fh.write(f" {idx}  1  {force_constant}  {force_constant}  {force_constant}\n")
                    fh.write("\n")

    _write(output_top, 10_000)
    _write(hardrest_top, 1_000_000)


def copy_molecules_section(src_top: str | Path, dst_top: str | Path) -> None:
    """Replace the ``[ molecules ]`` section of *dst_top* with the one from *src_top*.

    ``gmx solvate -p`` appends the water count (``SOL  N``) to the topology it
    is given.  The soft- and hard-restraint topologies are otherwise
    identical, so after solvating with one of them the other must receive the
    same ``[ molecules ]`` block or ``grompp`` will reject it for an atom-count
    mismatch.

    Parameters
    ----------
    src_top:
        Topology whose ``[ molecules ]`` section is authoritative.
    dst_top:
        Topology to update in place.
    """
    def _split(lines: list[str]) -> tuple[list[str], list[str]]:
        for i, line in enumerate(lines):
            if line.strip() == "[ molecules ]":
                return lines[:i], lines[i:]
        raise ValueError("no [ molecules ] section found")

    src_lines = Path(src_top).read_text().splitlines(keepends=True)
    dst_lines = Path(dst_top).read_text().splitlines(keepends=True)
    _, src_mol = _split(src_lines)
    dst_head, _ = _split(dst_lines)
    Path(dst_top).write_text("".join(dst_head + src_mol))


def set_molecule_count(top: str | Path, name: str, count: int) -> None:
    """Set the count of molecule *name* in the ``[ molecules ]`` section of *top*.

    Used after water molecules are deleted from the ``gmx solvate`` output so
    that ``SOL  N`` matches the coordinate file again.

    Raises
    ------
    ValueError
        If the section or the molecule line is missing.
    """
    path = Path(top)
    lines = path.read_text().splitlines(keepends=True)
    in_mol = False
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("["):
            in_mol = s == "[ molecules ]"
            continue
        if in_mol and s and not s.startswith(";") and s.split()[0] == name:
            lines[i] = f"{name:<16s} {count}\n"
            path.write_text("".join(lines))
            return
    raise ValueError(f"no '{name}' line in [ molecules ] of {path}")


def remove_molecule(top: str | Path, name: str) -> int:
    """Delete every ``name`` line from the ``[ molecules ]`` section of *top*.

    ``gmx solvate -p`` *appends* ``SOL  N``; when ``step_solvate`` is re-run
    on a topology that already carries a ``SOL`` line the counts would add
    up.  Calling this first makes the step idempotent.

    Returns
    -------
    int
        Number of lines removed (0 if none, no error).
    """
    path = Path(top)
    lines = path.read_text().splitlines(keepends=True)
    out, in_mol, removed = [], False, 0
    for line in lines:
        s = line.strip()
        if s.startswith("["):
            in_mol = s == "[ molecules ]"
        elif in_mol and s and not s.startswith(";") and s.split()[0] == name:
            removed += 1
            continue
        out.append(line)
    if removed:
        path.write_text("".join(out))
    return removed
