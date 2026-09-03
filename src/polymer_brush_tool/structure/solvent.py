"""Post-processing of ``gmx solvate`` output.

``gmx solvate`` knows nothing about the GROMACS wall at z = 0: it fills the
whole box, including the slab between the wall and the grafting atoms.  Water
there sits inside the substrate and in the repulsive core of the 10-4 wall
potential, so it is removed before the topologies are finalised.
"""

from __future__ import annotations

from pathlib import Path


def remove_water_below(
    gro_in: str | Path,
    gro_out: str | Path,
    z_min: float,
    resname: str = "SOL",
) -> tuple[int, int]:
    """Delete whole *resname* molecules whose first atom lies below ``z_min`` (Å).

    The first atom of each residue (``OW`` for TIP3P) decides.  All other
    residues are copied unchanged; atom serial numbers are renumbered so the
    file stays consecutive.  The title, atom count and box line are rewritten.

    Parameters
    ----------
    gro_in, gro_out:
        Input and output ``.gro`` paths (may be the same file).
    z_min:
        Height in Å.  Molecules with first-atom z < ``z_min`` are removed.
    resname:
        Residue name of the solvent (default ``"SOL"``).

    Returns
    -------
    (kept, removed)
        Number of *resname* molecules kept and removed.
    """
    lines = Path(gro_in).read_text().splitlines()
    title, natoms = lines[0], int(lines[1].strip())
    atom_lines = lines[2 : 2 + natoms]
    box_line = lines[2 + natoms]

    # group consecutive atom lines into residues by (resid, resname) columns
    residues: list[list[str]] = []
    last_key = None
    for line in atom_lines:
        key = line[0:10]
        if key != last_key:
            residues.append([])
            last_key = key
        residues[-1].append(line)

    kept_lines: list[str] = []
    kept = removed = 0
    for res in residues:
        if res[0][5:10].strip() == resname:
            z_nm = float(res[0][36:44])
            if z_nm * 10.0 < z_min:
                removed += 1
                continue
            kept += 1
        kept_lines.extend(res)

    out = [title, f"{len(kept_lines):5d}"]
    for i, line in enumerate(kept_lines, start=1):
        out.append(f"{line[:15]}{i % 100000:5d}{line[20:]}")
    out.append(box_line)
    Path(gro_out).write_text("\n".join(out) + "\n")
    return kept, removed
