"""Tests for structure generation helpers (no external tools required)."""

import pytest
from unittest.mock import MagicMock, patch
from pathlib import Path

import numpy as np

from polymer_brush_tool.structure.atoms import find_atom_index, find_linker_indices


# --------------------------------------------------------------------------
# find_atom_index
# --------------------------------------------------------------------------

def _make_atoms(atomtypes, residuenames):
    """Create a minimal mock ASE Atoms object."""
    atoms = MagicMock()
    atoms.arrays = {
        "atomtypes": atomtypes,
        "residuenames": residuenames,
    }
    return atoms


class TestFindAtomIndex:
    def test_found(self):
        atoms = _make_atoms(
            ["C1", "H1", "C2", "H1"],
            ["HMP", "HMP", "MMP", "MMP"],
        )
        assert find_atom_index(atoms, "H1", "HMP") == 1

    def test_returns_first_match(self):
        atoms = _make_atoms(
            ["H1", "H1"],
            ["HMP", "HMP"],
        )
        assert find_atom_index(atoms, "H1", "HMP") == 0

    def test_not_found_raises(self):
        atoms = _make_atoms(["C1"], ["HMP"])
        with pytest.raises(ValueError, match="No atom"):
            find_atom_index(atoms, "H99", "HMP")

    def test_resname_case_sensitive(self):
        atoms = _make_atoms(["H1"], ["hmp"])  # lower-case resname
        with pytest.raises(ValueError):
            find_atom_index(atoms, "H1", "HMP")  # upper-case expected


# --------------------------------------------------------------------------
# find_linker_indices
# --------------------------------------------------------------------------

class TestFindLinkerIndices:
    def test_single_linker(self):
        atoms = _make_atoms(
            ["C1", "H1", "C2"],
            ["HMP", "HMP", "MMP"],
        )
        indices = find_linker_indices(
            atoms, [{"resname": "HMP", "atomname": "H1"}]
        )
        assert indices == [1]

    def test_multiple_chains_same_residue(self):
        """Multiple chains in the grafted PDB each have an HMP/H1 atom."""
        atoms = _make_atoms(
            ["H1", "C1", "H1"],
            ["HMP", "MMP", "HMP"],
        )
        indices = find_linker_indices(
            atoms, [{"resname": "HMP", "atomname": "H1"}]
        )
        assert indices == [0, 2]

    def test_loop_topology_two_linkers(self):
        atoms = _make_atoms(
            ["H1", "C1", "H24"],
            ["HMP", "MMP", "TMP"],
        )
        indices = find_linker_indices(
            atoms,
            [
                {"resname": "HMP", "atomname": "H1"},
                {"resname": "TMP", "atomname": "H24"},
            ],
        )
        assert set(indices) == {0, 2}

    def test_no_match_returns_empty(self):
        atoms = _make_atoms(["C1"], ["HMP"])
        indices = find_linker_indices(
            atoms, [{"resname": "TMP", "atomname": "H99"}]
        )
        assert indices == []


# --------------------------------------------------------------------------
# graft_brush (PACKMOL file content, no actual PACKMOL execution)
# --------------------------------------------------------------------------

class TestGraftBrushInput:
    """Test that graft_brush writes the correct PACKMOL input file."""

    def test_packmol_input_written(self, tmp_path):
        # Create a minimal fake chain PDB
        chain_pdb = tmp_path / "aligned_chain.pdb"
        chain_pdb.write_text("ATOM      1  C1  HMP A   1       0.000   0.000   0.000\n")

        # Mock ASE read and PACKMOL execution
        mock_atoms = MagicMock()
        mock_atoms.get_positions.return_value = np.array([[0.0, 0.0, 0.0]])
        mock_atoms.cell = np.diag([30.0, 30.0, 50.0])

        with (
            patch("polymer_brush_tool.structure.graft.read", return_value=mock_atoms),
            patch("polymer_brush_tool.structure.graft.runner.run"),
        ):
            from polymer_brush_tool.structure.graft import graft_brush
            graft_brush(chain_pdb, 30.0, 30.0, 2, 2, work_dir=tmp_path)

        packmol_in = tmp_path / "graft_brush.in"
        assert packmol_in.exists()
        text = packmol_in.read_text()

        # Check for 4 structure blocks (2×2 grid); "end structure" closes each
        assert text.count("end structure") == 4
        assert "chain A" in text
        assert "chain B" in text
        assert "pbc 30.0000 30.0000" in text

    def test_chain_id_sequence(self, tmp_path):
        """Each chain should get a unique alphabetical chain ID."""
        chain_pdb = tmp_path / "aligned_chain.pdb"
        chain_pdb.write_text("ATOM      1  C1  HMP A   1       0.000   0.000   0.000\n")

        mock_atoms = MagicMock()
        mock_atoms.get_positions.return_value = np.array([[0.0, 0.0, 0.0]])
        mock_atoms.cell = np.diag([30.0, 30.0, 50.0])

        with (
            patch("polymer_brush_tool.structure.graft.read", return_value=mock_atoms),
            patch("polymer_brush_tool.structure.graft.runner.run"),
        ):
            from polymer_brush_tool.structure.graft import graft_brush
            graft_brush(chain_pdb, 30.0, 30.0, 2, 1, work_dir=tmp_path)

        text = (tmp_path / "graft_brush.in").read_text()
        assert "chain A" in text
        assert "chain B" in text
        assert "chain C" not in text


# --------------------------------------------------------------------------
# pin_linkers_to_substrate
# --------------------------------------------------------------------------

class TestPinLinkersToSubstrate:
    def _atoms(self):
        from ase import Atoms
        # two "chains": linker atoms 0 and 3 at slightly different heights,
        # a non-linker atom (2) is the lowest point of the system
        return Atoms("H4", positions=[[0, 0, 5.0], [0, 0, 8.0], [1, 1, 4.0], [3, 3, 5.5]])

    def test_default_height_is_bond_length(self):
        from polymer_brush_tool.structure.atoms import pin_linkers_to_substrate
        atoms = pin_linkers_to_substrate(self._atoms(), [0, 3])
        z = atoms.get_positions()[:, 2]
        assert z[0] == pytest.approx(1.5)
        assert z[3] == pytest.approx(1.5)

    def test_rigid_shift_preserves_other_atoms(self):
        from polymer_brush_tool.structure.atoms import pin_linkers_to_substrate
        atoms = pin_linkers_to_substrate(self._atoms(), [0, 3], height=1.5)
        z = atoms.get_positions()[:, 2]
        # min_z was 4.0 (atom 2): linkers moved onto that plane, whole system +(1.5 - 4.0)
        assert z[1] == pytest.approx(8.0 - 4.0 + 1.5)
        assert z[2] == pytest.approx(1.5)
        assert atoms.get_positions()[1, :2].tolist() == [0, 0]

    def test_height_zero_puts_linkers_on_wall(self):
        from polymer_brush_tool.structure.atoms import pin_linkers_to_substrate
        atoms = pin_linkers_to_substrate(self._atoms(), [0, 3], height=0.0)
        assert atoms.get_positions()[[0, 3], 2] == pytest.approx([0.0, 0.0])

    def test_rejects_bad_input(self):
        from polymer_brush_tool.structure.atoms import pin_linkers_to_substrate
        with pytest.raises(ValueError):
            pin_linkers_to_substrate(self._atoms(), [])
        with pytest.raises(ValueError):
            pin_linkers_to_substrate(self._atoms(), [0], height=-1.0)


# --------------------------------------------------------------------------
# remove_water_below
# --------------------------------------------------------------------------

_GRO = """\
brush
    9
    1HMP     C1    1   0.100   0.100   0.150
    1HMP     H1    2   0.100   0.100   0.250
    2SOL     OW    3   0.500   0.500   0.020
    2SOL    HW1    4   0.500   0.500   0.110
    2SOL    HW2    5   0.560   0.500   0.050
    3SOL     OW    6   1.000   1.000   0.600
    3SOL    HW1    7   1.000   1.000   0.690
    3SOL    HW2    8   1.060   1.000   0.630
    4SOL     OW    9   1.500   1.500   0.290
   2.98100   2.98100   6.30700
"""


class TestRemoveWaterBelow:
    def test_removes_low_water_and_renumbers(self, tmp_path):
        from polymer_brush_tool.structure.solvent import remove_water_below
        src = tmp_path / "in.gro"; dst = tmp_path / "out.gro"
        src.write_text(_GRO)
        kept, removed = remove_water_below(src, dst, z_min=3.0)
        assert (kept, removed) == (1, 2)          # OW at 0.2 Å and 2.9 Å removed, 6.0 Å kept
        lines = dst.read_text().splitlines()
        assert lines[0] == "brush"
        assert int(lines[1]) == 5
        atom_lines = lines[2:7]
        assert [l[5:10].strip() for l in atom_lines] == ["HMP", "HMP", "SOL", "SOL", "SOL"]
        assert [int(l[15:20]) for l in atom_lines] == [1, 2, 3, 4, 5]
        assert atom_lines[2][36:44].strip() == "0.600"
        assert lines[-1].strip() == "2.98100   2.98100   6.30700"

    def test_zero_threshold_keeps_everything(self, tmp_path):
        from polymer_brush_tool.structure.solvent import remove_water_below
        src = tmp_path / "in.gro"
        src.write_text(_GRO)
        kept, removed = remove_water_below(src, src, z_min=0.0)
        assert (kept, removed) == (3, 0)
        assert src.read_text().splitlines()[1].strip() == "9"

    def test_polymer_below_threshold_untouched(self, tmp_path):
        from polymer_brush_tool.structure.solvent import remove_water_below
        src = tmp_path / "in.gro"; dst = tmp_path / "out.gro"
        src.write_text(_GRO)
        remove_water_below(src, dst, z_min=100.0)
        assert "HMP" in dst.read_text() and "SOL" not in dst.read_text()
