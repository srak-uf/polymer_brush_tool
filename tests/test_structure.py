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
