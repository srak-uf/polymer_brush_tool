"""Tests for copy_molecules_section and the shipped templates."""

from importlib import resources
from pathlib import Path

import pytest

from polymer_brush_tool.ff.topology import copy_molecules_section

REPO = Path(__file__).resolve().parent.parent

_HARD = """\
[ moleculetype ]
mol 3

[ system ]
brush

[ molecules ]
; Compound  #mols
mol  4
SOL  1234
"""

_SOFT = """\
[ moleculetype ]
mol 3

[ system ]
brush

[ molecules ]
; Compound  #mols
mol  4
"""


class TestCopyMoleculesSection:
    def test_sol_line_propagates(self, tmp_path):
        hard = tmp_path / "hard.top"
        soft = tmp_path / "soft.top"
        hard.write_text(_HARD)
        soft.write_text(_SOFT)
        copy_molecules_section(hard, soft)
        text = soft.read_text()
        assert "SOL  1234" in text
        # Everything before [ molecules ] is untouched
        assert text.startswith("[ moleculetype ]\nmol 3\n")
        assert text.count("[ molecules ]") == 1

    def test_missing_section_raises(self, tmp_path):
        hard = tmp_path / "hard.top"
        soft = tmp_path / "soft.top"
        hard.write_text("[ system ]\nx\n")
        soft.write_text(_SOFT)
        with pytest.raises(ValueError, match="molecules"):
            copy_molecules_section(hard, soft)


class TestBundledTemplates:
    """The bundled config templates must match the commented examples in examples/."""

    @pytest.mark.parametrize(
        "pkg_name, repo_path",
        [
            ("linear_config.yaml", "examples/linear_config.yaml"),
            ("loop_config.yaml", "examples/loop_config.yaml"),
        ],
    )
    def test_template_matches_repo_file(self, pkg_name, repo_path):
        pkg = (resources.files("polymer_brush_tool") / "templates" / pkg_name).read_text()
        repo = (REPO / repo_path).read_text()
        assert pkg == repo, f"{pkg_name} differs from {repo_path}; re-copy it"


class TestSetMoleculeCount:
    def test_sol_count_rewritten(self, tmp_path):
        from polymer_brush_tool.ff.topology import set_molecule_count
        top = tmp_path / "hard.top"
        top.write_text(_HARD)
        set_molecule_count(top, "SOL", 1200)
        text = top.read_text()
        assert "SOL              1200" in text and "SOL  1234" not in text
        assert "mol  4" in text                     # other molecules untouched

    def test_missing_line_raises(self, tmp_path):
        from polymer_brush_tool.ff.topology import set_molecule_count
        top = tmp_path / "soft.top"
        top.write_text(_SOFT)
        with pytest.raises(ValueError, match="SOL"):
            set_molecule_count(top, "SOL", 1)


class TestRemoveMolecule:
    def test_removes_all_sol_lines(self, tmp_path):
        from polymer_brush_tool.ff.topology import remove_molecule
        top = tmp_path / "hard.top"
        top.write_text(_HARD + "SOL  99\n")
        assert remove_molecule(top, "SOL") == 2
        text = top.read_text()
        assert "SOL" not in text and "mol  4" in text
        assert text.count("[ moleculetype ]") == 1     # earlier sections untouched

    def test_no_sol_is_noop(self, tmp_path):
        from polymer_brush_tool.ff.topology import remove_molecule
        top = tmp_path / "soft.top"
        top.write_text(_SOFT)
        assert remove_molecule(top, "SOL") == 0
        assert top.read_text() == _SOFT
