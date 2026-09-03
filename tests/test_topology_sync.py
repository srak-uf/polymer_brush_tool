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
    """The package ships copies of files that also live in the repo; keep them in sync."""

    @pytest.mark.parametrize(
        "pkg_name, repo_path",
        [
            ("min_vac.mdp", "01_FF_template/min_vac.mdp"),
            ("nvt_vac.mdp", "01_FF_template/nvt_vac.mdp"),
            ("tip3p.itp", "02_MD_template/tip3p.itp"),
            ("linear_config.yaml", "examples/linear_config.yaml"),
            ("loop_config.yaml", "examples/loop_config.yaml"),
        ],
    )
    def test_template_matches_repo_file(self, pkg_name, repo_path):
        pkg = (resources.files("polymer_brush_tool") / "templates" / pkg_name).read_text()
        repo = (REPO / repo_path).read_text()
        assert pkg == repo, f"{pkg_name} differs from {repo_path}; re-copy it"
