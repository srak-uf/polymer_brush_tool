"""Tests for force-field preparation functions (no external tools required)."""

import pytest
from pathlib import Path

from polymer_brush_tool.ff.fragments import write_fragment, write_tleap, write_tleap_grafted
from polymer_brush_tool.ff.topology import insert_tip3p_top, insert_restraint_top


# --------------------------------------------------------------------------
# write_fragment
# --------------------------------------------------------------------------

class TestWriteFragment:
    def test_head_monomer(self, tmp_path):
        """HEAD monomer: no headname, only tail side."""
        p = write_fragment(
            "hmp",
            tailname="C2",
            omitnames=["H24"],
            post_tailtype="c3",
            work_dir=tmp_path,
        )
        text = p.read_text()
        assert "TAIL_NAME  C2" in text
        assert "POST_TAIL_TYPE  c3" in text
        assert "OMIT_NAME  H24" in text
        assert "CHARGE 0.000" in text
        assert "HEAD_NAME" not in text

    def test_mid_monomer(self, tmp_path):
        """MID monomer: both head and tail sides."""
        p = write_fragment(
            "mmp",
            headname="C11",
            tailname="C2",
            omitnames=["H23", "H24"],
            pre_headtype="c3",
            post_tailtype="c3",
            work_dir=tmp_path,
        )
        text = p.read_text()
        assert "HEAD_NAME  C11" in text
        assert "TAIL_NAME  C2" in text
        assert "OMIT_NAME  H23" in text
        assert "OMIT_NAME  H24" in text

    def test_tail_monomer(self, tmp_path):
        """TAIL monomer: only head side."""
        p = write_fragment(
            "tmp",
            headname="C11",
            omitnames=["H23"],
            pre_headtype="c3",
            work_dir=tmp_path,
        )
        text = p.read_text()
        assert "HEAD_NAME  C11" in text
        assert "TAIL_NAME" not in text

    def test_output_path(self, tmp_path):
        p = write_fragment("xyz", work_dir=tmp_path)
        assert p == tmp_path / "xyz.chain"
        assert p.exists()

    def test_missing_pre_headtype_raises(self, tmp_path):
        with pytest.raises(ValueError, match="pre_headtype"):
            write_fragment("mmp", headname="C11", omitnames=["H23"], work_dir=tmp_path)

    def test_missing_post_tailtype_raises(self, tmp_path):
        with pytest.raises(ValueError, match="post_tailtype"):
            write_fragment("hmp", tailname="C2", omitnames=["H24"], work_dir=tmp_path)

    def test_wrong_omitnames_length_raises(self, tmp_path):
        with pytest.raises(ValueError, match="Length of omitnames"):
            write_fragment(
                "mmp",
                headname="C11",
                tailname="C2",
                omitnames=["H23"],     # should have 2 entries
                pre_headtype="c3",
                post_tailtype="c3",
                work_dir=tmp_path,
            )

    def test_omitnames_string_accepted(self, tmp_path):
        """A plain string for omitnames should be wrapped in a list."""
        p = write_fragment(
            "hmp",
            tailname="C2",
            omitnames="H24",           # string, not list
            post_tailtype="c3",
            work_dir=tmp_path,
        )
        assert "OMIT_NAME  H24" in p.read_text()

    def test_omitnames_bad_type_raises(self, tmp_path):
        with pytest.raises(TypeError, match="omitnames must be"):
            write_fragment(
                "hmp",
                tailname="C2",
                omitnames=42,          # neither str nor list
                post_tailtype="c3",
                work_dir=tmp_path,
            )


# --------------------------------------------------------------------------
# write_tleap
# --------------------------------------------------------------------------

class TestWriteTleap:
    def test_output_path(self, tmp_path):
        p = write_tleap("hmp", "mmp", "tmp", 3, work_dir=tmp_path)
        assert p == tmp_path / "build_chain.tleap"
        assert p.exists()

    def test_sequence_content(self, tmp_path):
        p = write_tleap("hmp", "mmp", "tmp", 3, work_dir=tmp_path)
        text = p.read_text()
        assert "source leaprc.gaff2" in text
        assert "loadamberprep hmp.prepi" in text
        assert "loadamberprep mmp.prepi" in text
        assert "loadamberprep tmp.prepi" in text
        assert "HMP MMP MMP MMP  TMP" in text
        assert "saveAmberParm chain chain.prmtop chain.inpcrd" in text

    def test_zero_mid_units(self, tmp_path):
        p = write_tleap("hmp", "mmp", "tmp", 0, work_dir=tmp_path)
        text = p.read_text()
        # Sequence: HMP  TMP with no MMP
        assert "MMP" not in text


# --------------------------------------------------------------------------
# write_tleap_grafted
# --------------------------------------------------------------------------

class TestWriteTleapGrafted:
    def test_box_dimensions_written(self, tmp_path):
        p = write_tleap_grafted(
            "hmp", "mmp", "tmp",
            30.0, 30.0, 100.0,
            work_dir=tmp_path,
        )
        text = p.read_text()
        assert "set mol box {30.0 30.0 100.0}" in text

    def test_default_pdb_file(self, tmp_path):
        p = write_tleap_grafted("hmp", "mmp", "tmp", 1, 1, 1, work_dir=tmp_path)
        assert "loadpdb grafted_chain.pdb" in p.read_text()


# --------------------------------------------------------------------------
# insert_tip3p_top
# --------------------------------------------------------------------------

_SAMPLE_TOP = """\
[ defaults ]
1 2 yes 0.5 0.833

[ atomtypes ]
; ...

[ moleculetype ]
; mol 1
"""


class TestInsertTip3pTop:
    def test_include_inserted(self, tmp_path):
        src = tmp_path / "in.top"
        out = tmp_path / "out.top"
        src.write_text(_SAMPLE_TOP)
        insert_tip3p_top(src, out)
        text = out.read_text()
        assert '#include "tip3p.itp"' in text
        # The include should come BEFORE the [ moleculetype ] line
        idx_include = text.index('#include "tip3p.itp"')
        idx_mol = text.index("[ moleculetype ]")
        assert idx_include < idx_mol

    def test_no_moleculetype_unchanged(self, tmp_path):
        src = tmp_path / "in.top"
        out = tmp_path / "out.top"
        content = "[ defaults ]\n1 2 yes 0.5 0.833\n"
        src.write_text(content)
        insert_tip3p_top(src, out)
        assert out.read_text() == content


# --------------------------------------------------------------------------
# insert_restraint_top
# --------------------------------------------------------------------------

_SAMPLE_TOP_RESTRAINT = """\
[ moleculetype ]
mol 3

[ atoms ]
;   nr   type  resnr residue  atom   cgnr     charge       mass
     1   c3      1    HMP      C1      1    -0.0930   12.0100
     2   c3      1    HMP      C2      2    -0.0930   12.0100
     3   hc      1    HMP      H1      3     0.0310    1.0080

[ bonds ]
; ...

[ system ]
grafted brush system
"""


class TestInsertRestraintTop:
    def test_soft_restraints_written(self, tmp_path):
        src = tmp_path / "in.top"
        out = tmp_path / "out.top"
        src.write_text(_SAMPLE_TOP_RESTRAINT)
        insert_restraint_top(src, out, linker_indices=[1, 3])
        text = out.read_text()
        assert "[ position_restraints ]" in text
        assert "1  1  10000  10000  10000" in text
        assert "3  1  10000  10000  10000" in text

    def test_hard_restraint_file_created(self, tmp_path):
        src = tmp_path / "in.top"
        out = tmp_path / "out.top"
        src.write_text(_SAMPLE_TOP_RESTRAINT)
        insert_restraint_top(src, out, linker_indices=[1])
        hard_out = tmp_path / "hardrest_out.top"
        assert hard_out.exists()
        text = hard_out.read_text()
        assert "1000000" in text

    def test_invalid_index_skipped(self, tmp_path):
        """An index that is not in the [ atoms ] section should be silently skipped."""
        src = tmp_path / "in.top"
        out = tmp_path / "out.top"
        src.write_text(_SAMPLE_TOP_RESTRAINT)
        insert_restraint_top(src, out, linker_indices=[999])
        text = out.read_text()
        assert "999" not in text

    def test_restraints_before_system(self, tmp_path):
        src = tmp_path / "in.top"
        out = tmp_path / "out.top"
        src.write_text(_SAMPLE_TOP_RESTRAINT)
        insert_restraint_top(src, out, linker_indices=[1])
        text = out.read_text()
        idx_rest = text.index("[ position_restraints ]")
        idx_sys = text.index("[ system ]")
        assert idx_rest < idx_sys
