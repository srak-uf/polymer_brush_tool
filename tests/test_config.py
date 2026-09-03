"""Tests for BrushConfig loading and computed properties."""

import math
import textwrap
from pathlib import Path

import pytest
import yaml

from polymer_brush_tool.config import BrushConfig, MonomerSpec


# --------------------------------------------------------------------------
# BrushConfig: computed geometry helpers
# --------------------------------------------------------------------------

class TestBrushConfigGeometry:
    def test_box_x_square(self, linear_config):
        # rho=0.45, nx=ny=2 → box_x = sqrt(4/0.45)*10
        expected = math.sqrt(2 * 2 / 0.45) * 10
        assert math.isclose(linear_config.box_x(), expected, rel_tol=1e-9)

    def test_box_y_default_square(self, linear_config):
        assert math.isclose(linear_config.box_y(), linear_config.box_x(), rel_tol=1e-9)

    def test_box_y_independent_of_xyratio(self, linear_config):
        # xyratio only affects the PACKMOL alignment box; the simulation box
        # stays square so the graft density is exactly rho (legacy behaviour).
        linear_config.xyratio = 0.5
        assert math.isclose(linear_config.box_y(), linear_config.box_x(), rel_tol=1e-9)
        area_nm2 = linear_config.box_x() * linear_config.box_y() / 100
        assert math.isclose(linear_config.nx * linear_config.ny / area_nm2, linear_config.rho, rel_tol=1e-9)

    def test_n_cc_all(self, linear_config):
        # head_n_cc=1, mid_n_cc=1, tail_n_cc=1, n_mid=12
        # (1+1) + 12*(1+1) + (1+1) + 2 = 2 + 24 + 2 + 2 = 30
        assert linear_config.n_cc_all() == 30

    def test_polymer_length_auto(self, linear_config):
        expected = 1.54 * 30 * 0.8
        assert math.isclose(linear_config.polymer_length(), expected, rel_tol=1e-9)

    def test_polymer_length_explicit(self, linear_config):
        linear_config.d_polymer = 99.9
        assert linear_config.polymer_length() == 99.9

    def test_loop_height(self, linear_config):
        # Half of polymer_length
        assert math.isclose(
            linear_config.loop_height(),
            linear_config.polymer_length() / 2,
            rel_tol=1e-9,
        )


# --------------------------------------------------------------------------
# BrushConfig: YAML round-trip
# --------------------------------------------------------------------------

class TestBrushConfigFromYaml:
    def _write_yaml(self, tmp_path: Path, data: dict) -> Path:
        p = tmp_path / "config.yaml"
        p.write_text(yaml.dump(data))
        return p

    def test_basic_round_trip(self, tmp_path, linear_config):
        data = {
            "topology": "linear",
            "n_mid_repeat_units": 12,
            "rho": 0.45,
            "nx": 2,
            "ny": 2,
            "head": {
                "resname": "hmp",
                "ac_file": "mpc.ac",
                "omitnames": ["H24"],
                "tailname": "C2",
                "post_tailtype": "c3",
                "termname": "H1",
                "n_cc": 1,
            },
            "mid": {
                "resname": "mmp",
                "ac_file": "mpc.ac",
                "omitnames": ["H23", "H24"],
                "headname": "C11",
                "tailname": "C2",
                "pre_headtype": "c3",
                "post_tailtype": "c3",
                "n_cc": 1,
            },
            "tail": {
                "resname": "tmp",
                "ac_file": "mpc.ac",
                "omitnames": ["H23"],
                "headname": "C11",
                "pre_headtype": "c3",
                "termname": "H24",
                "n_cc": 1,
            },
        }
        yaml_path = self._write_yaml(tmp_path, data)
        cfg = BrushConfig.from_yaml(yaml_path)

        assert cfg.topology == "linear"
        assert cfg.n_mid_repeat_units == 12
        assert isinstance(cfg.head, MonomerSpec)
        assert cfg.head.resname == "hmp"
        assert cfg.mid.omitnames == ["H23", "H24"]
        assert cfg.tail.termname == "H24"

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            BrushConfig.from_yaml("/does/not/exist.yaml")

    def test_optional_bottom_atom_index(self, tmp_path):
        data = {
            "topology": "linear",
            "n_mid_repeat_units": 5,
            "rho": 0.45,
            "nx": 1,
            "ny": 1,
            "bottom_atom_index": 37,
            "head": {"resname": "hmp", "ac_file": "x.ac", "omitnames": ["H1"],
                     "termname": "H9", "tailname": "C2", "post_tailtype": "c3"},
            "mid": {"resname": "mmp", "ac_file": "x.ac", "omitnames": ["H1", "H2"],
                    "headname": "C1", "tailname": "C2",
                    "pre_headtype": "c3", "post_tailtype": "c3"},
            "tail": {"resname": "tmp", "ac_file": "x.ac", "omitnames": ["H1"],
                     "termname": "H8", "headname": "C1", "pre_headtype": "c3"},
        }
        yaml_path = self._write_yaml(tmp_path, data)
        cfg = BrushConfig.from_yaml(yaml_path)
        assert cfg.bottom_atom_index == 37
