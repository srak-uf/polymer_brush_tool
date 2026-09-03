"""Tests for BrushConfig.validate and YAML error reporting."""

import dataclasses

import pytest
import yaml

from polymer_brush_tool.config import BrushConfig, ConfigError


def test_valid_linear_passes(linear_config):
    linear_config.validate()


def test_loop_requires_d_polymer(linear_config):
    cfg = dataclasses.replace(linear_config, topology="loop", d_polymer=None)
    with pytest.raises(ConfigError, match="d_polymer"):
        cfg.validate()
    cfg.d_polymer = 14.9
    cfg.validate()


def test_bad_topology(linear_config):
    cfg = dataclasses.replace(linear_config, topology="ring")
    with pytest.raises(ConfigError, match="topology"):
        cfg.validate()


def test_missing_monomer(linear_config):
    cfg = dataclasses.replace(linear_config, tail=None)
    with pytest.raises(ConfigError, match="'tail'"):
        cfg.validate()


def test_missing_termname(linear_config):
    cfg = dataclasses.replace(linear_config, head=dataclasses.replace(linear_config.head, termname=None))
    with pytest.raises(ConfigError, match="head.termname"):
        cfg.validate()


def test_bad_linker_spec(linear_config):
    cfg = dataclasses.replace(linear_config, linker_atoms=[{"resname": "HMP"}])
    with pytest.raises(ConfigError, match="linker_atoms"):
        cfg.validate()


def test_unknown_yaml_key_reported(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text(yaml.dump({"topology": "linear", "rhoo": 0.4}))
    with pytest.raises(ConfigError, match="rhoo"):
        BrushConfig.from_yaml(p)


def test_example_configs_load_and_validate():
    from pathlib import Path
    repo = Path(__file__).resolve().parent.parent
    lin = BrushConfig.from_yaml(repo / "examples/linear_config.yaml")
    loop = BrushConfig.from_yaml(repo / "examples/loop_config.yaml")
    assert lin.topology == "linear" and lin.d_polymer is None
    assert loop.topology == "loop" and loop.d_polymer == 14.9
    assert (loop.nx, loop.ny, loop.xyratio) == (1, 2, 0.5)
