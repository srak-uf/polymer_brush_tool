"""Shared pytest fixtures."""

import pytest
from pathlib import Path

from polymer_brush_tool.config import BrushConfig, MonomerSpec


@pytest.fixture
def mpc_head():
    return MonomerSpec(
        resname="hmp",
        ac_file="mpc.ac",
        termname="H1",
        tailname="C2",
        omitnames=["H24"],
        post_tailtype="c3",
        n_cc=1,
    )


@pytest.fixture
def mpc_mid():
    return MonomerSpec(
        resname="mmp",
        ac_file="mpc.ac",
        headname="C11",
        tailname="C2",
        omitnames=["H23", "H24"],
        pre_headtype="c3",
        post_tailtype="c3",
        n_cc=1,
    )


@pytest.fixture
def mpc_tail():
    return MonomerSpec(
        resname="tmp",
        ac_file="mpc.ac",
        termname="H24",
        headname="C11",
        omitnames=["H23"],
        pre_headtype="c3",
        n_cc=1,
    )


@pytest.fixture
def linear_config(mpc_head, mpc_mid, mpc_tail):
    return BrushConfig(
        topology="linear",
        n_mid_repeat_units=12,
        rho=0.45,
        nx=2,
        ny=2,
        head=mpc_head,
        mid=mpc_mid,
        tail=mpc_tail,
    )
