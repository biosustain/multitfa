import sys

import numpy as np
import optlang
import pytest
from numpy.testing._private.utils import assert_approx_equal, assert_array_equal

from .load_test_model import build_test_model


@pytest.fixture
def tfa_rxn():
    tfa_model = build_test_model()
    return tfa_model.reactions.get_by_id("ATPS4r")


def test_thermo_vars(tfa_rxn):
    assert tfa_rxn.delG_forward
    assert tfa_rxn.delG_reverse
    assert tfa_rxn.indicator_forward
    assert tfa_rxn.indicator_reverse


def test_thermo_cons(tfa_rxn):
    assert tfa_rxn.delG_constraint
    assert tfa_rxn.delG_constraint[0].ub == tfa_rxn.delG_prime + tfa_rxn.delG_transport
    assert tfa_rxn.indicator_forward.type == "binary"
    assert tfa_rxn.indicator_reverse.type == "binary"


def test_stoichiometry(tfa_rxn):
    assert_array_equal(
        tfa_rxn.cal_stoichiometric_matrix().tolist(),
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            3.0,
            -4.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )


def test_delG_prime(tfa_rxn):
    assert_approx_equal(tfa_rxn.delG_prime, 32.396, significant=3)
