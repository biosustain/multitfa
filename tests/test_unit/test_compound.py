from numpy.testing._private.utils import assert_approx_equal
import pytest
import numpy as np
from multitfa.core import Thermo_met

from .load_test_model import build_test_model


@pytest.fixture
def tfa_met():
    tfa_model = build_test_model()
    return tfa_model.metabolites.get_by_id("atp_c")


def test_variables(tfa_met):
    assert tfa_met.concentration_variable
    assert tfa_met.delG_err_variable
    assert np.any(tfa_met.compound_vector)


def test_thermo_property(tfa_met):
    assert_approx_equal(tfa_met.delG_f, -2259.1882733696866, significant=3)
    assert not tfa_met.is_proton
    assert tfa_met.equilibrator_accession
    assert tfa_met.Kegg_id == "bigg.metabolite:atp"


def test_vars_bounds(tfa_met):
    assert np.log(tfa_met.concentration_min) == tfa_met.concentration_variable.lb
    assert np.log(tfa_met.concentration_max) == tfa_met.concentration_variable.ub