from cobra.util.solver import linear_reaction_coefficients
import pytest
import numpy as np
import sys

from .load_test_model import build_test_model


@pytest.fixture
def tfa_model():
    return build_test_model()


def test_no_cons(tfa_model):
    # tfa_model = build_core_model()
    num_cons = 2 * 3 * (
        len(tfa_model.reactions) - len(tfa_model.Exclude_reactions)
    ) + len(tfa_model.metabolites)

    num_vars = (
        2 * (len(tfa_model.metabolites))
        + 4 * (len(tfa_model.reactions) - len(tfa_model.Exclude_reactions))
        + 2 * len(tfa_model.reactions)
    )

    assert num_cons == len(tfa_model.constraints)
    assert num_vars == len(tfa_model.variables)