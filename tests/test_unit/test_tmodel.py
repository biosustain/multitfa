import sys

import numpy as np
import optlang
import pytest
from cobra.util.solver import linear_reaction_coefficients
from numpy.testing._private.utils import assert_almost_equal
from optlang.util import solve_with_glpsol

from .load_test_model import build_test_model


@pytest.fixture
def tfa_model():
    return build_test_model()


def test_num_cons_vars(tfa_model):
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


def test_solver_instances(tfa_model):
    if optlang.available_solvers["GUROBI"]:
        tfa_model.solver = "gurobi"
        assert tfa_model.gurobi_interface
    elif optlang.available_solvers["CPLEX"]:
        tfa_model.solver = "cplex"
        assert tfa_model.cplex_interface
    else:
        pass


def test_optimization(tfa_model):
    solution = tfa_model.optimize()
    assert_almost_equal(abs(solution.objective_value), 0.8739, decimal=3)
