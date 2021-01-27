from multitfa.analysis import generate_n_sphere_sample, preprocess_model
import numpy as np
import pytest


def test_sphere():
    sphere_sample = generate_n_sphere_sample(3)
    assert abs(np.sqrt(np.sum(np.square(sphere_sample))) - 1) < 1e-3


from .load_test_model import build_test_model


@pytest.fixture
def tfa_model():
    return build_test_model()


def test_preprocess(tfa_model):
    test_model = preprocess_model(tfa_model)

    check_vars = [
        var
        for var in test_model.variables
        if var.name.startswith("component_") or var.name.startswith("dG_err_")
    ]
    assert check_vars == []