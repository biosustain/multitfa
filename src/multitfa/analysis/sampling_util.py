import numpy as np
from optlang import Constraint
from scipy import stats

from ..util.constraints import *
from ..util.linalg_fun import *
from ..util.thermo_constants import *


def generate_n_sphere_sample(n_variables):
    """Generates unit n-sphere sample. Works by picking random sample from normal distribution and normalized by radius.

    :param n_variables: Dimension of the n-sphere.
    :type n_variables: float
    :return: unit n-sphere sample
    :rtype: np.ndarray
    """
    # n-sphere sample from N(0,1)
    random_sample = np.random.normal(loc=0, scale=1.0, size=(n_variables))
    circle_radius = np.sqrt(np.sum(np.square(random_sample)))
    normalized_sample = random_sample / circle_radius

    return normalized_sample


def generate_ellipsoid_sample(cholesky):
    """sampling on the surface of n-dimensional ellipsoid
    sample on n-ellipsoid  is linear transformation of unit n-sphere
    N(mu,var) = mu + A @ N(0,1)
    A@A' = var, A is cholesky matrix

    :param cholesky: cholesky matrix
    :type cholesky: np.ndarray
    :return: ellipsoid sample
    :rtype: np.ndarray
    """

    n_dimensions = len(cholesky)
    chi_crit_val = chi2.isf(q=0.05, df=n_dimensions)
    n_sphere_sample = generate_n_sphere_sample(n_dimensions)
    ellipsoid_sample = np.sqrt(chi_crit_val) * cholesky @ n_sphere_sample

    return ellipsoid_sample


def preprocess_model(model):
    # First remove the delG constraint and associated variables, we will add them later
    remove_vars = [
        var
        for var in model.variables
        if var.name.startswith("component_") or var.name.startswith("dG_err_")
    ]
    remove_cons = [
        cons
        for cons in model.constraints
        if cons.name.startswith("delG_") or cons.name.startswith("std_dev_")
    ]

    # Pick indices of components present in the current model
    model_component_indices = [
        i
        for i in range(model.compound_vector_matrix.shape[1])
        if np.any(model.compound_vector_matrix[:, i])
    ]

    # Reduced the compound_vector to contain only the non zero entries
    model_compound_vector = model.compound_vector_matrix[:, model_component_indices]

    # Now extract the sub covariance matrix containing only the components present in the model
    component_model_covariance = covariance[:, model_component_indices][
        model_component_indices, :
    ]

    # Now separate the compounds that have variance > 1000 and others to avoid numerical issues
    high_variance_indices = np.where(np.diag(component_model_covariance) > 1000)[0]
    low_variance_indices = np.where(np.diag(component_model_covariance) < 1000)[0]

    # Calculate cholesky matrix for two different covariance matrices
    if len(low_variance_indices) > 0:
        small_component_covariance = component_model_covariance[
            :, low_variance_indices
        ][low_variance_indices, :]

        cholesky_small_variance = matrix_decomposition(small_component_covariance)
        chi2_value_small = stats.chi2.isf(
            q=0.05, df=cholesky_small_variance.shape[1]
        )  # Chi-square value to map confidence interval

        sphere_s_vars = np.array(
            [
                model.problem.Variable("Sphere_s_{}".format(i), lb=-1, ub=1)
                for i in range(cholesky_small_variance.shape[1])
            ]
        )  # adding sphere variables for low variance compounds
        model.add_cons_vars(sphere_s_vars.tolist())

        for i in high_variance_indices:
            zeros_axis = np.zeros((cholesky_small_variance.shape[1],))
            cholesky_small_variance = np.insert(
                cholesky_small_variance, i, zeros_axis, axis=0
            )

        metabolite_sphere_small = (
            model_compound_vector @ cholesky_small_variance
        )  # This is a fixed term compound_vector @ cholesky

    if len(high_variance_indices) > 0:
        large_component_covariance = component_model_covariance[
            :, high_variance_indices
        ][
            high_variance_indices, :
        ]  # Covariance matrix for the high variance components

        cholesky_large_variance = matrix_decomposition(large_component_covariance)
        chi2_value_high = stats.chi2.isf(q=0.05, df=cholesky_large_variance.shape[1])

        sphere_l_vars = np.array(
            [
                model.problem.Variable("Sphere_l_{}".format(i), lb=-1, ub=1)
                for i in range(cholesky_large_variance.shape[1])
            ]
        )  # adding sphere variables for high variance compounds
        model.add_cons_vars(sphere_l_vars.tolist())

        # Insert empty rows for the low_variance_components
        for i in low_variance_indices:
            zeros_axis = np.zeros((cholesky_large_variance.shape[1],))
            cholesky_large_variance = np.insert(
                cholesky_large_variance, i, zeros_axis, axis=0
            )
        metabolite_sphere_large = (
            model_compound_vector @ cholesky_large_variance
        )  # This is a fixed term compound_vector @ cholesky

    small_sphere_vars = np.array(
        [var for var in model.variables if var.name.startswith("Sphere_s_")]
    )
    large_sphere_vars = np.array(
        [var for var in model.variables if var.name.startswith("Sphere_l_")]
    )

    delG_constraints = []
    for rxn in model.reactions:
        if rxn.id in model.Exclude_reactions:
            continue
        S_vector = rxn.cal_stoichiometric_matrix()
        concentration_term = sum(
            stoic * metabolite.concentration_variable
            for metabolite, stoic in iteritems(rxn.metabolites)
            if metabolite.equilibrator_accession.inchi_key != PROTON_INCHI_KEY
        )

        if len(high_variance_indices) > 0:
            coefficients_high_var = (
                np.sqrt(chi2_value_high) * S_vector @ metabolite_sphere_large
            )
            err_expression_large = (
                coefficients_high_var[np.nonzero(coefficients_high_var)]
                @ large_sphere_vars[np.nonzero(coefficients_high_var)]
            )
        else:
            err_expression_large = 0

        if len(low_variance_indices) > 0:
            coefficients_small_var = (
                np.sqrt(chi2_value_small) * S_vector @ metabolite_sphere_small
            )
            err_expression_small = (
                coefficients_small_var[np.nonzero(coefficients_small_var)]
                @ small_sphere_vars[np.nonzero(coefficients_small_var)]
            )
        else:
            err_expression_small = 0

        lhs_forward = (
            rxn.delG_forward
            - RT * concentration_term
            - err_expression_small
            - err_expression_large
        )
        lhs_reverse = (
            rxn.delG_reverse
            + RT * concentration_term
            + err_expression_small
            + err_expression_large
        )
        rhs = rxn.delG_prime + rxn.delG_transport

        delG_f = model.problem.Constraint(
            lhs_forward,
            lb=rhs,
            ub=rhs,
            name="delG_{}".format(rxn.forward_variable.name),
        )
        delG_r = model.problem.Constraint(
            lhs_reverse,
            lb=-rhs,
            ub=-rhs,
            name="delG_{}".format(rxn.reverse_variable.name),
        )
        delG_constraints.extend([delG_f, delG_r])
    model.add_cons_vars(delG_constraints)

    return model


def compare_dataframes(df1, df2):

    flags = []

    for i in range(len(df1)):
        range1 = df1["maximum"][i] - df1["minimum"][i]
        range2 = df2["maximum"][i] - df2["minimum"][i]

        if range2 > range1 * 1.05:
            flags.append("Y")
        else:
            flags.append("N")
    return flags


def extreme_value_distribution(data_set):
    """Fits the Gibbs free energy data to the Generalized extreme value distribution and predicts the extreme value at 95 % CI.
    Uses Scipy genextreme function

    Arguments:
        data_set [list] -- The max or min range of Gibbs free energy values

    Returns:
        tuple -- min or max value predicted from GEV at 99% confidence
    """
    c, loc, scale = stats.genextreme.fit(data_set)
    min_extreme, max_extreme = stats.genextreme.interval(0.99, c, loc, scale)

    return min_extreme, max_extreme
