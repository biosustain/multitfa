from copy import deepcopy

import numpy as np
from pandas import DataFrame, Series, concat

from .sampling_util import *
from .variability import variability
from copy import copy


def cutoff_sampling(
    model_variability,
    cutoff=100,
    variable_list=None,
    min_growth=False,
    fraction_of_optim=0.9,
):
    """Implements the quadratic constraint using repeated sampling on the surface of ellipsoid. Exits when 100 consecutive samples represent better solution. We fix the component sphere variables lb & ub to the sampled covariance and solve the problem.

    Parameters
    ----------
    model_variability : cobra model
        cobra model to which variability to apply
    cutoff : int, optional
        number of consecutive samples before we exit the sampling, by default 100
    variable_list : list, optional
        list of variable names to perform TVA on, by default None
    min_growth : bool, optional
        Boolean, to add minimum growth constraint or not, by default False
    fraction_of_optim : float, optional
        fraction of original growth/flux value, by default 0.9

    Returns
    -------
    tuple
        tuple of no.of samples taken to achieve optima and optimal ranges of variables (pd.Dataframe)

    Raises
    ------
    ValueError
        Initial check to see if model is feasible with given constraints
    """
    model = preprocess_model(model_variability)

    # Retrieve small and large sphere variables
    small_sphere_vars = [
        var for var in model.variables if var.name.startswith("Sphere_s_")
    ]
    large_sphere_vars = [
        var for var in model.variables if var.name.startswith("Sphere_l_")
    ]

    if variable_list == None:
        variables = [var.name for var in model.solver.variables]
    else:
        variables = [var for var in variable_list]

    if np.isnan(model.slim_optimize()):
        raise ValueError("model infeasible with given constraints")

    if min_growth:
        if model.solver.objective.direction == "max":
            fva_old_objective = model.problem.Variable(
                "fva_old_objective", lb=fraction_of_optim * model.solver.objective.value
            )
        else:
            fva_old_objective = model.problem.Variable(
                "fva_old_objective", ub=fraction_of_optim * model.solver.objective.value
            )
        # Add the minimal growth/production constraint
        fva_old_obj_constraint = model.problem.Constraint(
            model.solver.objective.expression - fva_old_objective,
            lb=0,
            ub=0,
            name="fva_old_objective_constraint",
        )

        model.add_cons_vars([fva_old_obj_constraint, fva_old_objective])

    representative_ranges = DataFrame(
        np.zeros((len(variables), 2)), columns=["minimum", "maximum"]
    )

    n_improvement, total_samples = (0, 0)
    while n_improvement < cutoff:
        total_samples = total_samples + 1

        # Sample for components energy covariance ellipsoid
        small_sphr_sample = generate_n_sphere_sample(len(small_sphere_vars))
        large_sphr_sample = generate_n_sphere_sample(len(large_sphere_vars))

        # Fix the component variable lb, ub to sampled formation energy
        if len(small_sphere_vars) > 0:
            for i in range(len(small_sphere_vars)):
                # To avoid lb > ub fix them to larger bounds
                small_sphere_vars[i].lb = -1000
                small_sphere_vars[i].ub = 1000

                small_sphere_vars[i].lb = small_sphr_sample[i]
                small_sphere_vars[i].ub = small_sphr_sample[i]

        if len(large_sphere_vars) > 0:
            for i in range(len(large_sphere_vars)):
                # To avoid lb > ub fix them to larger bounds
                large_sphere_vars[i].lb = -1000
                large_sphere_vars[i].ub = 1000

                large_sphere_vars[i].lb = large_sphr_sample[i]
                large_sphere_vars[i].ub = large_sphr_sample[i]

        tva_ranges = variability(model, variable_list=variables)
        if tva_ranges.empty or tva_ranges.isnull().all()["maximum"]:
            total_samples = total_samples - 1
            continue
        else:
            flags = compare_dataframes(representative_ranges, tva_ranges)
            Y_count = flags.count("Y")
            # print(tva_ranges)
            if Y_count > 0.05 * len(flags):
                n_improvement = 0
                representative_ranges = tva_ranges
            else:
                n_improvement = n_improvement + 1

        lastone = tva_ranges
    return representative_ranges, total_samples


def gev_sampling(
    model_variability,
    cutoff=1000,
    variable_list=None,
    min_growth=False,
    fraction_of_optim=0.9,
):
    """Implements the quadratic constraint using repeated sampling on the surface of ellipsoid. After sampling for fixed number of times, we use generalised extreme value distribution to predict the possible extremum of the distribution. We fix the component sphere variables lb & ub to the sampled covariance and solve the problem.

    Parameters
    ----------
    model_variability : cobra model
        cobra model to which variability to apply
    cutoff : int, optional
        number of samples to train extreme value distribution, by default 1000
    variable_list : list, optional
        list of variable names to perform TVA on, by default None
    min_growth : bool, optional
        Boolean, to add minimum growth constraint or not, by default False
    fraction_of_optim : float, optional
        fraction of original growth/flux value, by default 0.9

    Returns
    -------
    pd.DataFrame
        pd.DataFrame of ranges of values for the variables
    """
    model = preprocess_model(model_variability)

    # Retrieve small and large sphere variables
    small_sphere_vars = [
        var for var in model.variables if var.name.startswith("Sphere_s_")
    ]
    large_sphere_vars = [
        var for var in model.variables if var.name.startswith("Sphere_l_")
    ]

    if variable_list == None:
        variables = [var.name for var in model.solver.variables]
    else:
        variables = [var for var in variable_list]

    if min_growth:
        if model.solver.objective.direction == "max":
            fva_old_objective = model.problem.Variable(
                "fva_old_objective", lb=fraction_of_optim * model.solver.objective.value
            )
        else:
            fva_old_objective = model.problem.Variable(
                "fva_old_objective", ub=fraction_of_optim * model.solver.objective.value
            )
        # Add the minimal growth/production constraint
        fva_old_obj_constraint = model.problem.Constraint(
            model.solver.objective.expression - fva_old_objective,
            lb=0,
            ub=0,
            name="fva_old_objective_constraint",
        )

        model.add_cons_vars([fva_old_obj_constraint, fva_old_objective])

    Whole_ranges = DataFrame(
        {
            "minimum": Series(index=variables, data=np.empty(len(variables))),
            "maximum": Series(index=variables, data=np.empty(len(variables))),
        }
    )

    total_samples = 0
    mins = []
    maxs = []

    while total_samples < cutoff:
        total_samples = total_samples + 1

        # Sample for components energy covariance ellipsoid
        small_sphr_sample = generate_n_sphere_sample(len(small_sphere_vars))
        large_sphr_sample = generate_n_sphere_sample(len(large_sphere_vars))

        # Fix the component variable lb, ub to sampled formation energy
        if len(small_sphere_vars) > 0:
            for i in range(len(small_sphere_vars)):
                # To avoid lb > ub fix them to larger bounds
                small_sphere_vars[i].lb = -1000
                small_sphere_vars[i].ub = 1000

                small_sphere_vars[i].lb = small_sphr_sample[i]
                small_sphere_vars[i].ub = small_sphr_sample[i]

        if len(large_sphere_vars) > 0:
            for i in range(len(large_sphere_vars)):
                # To avoid lb > ub fix them to larger bounds
                large_sphere_vars[i].lb = -1000
                large_sphere_vars[i].ub = 1000

                large_sphere_vars[i].lb = large_sphr_sample[i]
                large_sphere_vars[i].ub = large_sphr_sample[i]

        tva_ranges = variability(model, variable_list=variables)

        if tva_ranges.empty or tva_ranges.isnull().all()["maximum"]:
            total_samples = total_samples - 1
            continue

        Whole_ranges = concat([Whole_ranges, tva_ranges], axis=1)

    for var in variables:
        min_delG, max_d = extreme_value_distribution(Whole_ranges.loc[var, "minimum"])
        min_d, max_delG = extreme_value_distribution(Whole_ranges.loc[var, "maximum"])
        mins.append(min_delG)
        maxs.append(max_delG)

    representative_ranges = DataFrame(
        {
            "minimum": Series(index=variables, data=mins),
            "maximum": Series(index=variables, data=maxs),
        }
    )

    return representative_ranges


def sampling(
    model,
    cutoff=100,
    variable_list=None,
    min_growth=False,
    fraction_of_optim=0.9,
    exit_strat="gev",
):

    if exit_strat == "gev":
        var_ranges = gev_sampling(
            model,
            cutoff=cutoff,
            variable_list=variable_list,
            min_growth=min_growth,
            fraction_of_optim=fraction_of_optim,
        )
    else:
        var_ranges = cutoff_sampling(
            model,
            cutoff=cutoff,
            variable_list=variable_list,
            min_growth=min_growth,
            fraction_of_optim=fraction_of_optim,
        )
    return var_ranges


def generate_valid_sample(model, min_growth=False, fraction_of_optim=0.9):

    model_copy = deepcopy(model)
    model_copy.slim_optimize()
    if np.isnan(model_copy.slim_optimize()):
        raise ValueError("infeasible model")

    if min_growth:
        if model_copy.solver.objective.direction == "max":
            fva_old_objective = model_copy.problem.Variable(
                "fva_old_objective",
                lb=fraction_of_optim * model_copy.solver.objective.value,
            )
        else:
            fva_old_objective = model_copy.problem.Variable(
                "fva_old_objective",
                ub=fraction_of_optim * model_copy.solver.objective.value,
            )
        # Add the minimal growth/production constraint
        fva_old_obj_constraint = model_copy.problem.Constraint(
            model_copy.solver.objective.expression - fva_old_objective,
            lb=0,
            ub=0,
            name="fva_old_objective_constraint",
        )

        model_copy.add_cons_vars([fva_old_obj_constraint, fva_old_objective])

    while True:

        # Sample for formation energy covariance ellipsoid
        formation_sample = generate_ellipsoid_sample(model_copy.cholskey_matrix)

        # Fix the formation energy variable lb, ub to sampled formation energy
        for metabolite in model_copy.metabolites:
            if metabolite.std_dev > 50:
                continue
            # To avoid lb > ub fix them to larger bounds
            metabolite.compound_variable.lb = -1e9
            metabolite.compound_variable.ub = 1e9

            metabolite.compound_variable.lb = formation_sample[
                model_copy.metabolites.index(metabolite)
            ]
            metabolite.compound_variable.ub = formation_sample[
                model_copy.metabolites.index(metabolite)
            ]
        print(model_copy.slim_optimize())
        warm_start = {}
        if not np.isnan(model_copy.slim_optimize()):
            break

    for var in model_copy.variables:
        warm_start[var.name] = model_copy.solver.primal_values[var.name]

    return warm_start
