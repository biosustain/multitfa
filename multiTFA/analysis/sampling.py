from .sampling_util import (
    generate_constraints_sample,
    generate_ellipsoid_sample,
    compare_dataframes,
    extreme_value_distribution,
)
from ..util.constraints import stddev_sampling_rhs
from .variability import variability
from pandas import DataFrame, Series, concat
import numpy as np
from warnings import warn


def cutoff_sampling(
    model, cutoff=100, variable_list=None, min_growth=False, fraction_of_optim=0.9
):
    """ Implements the quadratic constraint using repeated sampling on the surface of ellipsoid. Exits when 100 consecutive samples represent better solution. We fix the formation energy variable lb & ub to the sampled covariance and solve the problem.

    :param model: tmodel model with constraints
    :type model: core.model
    :param cutoff: number of consecutive samples with better solution, defaults to 100
    :type cutoff: int, optional
    :param variable_list: variable list run variability analysis, defaults to None
    :type variable_list: List, optional
    :param min_growth: if minimum growth constraint is needed, if yes, then a constraint is added to enforce cutoff percentage of objective is maintained, defaults to False
    :type min_growth: bool, optional
    :param fraction_of_optim: Percentage fraction of optimal growth , defaults to 0.9
    :type fraction_of_optim: float, optional
    :return: Ranges of variables
    :rtype: Pd.DataFrame
    """
    if variable_list == None:
        variables = [var.name for var in model.solver.variables]
    else:
        variables = [var for var in variable_list]

    if min_growth:
        if model.solver.objective.direction == "max":
            fva_old_objective = model.problem.variable(
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

        # Sample for formation energy covariance ellipsoid
        formation_sample = generate_ellipsoid_sample(model.cholskey_matrix)

        # Fix the formation energy variable lb, ub to sampled formation energy
        for metabolite in model.metabolites:
            # To avoid lb > ub fix them to larger bounds
            metabolite.compound_variable.lb = -1000
            metabolite.compound_variable.ub = 1000

            metabolite.compound_variable.lb = formation_sample[
                model.metabolits.index(metabolite)
            ]
            metabolite.compound_variable.ub = formation_sample[
                model.metabolites.index(metabolite)
            ]

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
    model, cutoff=1000, variable_list=None, min_growth=False, fraction_of_optim=0.9
):
    """[summary]

    :param model: [description]
    :type model: [type]
    :param cutoff: [description], defaults to 1000
    :type cutoff: int, optional
    :param variable_list: [description], defaults to None
    :type variable_list: [type], optional
    :param min_growth: [description], defaults to False
    :type min_growth: bool, optional
    :param fraction_of_optim: [description], defaults to 0.9
    :type fraction_of_optim: float, optional
    :return: [description]
    :rtype: [type]
    """
    if variable_list == None:
        variables = [var.name for var in model.solver.variables]
    else:
        variables = [var for var in variable_list]

    if min_growth:
        if model.solver.objective.direction == "max":
            fva_old_objective = model.problem.variable(
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

        # Sample for formation energy covariance ellipsoid
        formation_sample = generate_ellipsoid_sample(model.cholskey_matrix)

        # Fix the formation energy variable lb, ub to sampled formation energy
        for metabolite in model.metabolites:
            # To avoid lb > ub fix them to larger bounds
            metabolite.compound_variable.lb = -1000
            metabolite.compound_variable.ub = 1000

            metabolite.compound_variable.lb = formation_sample[
                model.metabolites.index(metabolite)
            ]
            metabolite.compound_variable.ub = formation_sample[
                model.metabolites.index(metabolite)
            ]

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
