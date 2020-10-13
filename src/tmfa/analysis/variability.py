from copy import copy

import numpy as np
from pandas import DataFrame, Series


def variability(
    model_variability, variable_list=None, min_growth=False, fraction_of_optim=0.9
):
    """Perform thermodynamic variability analysis.

    Determine the minimum and maximum values for the input variables (Flux, Gibbs free
    energies and metabolite concentrations). if min_growth constraint is applied then
    growth is maintained at given percentage of optimum.

    :param model: thermodynamic model
    :type model: core.model
    :param variable_list: list of variables, if None performs on all the variables in the model, defaults to None
    :type variable_list: list of str, optional
    :param min_growth: if minimum growth constraint is needed, defaults to False
    :type min_growth: bool, optional
    :param fraction_of_optim: fraction of optimal growth rate needs to be maintained, defaults to 0.9
    :type fraction_of_optim: float, optional
    :raises ValueError: [description]
    :return: minimum and max of variables
    :rtype: pd.Dataframe

    """

    model = copy(model_variability)
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

    fluxes_min = np.empty(len(variables))
    fluxes_max = np.empty(len(variables))
    rxn_name = list()

    rxn_ids = [rxn.id for rxn in model.reactions]

    for i in range(len(variables)):
        # For reaction flux, objective is forward - reverse variables
        if variables[i] in rxn_ids:
            rxn = model.reactions.get_by_id(variables[i])
            objective_exp = 1 * rxn.forward_variable - 1 * rxn.reverse_variable
        else:
            var = model.solver.variables[variables[i]]
            objective_exp = 1 * var

        rxn_name.append(variables[i])
        model.objective = objective_exp

        # minimization
        model.objective_direction = "min"
        _ = model.slim_optimize()
        objective_value = model.objective.value
        fluxes_min[i] = objective_value

        # maximiztion
        model.objective_direction = "max"
        _ = model.slim_optimize()
        objective_value = model.objective.value
        fluxes_max[i] = objective_value

    return DataFrame(
        {
            "minimum": Series(index=rxn_name, data=fluxes_min),
            "maximum": Series(index=rxn_name, data=fluxes_max),
        }
    )


from gurobipy import *


def variability_legacy_gurobi(
    model_variability,
    variable_list=None,
    min_growth=False,
    fraction_of_optim=0.9,
    warm_start={},
):
    model = copy(model_variability)
    if variable_list == None:
        variables = model.gurobi_interface.getVars()
    else:
        variables = [var for var in variable_list]

    model.gurobi_interface.optimize()
    # if warm start not provided start with box solution
    # if len(warm_start) == 0:
    #    model.slim_optimize()
    #    if model.solver.status == "optimal":
    #        for var in model.variables:
    #            if var.name.startswith("indicator_"):
    #                warm_start[var.name] = model.solver.primal_values[var.name]
    print("updated")
    if min_growth:
        if model.solver.objective.direction == "max":
            fva_old_objective = model.gurobi_interface.addVar(
                name="fva_old_objective",
                lb=fraction_of_optim * model.gurobi_interface.ObjVal,
            )
        else:
            fva_old_objective = model.gurobi_interface.addVar(
                name="fva_old_objective",
                ub=fraction_of_optim * model.gurobi_interface.ObjVal,
            )

        # Add the minimal growth/production constraint
        old_obj = model.gurobi_interface.getObjective()
        model.gurobi_interface.addConstr(
            old_obj - fva_old_objective,
            GRB.EQUAL,
            0,
            name="fva_old_objective_constraint",
        )

        model.gurobi_interface.update()

    fluxes_min = np.empty(len(variables))
    fluxes_max = np.empty(len(variables))
    rxn_name = list()

    rxn_ids = [rxn.id for rxn in model.reactions]

    for i in range(len(variables)):
        print(variables[i])
        # if the variable is reactions optimize for forward - reverse variables else optimize for the variable
        if variables[i] in rxn_ids:
            rxn = model.reactions.get_by_id(variables[i])
            for_var = model.gurobi_interface.getVarByName(rxn.forward_variable.name)
            rev_var = model.gurobi_interface.getVarByName(rxn.reverse_variable.name)
            obj_exp = for_var - rev_var
        else:
            obj_exp = model.gurobi_interface.getVarByName(variables[i])

        rxn_name.append(variables[i])

        # minimization
        if len(warm_start) != 0:
            for var in model.gurobi_interface.getVars():
                if var in warm_start:
                    var.Start = warm_start[var]

        model.gurobi_interface.setObjective(obj_exp, GRB.MINIMIZE)
        model.gurobi_interface.update()
        model.gurobi_interface.optimize()
        objective_value = model.gurobi_interface.ObjVal
        fluxes_min[i] = objective_value
        # print(rxn.id, "min", objective_value)

        warm_start = {}
        for var in model.gurobi_interface.getVars():
            if var.VarName.startswith("indicator_"):
                warm_start[var] = var.x

        # maximiztion
        for var in model.gurobi_interface.getVars():
            if var in warm_start:
                var.Start = warm_start[var]

        model.gurobi_interface.setObjective(obj_exp, GRB.MAXIMIZE)
        model.gurobi_interface.update()
        model.gurobi_interface.optimize()
        objective_value = model.gurobi_interface.ObjVal
        fluxes_max[i] = objective_value
        # print(rxn.id, "max", objective_value)

    return DataFrame(
        {
            "minimum": Series(index=rxn_name, data=fluxes_min),
            "maximum": Series(index=rxn_name, data=fluxes_max),
        }
    )


import os
import string
import tempfile
from random import choices


def variability_legacy_cplex(
    model,
    variable_list=None,
    min_growth=False,
    biomass_rxn=None,
    fraction_of_optim=0.9,
    warm_start={},
):
    # Instead of copying the whole model, just copy the cplex solver object by writing to a file and reading again.
    from cplex import Cplex, SparsePair

    tmp_dir = (
        os.path.normpath(os.path.dirname(os.path.abspath(__file__)))
        + os.sep
        + os.pardir
        + os.sep
        + "tmp"
    )
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # Instantiate Cplex model
    cplex_model = Cplex()
    rand_str = "".join(choices(string.ascii_lowercase + string.digits, k=6))

    # write cplex model to mps file and re read
    with tempfile.TemporaryDirectory() as td:
        temp_filename = os.path.join(td, rand_str + ".mps")
        model.cplex_interface.write(temp_filename)
        cplex_model.read(temp_filename)

    # Make shorts for sense
    max_sense = cplex_model.objective.sense.maximize
    min_sense = cplex_model.objective.sense.minimize

    if variable_list == None:
        variables = model.cplex_interface.variables.get_names()
    else:
        variables = [var for var in variable_list]

    vars_list_cplex = cplex_model.variables.get_names()
    if min_growth:

        if biomass_rxn == None:
            raise ValueError(
                "Please provide growth reaction name to add the minimum growth constraint"
            )
        # Make the biomass reaction objective

        for varname in vars_list_cplex:
            cplex_model.objective.set_linear(varname, 0)

        biomass_for_name = model.reactions.get_by_id(biomass_rxn).forward_variable.name
        biomass_rev_name = model.reactions.get_by_id(biomass_rxn).reverse_variable.name

        cplex_model.objective.set_linear(
            [(biomass_for_name, 1), (biomass_rev_name, -1)]
        )

        cplex_model.solve()

        if model.solver.objective.direction == "max":
            fva_old_objective = cplex_model.variables.add(names=["fva_old_objective"])
            cplex_model.variables.set_lower_bounds = (
                "fva_old_objective",
                fraction_of_optim * cplex_model.solution.get_objective_value(),
            )

        else:
            fva_old_objective = cplex_model.variables.add(names=["fva_old_objective"])
            cplex_model.variables.set_upper_bounds = (
                "fva_old_objective",
                fraction_of_optim * cplex_model.solution.get_objective_value(),
            )

        # Add the minimal growth/production constraint
        _ = cplex_model.linear_constraints.add(
            lin_expr=[
                SparsePair(
                    ind=[biomass_for_name, biomass_rev_name, fva_old_objective],
                    val=[1, -1, -1],
                )
            ],
            senses="E",
            rhs=[0],
            names="min_growth",
        )

    fluxes_min = np.empty(len(variables))
    fluxes_max = np.empty(len(variables))
    rxn_name = list()

    rxn_ids = [rxn.id for rxn in model.reactions]

    for i in range(len(variables)):
        # Reset objective vector for each iteration
        for varname in vars_list_cplex:
            cplex_model.objective.set_linear(varname, 0)

        print(variables[i])
        # if the variable is reactions optimize for forward - reverse variables else optimize for the variable
        if variables[i] in rxn_ids:
            rxn = model.reactions.get_by_id(variables[i])
            cplex_model.objective.set_linear(
                [(rxn.forward_variable.name, 1), (rxn.reverse_variable.name, -1)]
            )

        else:
            cplex_model.objective.set_linear(variables[i], 1)

        rxn_name.append(variables[i])

        # minimization
        cplex_model.objective.set_sense(min_sense)
        cplex_model.solve()
        objective_value = cplex_model.solution.get_objective_value()
        fluxes_min[i] = objective_value
        # print(rxn.id, "min", objective_value)

        # maximiztion
        cplex_model.objective.set_sense(max_sense)
        cplex_model.solve()
        objective_value = cplex_model.solution.get_objective_value()
        fluxes_max[i] = objective_value
        # print(rxn.id, "max", objective_value)

    return DataFrame(
        {
            "minimum": Series(index=rxn_name, data=fluxes_min),
            "maximum": Series(index=rxn_name, data=fluxes_max),
        }
    )
