from copy import copy

import numpy as np
from pandas import DataFrame, Series


def variability(model_variability, variable_list=None):
    """Perform thermodynamic variability analysis.

    Determine the minimum and maximum values for the input variables (Flux, Gibbs free
    energies and metabolite concentrations). if min_growth constraint is applied then
    growth is maintained at given percentage of optimum.

    Parameters
    ----------
    model_variability : multitfa.core.tmodel
        multitfa model after thermodynamic constraints are added
    variable_list : List, optional
        List of variables to perform TVA on, by default None


    Returns
    -------
    pd.DataFrame
        Dataframe of min max ranges of variables

    Raises
    ------
    ValueError
        If model is infeasible with initial constraints raises valueerror
    """

    model = copy(model_variability)
    if variable_list == None:
        variables = [var.name for var in model.solver.variables]
    else:
        variables = [var for var in variable_list]

    if np.isnan(model.slim_optimize()):
        raise ValueError("model infeasible with given constraints")

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


def variability_legacy_gurobi(
    model_variability,
    variable_list=None,
    warm_start={},
    params=False,
):
    """Custom function to perform TVA on MIQC problem using gurobi.

    Parameters
    ----------
    model_variability : multitfa.core.tmodel
        multitfa model after thermodynamic constraints are added
    variable_list : List, optional
        List of variables to perform TVA on, by default None
    warm_start : dict, optional
        Optionally can specify warm start to speed up problem. variable name and initial solution, by default {}
    params : Bool, optional
        If True sets the Timelimit option to 300 sec

    Returns
    -------
    pd.DataFrame
        Dataframe of min max ranges of variables
    """
    from gurobipy import GRB

    gurobi_interface = model_variability.gurobi_interface.copy()

    if variable_list == None:
        variables = gurobi_interface.getVars()
    else:
        variables = [var for var in variable_list]

    gurobi_interface.optimize()

    # Set the time limit searching for solution, useful for pathlogical variables taking long time
    if params:
        gurobi_interface.params.TimeLimit = 300

    fluxes_min = np.empty(len(variables))
    fluxes_max = np.empty(len(variables))
    rxn_name = list()

    rxn_ids = [rxn.id for rxn in model_variability.reactions]

    for i in range(len(variables)):
        # if the variable is reactions optimize for forward - reverse variables else optimize for the variable
        if variables[i] in rxn_ids:
            rxn = model_variability.reactions.get_by_id(variables[i])
            for_var = gurobi_interface.getVarByName(rxn.forward_variable.name)
            rev_var = gurobi_interface.getVarByName(rxn.reverse_variable.name)
            obj_exp = for_var - rev_var
        else:
            obj_exp = gurobi_interface.getVarByName(variables[i])

        rxn_name.append(variables[i])

        # minimization
        if len(warm_start) != 0:
            for var in gurobi_interface.getVars():
                if var in warm_start:
                    var.Start = warm_start[var]

        gurobi_interface.setObjective(obj_exp, GRB.MINIMIZE)
        gurobi_interface.update()
        gurobi_interface.optimize()
        objective_value = gurobi_interface.ObjVal
        fluxes_min[i] = objective_value
        # print(rxn.id, "min", objective_value)

        warm_start = {}
        for var in gurobi_interface.getVars():
            if var.VarName.startswith("indicator_"):
                warm_start[var] = var.x

        # maximiztion
        for var in gurobi_interface.getVars():
            if var in warm_start:
                var.Start = warm_start[var]

        gurobi_interface.setObjective(obj_exp, GRB.MAXIMIZE)
        gurobi_interface.update()
        gurobi_interface.optimize()
        objective_value = gurobi_interface.ObjVal
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
    params=False,
):
    """Custom function to perform TVA on MIQC problem using gurobi.

    Parameters
    ----------
    model : multitfa.core.tmodel
        multitfa model after thermodynamic constraints are added
    variable_list : List, optional
        List of variables to perform TVA on, by default None
    params : Bool, optional
        If True sets the Timelimit option to 300 sec and reduced the mip gap to 0.005

    Returns
    -------
    pd.DataFrame
        Dataframe of min max ranges of variables

    Raises
    ------
    ValueError
        [description]
    """
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

    cplex_model.set_log_stream(None)
    cplex_model.set_error_stream(None)
    cplex_model.set_warning_stream(None)
    cplex_model.set_results_stream(None)

    if params:
        # print("lol")
        cplex_model.parameters.mip.tolerances.mipgap = 0.005
        cplex_model.parameters.timelimit = 300
        cplex_model.parameters.mip.limits.probetime = 300

    # Make shorts for sense
    max_sense = cplex_model.objective.sense.maximize
    min_sense = cplex_model.objective.sense.minimize

    if variable_list == None:
        variables = model.cplex_interface.variables.get_names()
    else:
        variables = [var for var in variable_list]

    vars_list_cplex = cplex_model.variables.get_names()

    fluxes_min = np.empty(len(variables))
    fluxes_max = np.empty(len(variables))
    rxn_name = list()

    rxn_ids = [rxn.id for rxn in model.reactions]

    for i in range(len(variables)):
        # Reset objective vector for each iteration
        for varname in vars_list_cplex:
            cplex_model.objective.set_linear(varname, 0)

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

        # maximiztion
        cplex_model.objective.set_sense(max_sense)
        cplex_model.solve()
        objective_value = cplex_model.solution.get_objective_value()
        fluxes_max[i] = objective_value

    return DataFrame(
        {
            "minimum": Series(index=rxn_name, data=fluxes_min),
            "maximum": Series(index=rxn_name, data=fluxes_max),
        }
    )
