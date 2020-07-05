import numpy as np
from pandas import DataFrame, Series, option_context


def variability(model, variable_list=None, min_growth=False, fraction_of_optim=0.9):
    """Performs thermodynamic variability analysis on t_model. Determines the minimum and maximum values for the input variables (Flux, Gibbs free energies and metabolite concentrations). if min_growth constraint is applied then growth is maintained at given percentage of optimum. 

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
    model, variable_list=None, min_growth=False, fraction_of_optim=0.9, warm_start={}
):

    if variable_list == None:
        variables = model.gurobi_interface.getVars()
    else:
        variables = [var for var in variable_list]

    model.gurobi_interface.optimize()
    # if warm start not provided start with box solution
    if len(warm_start) == 0:
        model.slim_optimize()
        if model.solver.status == "optimal":
            for var in model.variables:
                if var.name.startswith("indicator_"):
                    warm_start[var.name] = model.solver.primal_values[var.name]

    if min_growth:
        if model.solver.objective.direction == "max":
            fva_old_objective = model.gurobi_interface.addVar(
                "fva_old_objective",
                lb=fraction_of_optim * model.gurobi_interface.ObjVal,
            )
        else:
            fva_old_objective = model.gurobi_interface.addVar(
                "fva_old_objective",
                ub=fraction_of_optim * model.gurobi_interface.ObjVal,
            )

        # Add the minimal growth/production constraint
        old_obj = model.gurobi_interface.getObjective()
        model.gurobi_interface.addConstr(
            old_obj - fva_old_objective,
            lb=0,
            ub=0,
            name="fva_old_objective_constraint",
        )

        model.gurobi_interface.update()

    fluxes_min = np.empty(len(variables))
    fluxes_max = np.empty(len(variables))
    rxn_name = list()

    rxn_ids = [rxn.id for rxn in model.reactions]

    for i in range(len(variables)):
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
