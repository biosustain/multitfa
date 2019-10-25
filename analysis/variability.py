from optlang import Variable, Constraint, Objective
from numpy import empty
from pandas import DataFrame, Series, option_context


def variability(model, fraction_of_optim = 0.9 , variable_list = None):
        # Please set objective before proceedding to apply the minimal growth/production constraint
        # Check if model is feasible with all the constraints vefore starting TVA
       
        if variable_list == None:
            variables = [var.name for var in model.solver.variables]
        else:
            variables = [var for var in variable_list]         
       
        feasibility = model.slim_optimize()

        if feasibility == 'infeasible':
            raise ValueError('model infeasible with given constraints')

        """if model.solver.objective.direction == 'max':
            fva_old_objective = Variable('fva_old_objective',
                                         lb = fraction_of_optim * model.solver.objective.value)
        else:
            fva_old_objective = Variable('fva_old_objective',
                                         ub = fraction_of_optim * model.solver.objective.value)
        
        # Add the minimal growth/production constraint
        fva_old_obj_constraint = Constraint(
            model.solver.objective.expression - fva_old_objective, lb=0, ub=0,
            name="fva_old_objective_constraint")

        model.add_cons_vars([fva_old_obj_constraint,fva_old_objective])
        """
        fluxes_min = empty(len(variables))
        fluxes_max = empty(len(variables))  
        rxn_name = list()  
        
        rxn_ids = [rxn.id for rxn in model.reactions]

        for i in range(len(variables)):
            # if the variable is reactions optimize for forward - reverse variables else optimize for the variable
            var = [var for var in model.solver.variables if variables[i] == var.name][0]
            if variables[i] in  rxn_ids:
                rxn_var = model.get_by_id(variables[i])
                objective_exp = 1* rxn_var.forward_variable - 1 * rxn_var.reverse_variable
            else:
                objective_exp = 1* var
            #print(objective_exp)
            #if len(var) > 1:
            #    forward_rxn_variable = [i for i in var if 'reverse' not in i.name][0]
            #    reverse_rxn_variable = [i for i in var if 'reverse' in i.name][0]
            #    objective_exp = 1* forward_rxn_variable -1 *reverse_rxn_variable
            rxn_name.append(variables[i])
            model.objective = Objective(objective_exp)
            #minimization
            model.objective.direction = 'min'
            _ = model.slim_optimize()
            objective_value = model.objective.value
            fluxes_min[i] = objective_value
            #print('min:',objective_value)

            # maximiztion
            model.objective.direction = 'max'
            _ = model.slim_optimize()
            objective_value = model.objective.value
            fluxes_max[i] = objective_value
            #print('max:',objective_value)
        return DataFrame({'minimum' : Series(index=rxn_name, data=fluxes_min),
                            'maximum': Series(index=rxn_name, data=fluxes_max)})

