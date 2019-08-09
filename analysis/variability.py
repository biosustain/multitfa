from optlang import Variable, Constraint, Objective
from numpy import zeros

def variability(model, fraction_of_optim, variable_list = None):
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

        ranges = zeros((len(variables),2))  # Initialize flux ranges object
        for i in range(len(variables)):
            var = [var for var in model.solver.variables if variables[i] in var.name]
            if len(var) > 1:
                forward_rxn_variable = [i for i in var if 'reverse' not in i.name][0]
                reverse_rxn_variable = [i for i in var if 'reverse' in i.name][0]
                objective_exp = 1* forward_rxn_variable -1 *reverse_rxn_variable
            print(variables[i])
            model.objective = Objective(objective_exp)
            #minimization
            model.objective.direction = 'min'
            _ = model.slim_optimize()
            objective_value = model.objective.value
            ranges[i,0] = objective_value
            print('min:',objective_value)

            # maximiztion
            model.objective.direction = 'max'
            _ = model.slim_optimize()
            objective_value = model.objective.value
            ranges[i,1] = objective_value
            print('max:',objective_value)
        return ranges    