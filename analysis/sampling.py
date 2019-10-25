from .sampling_util import generate_constraints_sample, generate_ellipsoid_sample, stddev_sampling_rhs, compare_dataframes
from analysis.variability import variability
from pandas import DataFrame
from numpy import zeros, isnan

def sampling(model, cutoff =100):
    
    met_ids = [i.id for i in model.metabolites]
    cons, vars_list, centroids = generate_constraints_sample(model)
    model.add_cons_vars(cons)
    model.add_cons_vars(vars_list)
    delG_var_list = [var.name for var in model.solver.variables if 'G_r_' in var.name]

    met_sample_dict = {}
    representative_ranges = DataFrame(zeros((len(delG_var_list),2)), columns=['minimum', 'maximum'])
    n_improvement = 0
    total_samples = 0
    while n_improvement < cutoff:
        total_samples = total_samples + 1
        formation_sample = generate_ellipsoid_sample(model.cholskey_matrix)
        
        for met in model.metabolites:
            met_sample_dict[met.id] = formation_sample[met_ids.index(met.id)]
        for rxn in model.reactions:
            if rxn.id in model.Exclude_reactions:
                continue
            delG_stddev = stddev_sampling_rhs(rxn, met_sample_dict) 
            rhs = centroids[rxn.id] + delG_stddev

            delG_for_name = 'delG_'+str(rxn.forward_variable_name)
            delG_rev_name = 'delG_'+str(rxn.reverse_variable_name)
            # Just to avoid lb > ub error
            model.constraints[delG_for_name].ub = 1000
            model.constraints[delG_for_name].lb = -1000
            model.constraints[delG_rev_name].ub = 1000
            model.constraints[delG_rev_name].lb = -1000            

            model.constraints[delG_for_name].ub = rhs
            model.constraints[delG_for_name].lb = rhs
            model.constraints[delG_rev_name].ub = -rhs
            model.constraints[delG_rev_name].lb = -rhs
        problems_const = []
        while isnan(model.slim_optimize()):
	        model.solver.problem.computeIIS()
	        for c in model.solver.problem.getConstrs():
		        if c.IISConstr:
			        problems_const.append(c)
			        model.solver.problem.remove(c)
        tva_ranges = variability(model, variable_list=delG_var_list)
        flags = compare_dataframes(representative_ranges, tva_ranges)

        if len(set(flags)) > 1:
            n_improvement = 0
            representative_ranges = tva_ranges
        else:
            n_improvement = n_improvement + 1

    return representative_ranges, total_samples, formation_sample