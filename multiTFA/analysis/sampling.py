from .sampling_util import generate_constraints_sample, generate_ellipsoid_sample, compare_dataframes, extreme_value_distribution, generate_box_sample, prepare_sampling_model
from ..util.constraints import add_constraints, stddev_sampling_rhs
from .variability import variability
from pandas import DataFrame, Series, concat
from numpy import zeros, isnan, empty,nan
from warnings import warn

def cutoff_sampling(model, cutoff = 100, sample_type = 0):
    """[summary]
    
    Arguments:
        model {[type]} -- [description]
    
    Keyword Arguments:
        cutoff {int} -- [description] (default: {100})
        sample_type {int} -- type of sample 0 for ellipsoid, anything else for box type (default: {0})
    
    Returns:
        [type] -- [description]
    """

    met_ids = [i.id for i in model.metabolites]

    prepare_sampling_model(model)


    delG_var_list = [var.name for var in model.solver.variables
                     if 'dG_' in var.name]

    
    representative_ranges = DataFrame(zeros((len(delG_var_list),2)), columns=['minimum', 'maximum'])
    
    met_sample_dict = {}
    n_improvement = 0
    total_samples = 0

    while n_improvement < cutoff:
        total_samples = total_samples + 1

        if sample_type == 0:
            formation_sample = generate_ellipsoid_sample(model.cholskey_matrix)
        else:
            formation_sample = generate_box_sample(model)
        
        for met in model.metabolites:
            met_sample_dict[met.id] = formation_sample[met_ids.index(met.id)]

        for rxn in model.reactions:
            if rxn.id in model.Exclude_reactions:
                continue
            delG_stddev = stddev_sampling_rhs(rxn, met_sample_dict) 
            rhs = rxn.delG_transform + delG_stddev

            delG_for_name = 'delG_'+str(rxn.forward_variable.name)
            delG_rev_name = 'delG_'+str(rxn.reverse_variable.name)
            
            # Just to avoid lb > ub error
            model.constraints[delG_for_name].ub = 1000
            model.constraints[delG_for_name].lb = -1000
            model.constraints[delG_rev_name].ub = 1000
            model.constraints[delG_rev_name].lb = -1000            

            model.constraints[delG_for_name].ub = rhs
            model.constraints[delG_for_name].lb = rhs
            model.constraints[delG_rev_name].ub = -rhs
            model.constraints[delG_rev_name].lb = -rhs
       

        tva_ranges = variability(model, variable_list=delG_var_list)
        if tva_ranges.empty or tva_ranges.isnull().all()['maximum']:
            total_samples = total_samples -1
            continue
        else:
            flags = compare_dataframes(representative_ranges, tva_ranges)
            Y_count = flags.count('Y')
            print(tva_ranges)
            if Y_count > 0.05 * len(flags):
                n_improvement = 0
                representative_ranges = tva_ranges
            else:
                n_improvement = n_improvement + 1

        lastone = tva_ranges
    return representative_ranges, total_samples


def gev_sampling(model, cutoff = 1000, sample_type = 0):
    """[summary]
    
    Arguments:
        model {[type]} -- [description]
    
    Keyword Arguments:
        cutoff {int} -- [description] (default: {1000})
        sample_type {int} -- [description] (default: {0})
    
    Returns:
        [type] -- [description]
    """

    met_ids = [i.id for i in model.metabolites]

    prepare_sampling_model(model)


    print("Thermodynamic constraints prepared")
    
    delG_var_list = [var.name for var in model.solver.variables
                     if 'dG_' in var.name]

    met_sample_dict = {}


    Whole_ranges = DataFrame({'minimum' : Series(index = delG_var_list, 
                                    data=empty(len(delG_var_list))),                             'maximum': Series(index=delG_var_list, data=empty(len(delG_var_list)))})

    total_samples = 0
    mins = empty(len(delG_var_list))
    maxs = empty(len(delG_var_list))
    
    while total_samples < cutoff:
        total_samples = total_samples + 1
        #print("sample number {}".format(total_samples))

        if sample_type == 0:
            formation_sample = generate_ellipsoid_sample(model.cholskey_matrix)
        else:
            formation_sample = generate_box_sample(model)
        
        
        for met in model.metabolites:
            met_sample_dict[met.id] = formation_sample[met_ids.index(met.id)]

        for rxn in model.reactions:
            if rxn.id in model.Exclude_reactions:
                continue
            delG_stddev = stddev_sampling_rhs(rxn, met_sample_dict) 
            rhs = rxn.delG_transform + delG_stddev
            #print(rxn.id, rhs)

            delG_for_name = 'delG_'+str(rxn.forward_variable.name)
            delG_rev_name = 'delG_'+str(rxn.reverse_variable.name)
            
            # Just to avoid lb > ub error
            model.constraints[delG_for_name].ub = 1000
            model.constraints[delG_for_name].lb = -1000
            model.constraints[delG_rev_name].ub = 1000
            model.constraints[delG_rev_name].lb = -1000            

            model.constraints[delG_for_name].ub = rhs
            model.constraints[delG_for_name].lb = rhs
            model.constraints[delG_rev_name].ub = -rhs
            model.constraints[delG_rev_name].lb = -rhs
       

        if isnan(model.slim_optimize()):
            total_samples = total_samples -1
            continue
            
        print(model.slim_optimize())
        tva_ranges = variability(model, variable_list=delG_var_list)
        print(tva_ranges)
        if tva_ranges.empty or tva_ranges.isnull().all()['maximum']:
            total_samples = total_samples -1
            continue

        #print(tva_ranges)

        Whole_ranges = concat([Whole_ranges, tva_ranges], axis = 1)
    
    for delG in delG_var_list:
        min_delG , max_d = extreme_value_distribution(Whole_ranges.loc[delG,'minimum'])
        min_d , max_delG = extreme_value_distribution(Whole_ranges.loc[delG,'maximum'])
        mins.append(min_delG)
        maxs.append(max_delG)

    representative_ranges = DataFrame({'minimum' : Series(index=delG_var_list, data=mins),
                            'maximum': Series(index=delG_var_list, data=maxs)})

    return representative_ranges

def sampling(model, sample_type = 'box', exit_strat = 'gev', num_samples = 100):
    
    if sample_type == 'box':
        if exit_strat == 'gev':
            delG_range = gev_sampling(model = model, sample_type = 1, cutoff= num_samples)
        if exit_strat == 'cut_off':
            delG_range = cutoff_sampling(model = model, sample_type= 1, cutoff= num_samples)
    elif sample_type == 'ellipse':
        if exit_strat == 'gev':
            delG_range = gev_sampling(model = model, sample_type = 0)
        if exit_strat == 'cut_off':
            delG_range = cutoff_sampling(model = model, sample_type= 0)
    
    return delG_range