from six import iteritems
from optlang import Constraint
from util.thermo_constants import Vmax, RT, K

def directionality(reaction):
    directionality_constraint_f = Constraint(reaction.forward_variable 
                                    - Vmax * reaction.indicator_forward, ub = 0,
                                    name = 'directionality_{}'.format(reaction.forward_variable_name))

    directionality_constraint_r = Constraint(reaction.reverse_variable 
                                    - Vmax * reaction.indicator_reverse, ub = 0,
                                    name = 'directionality_{}'.format(reaction.reverse_variable_name))
    
    return directionality_constraint_f, directionality_constraint_r

def delG_indicator(reaction):
    delG_indicator_constraint_f = Constraint(reaction.delG_forward -K + 
                                    K * reaction.indicator_forward, ub = 0,
                                    name = 'ind_{}'.format(reaction.forward_variable_name))

    delG_indicator_constraint_r = Constraint(reaction.delG_reverse -K + 
                                    K * reaction.indicator_reverse, ub = 0,
                                    name = 'ind_{}'.format(reaction.reverse_variable_name))
    
    return delG_indicator_constraint_f, delG_indicator_constraint_r

def massbalance_constraint(model):
    mass_balance = []
    for met in model.metabolites:
        mass_balance.append(met.constraint)
    return mass_balance

def concentration_exp(reaction, conc_variables, Kegg_map):
    conc_exp = sum(stoic * conc_variables[metabolite.id] 
                    for metabolite, stoic in iteritems(reaction.metabolites) 
                        if Kegg_map[metabolite.id] not in ['C00080','cpd00067'])
    return conc_exp

def Ci_exp(reaction, S_matrix, cholskey_matrix, z_f_variable, rxn_ind):
    z_f_exp = (S_matrix.T[rxn_ind,:] @ cholskey_matrix) .dot(z_f_variable)

    return z_f_exp

def delG_rhs(reaction):
    rxn_delG = reaction.calculate_rxn_delG()
    
    return rxn_delG + reaction.transport_delG

def stddev_sampling_rhs(reaction, met_sample_dict): # Have to define met_sample_dict, {metid:sampled_value}
    rxn_delG_stdev = 0
    for metabolite, stoic in iteritems(reaction.metabolites):
        rxn_delG_stdev += stoic * (float(met_sample_dict[metabolite.id]))

    return rxn_delG_stdev

def add_constraints(model, constraint):
    for cons in constraint:
        if (cons in model.constraints) or (cons in model.variables):
            continue
        model.add_cons_vars(cons)
    