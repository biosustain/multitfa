from six import iteritems
from numpy import log, isnan, where, delete, array, linalg, sqrt, diag
from  .thermo_constants import Vmax, RT, K
from scipy.stats import chi2
from .posdef import nearestPD, isPD


""" This is a supplementary script to create all thermodynamic varibales and consraints to add to the model.
"""

def reaction_variables(reaction):
    """Creates the reaction DelG variables (model.problem.Variable), No need to import optlang
    
    Arguments:
        reaction {model.reaction} -- reaction object
    
    Returns:
        [List] -- List of delG and indicator for forward and reverse reaction objects
    """
    if reaction.model is not None:
        delG_forward = reaction.model.problem.Variable(
                                    'dG_{}'.format(
                                    reaction.forward_variable.name),
                                    lb = -1000, ub = 1000)
        delG_reverse = reaction.model.problem.Variable(
                        'dG_{}'.format(reaction.reverse_variable.name),
                         lb = -1000, ub = 1000)
        indicator_forward = reaction.model.problem.Variable(
                            'indicator_{}'.format(reaction.forward_variable.name),
                             lb = 0, ub = 1, type = 'binary')
        indicator_reverse = reaction.model.problem.Variable(
                            'indicator_{}'.format(reaction.reverse_variable.name),
                             lb = 0, ub = 1, type = 'binary')
        return [delG_forward, delG_reverse, indicator_forward, indicator_reverse]
    else:
        return None

def metabolite_variables(metabolite):
    """ Metabolite concentration and Confidence interval variables. Ci is set between 2 S.D (-1.96 to 1.96)
    
    Arguments:
        metabolite {model.metabolite} -- core.compound object
    
    Returns:
        [tuple] -- Tuple of concentration variable and C.I variable
    """
    if metabolite.model is not None:
        conc_variable = metabolite.model.problem.Variable(
                            'lnc_{}'.format(metabolite.id),
                            lb = log(metabolite.concentration_min),
                            ub = log(metabolite.concentration_max)
                            )
        Ci_variable = metabolite.model.problem.Variable(
                            'Ci_{}'.format(metabolite.id),
                            lb = -1.96, ub = 1.96)
        
        met_variable = metabolite.model.problem.Variable(
                            'met_{}'.format(metabolite.id),
                            lb = -1000, ub = 1000)
        return conc_variable, Ci_variable, met_variable
    else:
        return None


def directionality(reaction):
    """ Reaction directionality constraint
    Vi - Vmax * Zi <= 0
    
    Arguments:
        reaction {core.reaction} -- reaction object
    
    Returns:
        [tuple] -- tuple of direactionality constraints
    """
    if reaction.model is not None:
        directionality_constraint_f = reaction.model.problem.Constraint(
                                    reaction.forward_variable - Vmax * reaction.indicator_forward,
                                     ub = 0, name = 'directionality_{}'.format(reaction.forward_variable.name))

        directionality_constraint_r = reaction.model.problem.Constraint(
                                    reaction.reverse_variable - Vmax * reaction.indicator_reverse,
                                     ub = 0, name = 'directionality_{}'.format(reaction.reverse_variable.name))
    
        return directionality_constraint_f, directionality_constraint_r
    else:
        return None

def delG_indicator(reaction):
    """ Indicator constraints to ensure delG < 0 always
    delG -K + K*Zi <= 0 
    Arguments:
        reaction {core.reaction} -- reaction object
    
    Returns:
        [tuple] -- tuple of indicator constraints
    """
    if reaction.model is not None:

        delG_indicator_constraint_f = reaction.model.problem.Constraint(reaction.delG_forward -K + 
                                    K * reaction.indicator_forward, ub = 0,
                                    name = 'ind_{}'.format(reaction.forward_variable.name))

        delG_indicator_constraint_r = reaction.model.problem.Constraint(reaction.delG_reverse -K + 
                                    K * reaction.indicator_reverse, ub = 0,
                                    name = 'ind_{}'.format(reaction.reverse_variable.name))
    
        return delG_indicator_constraint_f, delG_indicator_constraint_r
    else:
        return None

def massbalance_constraint(model):
    """metabolite mass balance constraints, copying from cobra model
    
    Arguments:
        metabolites {core.compound} -- compound object
    
    Returns:
        [List] -- List of mass balance consttaints
    """
    mass_balance = []
    for met in model.metabolites:
        mass_balance.append(met.constraint)
    return mass_balance

def concentration_exp(reaction):
    """ Concentration term for the delG constraint
    S.T @ ln(X) implemented as sum(stoic * met.conc_var for all mets in reaction)
    
    Arguments:
        reaction {core.reaction} -- reaction object
    
    Returns:
        [symbolic exp] -- Concentration expression
    """
    conc_exp = sum(stoic * metabolite.concentration_variable 
                    for metabolite, stoic in iteritems(reaction.metabolites) 
                        if metabolite.Kegg_id not in ['C00080','cpd00067'])
    return conc_exp

def Ci_exp(reaction, z_f_variable):
    """ Confidence interval term for the delG constraint
    S.T @ Cholesky @ Z_f variables 
    
    Arguments:
        reaction {core.reaction} -- reaction object
        z_f_variable {List} -- C.I variables 
    
    Returns:
        [symbolic exp] -- C.I expression
    """

    S_matrix = reaction.S_matrix 
    z_f_exp = (S_matrix.T[0, :] @ reaction.model.cholskey_matrix) .dot(z_f_variable)

    return z_f_exp

def met_exp_qp(reaction):
    
    return sum(stoic * metabolite.compound_variable 
                    for metabolite, stoic in iteritems(reaction.metabolites) 
                        if metabolite.Kegg_id not in ['C00080','cpd00067'])


def stddev_sampling_rhs(reaction, met_sample_dict): # Have to define met_sample_dict, {metid:sampled_value}
    """ Used for sampling. to calculate the rxn standard deviation from the sampled fromation energies
    
    Arguments:
        reaction {core.reaction} -- reaction object
        met_sample_dict {Dict} -- Formation sample of metabolites

    Returns:
        float -- total S.D of reaction sample (S.T @ Sample)
    """
    rxn_delG_stdev = 0
    for metabolite, stoic in iteritems(reaction.metabolites):
        rxn_delG_stdev += stoic * (float(met_sample_dict[metabolite.id]))

    return rxn_delG_stdev

def add_constraints(model, constraint):
    """Check for duplicate constraints in the model and add 
    
    Arguments:
        model {core.model} -- tmodel
        constraint {model.problem.Constraint} -- Constraint to add
    """
    for cons in constraint:
        if (cons in model.constraints) or (cons in model.variables):
            continue
        model.add_cons_vars(cons)
    

def MIQP(model):
    
    met_ids = [met.id for met in model.metabolites]
    

    gurobi_interface = model.solver.problem.copy()
    metid_grbvars_dict = {}
    for var in gurobi_interface.getVars():
        if var.VarName[4:] in met_ids:
            metid_grbvars_dict[var.VarName[4:]] = var
    #print(metid_grbvars_dict)


    #std_dg, cov_dg = model.covariance_matrix()

    problem_indices = []
    problem_mets = []
    for met in model.metabolites:
        if met.delG_f == 0 or list(set(isnan(model.std_dG)))[0] or sqrt(diag(model.cov_dG)[model.metabolites.index(met)]) > 100:
            problem_indices.append(model.metabolites.index(met))
            problem_mets.append(met)

    #problem_indices = sorted(list(where(isnan(std_dg))[0]) + list(where(std_dg == 0)[0]))
    #problem_indices = [model.metabolites.index(met) for met in model.problem_metabolites]
    #problem_mets = [model.metabolites[ind] for ind in problem_indices]
    cov_dg = model.cov_dG
    cov_dg = delete(cov_dg, problem_indices, axis=0)
    cov_dg = delete(cov_dg, problem_indices, axis=1)

    if not isPD(cov_dg):
        nearPD = nearestPD(cov_dg)
        inv_cov = linalg.inv(nearPD)
    else:
        inv_cov = linalg.inv(cov_dg)

    chi_crit_val = chi2.isf(q = 0.05, df = len(cov_dg))
    
    delG_mean = []
    met_var = []
    for met in model.metabolites:
        if met in problem_mets:
            continue
        delG_mean.append(met.delG_f)
        met_var.append(metid_grbvars_dict[met.id])

    delG_mean = array(delG_mean)
    met_var = array(met_var)


    mu_var = array([met_var])

    

    lhs = mu_var @ inv_cov @ mu_var.T    
    cons = lhs[0]
    gurobi_interface.addConstr(cons[0] <= chi_crit_val, "qp_constraint")
    gurobi_interface.update()
    gurobi_interface.write('test.lp')

    #gurobi_interface.optimize()

    #for v in gurobi_interface.getVars():
    #    print('%s %g' % (v.varName, v.x))

    return gurobi_interface