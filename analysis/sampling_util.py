from scipy.stats import chi2
from numpy import sqrt, random, sum, square
from util.constraints import *
from optlang import Constraint

def generate_n_sphere_sample(n_variables):
   
    chi_crit_val = chi2.isf(q = 0.05, df = n_variables)
    # n-sphere sample from N(0,1)
    random_sample = random.normal(loc=0,scale=1.0,size=(n_variables))
    circle_radius = sqrt(sum(square(random_sample)))/sqrt(chi_crit_val)
    normalized_sample = random_sample/circle_radius
    
    return normalized_sample

def generate_ellipsoid_sample(cholesky):
    """sampling on the surface of n-dimensional ellipsoid
    sample on n-ellipsoid  is linear transformation of n-sphere
    N(mu,var) = mu + A @ N(0,1)
    A@A' = var, A is cholesky matrix
    
    Arguments:
        cholesky {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """

    n_dimensions = len(cholesky)
    n_sphere_sample = generate_n_sphere_sample(n_dimensions)

    ellipsoid_sample = cholesky @ n_sphere_sample

    return ellipsoid_sample

def generate_constraints_sample(model):
    """[summary]
    
    Arguments:
        model {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """

    massbalance = massbalance_constraint(model)

    #formation_sample = generate_ellipsoid_sample(model.cholskey_matrix)
    met_ids = [i.id for i in model.metabolites]

    conc_variables = {}; delG_centroids = {}
    constraints_list = []; var_list = [] #met_sample_dict = {}
    flux_variables = []; delG_variables = []; indicator_varibales = []

    for met in model.metabolites:
        #met_ind = met_ids.index(met.id)
        conc_variables[met.id] = met.conc_variable
        #met_sample_dict[met.id] = formation_sample[met_ind]

    for rxn in model.reactions:
        if rxn.id in model.Exclude_reactions:
            continue

        flux_variables.extend([rxn.forward_variable,rxn.reverse_variable])
        delG_variables.extend([rxn.delG_forward,rxn.delG_reverse])
        indicator_varibales.extend([rxn.indicator_forward, rxn.indicator_reverse])

        directionality_constraint_f, directionality_constraint_r = directionality(rxn)
        delG_indicator_f,delG_indicator_r = delG_indicator(rxn)

        conc_exp = concentration_exp(rxn, conc_variables, model.Kegg_map)

        delG_centroid = delG_rhs(rxn)
        delG_centroids[rxn.id] = delG_centroid
        #delG_stddev = stddev_sampling_rhs(rxn, met_sample_dict) 

        lhs_forward = rxn.delG_forward - RT * conc_exp
        lhs_reverse = rxn.delG_reverse + RT * conc_exp
        rhs = delG_centroid

        delG_constraint_f = Constraint(lhs_forward, lb = rhs, ub = rhs,
                             name = 'delG_{}'.format(rxn.forward_variable_name))
        delG_constraint_r = Constraint(lhs_reverse, lb = -rhs, ub = -rhs,
                             name = 'delG_{}'.format(rxn.reverse_variable_name))
        
        constraints_list.extend([directionality_constraint_f, directionality_constraint_r,
                                    delG_indicator_f, delG_indicator_r,
                                    delG_constraint_f, delG_constraint_r])
        constraints_list.extend(massbalance)
        
        var_list = flux_variables + delG_variables + indicator_varibales

    return constraints_list, var_list, delG_centroids  



def update_model(model):
    cons, var_list, _ = generate_constraints_sample(model)
    model.add_cons_vars(cons)
    model.add_cons_vars(var_list)

def compare_dataframes(df1, df2):

    flags = []
    for i in range(len(df1)):
        range1 = df1['maximum'][i] - df1['minimum'][i]
        range2 = df2['maximum'][i] - df2['minimum'][i]

        if range2 > range1*1.05:
            flags.append('Y')
        else:
            flags.append('N')
    return flags    


