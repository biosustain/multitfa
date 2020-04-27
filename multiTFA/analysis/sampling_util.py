from scipy.stats import chi2, genextreme
from numpy import sqrt, random, sum, square, diag
from ..util.constraints import *
from optlang import Constraint
from ..util.thermo_constants import RT
from random import uniform

def generate_box_sample(model):
    
    met_stdev_dict = []

    covariance = model.covariance_matrix()
    stdev = sqrt(diag(covariance))
    for met in model.metabolites:
        ind = model.metabolites.index(met)
        rand_sample = uniform(-1.96*stdev[ind], 1.96*stdev[ind])
        met_stdev_dict.append(rand_sample)
        
    return met_stdev_dict

def generate_n_sphere_sample(n_variables):
   
    #chi_crit_val = chi2.isf(q = 0.05, df = n_variables)
    # n-sphere sample from N(0,1)
    random_sample = random.normal(loc=0,scale=1.0,size=(n_variables))
    circle_radius = sqrt(sum(square(random_sample)))
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
    chi_crit_val = chi2.isf(q = 0.05, df = n_dimensions)
    n_sphere_sample = generate_n_sphere_sample(n_dimensions)

    ellipsoid_sample = chi_crit_val * cholesky @ n_sphere_sample

    return ellipsoid_sample

def generate_constraints_sample(model):
    """[summary]
    
    Arguments:
        model {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """

    #massbalance = massbalance_constraint(model)

    #formation_sample = generate_ellipsoid_sample(model.cholskey_matrix)
    #met_ids = [i.id for i in model.metabolites]

    
    constraints_list = []
        

    for rxn in model.reactions:
        if rxn.id in model.Exclude_reactions:
            continue

        directionality_constraint_f, directionality_constraint_r = directionality(rxn)
        delG_indicator_f,delG_indicator_r = delG_indicator(rxn)

        conc_exp = concentration_exp(rxn)
        #delG_centroids[rxn.id] = rxn.delG_transform

        lhs_forward = rxn.delG_forward - RT * conc_exp
        lhs_reverse = rxn.delG_reverse + RT * conc_exp
        rhs = rxn.delG_transform

        delG_constraint_f = Constraint(lhs_forward, lb = rhs, ub = rhs,
                             name = 'delG_{}'.format(rxn.forward_variable.name))
        delG_constraint_r = Constraint(lhs_reverse, lb = -rhs, ub = -rhs,
                             name = 'delG_{}'.format(rxn.reverse_variable.name))
        
        constraints_list.extend([directionality_constraint_f, directionality_constraint_r,
                                    delG_indicator_f, delG_indicator_r,
                                    delG_constraint_f, delG_constraint_r])
        
        
    #var_list = flux_variables + delG_variables + indicator_varibales
    #constraints_list.extend(massbalance)

    return constraints_list  


def prepare_sampling_model(Thermo_model):

    if len(Thermo_model.constraints) <= len(Thermo_model.metabolites):
        Thermo_model.update()

    for const in Thermo_model.constraints:
        if 'delG_' in const.name:
            for var in const.variables:
                if 'Ci' in var.name:
                    const.set_linear_coefficients({var:0})



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

def extreme_value_distribution(data_set):
    """ Fits the Gibbs free energy data to the Generalized extreme value distribution and predicts the extreme value at 95 % CI.
    Uses Scipy genextreme function
    
    Arguments:
        data_set [list] -- The max or min range of Gibbs free energy values
    
    Returns:
        tuple -- min or max value predicted from GEV at 95% confidence
    """
    c, loc, scale = genextreme.fit(data_set)
    min_extreme, max_extreme = genextreme.interval(0.95,c,loc,scale)
    
    return min_extreme, max_extreme


