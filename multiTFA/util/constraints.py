from six import iteritems
import numpy as np
from .thermo_constants import Vmax, RT, K
from scipy.stats import chi2
from scipy import linalg
from .posdef import nearestPD, isPD
from .util_func import findcorrelatedmets


""" This is a supplementary script to create all thermodynamic variables and constraints to add to the model.
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
            "dG_{}".format(reaction.forward_variable.name), lb=-1e5, ub=1e5
        )
        delG_reverse = reaction.model.problem.Variable(
            "dG_{}".format(reaction.reverse_variable.name), lb=-1e5, ub=1e5
        )
        indicator_forward = reaction.model.problem.Variable(
            "indicator_{}".format(reaction.forward_variable.name),
            lb=0,
            ub=1,
            type="binary",
        )
        indicator_reverse = reaction.model.problem.Variable(
            "indicator_{}".format(reaction.reverse_variable.name),
            lb=0,
            ub=1,
            type="binary",
        )
        return [delG_forward, delG_reverse, indicator_forward, indicator_reverse]
    else:
        return None


def metabolite_variables(metabolite):
    """ Metabolite concentration and formation energy variables. formation energy is varied between mean +/- 2 S.D for box type constraints or calculated separately otherwise based on eigen values and vectors.
    
    Arguments:
        metabolite {model.metabolite} -- core.thermo_met object
    
    Returns:
        [tuple] -- Tuple of concentration variable and formation energy variable
    """
    if metabolite.model is not None:
        conc_variable = metabolite.model.problem.Variable(
            "lnc_{}".format(metabolite.id),
            lb=np.log(metabolite.concentration_min),
            ub=np.log(metabolite.concentration_max),
        )

        return conc_variable
    else:
        return None


def formation_variable(model):
    formation_variables = []
    if model is not None:
        for met2dbid in model.kegg2met_dict:
            formation_var = model.problem.Variable(
                "met_{}".format(met2dbid), lb=-100, ub=100
            )
            formation_variables.append(formation_var)
        return formation_variables
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
            ub=0,
            name="directionality_{}".format(reaction.forward_variable.name),
        )

        directionality_constraint_r = reaction.model.problem.Constraint(
            reaction.reverse_variable - Vmax * reaction.indicator_reverse,
            ub=0,
            name="directionality_{}".format(reaction.reverse_variable.name),
        )

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

        delG_indicator_constraint_f = reaction.model.problem.Constraint(
            reaction.delG_forward - K + K * reaction.indicator_forward,
            ub=0,
            name="ind_{}".format(reaction.forward_variable.name),
        )

        delG_indicator_constraint_r = reaction.model.problem.Constraint(
            reaction.delG_reverse - K + K * reaction.indicator_reverse,
            ub=0,
            name="ind_{}".format(reaction.reverse_variable.name),
        )

        return delG_indicator_constraint_f, delG_indicator_constraint_r
    else:
        return None


def concentration_exp(reaction):
    """ Concentration term for the delG constraint
    S.T @ ln(X) implemented as sum(stoic * met.conc_var for all mets in reaction)
    
    Arguments:
        reaction {core.reaction} -- reaction object
    
    Returns:
        [symbolic exp] -- Concentration expression
    """
    conc_exp = sum(
        stoic * metabolite.concentration_variable
        for metabolite, stoic in iteritems(reaction.metabolites)
        if metabolite.Kegg_id not in ["C00080", "cpd00067"]
    )
    return conc_exp


def formation_exp(reaction):

    return sum(
        stoic * metabolite.compound_variable
        for metabolite, stoic in iteritems(reaction.metabolites)
        if metabolite.Kegg_id not in ["C00080", "cpd00067"]
    )


def stddev_sampling_rhs(
    reaction, met_sample_dict
):  # Have to define met_sample_dict, {metid:sampled_value}
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


def quad_constraint(covar, mets, met_var_dict):
    """ generates lhs and rhs for qudratic constraints for the metabolites

    Arguments:
        covar {[type]} -- [description]
        mets {[type]} -- [description]
        met_var_dict {[type]} -- [description]

    Returns:
        [type] -- [description]
    """

    # Check if is pos-def
    if not isPD(covar):
        nearPD = nearestPD(covar)
    else:
        nearPD = covar
    inv_cov = linalg.inv(nearPD)

    if not isPD(inv_cov):
        inv_cov_pd = nearestPD(inv_cov)
    else:
        inv_cov_pd = inv_cov

    chi_crit_val = chi2.isf(q=0.05, df=len(covar))  # Chi square
    met_var = [met_var_dict[met.id] for met in mets]
    centroids = [met.delG_f for met in mets]
    lhs_vars = np.array(met_var) - np.array(centroids)
    lhs_vars = lhs_vars[:, np.newaxis]

    pre_lhs = lhs_vars.T @ inv_cov_pd @ lhs_vars
    lhs = pre_lhs[0]

    return lhs[0], chi_crit_val


def bounds_ellipsoid(covariance):
    """ Calculates the bounds of formation energy variables from the covariance matrix. It will help the solver to get to the solution space quick, for the quadratic constraint problem. 

    Half length of ellipse axis is sqrt(eigen_val * chisquare_val). We calculate unit eigen vectors and multiply my corresponding half length. reconstruct the matrix and pick the largest value of each row, that should define the upper bound for the metabolite formation energy.


    :param covariance: [description]
    :type covariance: [type]
    :return: [description]
    :rtype: [type]
    """

    # Check if covariance is positive definite otherwise we can encounter negative eigne values
    if not isPD(covariance):
        covariance_PD = nearestPD(covariance)
    else:
        covariance_PD = covariance

    # First calculate the half lengths of ellipsoid
    chi2_value = chi2.isf(q=0.05, df=len(covariance_PD))
    eig_val, eig_vec = linalg.eigh(covariance_PD)
    half_len = np.sqrt(chi2_value * eig_val)

    # Calculate unit eigen vectors and UB in various axis for formation energies
    bounds_mat = np.zeros((len(covariance), len(covariance)))

    for i in range(len(eig_vec)):
        scaling = np.sqrt(np.sum(np.square(eig_vec[:, i])))
        unit_vec = eig_vec[:, i] / scaling
        bounds_mat[:, i] = half_len[i] * unit_vec

    UB = []
    for i in range(len(bounds_mat)):
        UB.append(max(bounds_mat[i, :]))

    return UB
