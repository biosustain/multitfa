from six import iteritems
from numpy import (
    log,
    isnan,
    where,
    delete,
    array,
    linalg,
    sqrt,
    diag,
    shape,
    count_nonzero,
    square,
    zeros,
    newaxis,
)
from .thermo_constants import Vmax, RT, K
from scipy.stats import chi2
from .posdef import nearestPD, isPD
from .util_func import findcorrelatedmets


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
            "dG_{}".format(reaction.forward_variable.name), lb=-1000, ub=1000
        )
        delG_reverse = reaction.model.problem.Variable(
            "dG_{}".format(reaction.reverse_variable.name), lb=-1000, ub=1000
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
    """ Metabolite concentration and Confidence interval variables. Ci is set between 2 S.D (-1.96 to 1.96)
    
    Arguments:
        metabolite {model.metabolite} -- core.compound object
    
    Returns:
        [tuple] -- Tuple of concentration variable and C.I variable
    """
    if metabolite.model is not None:
        conc_variable = metabolite.model.problem.Variable(
            "lnc_{}".format(metabolite.id),
            lb=log(metabolite.concentration_min),
            ub=log(metabolite.concentration_max),
        )

        if metabolite.delG_f == 0 or isnan(metabolite.delG_f):
            lb_met = -100
            ub_met = 100
        else:
            lb_met = metabolite.delG_f - 1.96 * metabolite.std_dev
            ub_met = metabolite.delG_f + 1.96 * metabolite.std_dev
        met_variable = metabolite.model.problem.Variable(
            "met_{}".format(metabolite.id), lb=lb_met, ub=ub_met
        )
        return conc_variable, met_variable
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
    conc_exp = sum(
        stoic * metabolite.concentration_variable
        for metabolite, stoic in iteritems(reaction.metabolites)
        if metabolite.Kegg_id not in ["C00080", "cpd00067"]
    )
    return conc_exp


def met_exp_qp(reaction):

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
    chi_crit_val = chi2.isf(q=0.05, df=len(covar))

    met_var = [met_var_dict[met.id] for met in mets]
    centroids = [met.delG_f for met in mets]
    lhs_vars = array(met_var) - array(centroids)
    lhs_vars = lhs_vars[:, newaxis]

    pre_lhs = lhs_vars.T @ inv_cov @ lhs_vars
    lhs = pre_lhs[0]

    return lhs[0], chi_crit_val


def MIQP(model):

    if model.solver.__class__.__module__ == "optlang.gurobi_interface":

        solver_interface = model.gurobi_interface

        # Get metabolite variable from gurobi interface
        metid_vars_dict = {}
        for var in solver_interface.getVars():
            if var.VarName.startswith("met_"):
                metid_vars_dict[var.VarName[4:]] = var

    elif model.solver.__class__.__module__ == "optlang.cplex_interface":
        pass

    else:
        raise NotImplementedError(
            "Current solver does not support quadratic constraints, please use Gurobi or Cplex"
        )

    # Problem metabolites, if met.delGf == 0 or cholesky row is zeros then delete them
    delete_met, cov_mets, cov_met_inds = [], [], []

    for met in model.metabolites:
        if met.delG_f == 0 or isnan(met.delG_f):
            delete_met.append(met)
        else:
            cov_met_inds.append(model.metabolites.index(met))
            cov_mets.append(met)

    cov_dg = model.cov_dG
    # Pick indices of non zero non nan metabolites
    cov_dG = cov_dg[:, cov_met_inds]
    cov_dg = cov_dG[cov_met_inds, :]

    lhs, rhs = quad_constraint(
        cov_dg, cov_mets, metid_vars_dict
    )  # Calculate lhs, rhs for quadratic constraints

    # Calculate ellipsoid box bounds and set to variables
    bounds = bounds_ellipsoid(cov_dg)  # Check for posdef cov_dg
    for met in model.metabolites:
        if met in delete_met:
            continue
        metid_vars_dict[met.id].LB = met.delG_f - bounds[cov_mets.index(met)]
        metid_vars_dict[met.id].UB = met.delG_f + bounds[cov_mets.index(met)]

    # solver_interface.addQConstr(lhs <= rhs, "qp_constraint")
    solver_interface.update()

    solver_interface.write("QC_problem.lp")

    return solver_interface


def bounds_ellipsoid(covariance):
    """ Calculates the bounds of formation energy variables from the covariance matrix. It will help the solver to get to the solution space quick, for the quadratic constraint problem. 

    Half length of ellipse axis is sqrt(eigen_val * chisquare_val). We calculate unit eigen vectors and multiply my corresponding half length. reconstruct the matrix and pick the largest value of each row, that should define the upper bound for the metabolite formation energy.


    :param covariance: [description]
    :type covariance: [type]
    :return: [description]
    :rtype: [type]
    """

    # First calculate the half lengths of ellipsoid
    chi2_value = chi2.isf(q=0.05, df=len(covariance))
    eig_val, eig_vec = linalg.eig(covariance)
    half_len = sqrt(chi2_value * eig_val)

    # Calculate unit eigen vectors and UB in various axis for formation energies
    bounds_mat = zeros((len(covariance), len(covariance)))

    for i in range(len(eig_vec)):
        scaling = sqrt(sum(square(eig_vec[:, i])))
        unit_vec = eig_vec[:, i] / scaling
        bounds_mat[:, i] = half_len[i] * unit_vec

    UB = []
    for i in range(len(bounds_mat)):
        UB.append(max(bounds_mat[i, :]))

    return UB
