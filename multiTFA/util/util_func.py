import numpy as np
from .posdef import isPD, nearestPD


def cov2corr(covariance):
    """ Calculates correlation matrix from covariance matrix
    corr(i,j) = cov(i,j)/stdev(i) * stdev(j) 

    Arguments:
        covariance {np.ndarray} -- covariance matrix

    Returns:
        [np.ndarray] -- correlation matrix
    """

    stdev = np.sqrt(np.diag(covariance))
    outer_stdev = np.outer(stdev, stdev)
    correlation = covariance / outer_stdev
    correlation[covariance == 0] = 0

    return correlation


def findcorrelatedmets(covariance, metabolites):
    """[summary]

    Arguments:
        covariance {[type]} -- [description]
        metabolites {[type]} -- [description]

    Returns:
        [type] -- [description]
    """

    correlation_mat = cov2corr(covariance)

    # Check for nan in correlation matrix and keep track of them and the corresponding metabolites
    non_prob_ind, prob_ind, non_prob_mets, nan_mets = [], [], [], []
    for i in range(len(correlation_mat)):
        if not np.isnan(correlation_mat[:, i]).all():
            non_prob_ind.append(i)
            non_prob_mets.append(metabolites[i])
        else:
            prob_ind.append(i)
            nan_mets.append(metabolites[i])

    reduced_correlation = correlation_mat[:, non_prob_ind]
    reduced_correlation = reduced_correlation[non_prob_ind, :]

    reduced_cov = covariance[:, non_prob_ind]
    reduced_cov = reduced_cov[non_prob_ind, :]

    # removing metabolites with all zeros in correlation matrix
    zero_mets, final_mets, non_zero_ind = [], [], []
    for i in range(0, len(reduced_correlation)):
        if np.count_nonzero(reduced_correlation[:, i]) == 0:
            zero_mets.append(non_prob_mets[i])
        else:
            non_zero_ind.append(i)
            final_mets.append(non_prob_mets[i])

    final_correlation = reduced_correlation[:, non_zero_ind]
    final_correlation = final_correlation[non_zero_ind, :]

    final_cov = reduced_cov[:, non_zero_ind]
    final_cov = final_cov[non_zero_ind, :]

    # Find high variance metabolites
    ind_high_variances = list(np.where(np.sqrt(np.diag(final_cov)) > 10)[0])

    # Find indices that are highly correlated (corr > 0.7 | corr < -0.7) and check if they are same as high variance metabolites
    (
        new_ellipsoid_ind,
        old_ellipsoid_ind,
        no_ellipse_mets,
        new_ellipse_mets,
        old_ellipse_mets,
    ) = ([], [], [], [], [])
    for i in range(len(final_correlation)):
        if i in ind_high_variances:
            pos_corr = list(set(np.where(final_correlation[:, i] > 0.7)[0]))
            neg_corr = list(set(np.where(final_correlation[:, i] < -0.7)[0]))
            correlated_ind = pos_corr + neg_corr
            if len(correlated_ind) == 0:
                no_ellipse_mets.append(final_mets[i])
            if set(correlated_ind).intersection(set(ind_high_variances)) == set(
                correlated_ind
            ):
                new_ellipsoid_ind.append(i)
                new_ellipse_mets.append(final_mets[i])
            else:
                old_ellipsoid_ind.append(i)
                old_ellipse_mets.append(final_mets[i])
        else:
            old_ellipsoid_ind.append(i)
            old_ellipse_mets.append(final_mets[i])

    return (
        old_ellipse_mets,
        new_ellipse_mets,
    )
