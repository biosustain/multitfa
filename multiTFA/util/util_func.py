from functools import wraps 
""" Adopted from https://medium.com/@mgarod/dynamically-add-a-method-to-a-class-in-python-c49204b85bd6

Defining a decorator that accepts the method for a class

"""
def add_method(cls):
    def decorator(func):
        @wraps(func) 
        def wrapper(self, *args, **kwargs): 
            return func(*args, **kwargs)
        setattr(cls, func.__name__, wrapper)
        # Note we are not binding func, but wrapper which accepts self but does exactly the same as func
        return func # returning func means func can still be used normally
    return decorator

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


def findcorrelatedmets(model):
    
    covariance = model.cov_dG
    correlation_mat = cov2corr(covariance)

    # Check for nan in correlation matrix and keep track of them and the corresponding metabolites
    non_prob_ind = []
    prob_ind = []
    non_prob_mets = []
    nan_mets = []
    for i in range(len(correlation_mat)):
        if not np.isnan(correlation_mat[:,i]).all():
            non_prob_ind.append(i)
            non_prob_mets.append(model.metabolites[i].id)
        else:
            prob_ind.append(i)
            nan_mets.append(model.metabolites[i].id)

    reduced_correlation = correlation_mat[:, non_prob_ind]
    reduced_correlation = reduced_correlation[non_prob_ind, :]

    reduced_cov = covariance[:, non_prob_ind]
    reduced_cov = reduced_cov[non_prob_ind, :]

    #removing metabolites with all zeros in correlation matrix
    zero_mets = []
    final_mets = []
    non_zero_ind = []

    for i in range(0,len(reduced_correlation)):
        if np.count_nonzero(reduced_correlation[:,i]) == 0:
            zero_mets.append(non_prob_mets[i])
        else:
            non_zero_ind.append(i)
            final_mets.append(non_prob_mets[i])

    final_correlation = reduced_correlation[:, non_zero_ind]
    final_correlation = final_correlation[non_zero_ind, :]

    final_cov = reduced_cov[:, non_zero_ind]
    final_cov = final_cov[non_zero_ind, :]
    
    # Find high variance metabolites
    ind_high_variances = list(np.where(np.sqrt(np.diag(final_cov)) > 10 )[0]) 
    high_var_mets = [final_mets[i] for i in ind_high_variances]

    return final_correlation, final_mets, high_var_mets, ind_high_variances
    
