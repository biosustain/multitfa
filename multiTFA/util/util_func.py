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


def findcorrelatedmets(covariance):
    non_prob_ind = []
    prob_ind = []
    for i in range(len(covariance)):
        if not np.isnan(covariance[:,i]).all():
            non_prob_ind.append(i)
        else:
            prob_ind.append(i)

    reduced_covar = covariance[:, non_prob_ind]
    reduced_covar = reduced_covar[non_prob_ind, :]

    correlation_mat = cov2corr(covariance)

    return correlation_mat
    
