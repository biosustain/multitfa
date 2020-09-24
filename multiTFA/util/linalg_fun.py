import scipy as sp
import numpy as np


def matrix_decomposition(square_matrix):
    try:
        cholesky = sp.linalg.cholesky(square_matrix, lower=True)
    except:
        eig_val, eig_vec = sp.linalg.eigh(square_matrix)
        eig_val[eig_val < 1e-8] = 0
        outlier = len(np.where(eig_val == 0)[0])
        cholesky = eig_vec @ np.sqrt(np.diag(eig_val))
        cholesky = cholesky[:, outlier:]
    return cholesky
