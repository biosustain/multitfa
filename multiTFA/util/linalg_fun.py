import scipy as sp
import numpy as np


def matrix_decomposition(square_matrix):
    try:
        cholesky = sp.linalg.cholesky(square_matrix, lower=True)
    except:
        L, D, perm = sp.linalg.ldl(square_matrix)
        D[D < 1e-8] = 0
        independent_variables = np.where(np.diag(D) != 0)[0]
        cholesky = L @ np.sqrt(D)
        cholesky = cholesky[:, independent_variables]
    return cholesky
