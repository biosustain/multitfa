import numpy as np
from scipy import linalg


def matrix_decomposition(square_matrix):

    L, D, perm = linalg.ldl(square_matrix)
    D[D < 1e-8] = 0
    independent_variables = np.where(np.diag(D) != 0)[0]
    cholesky = L @ np.sqrt(D)
    cholesky = cholesky[:, independent_variables]

    return cholesky
