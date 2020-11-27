import os

import numpy as np

from ..comp_cache import comp_cache
from .posdef import isPD, nearestPD


data_dir = os.path.normpath(
    os.path.dirname(os.path.abspath(__file__)) + os.sep + os.pardir + os.sep + "Data"
)
cc_data = np.load(data_dir + os.sep + "compcont_cache.npz", allow_pickle=True)


def preparecompoundmatrices(Compounds):
    """Adopts component contribution method and preprocess covariance matrices

    Arguments:
        Compounds {core.compound} -- compound object

    Returns:
        np.ndarrays -- X, G arrays
    """
    training_cids = cc_data["cids"].tolist()
    Ng = len(cc_data["groups"])
    Ng_full_training = len(cc_data["v_g"])

    X = np.zeros((len(training_cids), len(Compounds)))
    G = np.zeros((Ng_full_training, len(Compounds)))

    for i, compound in enumerate(Compounds):
        if compound.Kegg_id in ["C00080", "cpd00067"]:
            continue
        comp = comp_cache(compound.Kegg_id)
        if comp.internal_id == -1:
            continue
        else:
            if comp.internal_id in training_cids:
                comp_index = training_cids.index(comp.internal_id)
                X[comp_index, i] = 1
            else:
                grp_vector = comp.group_vector
                if grp_vector == []:
                    g = np.zeros((Ng,))
                else:
                    g = np.array(grp_vector, dtype=float)
                G[:Ng, i] = g

    return X, G


def calculate_dGf(Compounds):
    """Calculates the standard delG and covarinace matrix for the input compounds

    Arguments:
        Compounds {core.compound} -- list of compound objects
        Kegg_map {Dict} -- kegg identifier dictionary of compound ids

    Returns:
        [type] -- [description]
    """
    X, G = preparecompoundmatrices(Compounds)
    C1 = cc_data["C1"]
    C2 = cc_data["C2"]
    C3 = cc_data["C3"]

    dG0_cc = cc_data["v_r"]
    dG0_gc = cc_data["v_g"]

    standard_dg = X.T @ dG0_cc + G.T @ dG0_gc
    cov_dg = X.T @ C1 @ X + X.T @ C2 @ G + G.T @ C2.T @ X + G.T @ C3 @ G
    return standard_dg, cov_dg


def cholesky_decomposition(dg, cov_dg):
    """Cholesky decomposition of the covarinace matrix. Uses numpy linalg function. First deletes
    the compound rows whose dG =0 or whole row of covariance is zero and calculate cholesky and reinsert the
    deleted rows. While calculating cholesky, if the covariance is not positive definite, we find the nearest
    positive definite matrix using the algorithm defined in the file.

    Arguments:
        dg {np.ndarray} -- delG of the compounds
        cov_dg {np.ndarray} -- covariance matrix

    Returns:
        [np.ndarray] -- cholesky matrix
    """

    exclude_indices = sorted(
        list(np.where(np.isnan(dg))[0]) + list(np.where(dg == 0)[0])
    )
    cov_dg = np.delete(cov_dg, exclude_indices, axis=0)
    cov_dg = np.delete(cov_dg, exclude_indices, axis=1)

    if not isPD(cov_dg):
        posdef_cov_dg = nearestPD(cov_dg)
        cholesky_matrix = np.linalg.cholesky(posdef_cov_dg)
    else:
        cholesky_matrix = np.linalg.cholesky(cov_dg)

    for i in exclude_indices:
        dim1, _dim2 = np.shape(cholesky_matrix)
        zeros_axis0 = np.zeros((dim1,))
        zeros_axis1 = np.zeros((dim1 + 1,))
        cholesky_matrix = np.insert(cholesky_matrix, i, zeros_axis0, axis=0)
        cholesky_matrix = np.insert(cholesky_matrix, i, zeros_axis1, axis=1)

    return cholesky_matrix
