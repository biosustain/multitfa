from util.posdef import nearestPD, isPD
from scipy.stats import chi2
from numpy import sin,cos, pi, linalg, sqrt

def generate_samples(model):
    pass

def calculate_axislen(covariance):

    if not isPD(covariance):
        cov_dg = nearestPD(covariance)
    else:
        cov_dg = covariance
    d_f = len(cov_dg)
    eig_value,eig_vec = linalg.eig(cov_dg)
    chisquare_critical = chi2.isf(q = 0.05, df = d_f)

    axis_len = [sqrt(eig*chisquare_critical) for eig in eig_value]

    return eig_value