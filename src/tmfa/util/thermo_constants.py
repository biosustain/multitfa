import numpy as np
import scipy as sp
import os
from component_contribution import CCModelParameters
import scipy.stats as stat

PROTON_INCHI_KEY = "GPRLSGONYQIRFK-UHFFFAOYSA-N"

data_dir = os.path.normpath(
    (os.path.dirname(os.path.abspath(__file__))) + os.sep + os.pardir + os.sep + "Data"
)

# covar_data = np.load(data_dir + os.sep + "covariance.npz")
covar_data = np.load(data_dir + os.sep + "component_data.npz")

covariance = covar_data["covariance"]
# cholesky_small = covar_data["cholesky_low"]
# cholesky_high = covar_data["cholesky_big"]
# chi2_value_small = stat.chi2.isf(q=0.05, df=cholesky_small.shape[1])
# chi2_value_high = stat.chi2.isf(q=0.05, df=cholesky_high.shape[1])

params = CCModelParameters.from_quilt()
rc_compound_ids = params.train_G.index.tolist()
MSE_rc = params.MSE.at["rc", "MSE"]
MSE_gc = params.MSE.at["gc", "MSE"]
MSE_inf = 1e10
G = params.train_G.values
Nc = params.dimensions.at["Nc", "number"]
Ng = params.dimensions.at["Ng", "number"]
mu = np.hstack([params.dG0_cc, params.dG0_gc[:Ng]])


R = 8.31e-3  # "kJ / mol / K"
FARADAY = 96.485  # "kJ / mol"
Vmax = 1000
K = 1e8
default_T = 298.15  # , "K")
default_I = 0.25  # , "M")
default_pH = 7.0
default_pMg = 10
RT = R * default_T
