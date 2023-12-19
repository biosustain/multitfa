from pathlib import Path

import numpy as np
from component_contribution import CCModelParameters


PROTON_INCHI_KEY = "GPRLSGONYQIRFK-UHFFFAOYSA-N"
DATA_DIR = Path(__file__).parent.parent / "data"


# covar_data = np.load(DATA_DIR / "covariance.npz")
covar_data = np.load(DATA_DIR / "component_data.npz")

covariance = covar_data["covariance"]


params = CCModelParameters.from_zenodo()
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
K = 1e6
default_T = 298.15  # , "K")
default_I = 0.25  # , "M")
default_pH = 7.0
default_pMg = 10
RT = R * default_T
