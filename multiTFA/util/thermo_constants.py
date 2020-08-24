import numpy as np
from component_contribution import CCModelParameters

params = CCModelParameters.from_quilt()
rc_compound_ids = params.train_G.index.tolist()
MSE_rc = params.MSE.at["rc", "MSE"]
MSE_gc = params.MSE.at["gc", "MSE"]
MSE_inf = 1e10
G = params.train_G.values
Nc = params.dimensions.at["Nc", "number"]
Ng = params.dimensions.at["Ng", "number"]

C1 = MSE_rc * params.V_rc + MSE_gc * params.V_gc + MSE_inf * params.V_inf
C2 = MSE_gc * params.P_N_rc @ G @ params.inv_GSWGS + MSE_inf * G @ params.P_N_gc
C3 = MSE_gc * params.inv_GSWGS + MSE_inf * params.P_N_gc

C = np.block([[C1, C2], [C2.T, C3]])
C = C[: (Nc + Ng), : (Nc + Ng)]
mu = np.hstack([params.dG0_cc, params.dG0_gc[:Ng]])

R = 8.31e-3  # "kJ / mol / K"
FARADAY = 96.485  # "kJ / mol"
Vmax = 1000
K = 10 * Vmax
default_T = 298.15  # , "K")
default_I = 0.25  # , "M")
default_pH = 7.0
default_pMg = 10
RT = R * default_T
