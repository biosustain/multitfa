from cobra import io

from multitfa.core import tmodel


model = io.load_matlab_model("e_coli_core.mat")


pH_I_T_dict = {
    "pH": {"c": 7.5, "e": 7, "p": 7},
    "I": {"c": 0.25, "e": 0, "p": 0},
    "T": {"c": 298.15, "e": 298.15, "p": 298.15},
}
del_psi_dict = {
    "c": {"c": 0, "e": 0, "p": 150},
    "e": {"c": 0, "e": 0, "p": 0},
    "p": {"c": -150, "e": 0, "p": 0},
}

import pandas as pd


del_psi = pd.DataFrame.from_dict(data=del_psi_dict)
comp_info = pd.DataFrame.from_dict(data=pH_I_T_dict)

Excl = [
    rxn.id
    for rxn in model.reactions
    if rxn.id.startswith("EX_") or rxn.id.startswith("DM_")
] + ["BIOMASS_Ecoli_core_w_GAM", "O2t", "H2Ot"]


tfa_model = tmodel(
    model, Exclude_list=Excl, compartment_info=comp_info, membrane_potential=del_psi
)
tfa_model.solver = "cplex"
from multitfa.analysis import variability_legacy_cplex


vars_analysis = [rxn.id for rxn in tfa_model.reactions if not rxn.id.startswith("DM_")]

ranges = variability_legacy_cplex(tfa_model, variable_list=vars_analysis)

f = open("ranges_ecoli_cplex.txt", "w")
for index, ele in ranges.iterrows():
    f.write("{}\t{}\t{}\n".format(index, ele["minimum"], ele["maximum"]))
