from cobra import io
from multiTFA.core import tmodel
import numpy as np
import time

model = io.load_matlab_model("multiTFA/model/small_ecoli.mat")

Kegg_map = {}
with open("multiTFA/model/ecoli_kegg_map.txt", "r") as f:
    for line in f:
        line = line.strip()
        line = line.split("\t")
        Kegg_map[line[0]] = line[1]

pH_I_T_dict = {
    "pH": {"c": 7.5, "e": 7, "p": 7},
    "I": {"c": 0.25, "e": 0, "p": 0},
    "T": {"c": 298.15, "e": 298.15, "p": 298.15},
}
del_psi_dict = {
    "c": {"e": 0, "p": 150},
    "e": {"c": 0, "p": 0},
    "p": {"c": -150, "e": 0},
}
Excl = [i.id for i in model.reactions if "DM_" in i.id]
conc_dict = {
    "min": {"atp_c": 1e-3, "adp_c": 4e-4, "amp_c": 2e-4},
    "max": {"atp_c": 1e-2, "adp_c": 7e-4, "amp_c": 3e-4},
}

t_model = tmodel.tmodel(
    model=model,
    Kegg_map=Kegg_map,
    pH_I_dict=pH_I_T_dict,
    del_psi_dict=del_psi_dict,
    Exclude_list=Excl,
    concentration_dict=conc_dict,
)

t_model.slim_optimize()

warm_start = {}
for var in t_model.variables:
    if var.name.startswith("indicator"):
        warm_start[var.name] = t_model.solver.primal_values[var.name]

for var in t_model.gurobi_interface.getVars():
    if var.VarName in warm_start:
        var.Start = warm_start[var.VarName]
t_model.gurobi_interface.params.OutputFlag = 1

t_model.gurobi_interface.optimize()

print(t_model.gurobi_interface.ObjVal)
