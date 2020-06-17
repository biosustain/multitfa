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
    pH_I_T_dict=pH_I_T_dict,
    del_psi_dict=del_psi_dict,
    Exclude_list=Excl,
    concentration_dict=conc_dict,
)

from multiTFA.util import constraints, util_func

solver_interface = constraints.MIQP(t_model)
# Generate warm start with expanded box bounds
solver_interface.optimize()
start_dict = {}
for var in solver_interface.getVars():
    if var.VarName.startswith("indicator_"):
        start_dict[var] = var.x
print("box optimization complete")

# Generate metid_vars dict
metid_vars_dict = {}
for var in solver_interface.getVars():
    if var.VarName.startswith("met_"):
        metid_vars_dict[var.VarName[4:]] = var

old_ellipse_mets, new_ellipse_mets, old_cov, new_cov = util_func.findcorrelatedmets(
    t_model.cov_dG, t_model.metabolites
)

old_lhs, old_rhs = constraints.quad_constraint(
    old_cov, old_ellipse_mets, metid_vars_dict
)
print("ellipse_constraints calculated")

for var in solver_interface.getVars():
    if var in start_dict:
        var.Start = start_dict[var]

# Now add old_ellipse_constraint
solver_interface.addQConstr(old_lhs <= old_rhs, "old_ellipse")
solver_interface.update()
solver_interface.write("codel.lp")
print("optimization started")
start = time.time()
solver_interface.optimize()
end = time.time()
print(end - start)
print(solver_interface.status)
print(solver_interface.ObjVal)

from multiTFA.analysis import variability_legacy

ranges = variability_legacy(t_model, solver_interface, 0.9)
print(ranges)
"""

from analysis import variability

rxn_ids = [i.id for i in t_model.reactions]
delgs = [var.name for var in t_model.solver.variables if "G_r_" in var.name]
conc_var = [var.name for var in t_model.solver.variables if "lnc_" in var.name]
vars_analysis = conc_var

ranges_vars = variability.variability(
    t_model, fraction_of_optim=0.9, variable_list=vars_analysis
)

for index, ele in ranges_vars.iterrows():
    print(index, ele["minimum"], ele["maximum"])
"""
