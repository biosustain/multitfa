from cobra import io
from core import tmodel
import numpy as np

model = io.load_matlab_model('model/cauto_working_model.mat')
w4 = model.reactions.get_by_id('W4')
w4.subtract_metabolites({model.metabolites.get_by_id("H2O[c]"): -1}) # unbalanced H2O

CO = model.reactions.get_by_id('EX_CO[c]')
CO.lower_bound = -0.03156
CO.upper_bound = -0.02957

H2 = model.reactions.get_by_id('EX_H2[c]')
H2.lower_bound = -0.0129
H2.upper_bound = -0.0131

#CO2 = model.reactions.get_by_id('EX_Co2[c]')
#CO2.lower_bound = 0.01223
#CO2.upper_bound = 0.01291

acetate = model.reactions.get_by_id('EX_acetate[e]')
acetate.lower_bound = 0.0039
acetate.upper_bound = 0.0042

ethanol = model.reactions.get_by_id('DM_Ethanol[c]')
ethanol.lower_bound = 0.00337
ethanol.upper_bound = 0.0036

bdoh = model.reactions.get_by_id('DM_butanediol[c]')
bdoh.lower_bound = 6.86e-5
bdoh.upper_bound = 0.0001
"""
solution = model.optimize()
for index, ele in solution.iterrows():
	print('{}\t{}\t{}'.format(index, ele['minimum'], ele['maximum']))
"""

with open('model/cauto_charge_info.txt','r') as f:
	for line in f:
		line = line.strip()
		line = line.split("\t")
		met = model.metabolites.get_by_id(line[0])
		met.charge = float(line[1])

Kegg_map = {}
with open('model/cauto_kegg.map','r') as f:
	for line in f:
		line = line.strip()
		line = line.split("\t")
		Kegg_map[line[0]] = line[1]


conc_dict = {'min':{},'max':{}}
with open('model/cauto_metconc_highbio.txt','r') as f:
	for line in f:
		line = line.strip()
		line = line.split("\t")
		conc_dict['min'][line[0]] = float(line[1])
		conc_dict['max'][line[0]] = float(line[2])

pH_I_T_dict = {'pH':{'c':6,'e':5},'I':{'c':0.1,'e':0.1},'T':{'c':298.15,'e':298.15}}
del_psi_dict = {'c':{'e':100},'e':{'c':-100}}
Excl = [i.id for i in model.reactions if 'EX_' in i.id or 'DM_' in i.id]



t_model = tmodel.tmodel(model = model, Kegg_map = Kegg_map, pH_I_T_dict = pH_I_T_dict, del_psi_dict = del_psi_dict, Exclude_list = Excl, concentration_dict = conc_dict)
t_model.update()
#aor_delg = [var for var in t_model.solver.variables if 'G_r_E2-AOR' == var.name][0] # just for debug only
#aor_delg.lb = 0
#e1 = [var for var in t_model.solver.variables if 'G_r_E1-AdhE' == var.name][0] # just for debug only
#e1.lb = 0



problems_const = []
while np.isnan(t_model.slim_optimize()):
	t_model.solver.problem.computeIIS()
	for c in t_model.solver.problem.getConstrs():
		if c.IISConstr:
			problems_const.append(c)
			t_model.solver.problem.remove(c)

from analysis import variability

rxn_ids = [i.id for i in t_model.reactions]
delgs = [var.name for var in t_model.solver.variables if 'dG_' in var.name]
conc_vars = [var.name for var in t_model.solver.variables if 'lnc_' in var.name]

vars_analysis = rxn_ids + delgs + conc_vars

ranges_vars = variability.variability(t_model,fraction_of_optim = 0.9,variable_list = vars_analysis)

for index, ele in ranges_vars.iterrows():
	print('{}\t{}\t{}'.format(index, ele['minimum'], ele['maximum']))
