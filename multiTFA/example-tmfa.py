from cobra import io
from core import tmodel
import numpy as np

model = io.load_matlab_model('model/small_ecoli.mat')

Kegg_map = {}
with open('model/ecoli_kegg_map.txt','r') as f:
	for line in f:
		line = line.strip()
		line = line.split("\t")
		Kegg_map[line[0]] = line[1]
		
pH_I_T_dict = {'pH':{'c':7.5,'e':7,'p':7},'I':{'c':0.25,'e':0,'p':0},'T':{'c':298.15,'e':298.15,'p':298.15}}
del_psi_dict = {'c':{'e':150,'p':150},'e':{'c':-150,'p':-150},'p':{'c':-150,'e':150}}
Excl = [i.id for i in model.reactions if 'DM_' in i.id]
conc_dict = {'min':{'atp_c':1e-3,'adp_c':4e-4,'amp_c':2e-4},'max':{'atp_c':1e-2,'adp_c':7e-4,'amp_c':3e-4}}




t_model = tmodel.tmodel(model = model, Kegg_map = Kegg_map, pH_I_T_dict = pH_I_T_dict, del_psi_dict = del_psi_dict, Exclude_list = Excl)
t_model.update()

problems_const = []
while np.isnan(t_model.slim_optimize()):
	t_model.solver.problem.computeIIS()
	for c in t_model.solver.problem.getConstrs():
		if c.IISConstr:
			problems_const.append(c)
			t_model.solver.problem.remove(c)

from analysis import variability

rxn_ids = [i.id for i in t_model.reactions]
delgs = [var.name for var in t_model.solver.variables if 'G_r_' in var.name]
conc_var = [var.name for var in t_model.solver.variables if 'lnc_' in var.name]
vars_analysis = conc_var

ranges_vars = variability.variability(t_model,fraction_of_optim = 0.9,variable_list = vars_analysis)

for index, ele in ranges_vars.iterrows():
	print(index, ele['minimum'], ele['maximum'])
