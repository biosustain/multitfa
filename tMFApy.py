from cobra import *
from optlang import Model, Variable, Constraint, Objective
from .util import reactiondGconstraint

def tMFA (cobra_model, Kegg_map, Exclude_list, concentration_dict, pH_dict, I_dict, T_dict):
	
	Exchanges = cobra_model.exchanges
	Excluded_rxns = Exchanges + Exclude_list
	
	Thermo_model = Model(name='Thermo model')
	
	massbalance_constraint = []
	for metabolite in cobra_model.metabolites:
		massbalance_constraint.append(metabolite.constraint)
	
	Thermo_model.add(massbalance_constraint)
	
	for rxn in cobra_model.reactions:
		if rxn not in Excluded_rxns:
			c1, c2, c1_G, c2_G, c_G_g = reactiondGconstraint(rxn, Keg_map, concentration_dict, pH_dict, I_dict, T_dict)
			Thermo_model.add([c1, c2, c1_G, c2_G, c_G_g])
			
	return Thermo_model
	
def solve_problem (Thermo_model):
	
	status = Thermo_model.optimize()
	
	if status != 'optimal':
		print("Error in constraints")
		exit()
