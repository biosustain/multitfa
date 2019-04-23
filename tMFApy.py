from cobra import *
from optlang import Model, Variable, Constraint, Objective
from util_functions import *
from constraints import *
import numpy as np


def tMFA (cobra_model, Kegg_map, Exclude_list, concentration_dict, pH_I_T_dict, cc_data):
	"""
	"""
		
	Thermo_model = Model(name='Thermo model')
	
	massbalance_constraint = massbalanceConstraint(cobra_model)
	
	
	delG_conc_constraint, delGindicatorconstraints, fluxcapacityconstraints = delGconstraints(cobra_model, Exclude_list, concentration_dict, Kegg_map, cc_data, pH_I_T_dict)
	
	
	Thermo_model.add(massbalance_constraint)
	Thermo_model.add(fluxcapacityconstraints)
	Thermo_model.add(delG_conc_constraint)
	Thermo_model.add(delGindicatorconstraints)
				
	return Thermo_model
	
def solve_problem (Thermo_model):
	
	status = Thermo_model.optimize()
	
	if status != 'optimal':
		print("Error in constraints")
		exit()
