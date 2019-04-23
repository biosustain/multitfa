from matrices import delGindicatorconstraint_lhs, delGindicatorconstraint_rhs
from cobra import *
import math
from optlang import *
from thermo_constants import *
import numpy as np


def massbalanceConstraint (cobra_model):
	"""
	Adds mass balance constraints for the metabolites in the model
	S.V = 0
	
	Inputs
	cobra_model: model - cobra.model
	
	Returns
	massbalance_constraint: list of metabolite mass balance constraints - list
	"""
	
	massbalance_constraint = []
	
	for metabolite in cobra_model.metabolites:
		massbalance_constraint.append(metabolite.constraint)
		
	return massbalance_constraint


def delGconstraints (cobra_model, Exclude_list, concentration_dict, Kegg_map, cccache_data, pH_I_T_dict):
	"""
	Generates the following constraints to the model
	vi - Vmax*zi < 0, 0<zi<1
	delGr -K + K*zi < 0
	DelGr - RT*xT*ln(x) - xT*map_rc*covf*bi - xT*map_gc*G*covg*bi = xT*map_rc*DelGf + xT*map_rc*G*DelGg + DelG_transport + DelG(pH,I,T)
	
	"""
	delGvar = []
	delGindicatorconstraints = []
	fluxcapacityconstraints = []
	
	core_rxn = [rxn.id for rxn in cobra_model.reactions if rxn.id not in Exclude_list]
	
	for reaction in core_rxn:
		rxn = cobra_model.reactions.get_by_id(reaction)
		G_r_indicator_for = Variable('G_r_{}'.format(rxn.forward_variable.name), lb =-1000, ub=1000);
		G_r_indicator_rev = Variable('G_r_{}'.format(rxn.reverse_variable.name), lb =-1000, ub=1000);
		indicator_for = Variable('indiactor_{}'.format(rxn.forward_variable.name), lb =0, ub=1, type = 'binary')
		indicator_rev = Variable('indiactor_{}'.format(rxn.reverse_variable.name), lb=0, ub=1, type='binary')
			
		c1_G = Constraint(G_r_indicator_for -K + K * indicator_for, ub = 0)
		c2_G = Constraint(G_r_indicator_rev -K + K * indicator_rev, ub = 0)
		c1 = Constraint(rxn.forward_variable - Vmax * indicator_for, ub = 0)
		c2 = Constraint(rxn.reverse_variable - Vmax * indicator_rev, ub = 0)
			
		fluxcapacityconstraints += [c1, c2]
		delGindicatorconstraints += [c1_G, c2_G]
		delGvar += [G_r_indicator_for,G_r_indicator_rev]
	
	lhs = delGindicatorconstraint_lhs(cobra_model, Exclude_list, Kegg_map, cccache_data)
	rhs, Ng, point_estimates, dGG_m = delGindicatorconstraint_rhs(cobra_model, Kegg_map, Exclude_list, cccache_data, pH_I_T_dict)
	
	group_vars = [Variable('group_{}'.format(i), lb =-1.96, ub=1.96) for i in range(0, Ng)]
	z_formation = [Variable('z_f_{}'.format(metabolite.id), lb =-1.96, ub=1.96) for metabolite in cobra_model.metabolites]
	z_inf = [Variable('z_inf_{}'.format(metabolite.id), lb =-1.96, ub=1.96) for metabolite in cobra_model.metabolites]
	
	Variables_conc = []
	
	for met in cobra_model.metabolites:
		
		if met.id in concentration_dict.keys():
			conc_var = Variable('lnc_{}'.format(met.id), lb = log(concentration_dict['lb_c'][met.id]), ub= log(concentration_dict['ub_c'][met.id])) # define dict
			Variables_conc.append(conc_var)
			
		else:
			lb = math.log(1e-5)
			ub = math.log(1e-2)
			conc_var = Variable('lnc_{}'.format(met.id), lb = lb, ub= ub) 
			Variables_conc.append(conc_var)
	
	
	delGvariables = np.array(delGvar + Variables_conc + z_formation + group_vars + z_inf)
	
	delG_conc_constraint = np.array([Constraint(row, ub=bound) for row, bound in zip(lhs.dot(delGvariables), rhs)]);
	
	return delG_conc_constraint, delGindicatorconstraints, fluxcapacityconstraints
