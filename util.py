import numpy as np
from optlang import Model, Variable, Constraint, Objective
from six import iteritems
from component_contribution import component_contribution_trainer, Reaction, Compound, compound_cache, inchi2gv
from math import log, sqrt
from .transport import find_transportMets, transport_constraints

Vmax = 1000
K = 1000*Vmax
ccache = compound_cache.CompoundCache()



def reactiondGconstraint(reaction, Kegg_map, concentration_dict, pH_dict, I_dict, T_dict):
	"""
	

	"""
	# Instantiate group decomposition
	groupdata = inchi2gv.init_groups_data()
	decomposer = inchi2gv.InChIDecomposer(groupdata)
	groupnames = groupdata.GetGroupNames()
	Ng = len(groupnames)
	
	compcont = component_contribution_trainer.ComponentContribution()
	
	
	indicator_forward = Variable('indiactor_{}'.format(reaction.forward_variable.name), lb =0, ub=1, type = 'binary')
	indicator_reverse = Variable('indiactor_{}'.format(reaction.reverse_variable.name), lb=0, ub=1, type='binary')
	G_r_indicator_for = Variable('G_r_{}'.format(reaction.forward_variable.name), lb =-1000, ub=1000);
	G_r_indicator_rev = Variable('G_r_{}'.format(reaction.reverse_variable.name), lb =-1000, ub=1000);
	c1 = Constraint(rxn.forward_variable - Vmax * indicator_for, lb = 0, ub = 0)
	c2 = Constraint(rxn.reverse_variable - Vmax * indicator_rev, lb = 0, ub = 0)
	c1_G = Constraint(G_r_indicator_for + K * indicator_for, lb = 1000, ub = 1000)
	c2_G = Constraint(G_r_indicator_rev + K * indicator_rev, lb = 1000, ub = 1000)

	data = np.load('component_contribution.npz')
	#stdev_formation = np.sqrt(np.diag(compcont_data['V_rc']))
	#stdev_group = np.sqrt(np.diag(compcont_data['V_gc']))
	
	Conc_min = 1e-5
	Conc_max = 1e-2 # Some default bounds for metabolite concentrations
	
	
	group_vars = []
	for i in range(0,Ng):
		dG_g_stdev = sqrt(np.diag(data['MSE_gc']*data['V_gc'])[i])
		lb_g = data['dG0_gc'][i,] - 2*dG_g_stdev
		ub_g = data['dG0_gc'][i,] + 2*dG_g_stdev
		group_variable = Variable('group_{}'.format(i), lb=lb_g, ub=ub_g) # Finish this
		group_vars.append(group_variable)
	#group_vars = [Variable('group_{}'.format(i), lb=0) for i in range(0, Ng)]
	
	conc_var_dict = {}
	formation_stoic = []
	formation_var_list = []
	group_vec_stoic = []
	dGG = 0
	
	for metabolite, stoic in iteritems(reaction._metabolites):
		# Instantiate for pH, I correction
		comp = ccache.get_compound(Kegg_map[metabolite.id])
		if Kegg_map[metabolite.id] in ['C00080', 'KEGG:C00080']:
			continue
		dGG += stoic * comp.transform_p_h_7(pH_dict[metabolite.compartment],I_dict[metabolite.compartment],T_dict[metabolite.compartment])
			
		if metabolite in concentration_dict:
			lb_c = concentration_dict[metabolite.id]['lb']
			ub_c = concentration_dict[metabolite.id]['ub']
		else:
			lb_c = Conc_min
			ub_c = Conc_max
		conc_var = Variable('lnc_{}'.format(metabolite.id), lb = lb_c, ub= ub_c)
		conc_var_list[conc_var] = stoic
		
		transport_mets = find_transportMets(cobra_model, reaction, Kegg_map)
		if len(transport_mets) = 0:
			transport_delG = 0
		else:
			transport_delG = transport_constraints(reaction)
		
		if metabolite is not in transport_mets:
			if Kegg_map[metabolite.id] != 'C00080':
				if Kegg_map[metabolite] in cids: # Add Kegg_map as input to function ###
					i = cids.index(Kegg_map[metabolite]) # for var bounds
					dG_f = data['dG0_cc'][i, 0]
					dG_f_stdev = sqrt(np.diag(data['MSE_rc']*data['V_rc'])[i])
					lb_f = dG_f - 2*dG_f_stdev
					ub_f = dG_f + 2*dG_f_stdev
					
					formation_var = Variable('G_f_{}'.format(metabolite.id), lb =lb_f, ub=ub_f)
					formation_var_list.append(formation_var)
					formation_stoic.append(stoic)
				else:
							
					try:
						group_vec = decomposer.smiles_to_groupvec(comp.smiles)
						g = group_vec.as_array()
						dG0_gc = data['dG0_gc'][0:self.Ng]
						stoic_g = [x * stoic for x in g]
					except inchi2gv.GroupDecompositionError:
						lb_f = -1000
						ub_f = 1000
						formation_var = Variable('G_f_{}'.format(metabolite.id), lb =lb_f, ub=ub_f)
						formation_var_list.append(formation_var)
						formation_stoic.append(stoic)
					
		c_G_g  = Constraint(G_r_indicator_for - sum(formation_stoic[i] * formation_var_list[i] for i in range(0,len(formation_stoic))) - a * b for a,b in zip(stoic_g, group_vars) - sum(R*T*conc_var_dict[key] * key for key in conc_var_dict.keys()) + transport_delG + ddG)			
	
	return c1, c2, c1_G, c2_G, c_G_g
