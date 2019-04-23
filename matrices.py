import numpy as np
from equilibrator_cache import *
from group_decompose import *
from thermo_constants import *
from six import iteritems
from util_functions import *

def stoichiometric_matrix(cobra_model, Exclude_list):
	"""
	stoichiometric matrix for the core reactions
	
	Inputs
	cobra_model: self explanatory
	Exclude_list: Reactions you want to exclude from thermodynamic analysis- List
	
	Returns
	reaction_names: Names of reactions in the order of stoichiometric matrix (includes reverse reactions) -List
	S: stoichiometric matrix - numpy.ndarray
	"""
	core_rxn = [rxn.id for rxn in cobra_model.reactions if rxn.id not in Exclude_list]
	
	n_reactions = len(core_rxn)
	n_metabolites =  len(cobra_model.metabolites)
	
	S_matrix = np.zeros((2*n_reactions,n_metabolites))
	
	reaction_index = 0
	reaction_names = [] 
	
	for rxn in core_rxn:
		reaction = cobra_model.reactions.get_by_id(rxn)
		for metabolite, stoic in iteritems(reaction._metabolites):
			S_matrix[reaction_index,cobra_model.metabolites.index(metabolite)] = stoic;
			S_matrix[reaction_index+1,cobra_model.metabolites.index(metabolite)] = -stoic;
		reaction_index=reaction_index + 2;
		reaction_names.append(reaction.forward_variable.name)
		reaction_names.append(reaction.reverse_variable.name)
	
	S = np.transpose(S_matrix)
	
	return reaction_names, S



def creategroupincidencematrix (cobra_model, Kegg_map):
	"""
	Creates group incidence matrix for the compounds in cobra model
	Inputs
	cobra_model: model 
	Kegg_map: map between cobra_model.metabolites and kegg ids
	
	Output
	num_groups: Number of groups - int
	groupincidencematrix: group incidence matrix (m*G) - numpy.ndarray
	"""
	
	num_groups = len(GroupDecomposer().groups_data.all_groups) # define groupdecomposer
	
	groupincidencelist = []
	
	for met in cobra_model.metabolites:
		try:
			groups_vector = decomposegroups(met, Kegg_map)
			groupincidencelist.append(groups_vector)
		except AttributeError:
			groups_vector = [i*0 for i in range(0,num_groups)]
			groupincidencelist.append(groups_vector)
	
	groupincidencematrix = np.array(groupincidencelist)
	
	return num_groups, groupincidencematrix


def mapping_matrices (cobra_model, cc_data, Kegg_map):
	
	"""
	mapping matrices for the reactions in the model to map between rc,gc and inf
	
	Inputs
	cobra_model: Self Explanatory -model
	cc_data: cache file from component contribution - numpy.ndarray
	Kegg_map: Dictionary of metabolites to Kegg names - Dictionary
	
	Output
	map_rc: matrix for rc map (met*met) - numpy.ndarray
	map_gc: matrix for rc (met*met) - numpy.ndarray
	map_inf: matrix for compounds not covered by rc or gc (met*met) - numpy.ndarray
	
	"""
	
	map_rc = np.zeros((len(cobra_model.metabolites),len(cobra_model.metabolites)))
	map_gc = np.zeros((len(cobra_model.metabolites),len(cobra_model.metabolites)))
	map_inf = np.zeros((len(cobra_model.metabolites),len(cobra_model.metabolites)))
	
	for metabolite in cobra_model.metabolites:
		met_id = 'KEGG:'+Kegg_map[metabolite.id]
		if met_id in cc_data['cids']:
			map_rc[cobra_model.metabolites.index(metabolite),cobra_model.metabolites.index(metabolite)] = 1
		else:
			try:
				grp_vector = decomposegroups(metabolite, Kegg_map)
				map_gc[cobra_model.metabolites.index(metabolite),cobra_model.metabolites.index(metabolite)] = 1
			except:
				map_inf[cobra_model.metabolites.index(metabolite),cobra_model.metabolites.index(metabolite)] = 1
			
	return map_rc, map_gc, map_inf

	
def extractcovariances (cobra_model, Kegg_map, cc_data): # add exception for redox pairs, can deal with them separately
	
	"""
	Extracts covarinace matrices of model metabolites from component contribution cache file
	
	Inputs
	cobra_model: model
	Kegg_map: Dictionary of metabolites to Kegg names - Dictionary
	cc_data: cache file from component contribution - numpy.ndarray
	
	Outputs
	cov_rc: covariance matrix extracted from cc_data for model metabolites: numpy.ndarray
	cov_gc: covariance for groups: numpy.ndarray
	cov_inf: covariance for neither: numpy.ndarray
	"""
	cids = cc_data['cids'].tolist()
	num_groups = len(GroupDecomposer().groups_data.all_groups)
	
	#V_rc = cc_data['MSE_rc']*cc_data['inv_SWS']
	V_rc = cc_data['cov_dG0']
	
	cov_rc = np.zeros((len(cobra_model.metabolites),len(cobra_model.metabolites)))
	cov_gc = cc_data['MSE_gc']*cc_data['inv_GSWGS']
	cov_inf = 1e10*np.eye(len(cobra_model.metabolites))
	
	cid_indices = []
	
	for met in cobra_model.metabolites:
		temp_check = 'KEGG:'+Kegg_map[met.id]
		if temp_check in cids:
			cid_indices.append(cids.index(temp_check))
		else:
			cid_indices.append('na')
	
	for i in range(len(cid_indices)):
		for j in range(len(cid_indices)):
			try:
				cov_rc[i,j] = V_rc[cid_indices[i],cid_indices[j]]
				cov_rc[j,i] = V_rc[cid_indices[j],cid_indices[i]]
			except:
				continue
	
	cov_gc = cov_gc[0:num_groups,0:num_groups]
	
	return cov_rc, cov_gc, cov_inf
	

def delGindicatorconstraint_lhs (cobra_model, Exclude_list, Kegg_map, cccache_data):
	
	"""
	Generates lhs matrix for delG constraint
	DelGr - RT*xT*ln(x) - xT*map_rc*covf*bi - xT*map_gc*G*covg*bi
	
	Inputs
	cobra_model: cobra model - cobra_model.model
	Exclude list: Reactions excluded from thermo analysis - list
	Kegg_map: Kegg ids of metabolites - dict
	cccache_data: cc cache data - numpy.ndarray
	
	Returns
	lhs: lhs matrix of delGconstraint - numpy.ndarray
	"""
	
	rxn_name, S = stoichiometric_matrix(cobra_model,Exclude_list)
	
	[m,n] = np.shape(S)
	
	S_T = np.transpose(S)
	
	Ng, G = creategroupincidencematrix(cobra_model, Kegg_map)
	
	cov_rc, cov_gc, cov_inf = extractcovariances(cobra_model, Kegg_map, cccache_data)
	
	map_rc, map_gc, map_inf = mapping_matrices(cobra_model, cccache_data, Kegg_map)
	
	lhs = np.concatenate((np.eye(n),-RT*S_T,-S_T.dot(map_rc).dot(cov_rc),-S_T.dot(map_gc).dot(G).dot(cov_gc),-S_T.dot(map_inf).dot(cov_inf)),axis = 1)
	
	debug_rc = -S_T.dot(map_rc).dot(cov_rc)
	debug_gc =  -S_T.dot(map_gc).dot(G).dot(cov_gc)
	
	return lhs, debug_rc, debug_gc
	
def delGindicatorconstraint_rhs (cobra_model, Kegg_map, Exclude_list, cccache_data, pH_I_T_dict ):

	"""
	Generates rhs matrix for delG constraint
	xT*map_rc*DelGf + xT*map_rc*G*DelGg + DelG_transport + DelG(pH,I,T)
	
	Inputs
	cobra_model: cobra model - cobra_model.model
	Exclude list: Reactions excluded from thermo analysis - list
	Kegg_map: Kegg ids of metabolites - dict
	cccache_data: cc cache data - numpy.ndarray
	pH_I_T_dict: pH, I, T dictionary of all compartments in model - nested dict
	
	Returns
	rhs: rhs matrix of delGconstraint - numpy.ndarray
	"""
	
	rxn_name, S = stoichiometric_matrix(cobra_model,Exclude_list)
	map_rc, map_gc, map_inf = mapping_matrices(cobra_model, cccache_data, Kegg_map)
	S_T = np.transpose(S)
	cids = cccache_data['cids'].tolist()
	dG0_rc = cccache_data['dG0_cc'] # changed from rc to cc
	dG0_gc = cccache_data['dG0_gc']	

	dGG_m = []
	
	for reaction_name in rxn_name:
		if "_reverse_" not in reaction_name:
			rxn = cobra_model.reactions.get_by_id(reaction_name)
			dGG_rxn = dGG(rxn, Kegg_map, pH_I_T_dict)
			dGG_m.append(dGG_rxn)
			dGG_m.append(-dGG_rxn)
	dGG_m = np.array(dGG_m)
	
	

	
	dG0f = []
	for met in cobra_model.metabolites:
		temp_met_id = 'KEGG:'+Kegg_map[met.id]
		if temp_met_id in cids:
			i = cids.index(temp_met_id)
			formation = dG0_rc[i,]
			dG0f.append(formation)
		else:
			formation = 0
			dG0f.append(formation)
	dG0f = np.array(dG0f)
	
	Ng, G = creategroupincidencematrix(cobra_model, Kegg_map)
	dG0g = dG0_gc[0:Ng]
	
	point_estimates = S_T.dot(map_rc).dot(dG0f) + S_T.dot(map_gc).dot(G).dot(dG0g)
	
	rhs = dGG_m + point_estimates #+transport
	
	return rhs, Ng, point_estimates, dGG_m
