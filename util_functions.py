from group_decompose import *
from equilibrator_cache import *
from component_contribution import *
from six import iteritems


def decomposegroups (metabolite, Kegg_map):
	"""
	decompose a compound in groups
	
	Inputs
	metabolite: metabolite object - cobra_model.metabolite
	Kegg_map: map between cobra_model.metabolites and kegg ids
	
	Returns
	groups_vector: list of group composition of compound - list
	"""
	num_groups = len(GroupDecomposer().groups_data.all_groups)
	met_id_cache = 'KEGG:'+Kegg_map[metabolite.id]
	try:
		cpd = ccache.get_compound(met_id_cache)
		mol = Molecule.FromInChI(cpd.inchi)
		g = GroupDecomposer().Decompose(mol)
		groups_vector = g.AsVector().as_array().tolist()
	except:
		groups_vector = [i*0 for i in range(0,num_groups)]
		
	return groups_vector

def dGG (reaction, Kegg_map, pH_I_T_dict):
	"""
	Calculates dGG of a reaction for given pH, I and T
	
	Inputs
	reaction: reaction object - cobra.reaction
	Kegg_map: Kegg ids of metabolites in model - dict
	pH_I_T_dict: pH, I, T dictionary of all compartments in model - nested dict 
	
	Returns:
	sum(dGG) : dGG of a reaction - float
	
	""" 
	dGG = []
	
	for metabolite, stoic in iteritems(reaction._metabolites):
		kegg_id = 'KEGG:'+Kegg_map[metabolite.id]
		
		if kegg_id in ['KEGG:C00080']: # add more exceptions like redox or something
			continue
		else:
			comp = CompoundCache().get_compound(kegg_id)
			try:
				dGG_comp = float(comp.transform(Q_(pH_I_T_dict['pH'][metabolite.compartment]),Q_(str(pH_I_T_dict['I'][metabolite.compartment])+'M'),Q_(str(pH_I_T_dict['T'][metabolite.compartment])+'K')))
				dGG_comp_m = stoic * dGG_comp
				dGG.append(dGG_comp_m)
			except:
				dGG_comp = 0
				dGG.append(dGG_comp)
	
	return float(sum(dGG))
