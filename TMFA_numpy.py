import numpy as np
from optlang import Model, Variable, Constraint, Objective
import find_transportRxns 
import sys
from Thermodata import *
try:
    from scipy.sparse import dok_matrix, lil_matrix
except ImportError:
    dok_matrix, lil_matrix = None, None

from six import iteritems
from numpy import array

	
def stoichiometricmarix_core (cobra_model, array_type = 'dense', dtype =None):
	
	""" To calculate stoichiometric matrices separately transport reactions and other reactions in the model. Taken from cobra.util.array
	"""
	transport_rxns = find_transportRxns(cobra_model);
	
	if array_type not in ('DataFrame', 'dense') and not dok_matrix:
		raise ValueError('Sparse matrices require scipy')
	if dtype is None:
		dtype = np.float64
	array_constructor = {'dense': np.zeros, 'dok': dok_matrix, 'lil': lil_matrix,'DataFrame': np.zeros}
	
	reactions_core = [r for r in cobra_model.reactions if r not in transport_rxns];
	
	n_metabolites = len(cobra_model.metabolites)
	n_reactions_core = len(reactions_core)
	n_reactions_transport = len(transport_rxns)
	
	# Define stoichiometric matrices for both trasnport and normal reactions
	array_core = array_constructor[array_type]((2*n_reactions_core, n_metabolites), dtype=dtype) # Stoichiometry for reverse reaction as well
	array_transport = array_constructor[array_type]((n_reactions_transport, n_metabolites), dtype=dtype)
	

	core_reactions=[];
	
	reaction_index_core=0; # Temp variables to add reverse reaction stoichiometry
	transport_rxns_index=0;
	
	for reaction in cobra_model.reactions:
		if reaction in transport_rxns:
			for metabolite, stoic in iteritems(reaction._metabolites):
				array_transport[transport_rxns.index(reaction),cobra_model.metabolites.index(metabolite)] = stoic;				
		else:
			for metabolite, stoic in iteritems(reaction._metabolites):
				array_core[reaction_index_core,cobra_model.metabolites.index(metabolite)] = stoic;
				array_core[reaction_index_core+1,cobra_model.metabolites.index(metabolite)] = -stoic; 
			reaction_index_core=reaction_index_core+2;
			core_reactions.append(reaction)
		#print (reaction.id,reaction_index_core)	
	
	return array_transport, array_core, transport_rxns, core_reactions


def ThermoModel_setupmatrices(cobra_model):
	
	""" Function to add constraints for modified TMFA with component contribution method with consistent estimation of errors
	Sv=0 ---- Done, Use forward and reverse flux variables from cobrapy
	Vmin < Vi < Vmax
	vi-zi*vmax <= 0, {i=1,...,r}
    dG_ri-K+K*zi<0, {i=1,...,r}
    Pr*ST*dG0_f+Pn*ST*G*dG0_g+RT*SUM[sij*ln(xj)]=dG_ri, {i=1,...,r}
    
	"""
	vmax=1000;
	K=100000;
	
	reactions = [r for r in cobra_model.reactions];
	metabolites = [m for m in cobra_model.metabolites];
	
	S_transport, S_core, transport_rxns, core_reactions = stoichiometricmarix_core(cobra_model);
	S_T = np.transpose(S_core);
	STG= np.multiply(S_T,G);
	
	n_met,n_core_rxns=np.shape(S_core);
	n_met,n_transport_rxn=np.shape(S_transport);
	n_groups=len(group_vector);
	
	# variable order: Vi, Zi, G_r, G_f, G_g, ln(x)
	mass_balance = np.concatenate(S_core, S_transport, np.zeros(n_met,2*n_met+n_groups+n_core_rxns));
	indicator = np.concatenate(np.eye(n_core_rxns),np.zeros(n_core_rxns,n_transport_rxn),-vmax*np.eye(n_core_rxns),np.zeros(n_core_rxns),np.zeros(n_core_rxns,n_met),np.zeros(n_core_rxns,n_groups),np.zeros(n_core_rxns,n_met)); 
	G_r_indicator = np.concatenate(np.zeros(n_core_rxns),np.zeros(n_core_rxns,n_transport_rxn),K*np.eye(n_core_rxns),np.eye(n_core_rxns),np.zeros(n_core_rxns,n_met),np.zeros(n_core_rxns,n_groups),np.zeros(n_core_rxns,n_met));
	G_r_concentration = np.concatenate(np.zeros(n_core_rxns),np.zeros(n_core_rxns,n_transport_rxn),np.zeros(n_core_rxns),np.eye(n_core_rxns),-np.multiply(P_r,S_T),-np.multiply(P_n,STG),-RT*S_T);
	
	lhs = np.concatenate(mass_balance,indicator,G_r_indicator,G_r_concentration); # LHS matrix
	
	# Declare variables now ###
	
	mass_var = np.array((rxn.forward_variable, rxn.reverse_variable) for rxn in core_reactions); # Flux variable
	
	# Vraiables for reaction indicator
	temp_ind = [];
	for rxn in core_reactions:
		temp = Variable('indiactor_{}'.format(rxn.forward_variable.name), lb =0, ub=1, type = 'binary');
		temp1 = Variable('indiactor_{}'.format(rxn.reverse_variable.name), lb =0, ub=1, type = 'binary');
		temp_ind.append(temp);
		temp_ind.append(temp1);
		
	indicator_var = array(temp_ind);
	
	# Vraiables for G_r
	temp_G_r = [];
	for rxn in core_reactions:
		temp = Variable('G_r_{}'.format(rxn.forward_variable.name), lb =-1000, ub=1000);
		temp1 = Variable('G_r_{}'.format(rxn.reverse_variable.name), lb =-1000, ub=1000);
		temp_G_r.append(temp);
		temp_G_r.append(temp1);
		
	G_r_var = array(temp_G_r);
	
	G_f_var = np.array([Variable('G_f_{}_{}'.format(met.id,met.compartment), lb=0) for met in cobra_model.metabolites]) # G_f_variable
	
	G_g_var = np.array([Variable('group_{}'.format(i), lb=0) for i in range(1, n_groups)]); # Group variables
	
	ln_conc_var = np.array([Variable('ln_{}_{}'.format(met.id,met.compartment), lb=0) for met in cobra_model.metabolites]) # ln(metabolite concentration) variable
	
	TMFA_variables = np.concatenate(mass_var, indicator_var, G_r_var, G_f_var, G_g_var, ln_conc_var);
	
	# Setup RHS now
	rhs = np.concatenate(np.zeros(len(cobra_model.metabolites),), np.zeros(len(cobra_model.metabolites),), K*np.ones(len(cobra_model.metabolites),), np.zeros(len(core_reactions),));
	
	# Define constraints from matrices
	c = np.array([Constraint(row, ub=bound) for row, bound in zip(lhs.dot(TMFA_variables), rhs)]);
	
	Thermo_model = Model(name='Thermo model');
	Thermo_model.objective = obj;
	Thermo_model.add(c);
	
	return Thermo_model	
	#array_transport, array_core, trasnport, core_reactions = ThermoModel_setupmatrices(cobra_model);
	
			


def setup_problem ():
	""" Setup MIP problem using the matrices generated  
	"""
	
	obj = Objective(w.dot(x), direction='max');
	c = np.array([Constraint(row, ub=bound) for row, bound in zip(A.dot(x), bounds)]
	model = Model(name='Numpy model')
	model.objective = obj
	model.add(c)


	status = model.optimize()
