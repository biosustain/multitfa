from optlang import Constraint, Objective, Variable, Model
#from cobra import Model
from .reaction import reaction
import numpy as np
from util.thermo_constants import Vmax, K, RT
from copy import deepcopy
from util.dGf_calculation import calculate_dGf, cholesky_decomposition
from .compound import compound
from numpy import array, dot
from six import iteritems
from collections import OrderedDict


class tmodel(Model):

    def __init__(self,model, Kegg_map = {}, Exclude_list = [], pH_I_T_dict ={}, concentration_dict = {'min':{},'max':{}},
                    del_psi_dict = {}, cholskey_matrix = [],problematic_rxns = [], debug = False ):
        """[summary]
        
        Keyword Arguments:
            model {[type]} -- [description] (default: {Model()})
            Kegg_map {dict} -- [description] (default: {{}})
            Exclude_list {list} -- [description] (default: {[]})
            pH_I_T_dict {dict} -- [description] (default: {{}})
            cholskey_matrix {list} -- [description] (default: {[]})
        """
        #super().__init__()
        
        self.concentration_dict = concentration_dict
        self.del_psi_dict = del_psi_dict
        self.Kegg_map = Kegg_map 
        self.pH_I_T_dict = pH_I_T_dict 
        self.metabolites = [compound(comp = met, Kegg_map = self.Kegg_map,
                            concentration_min = self.concentration_dict['min'][met.id], 
                            concentration_max = self.concentration_dict['max'][met.id])
                            if met.id in self.concentration_dict['min'].keys()
                            else compound(comp = met, Kegg_map = self.Kegg_map)
                            for met in model.metabolites]
        self.cholskey_matrix = self.cal_cholskey_matrix()       
        self.reactions = [reaction(rxn = rxn, Kegg_map=self.Kegg_map,
                            pH_I_T_dict= self.pH_I_T_dict, 
                            del_psi_dict=self.del_psi_dict) for rxn in model.reactions]
        
        
        self.problematic_rxns = self.check_problematic_rxns()
        self.Exclude_reactions = list(set(Exclude_list + self.problematic_rxns))
        self._metabolites = [metabolite.id for metabolite in self.metabolites]
        self.stoichiometric_matrix = self.calculate_S_matrix()
        self._solver = model._solver
        self.model = Model()
        self.debug = debug
        

    def cal_cholskey_matrix(self):
        """[summary]
        
        Arguments:
            rxn {[type]} -- [description]
        
        Returns:
            [type] -- [description]
        """
        mets = self.metabolites
        ids = [met.id for met in mets]
        std_dg,cov_dg = calculate_dGf(ids, self.Kegg_map)
        chol_matrix = cholesky_decomposition(std_dg,cov_dg)
                
        return chol_matrix
    
    def problem_metabolites(self):

        problematic_metabolites = []

        for met in self.metabolites:
            met_index = list(self._metabolites).index(met.id)
            if self.Kegg_map[met.id] in ['C00080','cpd00067']:
                continue
            if np.count_nonzero(self.cholskey_matrix[:,met_index]) == 0:
                problematic_metabolites.append(met)
        return problematic_metabolites
    
    def check_problematic_rxns(self):

        problematic_rxns = []
        for met in self.metabolites:
            met_index = list(self.metabolites).index(met)
            if self.Kegg_map[met.id] in ['C00080','cpd00067']:
                continue
            if np.count_nonzero(self.cholskey_matrix[:,met_index]) == 0:
                problematic_rxns.append(met.reactions)
        if len(problematic_rxns) > 0:
            problematic_rxns = frozenset.union(* problematic_rxns)
            problems = [i.id for i in problematic_rxns]
            return problems
        else:
            return []
    
    def get_by_id(self,rxn_id):
        reaction_t = [rxn for rxn in self.reactions if rxn.id == rxn_id]
        return reaction_t[0]
                
    def calculate_S_matrix(self):
        
        core_rxn = [rxn.id for rxn in self.reactions if rxn.id not in self.Exclude_reactions]
        n_reactions = len(core_rxn)
        n_metabolites =  len(self.metabolites)
        S_matrix = np.zeros((2*n_reactions,n_metabolites))

        reaction_index = 0
	 
	
        for rxn in core_rxn:
            reaction = self.get_by_id(rxn)
            for metabolite, stoic in iteritems(reaction.metabolites):
                S_matrix[reaction_index,self._metabolites.index(metabolite.id)] = stoic
                S_matrix[reaction_index+1,self._metabolites.index(metabolite.id)] = -stoic
            reaction_index=reaction_index + 2
            
	
        S = np.transpose(S_matrix)
	
        return  S

    def reaction_variables(self):
        
        delGr_variables = OrderedDict()
        indicator_variables = OrderedDict()
        flux_variables = OrderedDict()

        core_rxn = [rxn for rxn in self.reactions if rxn.id not in self.Exclude_reactions]
        for rxn in core_rxn:
            delG_var_f = Variable('G_r_{}'.format(rxn.forward_variable.name), lb =-1000, ub=1000)
            delG_var_r = Variable('G_r_{}'.format(rxn.reverse_variable.name), lb =-1000, ub=1000)
            delGr_variables[rxn.forward_variable.name] = delG_var_f
            delGr_variables[rxn.reverse_variable.name] = delG_var_r
            flux_variables[rxn.forward_variable.name] = rxn.forward_variable
            flux_variables[rxn.reverse_variable.name] = rxn.reverse_variable

            indicator_var_f = Variable('indiactor_{}'.format(rxn.forward_variable.name), lb =0, ub=1, type = 'binary')
            indicator_var_r = Variable('indiactor_{}'.format(rxn.reverse_variable.name), lb =0, ub=1, type = 'binary')
            indicator_variables[rxn.forward_variable.name] = indicator_var_f
            indicator_variables[rxn.reverse_variable.name] = indicator_var_r

        return delGr_variables, indicator_variables, flux_variables
    
    def metabolite_variables(self):
        
        concentration_var = OrderedDict()
        z_f_variable = OrderedDict()
        dG_f_prime = []
        
        for met in self.metabolites:
            concentration_var[met.id] = met.conc_variable
            z_f_variable[met.id] = met.Ci_variable
            dG_f_prime.append(met.delG_f + met.transform(self.pH_I_T_dict['pH'][met.compartment],self.pH_I_T_dict['I'][met.compartment],
                                self.pH_I_T_dict['T'][met.compartment])) 
        
        dG_f_prime = np.array(dG_f_prime)

        return concentration_var, z_f_variable, dG_f_prime
    
    def thermo_constraints_lhs(self):
        
        S = self.stoichiometric_matrix
        _m,n = np.shape(S)

        lhs = np.concatenate((np.eye(n),-RT*S.T,-S.T @ self.cholskey_matrix),axis=1)
    
        return lhs

    
    def thermo_constraints_rhs(self):
        
        rhs = []
        core_rxn = [rxn for rxn in self.reactions if rxn.id not in self.Exclude_reactions]
        for rxn in core_rxn:
            rhs_val = rxn.calculate_rxn_delG() + rxn.transport_delG
            rhs.extend([rhs_val,-rhs_val])
        
        return np.array(rhs)
    
    def debug_variables(self):
        
        core_rxn = [rxn for rxn in self.reactions if rxn.id not in self.Exclude_reactions]
        debug_var = []
        for rxn in core_rxn:
            dbg_var_for = Variable('debug_{}'.format(rxn.forward_variable.name), lb =0, ub=1, type = 'binary')
            dbg_var_rev = Variable('debug_{}'.format(rxn.reverse_variable.name), lb =0, ub=1, type = 'binary')
            debug_var.extend([dbg_var_for,dbg_var_rev])

        return debug_var

    def _constraints(self):
        
        """[summary]
        
        Returns:
            [type] -- [description]
        """
        
        delGr_variables, indicator_variables, flux_variables = self.reaction_variables()
        concentration_var, z_f_variable, _ = self.metabolite_variables()
        _,n = np.shape(self.stoichiometric_matrix)
        
        # delG constraint
        lhs_variables = np.array(list(delGr_variables.values()) + list(concentration_var.values()) + list(z_f_variable.values()))
        lhs = self.thermo_constraints_lhs()
        rhs = self.thermo_constraints_rhs()
        delG_constraint = list([Constraint(row, ub=bound) for row, bound in zip(lhs.dot(lhs_variables), rhs)])
        
        # Flux indicator constraint
        flux_indicator_lhs  = np.concatenate((np.eye(n),-Vmax*np.eye(n)),axis=1)
        flux_indicator_variables = np.array(list(flux_variables.values()) + list(indicator_variables.values()))
        flux_indicator_constraint = list([Constraint(row, ub=bound) for row, bound in zip(flux_indicator_lhs.dot(flux_indicator_variables),np.zeros((n,1)))])
        
        # delG indicator constraint
        if self.debug == False:
            delG_indicator_lhs = np.concatenate((np.eye(n), K * np.eye(n)), axis = 1)
            delG_indicator_variables = np.array(list(delGr_variables.values()) + list(indicator_variables.values()))
            delG_indicator_constraint = list([Constraint(row, ub=bound) for row, bound in zip(delG_indicator_lhs.dot(delG_indicator_variables),K * np.ones((n,1)))])

            constraints = delG_constraint + flux_indicator_constraint + delG_indicator_constraint
            model_variables = list(delGr_variables.values()) + list(indicator_variables.values()) + list(flux_variables.values()) + list(concentration_var.values()) + list(z_f_variable.values())

        else:
            debug_variables = self.debug_variables()
            delG_indicator_lhs = np.concatenate((np.eye(n), K * np.eye(n), -K * np.eye(n)), axis = 1)
            delG_indicator_variables = np.array(list(delGr_variables.values()) + list(indicator_variables.values()) + list(debug_variables))
            delG_indicator_constraint = list([Constraint(row, ub=bound) for row, bound in zip(delG_indicator_lhs.dot(delG_indicator_variables),K * np.ones((n,1)))])

            constraints = delG_constraint + flux_indicator_constraint + delG_indicator_constraint
            model_variables = list(delGr_variables.values()) + list(indicator_variables.values()) + list(flux_variables.values()) + list(concentration_var.values()) + list(z_f_variable.values())

            debug_obj = Objective(sum(1* dbg_var for dbg_var in debug_variables),direction = "min")
        
        if self.debug == False:
            return constraints, model_variables
        else:
            return constraints, model_variables, debug_obj
    
    def massbalance_constraint(self):
        mass_balance = []
        for met in self.metabolites:
            mass_balance.append(met.constraint)
        return mass_balance

    """def set_objective(self,model):
        self.model.objective = model.objective
        self.model.objective_direction = model.objective_direction"""

    def update(self):
        if self.debug == True:
            thermo_constraints, _, debug_obj = self._constraints()
            mass_balance = self.massbalance_constraint()
            self.model.add(thermo_constraints)
            self.model.add(mass_balance)
            self.model.objective = debug_obj
        else:        
            thermo_constraints, _ = self._constraints()
            mass_balance = self.massbalance_constraint()
            self.model.add(thermo_constraints)
            self.model.add(mass_balance)
    
    def variability(self):

        ranges = np.zeros((len(self.model.variables),2))

        for i in range(len(self.model.variables)):
            self.model.objective = Objective(self.model.variables[i])
            #minimization
            self.model.objective.direction = 'min'
            _ = self.model.optimize()
            objective_value = self.model.objective.value
            ranges[i,0] = objective_value

            # maximiztion
            self.model.objective.direction = 'max'
            _ = self.model.optimize()
            objective_value = self.model.objective.value
            ranges[i,1] = objective_value
        
        return ranges
    
    
    def lp_matrices_matlab(self):
        """[Variable order:
            flux,excluded reactions, binary, delG, concentration, significance]
        """
        S = self.stoichiometric_matrix
        core_rxn = [rxn for rxn in self.reactions if rxn.id not in self.Exclude_reactions]
        m,n = np.shape(S)
        
        S_Ex = np.zeros((len(self.Exclude_reactions),len(self.metabolites)))
        i=0
        for rxn in self.Exclude_reactions:
            reaction = self.get_by_id(rxn)
            for metabolite, stoic in iteritems(reaction.metabolites):
                S_Ex[i,self._metabolites.index(metabolite.id)] = stoic
            i = i+1
        S_Exclude = S_Ex.T

        _,y = np.shape(S_Exclude)
        mass_matrix = np.concatenate((S,S_Exclude,np.zeros((m,2*m+2*n))),axis=1)
        flux_indicator = np.concatenate((np.eye(n),np.zeros((n,y)),-1000*np.eye(n),np.zeros((n,2*m+n))),axis=1)
        delG_indicator_matrix = np.concatenate((np.zeros((n,n)),np.zeros((n,y)),K*np.eye(n),np.eye(n),np.zeros((n,2*m))),axis=1)
        delG_matrix = np.concatenate((np.zeros((n,n)),np.zeros((n,y)),np.zeros((n,n)),np.eye(n),-RT*S.T,-S.T@self.cholskey_matrix),axis=1)

        lhs_matrix = np.concatenate((mass_matrix,flux_indicator,delG_indicator_matrix,delG_matrix),axis=0)

        delG_rhs = []
        for rxn in core_rxn:
            delG_rhs.extend([rxn.calculate_rxn_delG(),-1*rxn.calculate_rxn_delG()])

        rhs_matrix = []
        
        rhs_matrix.extend([0]*(m+n))
        rhs_matrix.extend([K]*n)
        rhs_matrix.extend(delG_rhs)  
        rhs_matrix = np.array(rhs_matrix)

        core_rxn_ids = [rxn.id for rxn in core_rxn]
        
        binaries = []; delGs = []; flux_var =[]
        for rxn_id in core_rxn_ids:
            flux_var.append(rxn_id)
            binaries.append(rxn_id+'_bi')
            delGs.append('delG_'+rxn_id)
            flux_var.append(rxn_id+'_rev')
            binaries.append(rxn_id+'_rev_bi')
            delGs.append('delG_'+rxn_id+'_rev')
        
        log_var =[]; sig_var = []
        for metid in self._metabolites:
            log_var.append('log_'+metid)
            sig_var.append('z_f_'+metid)
        var_list = flux_var + self.Exclude_reactions + binaries + delGs + log_var + sig_var
        
        core_flux_bounds_lb = []; core_flux_bounds_ub = []
        for rxn in core_rxn:
            core_flux_bounds_lb.extend([rxn.lower_bound,-rxn.upper_bound])
            core_flux_bounds_ub.extend([rxn.upper_bound,-rxn.lower_bound])
        
        exclude_lb = []; exclude_ub = []
        for rxn_id in self.Exclude_reactions:
            rxn = self.get_by_id(rxn_id)
            exclude_lb.append(rxn.lower_bound)
            exclude_ub.append(rxn.upper_bound)
        
        lower_bounds = core_flux_bounds_lb + exclude_lb + [0]*n + [-1000]*n + [1e-5]*m + [-1.96]*m
        upper_bounds = core_flux_bounds_ub + exclude_ub + [1]*n + [1000]*n + [0.01]*m + [1.96]*m
        
        return lhs_matrix, rhs_matrix, var_list,np.array(lower_bounds),np.array(upper_bounds)        