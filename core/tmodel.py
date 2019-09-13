from optlang import Constraint, Objective, Variable#, Model
from .reaction import reaction
import numpy as np
from util.thermo_constants import Vmax, K, RT
from copy import deepcopy
from util.dGf_calculation import calculate_dGf, cholesky_decomposition
from .compound import compound
from numpy import array, dot
from warnings import warn
from six import iteritems
from cobra import Model
from .solution import get_solution


class tmodel(Model):

    def __init__(self,model, Kegg_map = {}, Exclude_list = [], pH_I_T_dict ={}, 
                cholskey_matrix = [],problematic_rxns = [], concentration_dict = {'min':{},'max':{}},
                del_psi_dict = {}, objective_var = None, objective_dir = None, debug = False):

        """[summary]
        
        Keyword Arguments:
            model {[type]} -- [description] (default: {Model()})
            Kegg_map {dict} -- [description] (default: {{}})
            Exclude_list {list} -- [description] (default: {[]})
            pH_I_T_dict {dict} -- [description] (default: {{}})
            cholskey_matrix {list} -- [description] (default: {[]})
        """
        super().__init__()
        #self.model = Model(name = 'Thermo_model')
        self.reactions = model.reactions
        self.Kegg_map = Kegg_map
        self.Exclude_list = Exclude_list 
        self.pH_I_T_dict = pH_I_T_dict
        self.concentration_dict = concentration_dict
        self.del_psi_dict = del_psi_dict
        self.metabolites = [compound(comp = met, Kegg_map = self.Kegg_map,
                            concentration_min = self.concentration_dict['min'][met.id], 
                            concentration_max = self.concentration_dict['max'][met.id])
                            if met.id in self.concentration_dict['min'].keys()
                            else compound(comp = met, Kegg_map = self.Kegg_map)
                            for met in model.metabolites]

        self.reactions = [reaction(rxn = rxn, Kegg_map=self.Kegg_map,
                            pH_I_T_dict= self.pH_I_T_dict, 
                            del_psi_dict=self.del_psi_dict) for rxn in model.reactions]

        self.cholskey_matrix = self.cal_cholskey_matrix()
        self.problematic_rxns = self.check_problematic_rxns()
        self._objective_var = objective_var
        self._objective_dir = objective_dir
        self.Exclude_reactions = list(set(Exclude_list + self.problematic_rxns))
        self._metabolites = [metabolite.id for metabolite in self.metabolites]
        #self.model = Model()
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
	 
        rxn_order = []
        for rxn in core_rxn:
            reaction = self.get_by_id(rxn)
            rxn_order.extend([reaction.forward_variable.name, reaction.reverse_variable.name])
            for metabolite, stoic in iteritems(reaction.metabolites):
                S_matrix[reaction_index,self._metabolites.index(metabolite.id)] = stoic
                S_matrix[reaction_index+1,self._metabolites.index(metabolite.id)] = -stoic
            reaction_index=reaction_index + 2
            
        S = np.transpose(S_matrix)
	
        return rxn_order, S
 
    def _generate_constraints(self):
        
        z_f_variable = []; conc_variables = {}
        for met in self.metabolites:
            z_f_variable.append(met.Ci_variable)
            conc_variables[met.id] = met.conc_variable

        constraints_list = []
        flux_variables = []
        delG_variables = []
        indicator_vraibales = []
        
        for rxn in self.reactions:
            if rxn.id in self.Exclude_reactions:
                continue
            rxn_order, S_matrix = self.calculate_S_matrix()
            rxn_ind = rxn_order.index(rxn.id)
            rxn_delG = rxn.calculate_rxn_delG()

            flux_variables.extend([rxn.forward_variable,rxn.reverse_variable])
            delG_variables.extend([rxn.delG_forward,rxn.delG_reverse])
            indicator_vraibales.extend([rxn.indicator_forward, rxn.indicator_reverse])

            directionality_constraint_f = Constraint(rxn.forward_variable - Vmax * rxn.indicator_forward, ub = 0,
                                             name = 'directionality_{}'.format(rxn.forward_variable_name))
            directionality_constraint_r = Constraint(rxn.reverse_variable - Vmax * rxn.indicator_reverse, ub = 0,
                                             name = 'directionality_{}'.format(rxn.reverse_variable_name))
            delG_indicator_constraint_f = Constraint(rxn.delG_forward -K + K * rxn.indicator_forward, ub = 0,
                                                     name = 'delG_ind_{}'.format(rxn.forward_variable_name))
            delG_indicator_constraint_r = Constraint(rxn.delG_reverse -K + K * rxn.indicator_reverse, ub = 0,
                                                     name = 'delG_ind_{}'.format(rxn.reverse_variable_name))

            rhs = rxn_delG + rxn.transport_delG

            conc_exp = sum(stoic * conc_variables[metabolite.id] for metabolite, stoic in iteritems(rxn.metabolites) if self.Kegg_map[metabolite.id] not in ['C00080','cpd00067'])
            z_f_exp = (S_matrix.T[rxn_ind,:] @ self.cholskey_matrix) .dot(z_f_variable)


            lhs_forward = rxn.delG_forward - RT * conc_exp - z_f_exp                
            lhs_reverse = rxn.delG_reverse + RT * conc_exp + z_f_exp
            
            delG_constraint_f = Constraint(lhs_forward, lb = rhs, ub = rhs, name = 'delG_{}'.format(rxn.forward_variable_name))
            delG_constraint_r = Constraint(lhs_reverse, lb  = -rhs, ub = -rhs, name = 'delG_{}'.format(rxn.reverse_variable_name))

            constraints_list.extend([directionality_constraint_f, delG_indicator_constraint_f, directionality_constraint_r,\
                                     delG_indicator_constraint_r,delG_constraint_f,delG_constraint_r])
            
            model_variables = z_f_variable + flux_variables + delG_variables + indicator_vraibales
        
        return constraints_list, model_variables
        
    
    def massbalance_constraint(self):
        mass_balance = []
        for met in self.metabolites:
            mass_balance.append(met.constraint)
        return mass_balance

    def set_objective(self):

        if len(self.solver.variables) == 0:
            self.update()

        objective_variables = [var for var in self.solver.variables if self._objective_var == var.name]
        if len(objective_variables) == 0:
            raise ValueError('objective {} supplied not present in the model'.format(self._objective_var))
        elif len(objective_variables) > 1:
            forward_rxn_variable = [i for i in objective_variables if 'reverse' not in i.name][0]
            reverse_rxn_variable = [i for i in objective_variables if 'reverse' in i.name][0]
            objective_exp = 1* forward_rxn_variable -1 *reverse_rxn_variable
        elif len(objective_variables) == 1:
            objective_exp = 1 * objective_variables[0]
        else:
            warn('objective value not set')
        
        self.objective = Objective(objective_exp, direction = self._objective_dir)
        

    def update(self):

        mass_balance = self.massbalance_constraint()
        constraints, variables = self._generate_constraints()
        self.add_cons_vars(constraints)
        self.add_cons_vars(mass_balance)
        self.add_cons_vars(variables)
    
    
    def concentration_ratio_constraints(self,ratio_metabolites,ratio_lb, ratio_ub):
        for i in range(len(ratio_metabolites)):
            ratio_met1 = 'lnc_{}'.format(ratio_metabolites[i][0])
            var1 = [var for var in self.solver.variables if var.name == ratio_met1][0]
            ratio_met2 = 'lnc_{}'.format(ratio_metabolites[i][1])
            var2 = [var for var in self.solver.variables if var.name == ratio_met2][0]
            ratio_constraint = Constraint(1 * var1 - 1* var2, lb = ratio_lb[i], ub = ratio_ub[i])
            self.add_cons_vars(ratio_constraint)

    def optimize(self, raise_error = False):
        solution = get_solution(self, raise_error=raise_error)
        return solution
        
    def lp_matrices_matlab(self):
        """[Variable order:
            flux,excluded reactions, binary, delG, concentration, significance]
        """
        _, S = self.calculate_S_matrix()
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
      
