from cobra import Reaction
from six import iteritems
from numpy import zeros, dot
from .compound import compound
from util.thermo_constants import FARADAY
from optlang import Variable, Constraint
from copy import deepcopy
from util.thermo_constants import K, Vmax, RT

class reaction():
 
    def __init__(self,rxn, id = None, Kegg_map = {}, pH_I_T_dict ={}, del_psi_dict = {}, compartments = [],
                    isTrans = False, rxn_delG = 0 , transport_delG = 0, forward_variable = None, 
                    reverse_variable = None, delG_forward = None, delG_reverse = None, 
                    indicator_forward = None, indicator_reverse = None, concentration_dict = {'min':{},'max':{}},
                    stoichiometric_matrix = [], cholskey = []):
        
        #for k,v in rxn.__dict__.items():
        #    self.__dict__[k] = deepcopy(v)
        self.id = rxn.id
        self.Kegg_map = Kegg_map
        self.pH_I_T_dict = pH_I_T_dict
        self.concentration_dict = concentration_dict
        self.metabolites = {compound(comp = met, Kegg_map = self.Kegg_map,
                            concentration_min = self.concentration_dict['min'][met.id], 
                            concentration_max = self.concentration_dict['max'][met.id])
                            if met.id in self.concentration_dict['min'].keys() and met.id in self.concentration_dict['max'].keys()
                            else compound(comp = met, Kegg_map = self.Kegg_map) : value
                            for met,value in iteritems(rxn.metabolites)}
        self.reactants = [met for met, value in iteritems(self.metabolites) if value < 0]
        self.products = [met for met, value in iteritems(self.metabolites) if value > 0]
        self.model = rxn.model
        self.forward_variable = rxn.forward_variable
        self.reverse_variable = rxn.reverse_variable
        self.transport_metabolites = self.find_transportMets()
        self.del_psi_dict = del_psi_dict
        self.compartments = list(set(met.compartment for met in self.metabolites))
        self.isTrans = self.isTransport()
        #self.is_complex = self._iscomplex_transport()
        self.transport_delG = self.delG_transport()
        self.forward_variable_name = rxn.forward_variable.name
        self.reverse_variable_name = rxn.reverse_variable.name
        self.forward_variable = rxn.forward_variable
        self.reverse_variable = rxn.reverse_variable
        self.id = rxn.id
        self.delG_forward = Variable('G_r_{}'.format(self.forward_variable_name), lb =-1000, ub=1000)
        self.delG_reverse = Variable('G_r_{}'.format(self.reverse_variable_name), lb =-1000, ub=1000)
        self.indicator_forward = Variable('indiactor_{}'.format(self.forward_variable_name), lb =0, ub=1, type = 'binary')
        self.indicator_reverse = Variable('indiactor_{}'.format(self.reverse_variable_name), lb =0, ub=1, type = 'binary')
        self.debug_forward = Variable('debug_{}'.format(self.forward_variable_name), lb =0, ub=1, type = 'binary')
        self.debug_reverse = Variable('debug_{}'.format(self.reverse_variable_name), lb =0, ub=1, type = 'binary')
        self.reaction = rxn.reaction
        self.lower_bound = rxn.lower_bound
        self.upper_bound = rxn.upper_bound
        self.S_matrix = self.cal_stoichiometric_matrix()
        self.cholskey = cholskey

    def cal_stoichiometric_matrix(self):
        """[summary]
        
        Arguments:
            rxn {[type]} -- [description]
        
        Returns:
            [type] -- [description]
        """
        mets = [i.id for i in self.model.metabolites]
        S_matrix = zeros((2,len(mets)))

        for metabolite, stoic in iteritems(self.metabolites):
            S_matrix[0,mets.index(metabolite.id)] = stoic
            S_matrix[1,mets.index(metabolite.id)] = -stoic
        
        return S_matrix.T
    
    """def _generate_constraints(self):
        
        # vi - vmax * zi <=0
        flux_ind_f_cons = Constraint(self.forward_variable - Vmax * self.indicator_forward, ub = 0)
        flux_ind_r_cons = Constraint(self.reverse_variable - Vmax * self.indicator_reverse, ub = 0)

        # delGr - K + K * zi <= 0
        delG_ind_f_cons = Constraint(self.delG_forward + K * self.indicator_forward, ub = K)
        delG_ind_r_cons = Constraint(self.delG_reverse + K * self.indicator_reverse, ub = K)

        # delGr - RT@S.T@ln(x) - S.T@cholskey@significance_var = delG_transform
        met_conc_var = [met.conc_variable for met in self.metabolites]
        met_zf_var = [met.Ci_variable for met in self.metabolites]

        delG_f_cons = Constraint(self.delG_forward - RT*(self.S_matrix.T .dot(met_conc_var)) - ((self.S_matrix.T @ self.cholskey).dot(met_zf_var)), 
                        ub = self.delG_transport + self.calculate_rxn_delG())
        delG_r_cons = Constraint(self.delG_reverse + RT*(self.S_matrix.T .dot(met_conc_var)) + ((self.S_matrix.T @ self.cholskey).dot(met_zf_var)), 
                        ub = -1 * self.delG_transport - self.calculate_rxn_delG())
        return flux_ind_f_cons, flux_ind_r_cons, delG_ind_f_cons, delG_ind_r_cons, delG_f_cons, delG_r_cons
    """

    def get_coefficient(self, metabolite_id):
        _id_to_metabolites = {m.id: m for m in self.metabolites}
        return self.metabolites[_id_to_metabolites[metabolite_id]]    
    
    def calculate_rxn_delG(self):
        """[summary]
        
        Arguments:
            rxn {[type]} -- [description]
        
        Returns:
            [type] -- [description]
        """
        rxn_delG = 0
        for metabolite, stoic in iteritems(self.metabolites):
            pH = self.pH_I_T_dict['pH'][metabolite.compartment]
            ionic_strength = self.pH_I_T_dict['I'][metabolite.compartment]
            temperature = self.pH_I_T_dict['T'][metabolite.compartment]
            comp = compound(metabolite,Kegg_map = self.Kegg_map)
            transform = comp.transform(pH,ionic_strength,temperature)
            rxn_delG += stoic * (float(comp.delG_f) + transform)
        return rxn_delG
    
    def find_transportMets(self):
        """[summary]
        
        Arguments:
            rxn {[type]} -- [description]
        
        Returns:
            [type] -- [description]
        """
        products_KEGG = {}
        
        for i in self.products:
            products_KEGG[self.Kegg_map[i.id]] = i
        
        transport_mets = {}
        for i in self.reactants:
            if self.Kegg_map[i.id] in products_KEGG.keys():
                transport_mets[i] = products_KEGG[self.Kegg_map[i.id]]
        
        return transport_mets

    def isTransport(self):
        
        if len(self.compartments) > 1:
            self.isTrans = True
        else:
            self.isTrans = False
        return self.isTrans
    
    
    def _proton_balance(self):
        
        compartments = [met.compartment for met in self.metabolites]
        compartments = list(set(compartments))
        comp1 = []
        comp2 = []

        for metabolite in self.metabolites:
            if 'H' in metabolite.elements:
                num_H = metabolite.elements['H']
            else:
                num_H = 0
            
            if metabolite.compartment == self.compartments[0]:
                coeff = self.get_coefficient(metabolite.id)
                comp1.append(coeff * num_H)
                
            elif metabolite.compartment == self.compartments[1]:
                coeff = self.get_coefficient(metabolite.id)
                comp2.append(coeff * num_H)
            else:
                pass
            #assert abs(sum(comp1)) == abs(sum(comp2))
        
        return abs(sum(comp1))

    def _iscomplex_transport(self):
        is_complex = False
        proton_compounds = [met for met in self.metabolites if self.Kegg_map[met.id] in ['C00080','cpd00067']]
        in_coeff = self.get_coefficient(proton_compounds[0].id)
        out_coeff = self.get_coefficient(proton_compounds[1].id)

        if abs(in_coeff) != abs(out_coeff):
            is_complex = True
        return is_complex
    
    def _simple_proton(self):
        if not self._iscomplex_transport():
            proton_compounds = [met for met in self.metabolites if self.Kegg_map[met.id] in ['C00080','cpd00067']]
            in_coeff = self.get_coefficient(proton_compounds[0].id)
            out_coeff = self.get_coefficient(proton_compounds[1].id)
            assert(abs(in_coeff) == abs(out_coeff))
    
        return abs(in_coeff)

    def _proton_transport(self):
        if self._iscomplex_transport():
            h_transport = self._proton_balance()
        else:
            h_transport = self._simple_proton()
        return h_transport
                

    def delG_transport(self):
        """[summary]
        
        Arguments:
            rxn {[type]} -- [description]
        
        Returns:
            [type] -- [description]
        """
        if self.isTrans == True:
            if len(self.transport_metabolites) > 0:
                #print(self.id)
                reference_comp = list(self.transport_metabolites.keys())[0].compartment
                transported_comp = self.transport_metabolites[list(self.transport_metabolites.keys())[0]].compartment
                
                
                net_charge = []

                for transported in self.transport_metabolites:
                    if self.Kegg_map[transported.id] in ['C00080','cpd00067']:
                        coeff = self._proton_transport()
                    else:
                        coeff = self.get_coefficient(transported.id)

                    if transported.compartment == reference_comp:
                        charge_transported = abs(coeff) * transported.charge
                        net_charge.append(charge_transported)
                    else:
                        charge_transported = -abs(coeff) * transported.charge
                        net_charge.append(charge_transported)
                    
                del_psi = self.del_psi_dict[reference_comp][transported_comp]
                del_psi_G = sum(net_charge) * FARADAY * del_psi * 1e-3
                    
                transport_delG = del_psi_G
            else:
                transport_delG = 0    
        else:
            transport_delG = 0
        return transport_delG
