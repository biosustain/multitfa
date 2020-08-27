from cobra import Reaction
from cobra.core.dictlist import DictList
from six import iteritems
from numpy import zeros, dot
from .compound import Thermo_met
from ..util.thermo_constants import FARADAY
from copy import deepcopy
from ..util.thermo_constants import K, Vmax, RT, default_T
from copy import copy, deepcopy
from numpy import transpose


class thermo_reaction(Reaction):

    """ Class representation of thermo reaction Object

    Attributes:
        All cobra.Reaction attributes
        transport_metabolites: Metabolites involved in transport
        transport_delG: Gibbs free energy of transport
        delG_transform: transformed Gibbs free energy of reaction
        S_matrix: stoichiometry of reaction (matrix form)
    """

    def __init__(
        self, cobra_rxn, updated_model=None,
    ):

        """
        Arguments:
            cobra_rxn {Cobra reaction object} -- Cobra model reaction object
        
        Keyword Arguments:
            updated_model {tmodel} -- thermomodel object (default: {None})
            Kegg_map {dict} -- Kegg map of metabolites (default: {{}})
            pH_I_T_dict {dict} -- Dictionary of pH, I of compartments (default: {{}})
            del_psi_dict {dict} -- Dictionary of membrane potential of compartments (default: {{}})
            isTrans {bool} -- is it a trsnaport reaction (default: {False})
            rxn_delG {int} -- standard delG of reaction (default: {0})
            transport_delG {int} -- transport delG (default: {0})
            concentration_dict {dict} -- Dictionary of concentration of metabolites (default: {{'min':{},'max':{}}})
            stoichiometric_matrix {list} --  (default: {[]})
        """

        self._model = updated_model
        do_not_copy_by_ref = {"_model", "_metabolites", "_genes"}
        for attr, value in iteritems(cobra_rxn.__dict__):
            if attr not in do_not_copy_by_ref:
                self.__dict__[attr] = copy(value)

        self._metabolites = {}
        for met, stoic in iteritems(cobra_rxn._metabolites):
            new_met = self.model.metabolites.get_by_id(met.id)
            self._metabolites[new_met] = stoic
            new_met._reaction.add(self)
        self._genes = set()
        for gene in cobra_rxn._genes:
            new_gene = self.model.genes.get_by_id(gene.id)
            self._genes.add(new_gene)
            new_gene._reaction.add(self)

    @property
    def delG_forward(self):
        var_name = "dG_{}".format(self.forward_variable.name)
        if self.model is not None:
            if self.id not in self.model.Exclude_reactions:
                return self.model.variables[var_name]
            else:
                return None
        else:
            return None

    @property
    def delG_reverse(self):
        var_name = "dG_{}".format(self.reverse_variable.name)
        if self.model is not None:
            if self.id not in self.model.Exclude_reactions:
                return self.model.variables[var_name]
            else:
                return None
        else:
            return None

    @property
    def indicator_forward(self):
        var_name = "indicator_{}".format(self.forward_variable.name)
        if self.model is not None:
            if self.id not in self.model.Exclude_reactions:
                return self.model.variables[var_name]
            else:
                return None
        else:
            return None

    @property
    def indicator_reverse(self):
        var_name = "indicator_{}".format(self.reverse_variable.name)
        if self.model is not None:
            if self.id not in self.model.Exclude_reactions:
                return self.model.variables[var_name]
            else:
                return None
        else:
            return None

    @property
    def directionality_constraint(self):
        forward_cons = "directionality_{}".format(self.forward_variable.name)
        reverse_cons = "directionality_{}".format(self.reverse_variable.name)

        if self.model is not None:
            if self.id not in self.model.Exclude_reactions:
                return (
                    self.model.constraints[forward_cons],
                    self.model.constraints[reverse_cons],
                )
            else:
                return None
        else:
            return None

    @property
    def indicator_constraint(self):
        forward_cons = "ind_{}".format(self.forward_variable.name)
        reverse_cons = "ind_{}".format(self.reverse_variable.name)

        if self.model is not None:
            if self.id not in self.model.Exclude_reactions:
                return (
                    self.model.constraints[forward_cons],
                    self.model.constraints[reverse_cons],
                )
            else:
                return None
        else:
            return None

    @property
    def delG_constraint(self):
        forward_cons = "delG_{}".format(self.forward_variable.name)
        reverse_cons = "delG_{}".format(self.reverse_variable.name)

        if self.model is not None:
            if self.id not in self.model.Exclude_reactions:
                return (
                    self.model.constraints[forward_cons],
                    self.model.constraints[reverse_cons],
                )
            else:
                return None
        else:
            return None

    @property
    def delG_prime(self):
        try:
            return self._delG_prime
        except AttributeError:
            self._delG_prime = self.calculate_delG_prime()
            return self._delG_prime

    @delG_prime.setter
    def delG_prime(self, value):
        self._delG_prime = value  # Update the value in model constraints

    @property
    def net_charge_transport(self):
        try:
            return self._net_charge_transport
        except AttributeError:
            self._net_charge_transport = self.cal_net_charge()  # Adjust this
            return self._net_charge_transport

    @net_charge_transport.setter
    def net_charge_transport(self, value):
        self._net_charge_transport = value  # Same again

    @property
    def is_Trans(self):
        try:
            return self._is_Trans
        except AttributeError:
            if len(self.transport_metabolites) > 0:
                self._is_Trans = True
            else:
                self._is_Trans = False

            return self._is_Trans

    @property
    def transport_metabolites(self):
        try:
            return self._transport_metabolites
        except AttributeError:
            self._transport_metabolites = self.find_transportMets()
            return self._transport_metabolites

    @property
    def delG_transport(self):
        try:
            return self._delG_transport
        except AttributeError:
            if self.is_Trans:
                net_charge, direction = self.net_charge_transport
                # print(self.id, net_charge)
                ref_comp, target_comp = next(iter(direction.items()))
                # print(ref_comp, target_comp)
                # print(net_charge, direction)
                self._delG_transport = (
                    net_charge
                    * FARADAY
                    * self.model.membrane_potential[ref_comp][target_comp]
                    * 1e-3
                )
            else:
                self._delG_transport = 0

        return self._delG_transport

    @delG_transport.setter
    def delG_transport(self, value):
        self._delG_transport = value  # Update the consraints

    def calculate_transport_charge(self):
        pass

    def cal_stoichiometric_matrix(self):
        """ Reaction stoichiometry, two columns to represent forward and reverse
        
        Returns:
            np.ndarray -- S matrix of reaction
        """

        S_matrix = zeros((1, len(self.model.metabolites)))
        for metabolite, stoic in iteritems(self.metabolites):
            S_matrix[0, self.model.metabolites.index(metabolite)] = stoic
        stoichio = transpose(S_matrix)

        return stoichio

    def get_coefficient(self, metabolite_id):
        _id_to_metabolites = {m.id: m for m in self.metabolites}
        return self.metabolites[_id_to_metabolites[metabolite_id]]

    def calculate_delG_prime(self):

        """ Calculates standard Gibbs free energy of reactions (Needs to be transformed to pH and I to get transformed Gibbs free energy of reaction
        
        Returns:
            float -- standard Gibbs free energy of reaction
        """

        rxn_delG = 0
        for metabolite, stoic in iteritems(self.metabolites):
            rxn_delG += stoic * float(metabolite.delG_f)

        return rxn_delG

    def find_transportMets(self):
        """ Find the pairs of transport metabolites in a reaction. If a metabolite is present in more than one compartment, then the metabolite is considered as transport metabolite.

        Returns:
            Dict: key , value pair of transport metabolites in different compartments

            Example
            atp_c + h_c --> atp_e + h_c
            
            transport_pair = {atp_c : atp_e}
        """

        # First create a Kegg_id:metabolite dict for prodcucts
        Kegg_prod_map = {met.Kegg_id: met for met in self.products}

        # Now look for the products Kegg id's in reactants and make sure they are in different compartments
        transport_pair = {}
        for reactant in self.reactants:
            if reactant.Kegg_id in Kegg_prod_map:
                if reactant.compartment != Kegg_prod_map[reactant.Kegg_id].compartment:
                    transport_pair[reactant] = Kegg_prod_map[reactant.Kegg_id]

        return transport_pair

    def _proton_balance(self):

        """ Function to calculate if proton balance for some weird reactions like ATP synthase where protons are moved across membrane and some of them are consumed in the reaction
        """

        compartments = list(set([met.compartment for met in self.metabolites]))

        comp1, comp2 = ([], [])
        for metabolite in self.metabolites:
            num_H = 0
            if "H" in metabolite.elements:
                num_H = metabolite.elements["H"]

            if metabolite.compartment == compartments[0]:
                coeff = self.get_coefficient(metabolite.id)
                comp1.append(coeff * num_H)

            elif metabolite.compartment == compartments[1]:
                coeff = self.get_coefficient(metabolite.id)
                comp2.append(coeff * num_H)
            else:
                pass
            # assert abs(sum(comp1)) == abs(sum(comp2))

        return abs(sum(comp1))

    def _iscomplex_transport(self):

        is_complex = False
        proton_compounds = [
            met for met in self.metabolites if met.Kegg_id in ["C00080", "cpd00067"]
        ]
        in_coeff = self.get_coefficient(proton_compounds[0].id)
        out_coeff = self.get_coefficient(proton_compounds[1].id)

        if abs(in_coeff) != abs(out_coeff):
            is_complex = True

        return is_complex

    def _simple_proton(self):
        if not self._iscomplex_transport():
            proton_compounds = [
                met for met in self.metabolites if met.Kegg_id in ["C00080", "cpd00067"]
            ]
            in_coeff = self.get_coefficient(proton_compounds[0].id)
            out_coeff = self.get_coefficient(proton_compounds[1].id)
            assert abs(in_coeff) == abs(out_coeff)

        return abs(in_coeff)

    def _proton_transport(self):
        if self._iscomplex_transport():
            h_transport = self._proton_balance()
        else:
            h_transport = self._simple_proton()
        return h_transport

    def cal_net_charge(self):

        """ Calculates the transport Gibbs free energy of reactions
        delG_transport = net_charge * FARADAY * del_psi 

        Returns:
            Float -- transport Gibbs free energy of reaction
        """
        net_charge_transported = 0
        reference_comp, transported_comp = ("", "")

        if self.is_Trans:
            if len(self.transport_metabolites) > 0:
                for metabolite in self.transport_metabolites:
                    reference_comp = metabolite.compartment
                    transported_comp = self.transport_metabolites[
                        metabolite
                    ].compartment

                    break

                net_charge = []

                for transported in self.transport_metabolites:
                    if transported.Kegg_id in ["C00080", "cpd00067"]:
                        coeff = self._proton_transport()
                    else:
                        coeff = self.get_coefficient(transported.id)

                    if transported.compartment == reference_comp:
                        charge_transported = abs(coeff) * transported.charge
                        net_charge.append(charge_transported)
                    else:
                        charge_transported = -abs(coeff) * transported.charge
                        net_charge.append(charge_transported)
                net_charge_transported = sum(net_charge)

        return (net_charge_transported, {reference_comp: transported_comp})

