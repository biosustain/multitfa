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
        self,
        cobra_rxn,
        updated_model=None,
        Kegg_map={},
        pH_I_dict={},
        del_psi_dict={},
        isTrans=False,
        rxn_delG=0,
        transport_delG=0,
        concentration_dict={"min": {}, "max": {}},
        stoichiometric_matrix=[],
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

        self.Kegg_map = Kegg_map
        self.pH_I_dict = pH_I_dict
        self.concentration_dict = concentration_dict
        self.transport_metabolites = self.find_transportMets()
        self.del_psi_dict = del_psi_dict
        self.isTrans = self.isTransport()
        self.transport_delG = self.delG_transport()
        self.delG_transform = self.calculate_rxn_delG() + self.transport_delG
        self.S_matrix = self.cal_stoichiometric_matrix()
        self.transform = self.transform_adjustment()

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

    def calculate_rxn_delG(self):

        """ Calculates standard Gibbs free energy of reactions (Needs to be transformed to pH and I to get transformed Gibbs free energy of reaction
        
        Returns:
            float -- standard Gibbs free energy of reaction
        """

        rxn_delG = 0
        for metabolite, stoic in iteritems(self.metabolites):
            pH = self.pH_I_dict["pH"][metabolite.compartment]
            ionic_strength = self.pH_I_dict["I"][metabolite.compartment]
            transform = metabolite.transform(pH, ionic_strength)
            rxn_delG += stoic * (float(metabolite.delG_f) + transform)

        return rxn_delG

    def transform_adjustment(self):

        transform_adjust = 0
        for metabolite, stoic in iteritems(self.metabolites):
            pH = self.pH_I_dict["pH"][metabolite.compartment]
            ionic_strength = self.pH_I_dict["I"][metabolite.compartment]
            transform = metabolite.transform(pH, ionic_strength)
            transform_adjust += stoic * transform
        return transform_adjust

    def find_transportMets(self):

        """ Find transported species (metabolites) in the reaction

        if a compound is preset in both reactant and products it is considered as transported species
        
        Returns:
            Dict -- 
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

        """ Boolean -- Check if the reaction is transport reaction
        if the reaction is happening in more than one compratments then its transport reaction

        Returns:
            bool -- Boolean of transport or not
        """

        if len(self.compartments) > 1:
            self.isTrans = True
        else:
            self.isTrans = False
        return self.isTrans

    def _proton_balance(self):

        """ Function to calculate if proton balance for some weird reactions like ATP synthase where protons are moved across membrane and some of them are consumed in the reaction
        """

        compartments = [met.compartment for met in self.metabolites]
        compartments = list(set(compartments))
        comp1 = []
        comp2 = []

        for metabolite in self.metabolites:
            if "H" in metabolite.elements:
                num_H = metabolite.elements["H"]
            else:
                num_H = 0

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
            met
            for met in self.metabolites
            if self.Kegg_map[met.id] in ["C00080", "cpd00067"]
        ]
        in_coeff = self.get_coefficient(proton_compounds[0].id)
        out_coeff = self.get_coefficient(proton_compounds[1].id)

        if abs(in_coeff) != abs(out_coeff):
            is_complex = True
        return is_complex

    def _simple_proton(self):
        if not self._iscomplex_transport():
            proton_compounds = [
                met
                for met in self.metabolites
                if self.Kegg_map[met.id] in ["C00080", "cpd00067"]
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

    def delG_transport(self):

        """ Calculates the transport Gibbs free energy of reactions
        delG_transport = net_charge * FARADAY * del_psi 

        Returns:
            Float -- transport Gibbs free energy of reaction
        """

        if self.isTrans == True:
            if len(self.transport_metabolites) > 0:
                # print(self.id)
                reference_comp = list(self.transport_metabolites.keys())[0].compartment
                transported_comp = self.transport_metabolites[
                    list(self.transport_metabolites.keys())[0]
                ].compartment

                net_charge = []

                for transported in self.transport_metabolites:
                    if self.Kegg_map[transported.id] in ["C00080", "cpd00067"]:
                        coeff = self._proton_transport()
                    else:
                        coeff = self.get_coefficient(transported.id)

                    if transported.compartment == reference_comp:
                        # print(type(coeff), type(transported.charge))
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
