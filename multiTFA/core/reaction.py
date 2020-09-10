from cobra import Reaction
from cobra.core.dictlist import DictList
from six import iteritems
from numpy import zeros, dot
from .compound import Thermo_met
from ..util.thermo_constants import FARADAY
from copy import deepcopy
from ..util.thermo_constants import K, Vmax, RT, default_T
from copy import copy, deepcopy
import numpy as np


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
    def delG_transport(self):
        try:
            return self._delG_transport
        except AttributeError:
            self._delG_transport = self.calculate_delG_transport()
            return self._delG_transport

    @delG_transport.setter
    def delG_transport(self, value):
        self._delG_transport = value  # Update the consraints

    def calculate_transport_charge(self):

        n_charge = {}
        n_proton = {}
        # print("lol1")
        if len(self.compartments) > 1 and len(self.compartments) < 3:
            for metabolite, stoic in iteritems(self.metabolites):
                if metabolite.compartment in n_charge:
                    n_charge[metabolite.compartment].append(stoic * metabolite.charge)
                else:
                    n_charge[metabolite.compartment] = [stoic * metabolite.charge]

                if metabolite.compartment in n_proton:
                    n_proton[metabolite.compartment].append(
                        stoic * metabolite.elements.get("H", 0)
                    )
                else:
                    n_proton[metabolite.compartment] = [
                        stoic * metabolite.elements.get("H", 0)
                    ]
        return (
            {key: sum(val) for key, val in n_charge.items()},
            {key: sum(val) for key, val in n_proton.items()},
        )

    def calculate_delG_transport(self):
        if len(self.compartments) != 2:
            return 0
        else:
            charge_dict, proton_dict = self.calculate_transport_charge()
            comps = list(charge_dict.keys())

            electro_static_delG = (
                charge_dict[comps[0]]
                * FARADAY
                * self.model.membrane_potential[comps[1]][comps[0]]
                * 1e-3
            )

            proton_potential_adjustment = (
                proton_dict[comps[0]]
                * np.log(10)
                * (
                    self.model.compartment_info["pH"][comps[1]]
                    - self.model.compartment_info["pH"][comps[0]]
                )
            )

            return electro_static_delG + proton_potential_adjustment

    def cal_stoichiometric_matrix(self):
        """ Reaction stoichiometry, two columns to represent forward and reverse
        
        Returns:
            np.ndarray -- S matrix of reaction
        """

        S_matrix = zeros((len(self.model.metabolites)))
        for metabolite, stoic in iteritems(self.metabolites):
            S_matrix[self.model.metabolites.index(metabolite)] = stoic
        return S_matrix

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

