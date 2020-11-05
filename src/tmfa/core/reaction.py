from copy import copy

import numpy as np
from cobra import Reaction
from numpy import zeros
from six import iteritems

from ..util.thermo_constants import FARADAY


class thermo_reaction(Reaction):
    """
    Class representation of thermo reaction Object. We calculate the required thermodynamic constraints for performing tMFA. To do the constraints, we need Gibbs energy of reaction and transport.

    Parameters
    ----------
    cobra_rxn : cobra.core.Reaction
        Cobra reaction object, to copy the attributes from. We copy metabolites and genes.
    updated_model : core.tmodel, optional
        tmodel object, with updated thermo properties, by default None

    """

    def __init__(
        self,
        cobra_rxn,
        updated_model=None,
    ):
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
        """An optlang variable representing the Gibbs energy of the forward reaction.

        Returns
        -------
        optlang.interface.variable
            An optlang variable for the Gibbs energy of the forward half reaction or None if the reaction is not associated with the model or reaction is in the list of thermodynamic excluded reactions.
        """
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
        """An optlang variable representing the Gibbs energy of the reverse reaction.

        Returns
        -------
        optlang.interface.variable
            An optlang variable for the Gibbs energy of the reverse half reaction or None if the reaction is not associated with the model or reaction is in the list of thermodynamic excluded reactions.
        """
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
        """An optlang binary variable to dictate the flux through forward half reaction based on Gibbs energy.

        Returns
        -------
        optlang.interface.variable
            An optlang binary variable of forward half reaction or None if the reaction is not associated with the model or eaction is in the list of thermodynamic excluded reactions.
        """
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
        """An optlang binary variable to dictate the flux through reverse half reaction based on Gibbs energy.

        Returns
        -------
        optlang.interface.variable
            An optlang binary variable of reverse half reaction or None if the reaction is not associated with the model or eaction is in the list of thermodynamic excluded reactions.
        """
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
        """optlang constraints to determine reaction directionlity based on the binary variable. The reaction is split is forward and reverse variables. Either forward or reverse variables can carry positive flux only if the corresponding binary variable is '1' which will be determined by the Gibbs energy. If binary variable is '0', corresponding flux variable will be '0'.

        vi -Vmax * Zi <=0

        Returns
        -------
        tuple
            tuple of optlang.interface.constraint of forward and reverse reaction directionalities. if the reaction is not associated with model or reaction present in list of thermodynamic excluded reactions, then None
        """
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
        """optlang constraint to ensure that reaction Gibbs energy is always < 0.

        delG -K + K*Zi <= 0

        Returns
        -------
        tuple
            tuple of optlang.interface.constraint of forward and reverse reaction indicator constraints. if the reaction is not associated with model or reaction present in list of thermodynamic excluded reactions, then None
        """
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
        """optlang constraints to calculate the Gibbs energy of reactions. Forward and reverse reaction variables are represented separately.

        delGr - RT * S.T * ln(x) - S.T * delGferr - S.T @ delGf0 - delGtransport = 0

        Gibbs energy of a reaction is a linear combination of Gibbs energy of formation of metabolites at the compartment conditions and metabolite concentrations. There's an uncertainity associated with formation energy estimations. This uncertainity is modelled as error variable and is allowed to vary 2 standard devaitions from mean. Since the error is associated with formation energies they cancel out common error in a reaction.

        Returns
        -------
        tuple
            tuple of optlang.interface.constraint of forward and reverse reaction Gibbs energy constraints. if the reaction is not associated with model or reaction present in list of thermodynamic excluded reactions, then None
        """
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
        """Transformed Gibbs energy of a reaction. Standard Gibbs energy adjusted to the compartment conditions.

        Returns
        -------
        float
            transformed Gibbs energy of the reaction
        """
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
        """Gibbs energy of transport of the reaction.

        Returns
        -------
        float
            Gibbs energy of transport of the reaction
        """
        try:
            return self._delG_transport
        except AttributeError:
            self._delG_transport = self.calculate_delG_transport()
            return self._delG_transport

    @delG_transport.setter
    def delG_transport(self, value):
        self._delG_transport = value  # Update the consraints

    def calculate_transport_charge(self):
        """calculates the net charge and protons transported across the membrane. Net charge and protons transported for each compartment should be same. This function assumes that transport reactions involve only two compartments. Any reaction involving more than 2 compartments are not necessarily true transport reactions like biomass or lumped reaction, which are not really useful for transport considerations.

        Returns
        -------
        tuple
            tuple of two dictionaries. Dictionaries of transport charge and protons for each compartment.
        """
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
        """Gibbs energy of transport has two components, 1) proton transport across the membrane, This needs to be adjusted for compartment specific pH 2) charge transport, which is affected by membrane potential.

        Returns
        -------
        float
            Gibbs energy of transport
        """
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
        """stoichiometric vector of forward reaction. Vector has length of number of metabolites in the model.

        Returns
        -------
        np.ndarray
            numpy vector of reaction stoichiometry
        """
        S_matrix = zeros((len(self.model.metabolites)))
        for metabolite, stoic in iteritems(self.metabolites):
            S_matrix[self.model.metabolites.index(metabolite)] = stoic
        return S_matrix

    def get_coefficient(self, metabolite_id):
        _id_to_metabolites = {m.id: m for m in self.metabolites}
        return self.metabolites[_id_to_metabolites[metabolite_id]]

    def calculate_delG_prime(self):
        """Function to calculate the transformed Gibbs energy of reaction. Its stoichiometric combination of transformed Gibbs energy of formation of metabolites

        Returns
        -------
        float
            transformed reaction Gibbs energy
        """
        rxn_delG = 0
        for metabolite, stoic in iteritems(self.metabolites):
            rxn_delG += stoic * float(metabolite.delG_f)

        return rxn_delG
