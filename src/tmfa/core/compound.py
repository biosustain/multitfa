from copy import copy

from cobra import Metabolite
from equilibrator_api import Q_
from six import iteritems

from ..util.thermo_constants import *


class Thermo_met(Metabolite):
    """
    Class representation of thermodynamic metabolite object.

    It takes cobra Metabolite object as input and creates the Thermo_met object. Object
    is populated with properties required for thermodynamic calculations. We use
    equilibrator-api to do the Gibbs energy calculations.

    Parameters
    ----------
    metabolite : cobra.core.Metabolite
        cobra metabolite object
    updated_model : cobra.core.Model, optional
        cobra model, by default None

    """

    def __init__(
        self,
        metabolite,
        updated_model=None,
    ):
        self._model = updated_model
        self._reaction = set()
        do_not_copy_by_ref = {"_reaction", "_model"}
        for attr, value in iteritems(metabolite.__dict__):
            if attr not in do_not_copy_by_ref:
                self.__dict__[attr] = copy(value) if attr == "formula" else value

    @property
    def concentration_variable(self):
        """optlang variable for concentration of metabolite. This variable varies between lower and upper bounds of specified concentration. if not specified, uses generic concentration bounds.

        Returns
        -------
        optlang.interface.variable
            optlang interface variable for concentration of metabolite. if metabolite is not linked to model, then None
        """
        if self.model is not None:
            conc_var = "lnc_{}".format(self.id)
            return self.model.variables[conc_var]
        else:
            return None

    @property
    def concentration_min(self):
        """lower bound of metabolite concentration. Units should be in Moles. If not specified, it either assumes compartment specific information (metabolite.model.compartment_info["c_min]) or uses 1e-5 M.

        Returns
        -------
        float
            metabolite concentration lower bound. returns input value or compartment specific value or 1e-5
        """
        try:
            return self._concentration_min
        except AttributeError:
            if self.model.compartment_info is not None:
                try:
                    self._concentration_min = float(
                        self.model.compartment_info["c_min"][self.compartment]
                    )
                except KeyError:
                    self._concentration_min = 1e-5
                return self._concentration_min
            else:
                self._concentration_min = float(1e-5)
                return self._concentration_min

    @concentration_min.setter
    def concentration_min(self, value):
        self._concentration_min = value
        self.concentration_variable.set_bounds(
            lb=np.log(value), ub=np.log(self.concentration_max)
        )
        # Update the Gurobi and Cplex interface bounds
        if self.model.gurobi_interface is not None:
            self.model.gurobi_interface.getVarByName(
                self.concentration_variable.name
            ).LB = np.log(value)
            self.model.gurobi_interface.update()

        if self.model.cplex_interface is not None:
            self.model.cplex_interface.variables.set_lower_bounds(
                self.concentration_variable.name, np.log(value)
            )

    @property
    def concentration_max(self):
        """Upper bound of metabolite concentration. Units should be in Moles. If not specified, it either assumes compartment specific information (metabolite.model.compartment_info["c_max]) or uses 2e-2 M.

        Returns
        -------
        float
            metabolite concentration upper bound. returns input value or compartment specific value or 2e-2
        """
        try:
            return self._concentration_max
        except AttributeError:
            if self.model.compartment_info is not None:
                try:
                    self._concentration_max = float(
                        self.model.compartment_info["c_max"][self.compartment]
                    )
                except KeyError:
                    self._concentration_max = 2e-2
                return self._concentration_max
            else:
                self._concentration_max = float(2e-2)
                return self._concentration_max

    @concentration_max.setter
    def concentration_max(self, value):
        self._concentration_max = value
        self.concentration_variable.set_bounds(
            lb=np.log(self.concentration_min), ub=np.log(value)
        )
        # Update the Gurobi and cplex interface bounds
        if self.model.gurobi_interface is not None:
            self.model.gurobi_interface.getVarByName(
                self.concentration_variable.name
            ).UB = np.log(value)
            self.model.gurobi_interface.update()

        if self.model.cplex_interface is not None:
            self.model.cplex_interface.variables.set_upper_bounds(
                self.concentration_variable.name, np.log(value)
            )

    @property
    def Kegg_id(self):
        """External database identifier to match and retrieve thermodynamic properties from equilibrator-api. Will accept any identifier that equilibrator uses. Please refer to equilibrator-api documentation. Should use their corresponding format.

        for example for atp
        E.g:
        if bigg databse id: bigg.metabolite:atp
            Kegg: kegg:C00002

        Returns
        -------
        str
            External database id for metabolite
        """
        return self._Kegg_id

    @Kegg_id.setter
    def Kegg_id(self, value):
        self._Kegg_id = value
        # self._equilibrator_accession = api.get_compound(self.Kegg_id)

    @property
    def delG_err_variable(self):
        """optlang variable to represent the error on the formation energy of metabolite. Allowed to vary between 2 S.D from mean.

        Returns
        -------
        optlang.interface.variable
            metabolite Gibbs energy error variable, if metabolite is not associated with model, then None
        """
        if self.model is not None:
            err_var = "dG_err_{}".format(self.id)
            return self.model.variables[err_var]
        else:
            return None

    @property
    def compound_vector(self):
        """Metabolite component decomposition vector. Its useful for represent what compinents are present in the compound.

        Returns
        -------
        np.ndarray
            metabolite component decomposition vector
        """
        try:
            return self._compound_vector
        except AttributeError:
            self._compound_vector = self.get_compound_vector()
            return self._compound_vector

    @property
    def delG_f(self):
        """Gibbs energy of formation of the metabolite in kJ/mol.

        Returns
        -------
        float
            Formation energy of metabolite in kJ/mol
        """
        try:
            return self._delG_f
        except AttributeError:
            self._delG_f = self.calculate_delG_f()
            return self._delG_f

    @delG_f.setter
    def delG_f(self, value):
        self._delG_f = value

    @property
    def std_dev(self):
        """standard deviation of compound Gibbs formation energy.

        Returns
        -------
        float
            compound standard deviation of formation
        """
        try:
            return self._std_dev
        except AttributeError:
            variance = self.compound_vector @ covariance @ self.compound_vector.T
            self._std_dev = np.sqrt(variance[0][0])
            return self._std_dev

    @property
    def major_ms(self, p_mg=14):
        try:
            return self._major_ms
        except AttributeError:
            self._major_ms = self.abundant_ms(
                pH=Q_(self.model.compartment_info["pH"][self.compartment]),
                I=Q_(str(self.model.compartment_info["I"][self.compartment]) + " M"),
                temperature=Q_(str(default_T) + " K"),
                pMg=Q_(self.model.compartment_info["pMg"][self.compartment]),
            )
            return self._major_ms

    @property
    def is_exclude(self):
        """Function to check if metabolite needs to be excluded from thermodynamic analysis. if the compound component decomposition vector is all zeros then compound is excluded from thermodynamic analysis.

        Returns
        -------
        Boolean
            True if component decomposition vector is zeros, False otherwise.
        """
        if self.is_proton:
            return False
        elif ~self.compound_vector.any():
            return True
        else:
            return False

    @property
    def is_proton(self):
        """Checks if the compound is proton, matching against inchikey

        Returns
        -------
        Boolean
            True if proton, False otherwise
        """
        if self.equilibrator_accession:
            if self.equilibrator_accession.inchi_key == PROTON_INCHI_KEY:
                return True
        else:
            return False

    @property
    def equilibrator_accession(self):
        """retrieves metabolite accession from equilibrator to get the termodynamic properties. We do this in batch in self.model,metabolite_equilibrator_accessions.

        Returns
        -------
        equilibrator.compound
            equilibrator compound object
        """
        try:
            return self._equilibrator_accession
        except AttributeError:
            self._equilibrator_accession = (
                self.model.metabolite_equilibrator_accessions[self.id]
            )
            return self._equilibrator_accession

    @property
    def std_dev_expression(self):
        """optlang expression to represent the metabolite error in terms of components instead of formation energies.

        expression = compound_vector @ variables

        Returns
        -------
        sympy expression
            error expression in terms of component variables
        """
        if self.model._var_update:
            return (
                self.compound_vector[np.nonzero(self.compound_vector)]
                @ self.model.component_variables[np.nonzero(self.compound_vector)[1]]
            )

    def abundant_ms(self, pH, I, temperature, pMg):
        ddg_over_rts = [
            (
                ms.transform(
                    pH=pH,
                    ionic_strength=I,
                    T_in_K=temperature,
                    pMg=pMg,
                ),
                ms,
            )
            for ms in self.equilibrator_accession.microspecies
        ]
        min_ddg, min_ms = min(ddg_over_rts, key=lambda x: x[0])
        return min_ms

    def get_compound_vector(self):
        """This is the implementation of compound vector from component contribution. Checks if the compound is covered by group contribution, reactant contribution or neither

        Returns:
            comp_vector  np.array or None
                        Compound vector with index of the compound in training data (component contribution) or None

        """
        if self.equilibrator_accession:
            try:
                rc_index = rc_compound_ids.index(self.equilibrator_accession.id)
                comp_vector = np.zeros(Nc + Ng, dtype=float)
                comp_vector[rc_index] = 1
                return comp_vector[np.newaxis, :]
            except ValueError:
                if self.equilibrator_accession.group_vector:
                    comp_vector = np.hstack(
                        [
                            np.zeros(Nc, dtype=float),
                            self.equilibrator_accession.group_vector,
                        ]
                    )
                    return comp_vector[np.newaxis, :]
                else:
                    comp_vector = np.zeros(Nc + Ng, dtype=float)
                    return comp_vector[np.newaxis, :]
        else:
            comp_vector = np.zeros(Nc + Ng, dtype=float)
            return comp_vector[np.newaxis, :]

    def calculate_delG_f(self):
        """Calculates the standard transformed Gibbs formation energy of compound using component contribution method. pH, Ionic strength values are taken from model's compartment_info attribute

        Returns:
            float -- Transformed Gibbs energy of formation adjusted to pH, ionic strength of metabolite
        """

        std_dG_f = self.compound_vector @ mu
        if self.compound_vector.any():
            transform = self.equilibrator_accession.transform(
                p_h=Q_(self.model.compartment_info["pH"][self.compartment]),
                ionic_strength=Q_(
                    str(self.model.compartment_info["I"][self.compartment]) + " M"
                ),
                temperature=Q_(str(default_T) + " K"),
            )
            return std_dG_f[0] + transform.to_base_units().magnitude * 1e-3
        else:
            return std_dG_f[0]
