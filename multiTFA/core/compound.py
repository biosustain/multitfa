from cobra import Metabolite
import numpy as np
from ..util.thermo_constants import *
from six import iteritems
from copy import copy, deepcopy
from equilibrator_api import ComponentContribution, Q_
from equilibrator_cache import Compound
import time

api = ComponentContribution()


class Thermo_met(Metabolite):
    """ Class representation of thermodynamic metabolite object
    
    Arguments:
        Metabolite {cobra.core.Metabolite} -- Cobra metabolite object
    
    Attributes:
        All the cobra metabolite attributes
        Kegg_id: Kegg/Seed id of the metabolite
        delG_f: standard formation energy of the compound
        concentration_min: lb concentration
        concentration_max: ub concentration
    """

    def __init__(
        self, metabolite, updated_model=None,
    ):

        self._model = updated_model
        self._reaction = set()
        do_not_copy_by_ref = {"_reaction", "_model"}
        for attr, value in iteritems(metabolite.__dict__):
            if attr not in do_not_copy_by_ref:
                self.__dict__[attr] = copy(value) if attr == "formula" else value

    @property
    def concentration_variable(self):
        if self.model is not None:
            conc_var = "lnc_{}".format(self.id)
            return self.model.variables[conc_var]
        else:
            return None

    @property
    def concentration_min(self):
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
        return self._Kegg_id

    @Kegg_id.setter
    def Kegg_id(self, value):
        self._Kegg_id = value
        # self._equilibrator_accession = api.get_compound(self.Kegg_id)

    @property
    def compound_variable(self):
        if self.model is not None:
            conc_var = "met_{}".format(self.Kegg_id)
            return self.model.variables[conc_var]
        else:
            return None

    @property
    def compound_vector(self):
        try:
            return self._compound_vector
        except AttributeError:
            self._compound_vector = self.get_compound_vector()
            return self._compound_vector

    @property
    def delG_f(self):
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
        if self.is_proton:
            return False
        elif ~self.compound_vector.any():
            return True
        else:
            return False

    @property
    def is_proton(self):

        if self.equilibrator_accession:
            if self.equilibrator_accession.inchi_key == PROTON_INCHI_KEY:
                return True
        else:
            return False

    @property
    def sphere_var_expression(self):
        return (
            np.sqrt(chi2_value)
            * self.compound_vector
            @ cholesky
            @ self.model.sphere_variables
        )

    @property
    def equilibrator_accession(self):
        try:
            return self._equilibrator_accession
            # print(self._equilibrator_accession.id)
        except AttributeError:
            self._equilibrator_accession = api.get_compound(self.Kegg_id)
            return self._equilibrator_accession

    def abundant_ms(self, pH, I, temperature, pMg):
        ddg_over_rts = [
            (ms.transform(pH=pH, ionic_strength=I, T_in_K=temperature, pMg=pMg,), ms,)
            for ms in self.equilibrator_accession.microspecies
        ]
        min_ddg, min_ms = min(ddg_over_rts, key=lambda x: x[0])
        return min_ms

    def get_compound_vector(self):
        """ This is the implementation of compound vector from component contribution. Checks if the compound is covered by group contribution, reactant contribution or neither

        Returns:
            comp_vector  np.array or None
                        Compound vector with index of the compound in training data (component contribution) or None
            
        """
        # self._equilibrator_accession = api.get_compound(self.Kegg_id)
        # print(self._equilibrator_accession.id)
        # time.sleep(0.5)
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
        """ Calculates the standard transformed Gibbs formation energy of compound using component contribution method. pH, Ionic strength values are taken from model's compartment_info attribute
        
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

