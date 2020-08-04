from cobra import Metabolite
from ..comp_cache import comp_cache
from numpy import logaddexp, log, diag, sqrt
from ..util.dGf_calculation import calculate_dGf
from ..util.thermo_constants import RT, default_T
from six import iteritems
from copy import copy, deepcopy


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
                self._concentration_min = float(
                    self.model.compartment_info["c_min"][self.compartment]
                )
                return self._concentration_min
            self._concentration_min = float(1e-5)
            return self._concentration_min

    @concentration_min.setter
    def concentration_min(self, value):
        self._concentration_min = value
        self.concentration_variable.set_bounds(
            lb=log(value), ub=log(self.concentration_max)
        )
        # Update the Gurobi and Cplex interface bounds
        if self.model.gurobi_interface is not None:
            self.model.gurobi_interface.getVarByName(
                self.concentration_variable.name
            ).LB = log(value)
            self.model.gurobi_interface.update()

        if self.model.cplex_interface is not None:
            self.model.cplex_interface.variables.set_lower_bounds(
                self.concentration_variable.name, log(value)
            )

    @property
    def concentration_max(self):
        try:
            return self._concentration_max
        except AttributeError:
            if self.model.compartment_info is not None:
                self._concentration_max = float(
                    self.model.compartment_info["c_max"][self.compartment]
                )
                return self._concentration_max
            self._concentration_max = float(2e-2)
            return self._concentration_max

    @concentration_max.setter
    def concentration_max(self, value):
        self._concentration_max = value
        self.concentration_variable.set_bounds(
            lb=log(self.concentration_min), ub=log(value)
        )
        # Update the Gurobi and cplex interface bounds
        if self.model.gurobi_interface is not None:
            self.model.gurobi_interface.getVarByName(
                self.concentration_variable.name
            ).UB = log(value)
            self.model.gurobi_interface.update()

        if self.model.cplex_interface is not None:
            self.model.cplex_interface.variables.set_upper_bounds(
                self.concentration_variable.name, log(value)
            )

    @property
    def Kegg_id(self):
        return self._Kegg_id

    @Kegg_id.setter
    def Kegg_id(self, value):
        self._Kegg_id = value

    @property
    def compound_variable(self):
        if self.model is not None:
            conc_var = "met_{}".format(self.Kegg_id)
            return self.model.variables[conc_var]
        else:
            return None

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
        if self.model is not None:
            std_dev = sqrt(
                diag(self.model.cov_dG)[self.model.metabolites.index(self.id)]
            )

            return std_dev
        else:
            return 0

    def transform(self, pH, ionic_strength):
        """ Transform the delGf with respect to pH, ionic strength based on Alberty calculations
        
        Arguments:
            pH {float} -- pH
            ionic_strength {float} -- ionic strength in Molar
            temperature {float} -- temperature in Kelvin
        
        Returns:
            [type] -- [description]
        """
        ccache = comp_cache(self.Kegg_id)
        if ccache.microspecies != []:
            ddg_over_rt = sorted(
                -1 * ms.transform(pH, ionic_strength, default_T)
                for ms in ccache.microspecies
            )

            total_ddg_over_rt = ddg_over_rt[0]
            for x in ddg_over_rt[1:]:
                total_ddg_over_rt = logaddexp(total_ddg_over_rt, x)
        else:
            total_ddg_over_rt = 0
        return -RT * total_ddg_over_rt

    def calculate_delG_f(self):
        """ Calculates the standard formation energy of compound using component contribution method
        
        Returns:
            float -- delG_f
        """

        dG0f = calculate_dGf([self])
        pH = self.model.compartment_info["pH"][self.compartment]
        ionic_strength = self.model.compartment_info["I"][self.compartment]
        transform_value = self.transform(pH, ionic_strength)

        return float(dG0f[0]) + transform_value
