from cobra import Metabolite
from ..comp_cache import comp_cache
from numpy import logaddexp, log, diag, sqrt
from ..util.dGf_calculation import calculate_dGf
from ..util.thermo_constants import RT
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
        self,
        metabolite,
        updated_model=None,
        Kegg_map={},
        concentration_dict={"min": {}, "max": {}},
        delG_f=0,
        pH=None,
        ionic_strength=None,
        temperature=None,
    ):

        self._model = updated_model
        self._reaction = set()
        do_not_copy_by_ref = {"_reaction", "_model"}
        for attr, value in iteritems(metabolite.__dict__):
            if attr not in do_not_copy_by_ref:
                self.__dict__[attr] = copy(value) if attr == "formula" else value
        self.Kegg_id = Kegg_map[self.id]
        self.Kegg_map = Kegg_map
        self.delG_f = self.calculate_delG_f()
        # self.std_dev = sqrt(diag(self.model.cov_dG)[
        #                        self.model.metabolites.index(self)])

        if self.id in concentration_dict["min"].keys():
            self._concentration_min = concentration_dict["min"][self.id]
        else:
            self._concentration_min = float(1e-5)
        if self.id in concentration_dict["max"].keys():
            self._concentration_max = concentration_dict["max"][self.id]
        else:
            self._concentration_max = float(2e-2)

    @property
    def concentration_variable(self):
        if self.model is not None:
            conc_var = "lnc_{}".format(self.id)
            return self.model.variables[conc_var]
        else:
            return None

    @property
    def concentration_min(self):
        return self._concentration_min

    @concentration_min.setter
    def concentration_min(self, value):
        self._concentration_min = value
        self.concentration_variable.set_bounds(
            lb=log(value), ub=log(self.concentration_max)
        )

    @property
    def concentration_max(self):
        return self._concentration_max

    @concentration_max.setter
    def concentration_max(self, value):
        self._concentration_max = value
        self.concentration_variable.set_bounds(
            lb=log(self.concentration_min), ub=log(value)
        )

    @property
    def compound_variable(self):
        if self.model is not None:
            conc_var = "met_{}".format(self.id)
            return self.model.variables[conc_var]
        else:
            return None

    def transform(self, pH, ionic_strength, temperature):
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
                -1 * ms.transform(pH, ionic_strength, temperature)
                for ms in ccache.microspecies
            )

            total_ddg_over_rt = ddg_over_rt[0]
            for x in ddg_over_rt[1:]:
                total_ddg_over_rt = logaddexp(total_ddg_over_rt, x)
        else:
            total_ddg_over_rt = 0
        return -RT * total_ddg_over_rt

    def calculate_delG_f(self):
        """ Calculates the standars formation energy of compound using component contribution method
        
        Returns:
            float -- delG_f
        """

        dG0f = calculate_dGf([self], self.Kegg_map)

        return float(dG0f[0])
