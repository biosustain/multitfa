from cobra import Metabolite
from comp_cache import comp_cache
from numpy import logaddexp, log
from util.dGf_calculation import calculate_dGf
from util.thermo_constants import RT
from optlang import Variable


class compound():

    def __init__(self,comp, id = None, Kegg_map = {}, concentration_min = 1e-5, concentration_max = 1e-2,
     delG_f = 0, pH = None, ionic_strength = None, temperature = None, conc_variable = None, Ci_variable = None):
        super().__init__()
        self.comp = comp
        self.id = comp.id
        self.Kegg_map = Kegg_map
        self.delG_f = self.calculate_delG_f()
        self.concentration_min = concentration_min
        self.concentration_max = concentration_max
        self.conc_variable = Variable('lnc_{}'.format(comp.id), lb = log(self.concentration_min),ub= log(self.concentration_max))
        self.Ci_variable = Variable('z_f_{}'.format(self.id), lb =-1.96, ub=1.96)
        self.constraint = comp.constraint
        self.compartment = comp.compartment   
        self.reactions = comp.reactions  
        self.elements = comp.elements 
        self.charge = comp.charge

    def transform(self,pH,ionic_strength,temperature):
        ccache = comp_cache(self.Kegg_map[self.id])
        if ccache.microspecies != []:
            ddg_over_rt = sorted(-1 * ms.transform(pH,ionic_strength,temperature) for ms in ccache.microspecies)

            total_ddg_over_rt = ddg_over_rt[0]
            for x in ddg_over_rt[1:]:
                total_ddg_over_rt = logaddexp(total_ddg_over_rt,x)
        else:
            total_ddg_over_rt = 0
        return - RT * total_ddg_over_rt

    def calculate_delG_f(self):
        dG0f = calculate_dGf([self.id],self.Kegg_map)
       
        return float(dG0f[0])