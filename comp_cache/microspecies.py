from pickle import load
from numpy import log
import os
filepath = os.path.normpath(os.path.dirname(os.path.abspath(__file__)) + os.sep + os.pardir + os.sep + 'Data')
microspecies_data = load(open(filepath + os.sep + 'microspecies.pickle','rb'))
compound_microspecies = microspecies_data[0]


class microspecies():
    def __init__(self,ms,num_protons = 0, charge = 0, ddG_over_RT = 0):
        self.ms = ms
        self.num_protons = float(compound_microspecies['num_protons'][ms])
        self.charge = float(compound_microspecies['charge'][ms])
        self.ddG_over_RT = float(compound_microspecies['ddg_over_rt'][ms])

    @staticmethod
    def debye_hueckel(ionic_strength, temperature):
        a1 = 1.108 # (1/M**0.5)
        a2 = 1.546e-3 # (1/M**0.5/K)
        a3 = 5.959e-6 # (1 / M**0.5 / K**2)
        alpha = a1 - a2 * temperature + a3 * temperature ** 2
        B = 1.6 # (1/M**0.5)
        return alpha * ionic_strength ** 0.5 / (1.0 + B * ionic_strength ** 0.5)
    
    def legendre_transform(self,pH, ionic_strength, temperature):
        deb_hueck = self.debye_hueckel(ionic_strength, temperature)
        return self.num_protons * log(10) * float(pH)  + (self.num_protons - self.charge ** 2) * deb_hueck

    def transform(self,pH, ionic_strength,temperature):
        return self.ddG_over_RT + self.legendre_transform(pH,ionic_strength,temperature)
