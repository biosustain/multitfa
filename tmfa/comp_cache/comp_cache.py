from pickle import load
from .microspecies import microspecies
import os

data_dir = os.path.normpath(os.path.dirname(os.path.abspath(__file__)) + os.sep + os.pardir + os.sep + 'Data')
compound_data = load(open(data_dir + os.sep + 'compound_data.pickle','rb'))
Kegg2id = load(open(data_dir + os.sep + 'db2internal.pickle','rb'))
microspecies_data = load(open(data_dir + os.sep + 'microspecies.pickle','rb'))
microspecies_dict = microspecies_data[1]

class comp_cache():
    """[summary]
    
    Returns:
        [type] -- [description]
    """
    def __init__(self,comp,internal_id = None, mass = 0, atom_bag = {}, dissociation_constants = [], group_vector = [], inchikey = None, microspecies = []):
        self.comp = comp
        self.internal_id = self.get_internal_id(comp)
        self.atom_bag = self.get_comp_property(comp,'atom_bag')
        self.dissociation_constants = self.get_comp_property(comp,'dissociation_const')
        self.mass = self.get_comp_property(comp,'mass')
        self.group_vector = self.get_comp_property(comp,'group_vector')
        self.inchikey = self.get_comp_property(comp,'inchikey')
        self.microspecies = self.get_microspecies(comp)

    def get_internal_id(self,comp):
        if comp in Kegg2id.keys():
            return Kegg2id[comp]
        else:
            return -1
    def get_comp_property(self,comp,comp_property):
        internal_id = self.get_internal_id(comp)
        if internal_id != -1:
            return compound_data[comp_property][internal_id]
    def get_microspecies(self,comp):
        internal_id = self.get_internal_id(comp)
        if internal_id == -1:
            ms = []
        else:
            if str(internal_id) in microspecies_dict:
                msid = microspecies_dict[str(internal_id)]
                ms = [microspecies(ms) for ms in msid]
            else:
                ms = []
        return ms