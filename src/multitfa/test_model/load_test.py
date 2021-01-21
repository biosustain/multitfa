import pandas as pd
from ..core import tmodel
from cobra import io
from os.path import abspath, dirname, join

core_model_name = join(dirname(abspath(__file__)), "model", "e_coli_core.mat")


def load_test_data():

    model = io.load_matlab_model(core_model_name)
    pH_I_T_dict = {
        "pH": {"c": 7.5, "e": 7, "p": 7},
        "I": {"c": 0.25, "e": 0, "p": 0},
        "T": {"c": 298.15, "e": 298.15, "p": 298.15},
    }
    del_psi_dict = {
        "c": {"c": 0, "e": 0, "p": 150},
        "e": {"c": 0, "e": 0, "p": 0},
        "p": {"c": -150, "e": 0, "p": 0},
    }

    del_psi = pd.DataFrame.from_dict(data=del_psi_dict)
    comp_info = pd.DataFrame.from_dict(data=pH_I_T_dict)

    Excl = [
        rxn.id
        for rxn in model.reactions
        if rxn.id.startswith("EX_") or rxn.id.startswith("DM_")
    ] + ["BIOMASS_Ecoli_core_w_GAM", "O2t", "H2Ot"]

    tfa_model = tmodel(
        model, Exclude_list=Excl, compartment_info=comp_info, membrane_potential=del_psi
    )
    for met in tfa_model.metabolites:
        kegg_id = "bigg.metabolite:" + met.id[:-2]
        met.Kegg_id = kegg_id

    tfa_model.update()

    return tfa_model