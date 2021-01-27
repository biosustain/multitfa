"""Provide an example that generates the data for the paper figure."""


import pandas as pd
from cobra import Configuration
from cobra.io import load_model

# For Gurobi, use the other function `variability_legacy_gurobi`.
from multitfa.analysis import variability_legacy_cplex
from multitfa.core import tmodel


def build_core_model():
    model = load_model("e_coli_core")
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

    Excl = [rxn.id for rxn in model.boundary] + [
        "BIOMASS_Ecoli_core_w_GAM",
        "O2t",
        "H2Ot",
    ]

    tfa_model = tmodel(
        model, Exclude_list=Excl, compartment_info=comp_info, membrane_potential=del_psi
    )
    for met in tfa_model.metabolites:
        kegg_id = "bigg.metabolite:" + met.id[:-2]
        met.Kegg_id = kegg_id

    tfa_model.update()

    return tfa_model


def main():
    config = Configuration()
    config.solver = "cplex"
    tfa_model = build_core_model()
    vars_analysis = [
        rxn.id for rxn in tfa_model.reactions if not rxn.id.startswith("DM_")
    ]
    ranges = variability_legacy_cplex(tfa_model, variable_list=vars_analysis)

    ranges.to_csv("ranges_ecoli_cplex.tsv", header=True, index=True, sep="\t")


if __name__ == "__main__":
    main()
