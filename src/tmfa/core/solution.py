from numpy import empty, nan
from pandas import DataFrame, Series


class Solution:
    def __init__(
        self,
        objective_value,
        status,
        fluxes,
        reduced_costs=None,
        shadow_prices=None,
        Gibbs_energies=None,
        metabolite_concentrations=None,
        solver=None,
    ):

        self.solver = solver
        self.objective_value = objective_value
        self.status = status
        self.fluxes = fluxes
        self.reduced_costs = reduced_costs
        self.Gibbs_energies = Gibbs_energies
        self.metabolite_concentrations = metabolite_concentrations

    def __repr__(self):
        """String representation of the solution instance."""
        if self.status != "optimal":
            return "<Solution {0:s} at 0x{1:x}>".format(self.status, id(self))
        return "<Solution {0:.3f} at 0x{1:x}>".format(self.objective_value, id(self))

    def __getitem__(self, reaction_id):
        """
        Return the flux of a reaction.

        Parameters
        ----------
        reaction : str
            A model reaction ID.
        """
        return self.fluxes[reaction_id]

    get_primal_by_id = __getitem__

    def to_frame(self):
        """Return the fluxes and reduced costs as a data frame"""
        return DataFrame({"fluxes": self.fluxes, "reduced_costs": self.reduced_costs})


def get_solution(
    model,
    reactions=None,
    metabolites=None,
    Gibbs_energy=None,
    met_concentrations=None,
    raise_error=False,
):

    if model.solver.status != "optimal":
        raise ValueError("model status not optimal")

    if reactions is None:
        reactions = model.reactions
    if metabolites is None:
        metabolites = model.metabolites
    if Gibbs_energy is None:
        Gibbs_energy = [
            var for var in model.solver.variables if var.name.startswith("dG_")
        ]
    if met_concentrations is None:
        met_concentrations = [
            var for var in model.solver.variables if var.name.startswith("lnc_")
        ]

    rxn_index = list()
    delG_index = list()
    met_conc_index = list()
    fluxes = empty(len(reactions))
    Gibbs_energies = empty(len(Gibbs_energy))
    metabolite_concentrations = empty(len(met_concentrations))
    reduced = empty(len(reactions))
    var_primals = model.solver.primal_values
    shadow = empty(len(metabolites))

    if model.solver.is_integer:
        reduced.fill(nan)
        shadow.fill(nan)

        for (i, rxn) in enumerate(reactions):
            rxn_index.append(rxn.id)
            fluxes[i] = var_primals[rxn.id] - var_primals[rxn.reverse_id]

        for (i, delG) in enumerate(Gibbs_energy):
            delG_index.append(delG.name)
            Gibbs_energies[i] = var_primals[delG.name]

        for (i, conc) in enumerate(met_concentrations):
            met_conc_index.append(conc.name)
            metabolite_concentrations[i] = var_primals[conc.name]

        met_index = [met.id for met in metabolites]
    else:
        var_duals = model.solver.reduced_costs
        for (i, rxn) in enumerate(reactions):
            forward = rxn.id
            reverse = rxn.reverse_id
            rxn_index.append(forward)
            fluxes[i] = var_primals[forward] - var_primals[reverse]
            reduced[i] = var_duals[forward] - var_duals[reverse]

        for (i, delG) in enumerate(Gibbs_energy):
            delG_index.append(delG.name)
            Gibbs_energies[i] = var_primals[delG.name]

        for (i, conc) in enumerate(met_concentrations):
            met_conc_index.append(conc.name)
            metabolite_concentrations[i] = var_primals[conc.name]

        met_index = list()
        constr_duals = model.solver.shadow_prices
        for (i, met) in enumerate(metabolites):
            met_index.append(met.id)
            shadow[i] = constr_duals[met.id]

    return Solution(
        model.solver.objective.value,
        model.solver.status,
        Series(index=rxn_index, data=fluxes, name="fluxes"),
        Series(index=rxn_index, data=reduced, name="reduced_costs"),
        Series(index=met_index, data=shadow, name="shadow_prices"),
        Series(index=delG_index, data=Gibbs_energies, name="Gibbs_energies"),
        Series(
            index=met_conc_index,
            data=metabolite_concentrations,
            name="metabolite_concentrations",
        ),
    )


def get_legacy_solution(
    model,
    solver="gurobi",
    reactions=None,
    metabolites=None,
    Gibbs_energy=None,
    met_concentrations=None,
    raise_error=False,
):
    if solver == "gurobi":
        if model.gurobi_interface.status == 2:
            solver_status = "optimal"
    elif solver == "cplex":
        solver_status = model.cplex_interface.solve()
    else:
        pass

    if solver_status != "optimal":
        raise ValueError("model status not optimal")

    if reactions is None:
        reactions = model.reactions
    if metabolites is None:
        metabolites = model.metabolites
    if Gibbs_energy is None:
        Gibbs_energy = [
            var for var in model.solver.variables if var.name.startswith("dG_")
        ]
    if met_concentrations is None:
        met_concentrations = [
            var for var in model.solver.variables if var.name.startswith("lnc_")
        ]

    rxn_index = list()
    delG_index = list()
    met_conc_index = list()
    fluxes = empty(len(reactions))
    Gibbs_energies = empty(len(Gibbs_energy))
    metabolite_concentrations = empty(len(met_concentrations))
    reduced = empty(len(reactions))
    shadow = empty(len(metabolites))

    if solver == "gurobi":
        solver_interface = model.gurobi_interface
        reduced.fill(nan)
        shadow.fill(nan)

        for (i, rxn) in enumerate(reactions):
            rxn_index.append(rxn.id)
            fluxes[i] = (
                solver_interface.getVarByName(rxn.id).x
                - solver_interface.getVarByName(rxn.reverse_id).x
            )

        for (i, delG) in enumerate(Gibbs_energy):
            delG_index.append(delG.name)
            Gibbs_energies[i] = solver_interface.getVarByName(delG.name).x

        for (i, conc) in enumerate(met_concentrations):
            met_conc_index.append(conc.name)
            metabolite_concentrations[i] = solver_interface.getVarByName(conc.name).x

        met_index = [met.id for met in metabolites]

    elif solver == "cplex":
        solver_interface = model.cplex_interface
        reduced.fill(nan)
        shadow.fill(nan)

        for (i, rxn) in enumerate(reactions):
            rxn_index.append(rxn.id)
            fluxes[i] = (
                solver_interface.integer_var(name=rxn.id).solution_value
                - solver_interface.integer_var(name=rxn.reverse_id).solution_value
            )

        for (i, delG) in enumerate(Gibbs_energy):
            delG_index.append(delG.name)
            Gibbs_energies[i] = solver_interface.integer_var(
                name=delG.name
            ).solution_value

        for (i, conc) in enumerate(met_concentrations):
            met_conc_index.append(conc.name)
            metabolite_concentrations[i] = solver_interface.integer_var(
                conc.name
            ).solution_value

        met_index = [met.id for met in metabolites]
    return Solution(
        model.solver.objective.value,
        model.solver.status,
        Series(index=rxn_index, data=fluxes, name="fluxes"),
        Series(index=rxn_index, data=reduced, name="reduced_costs"),
        Series(index=met_index, data=shadow, name="shadow_prices"),
        Series(index=delG_index, data=Gibbs_energies, name="Gibbs_energies"),
        Series(
            index=met_conc_index,
            data=metabolite_concentrations,
            name="metabolite_concentrations",
        ),
    )
