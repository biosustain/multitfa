import os
import pickle
import string
import tempfile
from copy import copy, deepcopy
from random import choices

import optlang
from cobra import Model
from cobra.core.dictlist import DictList
from equilibrator_api import ComponentContribution
from scipy import stats
from six import iteritems

from ..util.constraints import delG_indicator, directionality
from ..util.linalg_fun import *
from ..util.thermo_constants import *
from .compound import Thermo_met
from .reaction import thermo_reaction
from .solution import get_legacy_solution, get_solution


api = ComponentContribution()

# present_dir = os.path.normpath(os.path.dirname(os.path.abspath(__file__)))
cwd = os.getcwd()

import logging


logs_dir = cwd + os.sep + "logs"
if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)

logging.basicConfig(
    filename=logs_dir + os.sep + "thermo_model.log", level=logging.DEBUG
)

cache_file = os.path.normpath(
    os.path.dirname(os.path.abspath(__file__))
    + os.sep
    + os.pardir
    + os.sep
    + os.pardir
    + os.sep
    + "cache"
    + os.sep
    + "compounds_cache.pickle"
)


class tmodel(Model):
    """tmodel is Class representation of thermodynamic metabolic flux analysis model. This class adds attributes and methods required for thermodynamic analysis of a COBRA model.

    Parameters
    ----------
    Model : cobra.core.model
        A Cobra model class
    model : instance of cobra.core.model
        tmodel requires a cobra model instance as input. Thermodynamic properties are added based on the stoichiometry of the underlying Cobra model
    Exclude_list : list, optional
        List of reactions user wants to exclude from thermodynamic analysis, For example, Exchange/Demand reactions,  by default []
    tolerance_integral : float, optional
        integrality tolerance of for the model , by default 1e-9
    compartment_info : pd.Dataframe, optional
        a pandas Dataframe containing the compartment information like pH, ionic strength, magnesium concentration etc. Row indices should be the compartment symbol and column indices should be property, by default None
    membrane_potential : pd.Dataframe, optional
        a pandas Dataframe containing membrane electrostatic potential information to calculate the delG of muli compartment transport. Values are read in a sequence that column represent the first compartment and row represent the compartment being transported to. Row & column indices should be the compartment symbols , by default None
    """

    def __init__(
        self,
        model,
        Exclude_list=[],
        tolerance_integral=1e-9,
        compartment_info=None,
        membrane_potential=None,
    ):

        self.compartment_info = compartment_info
        self.membrane_potential = membrane_potential

        do_not_copy_by_ref = {
            "metabolites",
            "reactions",
            "genes",
            "notes",
            "annotation",
        }
        for attr in model.__dict__:
            if attr not in do_not_copy_by_ref:
                self.__dict__[attr] = model.__dict__[attr]

        self.metabolites = DictList()
        do_not_copy_by_ref = {"_reaction", "_model"}
        for metabolite in model.metabolites:
            new_met = Thermo_met(
                metabolite=metabolite,
                updated_model=self,
            )
            self.metabolites.append(new_met)

        self.genes = DictList()
        for gene in model.genes:
            new_gene = gene.__class__(None)
            for attr, value in iteritems(gene.__dict__):
                if attr not in do_not_copy_by_ref:
                    new_gene.__dict__[attr] = (
                        copy(value) if attr == "formula" else value
                    )
            new_gene._model = self
            self.genes.append(new_gene)

        self.reactions = DictList()
        do_not_copy_by_ref = {"_model", "_metabolites", "_genes"}
        for reaction in model.reactions:
            new_reaction = thermo_reaction(
                cobra_rxn=reaction,
                updated_model=self,
            )
            self.reactions.append(new_reaction)

        try:
            self._solver = deepcopy(model.solver)
            # Cplex has an issue with deep copies
        except Exception:  # pragma: no cover
            self._solver = copy(model.solver)  # pragma: no cover

        self.Exclude_list = Exclude_list
        self.solver.configuration.tolerances.integrality = tolerance_integral
        self._var_update = False

    @property
    def gurobi_interface(self):
        """multiTFA at the moment supports two solvers Gurobi/Cplex for solving quadratic constraint problems. Optlang doesn't support adding QC, so we chose to add two separate solver interafaces to tmodel. This is gurobi solver interface. In addition to the linear constraints, this interface contain one extra constraint to represent sphere

        Returns
        -------
        gurobi model object
            Gurobi model containing the QC
        """
        try:
            return self._gurobi_interface
        except AttributeError:
            if self.solver.__class__.__module__ == "optlang.gurobi_interface":
                # self._gurobi_interface = self.solver.problem.copy()
                self._gurobi_interface = self.Quadratic_constraint()
                return self._gurobi_interface

    @property
    def cplex_interface(self):
        """Cplex interface to support QC

        Returns
        -------
        Cplex model
            Cplex model containing QC
        """
        try:
            return self._cplex_interface
        except AttributeError:
            if self.solver.__class__.__module__ == "optlang.cplex_interface":
                self._cplex_interface = self.Quadratic_constraint()
                return self._cplex_interface

    @property
    def problem_metabolites(self):
        """Metabolites for which we can't calculate the Gibbs free energy of formation using component contribution method. If the metabolite is not covered by reactant/group contribution method, then we have to write the metabolite and corresponding reactions from tMFA analysis or have to lump the reactions.

        Returns:
            List
                List of metabolite not covered by component contribution
        """
        problematic_metabolites = []
        for met in self.metabolites:
            if met.is_proton:
                continue
            if ~met.compound_vector.any():
                problematic_metabolites.append(met)
        return problematic_metabolites

    @property
    def Exclude_reactions(self):
        """Reactions that needs to be excluded from tMFA. This list includes both user excluded reactions and reactions involving non-coverage metabolites

        Returns
        -------
        List
            List of reactions Excluded from thermo analysis
        """
        try:
            return self._Exclude_reactions
        except AttributeError:
            self._Exclude_reactions = list(
                set(self.Exclude_list + self.problematic_rxns)
            )
            return self._Exclude_reactions

    @property
    def problematic_rxns(self):
        """List of reactions containing non-covered metabolites. These can either be written out or lumped

        Returns
        -------
        List
            List of non-covered reactions
        """
        try:
            return self._problematic_rxns
        except AttributeError:
            self._problematic_rxns = self.cal_problematic_rxns()
            return self._problematic_rxns

    @property
    def component_variables(self):
        return np.array(
            [var for var in self.variables if var.name.startswith("component_")]
        )

    @property
    def metabolite_equilibrator_accessions(self):
        try:
            return self._metabolite_equilibrator_accessions
        except AttributeError:
            self._metabolite_equilibrator_accessions = (
                self.populate_metabolite_properties()
            )
            return self._metabolite_equilibrator_accessions

    def populate_metabolite_properties(self):
        """Local cache file for equilibrator-api data. This is a temporary fix till equilibrator's cache

        Returns
        -------
        dict
            Dictionary of metabolite id to corresponding equilibrator compound object
        """
        if os.path.isfile(cache_file):
            with open(cache_file, "rb") as handle:
                metabolite_accessions, microspecies, mg_dissociation_data = pickle.load(
                    handle
                )
        else:
            metabolite_accessions = {}

        accessions = {}
        for metabolite in self.metabolites:
            if metabolite.Kegg_id in metabolite_accessions:
                accessions[metabolite.id] = metabolite_accessions[metabolite.Kegg_id]
            else:
                eq_accession = api.get_compound(metabolite.Kegg_id)
                accessions[metabolite.id] = eq_accession
                # update the cache file
                # if eq_accession is not None:
                #    metabolite_accessions[metabolite.Kegg_id] = eq_accession
                #    microspecies[metabolite.Kegg_id] = eq_accession.microspecies
                #    mg_dissociation_data[
                #        metabolite.Kegg_id
                #    ] = eq_accession.magnesium_dissociation_constants

        # Re-write the cache file with updated values
        # with open(cache_file, "wb") as handle:
        #    pickle.dump(
        #        [metabolite_accessions, microspecies, mg_dissociation_data], handle
        #    )

        return accessions

    @property
    def compound_vector_matrix(self):
        try:
            return self._compound_vector_matrix
        except AttributeError:
            # Initialize the matrix with zeros
            comp_vector = np.zeros((len(self.metabolites), Nc + Ng))
            for metabolite in self.metabolites:
                met_index = self.metabolites.index(metabolite)
                comp_vector[met_index, :] = metabolite.compound_vector
            self._compound_vector_matrix = comp_vector
            return self._compound_vector_matrix

    def core_stoichiometry(self):
        n_core_rxn = len(self.reactions) - len(self.Exclude_reactions)
        stoichiometry_core = np.zeros((n_core_rxn, len(self.metabolites)))
        i = 0
        rxn_var_name = []
        for reaction in self.reactions:
            if reaction.id in self.Exclude_reactions:
                continue
            rxn_stoichiometry = reaction.cal_stoichiometric_matrix()
            stoichiometry_core[i, :] = rxn_stoichiometry

            i = i + 2
            rxn_var_name.extend(
                [reaction.forward_variable.name, reaction.reverse_variable.name]
            )
        return (rxn_var_name, stoichiometry_core)

    def update_thermo_variables(self):
        """Generates reaction and metabolite variables required for thermodynamic analysis and adds to the model. We use two different methods to solve the tMFA problem. Traditional 'box' method employs MILP problem where components are allowed to vary between some s.d from mean. The other method uses MIQCP structure to use covariance matrix to capture covariance. Two methods share some common variables, where as MIQCP method requires independent variables to sample from original solution space.

        Common Variables: delG_reaction, indicator_reaction, concentration_metabolite
        MILP variables: metabolite error
        MIQCP variables: independent component variables

        """
        self._var_update = False

        # Add metabolite concentration variable and error variable for the metabolite
        conc_variables, dG_err_vars = ([], [])
        for metabolite in self.metabolites:
            conc_variable = self.problem.Variable(
                "lnc_{}".format(metabolite.id),
                lb=np.log(metabolite.concentration_min),
                ub=np.log(metabolite.concentration_max),
            )

            delG_err_variable = self.problem.Variable(
                "dG_err_{}".format(metabolite.id),
                lb=-1.96 * np.sqrt(metabolite.std_dev),
                ub=1.96 * np.sqrt(metabolite.std_dev),
            )
            conc_variables.append(conc_variable)
            dG_err_vars.append(delG_err_variable)

        self.add_cons_vars(conc_variables + dG_err_vars)

        # Adding the thermo variables for reactions, delG_reaction and indicator (binary)
        rxn_variables = []
        for rxn in self.reactions:
            if rxn.id in self.Exclude_reactions:
                continue
            delG_forward = self.problem.Variable(
                "dG_{}".format(rxn.forward_variable.name), lb=-1e6, ub=1e5
            )

            delG_reverse = self.problem.Variable(
                "dG_{}".format(rxn.reverse_variable.name), lb=-1e6, ub=1e5
            )

            indicator_forward = self.problem.Variable(
                "indicator_{}".format(rxn.forward_variable.name),
                lb=0,
                ub=1,
                type="binary",
            )

            indicator_reverse = self.problem.Variable(
                "indicator_{}".format(rxn.reverse_variable.name),
                lb=0,
                ub=1,
                type="binary",
            )
            rxn_variables.extend(
                [delG_forward, delG_reverse, indicator_forward, indicator_reverse]
            )
        self.add_cons_vars(rxn_variables)

        self._var_update = True

    def _generate_constraints(self):
        """Generates thermodynamic constraints for the model. See util/constraints.py for detailed explanation of constraints

        Vi - Vmax * Zi <= 0
        delGr - K + K * Zi <= 0
        delGr - RT * S.T * ln(x) - S.T @ delGf - delGtransport = 0


        Returns:
            List -- List of themrodynamic constraints
        """
        # First check if thermovariables are added to the model
        if not self._var_update:
            self.update_thermo_variables()

        rxn_constraints = []
        # Now add reaction variables and generate remaining constraints
        for rxn in self.reactions:
            if rxn.id in self.Exclude_reactions:
                logging.debug(
                    "Reaction {} is excluded from thermodyanmic analysis".format(rxn.id)
                )
                continue

            # Directionality constraint
            dir_f, dir_r = directionality(rxn)
            ind_f, ind_r = delG_indicator(rxn)

            rxn_constraints.extend([dir_f, dir_r, ind_f, ind_r])

            # Create two different constraints for box method and MIQC method

            # delG constraint for box
            concentration_term = sum(
                stoic * metabolite.concentration_variable
                for metabolite, stoic in iteritems(rxn.metabolites)
                if metabolite.equilibrator_accession.inchi_key != PROTON_INCHI_KEY
            )

            err_term = sum(
                stoic * metabolite.delG_err_variable
                for metabolite, stoic in iteritems(rxn.metabolites)
                if metabolite.equilibrator_accession.inchi_key != PROTON_INCHI_KEY
            )

            lhs_forward = rxn.delG_forward - RT * concentration_term - err_term
            lhs_reverse = rxn.delG_reverse + RT * concentration_term + err_term
            rhs = rxn.delG_prime + rxn.delG_transport

            delG_f = self.problem.Constraint(
                lhs_forward,
                lb=rhs,
                ub=rhs,
                name="delG_{}".format(rxn.forward_variable.name),
            )

            delG_r = self.problem.Constraint(
                lhs_reverse,
                lb=-rhs,
                ub=-rhs,
                name="delG_{}".format(rxn.reverse_variable.name),
            )
            rxn_constraints.extend([delG_f, delG_r])

        return rxn_constraints

    def update(self):
        """ Adds the generated thermo constaints to  model. Checks for duplication"""
        thermo_constraints = self._generate_constraints()

        for cons in thermo_constraints:
            if cons.name not in self.constraints:
                self.add_cons_vars([cons])
                logging.debug("Constraint {} added to the model".format(cons.name))
            else:
                logging.warning(
                    "Constraint {} already in the model, removing previous entry".format(
                        cons.name
                    )
                )
                self.solver.remove(cons.name)
                self.add_cons_vars([cons])

    def optimize(self, solve_method="QC", raise_error=False):
        """solves the model with given constraints. By default, we try to solve the model with quadratic constraints. Note: Quadratic constraints are supported by Gurobi/Cplex currently. if either of two solvers are not found, one can solve 'box' type MILP problem.

        :param solve_method: Method to solve the problem, defaults to "MIQC", any other string input leades to solving with box MILP
        :type solve_method: str, optional
        :param raise_error: , defaults to False
        :type raise_error: bool, optional
        :return: returns solution object
        :rtype: solution object (refer to Solution class)
        """

        if solve_method.lower() == "qc":
            if not (
                optlang.available_solvers["GUROBI"]
                or optlang.available_solvers["CPLEX"]
            ):
                logging.warning(
                    "GUROBI/CPLEX not available, Quadratic constraints are not supported by current solver"
                )
                print(
                    "GUROBI/CPLEX not available, Quadratic constraints are not supported by current solver, solving MIP problem instead."
                )
                self.slim_optimize()
                solution = get_solution(self, raise_error=raise_error)
                return solution

            if self.solver.__class__.__module__ == "optlang.gurobi_interface":
                self.gurobi_interface.optimize()
                solution = get_legacy_solution(self, solver="gurobi")

                return solution

            elif self.solver.__class__.__module__ == "optlang.cplex_interface":
                solution = self.cplex_interface.solve()
                solution = get_legacy_solution(self, solver="cplex")

                return solution

        elif solve_method.lower() == "mip":
            self.slim_optimize()
            solution = get_solution(self, raise_error=raise_error)

            return solution
        else:
            raise ValueError("Solver not understood")

    def Quadratic_constraint(self):
        """Adds Quadratic constraint to the model's Gurobi/Cplex Interface.
        (x-mu).T @ inv(cov) @ (x-mu) <= chi-square
        Note: This one creates one ellipsoidal constraint for all the metabolites that has non zero or non 'nan' formation energy, irrespective of the magnitude of variance. if the model is infeasible after adding this constraint, refer to util_func.py, find_correlated metabolites to add different ellipsoidal constraints to high variance and normal compounds to avoid possible numerical issues.

        Unable to retrieve quadratic constraints in Gurobi model, can see the QC when printed.

        :raises NotImplementedError: Implemented only for Gurobi/Cplex interfaces.
        :return: [description]
        :rtype: [type]
        """

        # Pick indices of components present in the current model
        model_component_indices = [
            i
            for i in range(self.compound_vector_matrix.shape[1])
            if np.any(self.compound_vector_matrix[:, i])
        ]

        # Reduced the compound_vector to contain only the non zero entries
        model_compound_vector = self.compound_vector_matrix[:, model_component_indices]

        # Now extract the sub covariance matrix containing only the components present in the model
        component_model_covariance = covariance[:, model_component_indices][
            model_component_indices, :
        ]

        # Now separate the compounds that have variance > 1000 and others to avoid numerical issues
        high_variance_indices = np.where(np.diag(component_model_covariance) > 1000)[0]
        low_variance_indices = np.where(np.diag(component_model_covariance) < 1000)[0]

        # Calculate cholesky matrix for two different covariance matrices
        if len(low_variance_indices) > 0:
            small_component_covariance = component_model_covariance[
                :, low_variance_indices
            ][low_variance_indices, :]
            cholesky_small_variance = matrix_decomposition(small_component_covariance)
            chi2_value_small = stats.chi2.isf(
                q=0.05, df=cholesky_small_variance.shape[1]
            )  # Chi-square value to map confidence interval

            for i in high_variance_indices:
                zeros_axis = np.zeros((cholesky_small_variance.shape[1],))
                cholesky_small_variance = np.insert(
                    cholesky_small_variance, i, zeros_axis, axis=0
                )

            metabolite_sphere_small = (
                model_compound_vector @ cholesky_small_variance
            )  # This is a fixed term compound_vector @ cholesky

        if len(high_variance_indices) > 0:
            large_component_covariance = component_model_covariance[
                :, high_variance_indices
            ][
                high_variance_indices, :
            ]  # Covariance matrix for the high variance components

            cholesky_large_variance = matrix_decomposition(large_component_covariance)
            chi2_value_high = stats.chi2.isf(
                q=0.05, df=cholesky_large_variance.shape[1]
            )

            # Insert empty rows for the low_variance_components
            for i in low_variance_indices:
                zeros_axis = np.zeros((cholesky_large_variance.shape[1],))
                cholesky_large_variance = np.insert(
                    cholesky_large_variance, i, zeros_axis, axis=0
                )
            metabolite_sphere_large = (
                model_compound_vector @ cholesky_large_variance
            )  # This is a fixed term compound_vector @ cholesky

        proton_indices = [
            self.metabolites.index(metabolite)
            for metabolite in self.metabolites
            if metabolite.equilibrator_accession is not None
            if metabolite.equilibrator_accession.inchi_key == PROTON_INCHI_KEY
        ]  # Get indices of protons in metabolite list to avoid double correcting them for concentrations

        if self.solver.__class__.__module__ == "optlang.cplex_interface":

            from cplex import Cplex, SparsePair, SparseTriple

            # Instantiate Cplex model
            cplex_model = Cplex()

            rand_str = "".join(choices(string.ascii_lowercase + string.digits, k=6))
            # write cplex model to mps file in random directory and re read
            with tempfile.TemporaryDirectory() as td:
                temp_filename = os.path.join(td, rand_str + ".mps")
                self.solver.problem.write(temp_filename)
                cplex_model.read(temp_filename)

            # Remove the unnecessary variables and constraints
            remove_vars = [
                var
                for var in cplex_model.variables.get_names()
                if var.startswith("component_") or var.startswith("dG_err_")
            ]  # Remove error variables

            remove_constrs = [
                cons
                for cons in cplex_model.linear_constraints.get_names()
                if cons.startswith("delG_") or cons.startswith("std_dev_")
            ]  # Remove delG constraint and re-add with component variables

            cplex_model.linear_constraints.delete(remove_constrs)  # Removing constr
            cplex_model.variables.delete(remove_vars)  # Removing Vars

            # QC for small variance components
            if len(low_variance_indices) > 0:
                indices_sphere1 = cplex_model.variables.add(
                    names=[
                        "Sphere1_{}".format(i)
                        for i in range(cholesky_small_variance.shape[1])
                    ],
                    lb=[-1] * cholesky_small_variance.shape[1],
                    ub=[1] * cholesky_small_variance.shape[1],
                )  # Adding independent component variables to the model, store the variable indices

                # Add the Sphere constraint
                cplex_model.quadratic_constraints.add(
                    quad_expr=SparseTriple(
                        ind1=indices_sphere1,
                        ind2=indices_sphere1,
                        val=len(indices_sphere1) * [1],
                    ),
                    sense="L",
                    rhs=1,
                    name="unit_normal_small_variance",
                )
            else:
                indices_sphere1 = []  # Just to adjust the matrix dimensions later

            # QC for large variance components
            if len(high_variance_indices) > 0:
                indices_sphere2 = cplex_model.variables.add(
                    names=[
                        "Sphere2_{}".format(i)
                        for i in range(cholesky_large_variance.shape[1])
                    ],
                    lb=[-1] * cholesky_large_variance.shape[1],
                    ub=[1] * cholesky_large_variance.shape[1],
                )  # Independent large variance components

                cplex_model.quadratic_constraints.add(
                    quad_expr=SparseTriple(
                        ind1=indices_sphere2,
                        ind2=indices_sphere2,
                        val=len(indices_sphere2) * [1],
                    ),
                    rhs=1,
                    sense="L",
                    name="unit_normal_high_variance",
                )
            else:
                indices_sphere2 = []  # Balancing matrix dimensions

            concentration_variables = [
                "lnc_{}".format(metabolite.id) for metabolite in self.metabolites
            ]

            # Add the delG constraints
            for reaction in self.reactions:
                if reaction.id in self.Exclude_reactions:
                    continue
                rxn_stoichiometry = reaction.cal_stoichiometric_matrix()
                rxn_stoichiometry = rxn_stoichiometry[np.newaxis, :]

                if len(low_variance_indices) > 0:
                    coefficient_matrix_small_variance = (
                        np.sqrt(chi2_value_small)
                        * rxn_stoichiometry
                        @ metabolite_sphere_small
                    )  # Coefficient array for small variance ellipsoid
                else:
                    coefficient_matrix_small_variance = np.array(())

                if len(high_variance_indices) > 0:
                    coefficient_matrix_large_variance = (
                        np.sqrt(chi2_value_high)
                        * rxn_stoichiometry
                        @ metabolite_sphere_large
                    )  # Coefficient array for large variance ellipsoid
                else:
                    coefficient_matrix_large_variance = np.array(())

                concentration_coefficients = RT * rxn_stoichiometry
                concentration_coefficients[0, proton_indices] = 0

                coefficients_forward = np.hstack(
                    (
                        np.array((1)),
                        -1 * concentration_coefficients.flatten(),
                        -1 * coefficient_matrix_small_variance.flatten(),
                        -1 * coefficient_matrix_large_variance.flatten(),
                    )
                )

                coefficients_reverse = np.hstack(
                    (
                        np.array((1)),
                        concentration_coefficients.flatten(),
                        coefficient_matrix_small_variance.flatten(),
                        coefficient_matrix_large_variance.flatten(),
                    )
                )

                variable_order_forward = (
                    ["dG_{}".format(reaction.forward_variable.name)]
                    + concentration_variables
                    + list(indices_sphere1)
                    + list(indices_sphere2)
                )
                variable_order_reverse = (
                    ["dG_{}".format(reaction.reverse_variable.name)]
                    + concentration_variables
                    + list(indices_sphere1)
                    + list(indices_sphere2)
                )

                rhs = reaction.delG_prime + reaction.delG_transport

                cplex_model.linear_constraints.add(
                    lin_expr=[
                        SparsePair(
                            ind=variable_order_forward,
                            val=coefficients_forward.tolist(),
                        )
                    ],
                    senses=["E"],
                    rhs=[rhs],
                    names=["delG_{}".format(reaction.forward_variable.name)],
                )  # delG constraint for forward reaction

                cplex_model.linear_constraints.add(
                    lin_expr=[
                        SparsePair(
                            ind=variable_order_reverse,
                            val=coefficients_reverse.tolist(),
                        )
                    ],
                    senses=["E"],
                    rhs=[-rhs],
                    names=["delG_{}".format(reaction.reverse_variable.name)],
                )  # delG constraint for reverse reaction

            return cplex_model

        elif self.solver.__class__.__module__ == "optlang.gurobi_interface":
            from gurobipy import GRB, LinExpr

            gurobi_model = self.solver.problem.copy()

            # Remove unnecessary variables and constraints and rebuild  appropriate ones
            remove_vars = [
                var
                for var in gurobi_model.getVars()
                if var.VarName.startswith("component_")
                or var.VarName.startswith("dG_err_")
            ]

            remove_constrs = [
                cons
                for cons in gurobi_model.getConstrs()
                if cons.ConstrName.startswith("delG_")
                or cons.ConstrName.startswith("std_dev_")
            ]

            gurobi_model.remove(remove_constrs + remove_vars)

            # Add sphere variables for smaller set and larger set separately
            if len(low_variance_indices) > 0:
                for i in range(cholesky_small_variance.shape[1]):
                    gurobi_model.addVar(lb=-1, ub=1, name="Sphere1_{}".format(i))

                gurobi_model.update()
                sphere1_variables = [
                    var
                    for var in gurobi_model.getVars()
                    if var.VarName.startswith("Sphere1_")
                ]

                gurobi_model.addQConstr(
                    np.sum(np.square(np.array(sphere1_variables))) <= 1,
                    name="unit_normal_small_variance",
                )
                gurobi_model.update()

            # QC for large variance components
            if len(high_variance_indices) > 0:
                for i in range(cholesky_large_variance.shape[1]):
                    gurobi_model.addVar(lb=-1, ub=1, name="Sphere2_{}".format(i))

                gurobi_model.update()
                sphere2_variables = [
                    var
                    for var in gurobi_model.getVars()
                    if var.VarName.startswith("Sphere2_")
                ]

                gurobi_model.addQConstr(
                    np.sum(np.square(np.array(sphere2_variables))) <= 1,
                    name="unit_normal_high_variance",
                )
                gurobi_model.update()

            # Create a list of metabolite concentration variables
            concentration_variables = []
            for metabolite in self.metabolites:
                varname = "lnc_{}".format(metabolite.id)
                conc_var = gurobi_model.getVarByName(varname)
                concentration_variables.append(conc_var)

            # Add the delG constraints
            for reaction in self.reactions:
                if reaction.id in self.Exclude_reactions:
                    continue
                rxn_stoichiometry = reaction.cal_stoichiometric_matrix()
                rxn_stoichiometry = rxn_stoichiometry[np.newaxis, :]

                if len(low_variance_indices) > 0:
                    coefficient_matrix_small_variance = (
                        np.sqrt(chi2_value_small)
                        * rxn_stoichiometry
                        @ metabolite_sphere_small
                    )  # Coefficient array for small variance ellipsoid
                else:
                    coefficient_matrix_small_variance = np.array(())

                if len(high_variance_indices) > 0:
                    coefficient_matrix_large_variance = (
                        np.sqrt(chi2_value_high)
                        * rxn_stoichiometry
                        @ metabolite_sphere_large
                    )  # Coefficient array for large variance ellipsoid
                else:
                    coefficient_matrix_large_variance = np.array(())

                concentration_coefficients = RT * rxn_stoichiometry
                concentration_coefficients[0, proton_indices] = 0

                coefficients_forward = np.hstack(
                    (
                        -1 * concentration_coefficients.flatten(),
                        -1 * coefficient_matrix_small_variance.flatten(),
                        -1 * coefficient_matrix_large_variance.flatten(),
                    )
                )

                coefficients_reverse = np.hstack(
                    (
                        concentration_coefficients.flatten(),
                        coefficient_matrix_small_variance.flatten(),
                        coefficient_matrix_large_variance.flatten(),
                    )
                )

                variable_order = (
                    concentration_variables + sphere1_variables + sphere2_variables
                )

                delG_err_forward = LinExpr(
                    coefficients_forward.tolist(), variable_order
                )
                delG_err_reverse = LinExpr(
                    coefficients_reverse.tolist(), variable_order
                )

                delG_for_var = gurobi_model.getVarByName(
                    "dG_{}".format(reaction.forward_variable.name)
                )
                delG_rev_var = gurobi_model.getVarByName(
                    "dG_{}".format(reaction.reverse_variable.name)
                )
                rhs = reaction.delG_prime + reaction.delG_transport

                gurobi_model.addConstr(
                    delG_for_var + delG_err_forward,
                    GRB.EQUAL,
                    rhs,
                    name="delG_{}".format(reaction.forward_variable.name),
                )

                gurobi_model.addConstr(
                    delG_rev_var + delG_err_reverse,
                    GRB.EQUAL,
                    -rhs,
                    name="delG_{}".format(reaction.reverse_variable.name),
                )

            gurobi_model.update()

            return gurobi_model

        else:
            raise NotImplementedError("Current solver doesn't support QC")
            logging.error("Current solver doesnt support problesm of type MIQC")

    def calculate_S_matrix(self):
        """Calculates the stoichiometric matrix (metabolites * Reactions)

        Returns:
            Tuple  -- Tuple of reaction order, np.ndarray of stoichiometric matrix
        """

        n_reactions = len(self.reactions)
        n_metabolites = len(self.metabolites)
        S_matrix = np.zeros((2 * n_reactions, n_metabolites))

        reaction_index = 0
        rxn_order = []
        for reaction in self.reactions:
            rxn_order.append(reaction.forward_variable.name)
            rxn_order.append(reaction.reverse_variable.name)
            for metabolite, stoic in iteritems(reaction.metabolites):
                S_matrix[reaction_index, self.metabolites.index(metabolite)] = stoic
                S_matrix[
                    reaction_index + 1, self.metabolites.index(metabolite)
                ] = -stoic
            reaction_index = reaction_index + 2

        S = np.transpose(S_matrix)

        return rxn_order, S

    def concentration_ratio_constraints(self, ratio_metabolites, ratio_lb, ratio_ub):
        """Function to add metabolite concentration ratio constraints to the model. E.g. ratio of redox pairs

        Arguments:
            ratio_metabolites {Tuple} -- Tuple of metabolite names to which we have concentration ratios in the same order

            ratio_lb {List} -- Lower bound of the ratio
            ratio_ub {List} -- Upper bound on the ratio
        """

        for i in range(len(ratio_metabolites)):
            ratio_met1 = self.metabolites.get_by_id(ratio_metabolites[i][0])

            ratio_met2 = self.metabolites.get_by_id(ratio_metabolites[i][1])

            ratio_constraint = self.problem.Constraint(
                1 * ratio_met1.concentration_variable
                - 1 * ratio_met2.concentration_variable,
                lb=ratio_lb[i],
                ub=ratio_ub[i],
            )

            self.add_cons_vars(ratio_constraint)

    def cal_problematic_rxns(self):
        """Reactions which can't be included in thermodynamic analysis
            reactions involving problematic metabolites

        Returns:
            List -- List of reactions excluded, this combined with 'Exclude_list' gives us 'Exclude_reactions'
        """

        problematic_rxns = []
        for met in self.metabolites:
            if met.is_exclude:
                problematic_rxns.append(met.reactions)

        if len(problematic_rxns) > 0:
            problematic_rxns = frozenset.union(*problematic_rxns)
            problems = [i.id for i in problematic_rxns]
            return problems
        else:
            return []

    def export_MIP_matrix(self):
        """Creates matrices structure of the MILP problem. Quadratic constraint is not exported.

        :return: lhs- lhs matrix representing all constraints
                rhs - rhs matrix
                var_names - variable name
                lb, ub- lower, upper bounds of variables
                cons_sense - constraint sense (eg: equal, less equal etc)

        :rtype: Tuple
        """

        rxn_var, S = self.calculate_S_matrix()
        n_mets, rxn_tot = np.shape(S)
        n_rxn_Excl = len(self.Exclude_reactions)
        n_core_rxns = rxn_tot - n_rxn_Excl

        # Create S.v = 0 and expand for other constraints
        mass_balance = np.concatenate(
            (S, np.zeros((n_mets, 2 * n_core_rxns + 2 * n_mets))), axis=1
        )
        rhs_mass_bal = [0] * n_mets
        sense_mass_bal = ["E"] * n_mets

        indicators, delGr, concentration, formation, rhs_delG, sense = (
            [],
            [],
            [],
            [],
            [],
            [],
        )
        # Intialize 6 constraints
        delG_cons_matrix = np.zeros((2 * n_core_rxns, 3 * n_core_rxns + 2 * n_mets))
        for reaction in self.reactions:
            if reaction.id in self.Exclude_reactions:
                continue
            S_rxn = reaction.S_matrix
            rxn_index = rxn_var.index(reaction.id)
            # vi-vmax*zi <= 0
            delG_cons_matrix[rxn_index, rxn_index] = 1
            delG_cons_matrix[rxn_index + 1, rxn_index] = 1
            delG_cons_matrix[rxn_index, rxn_index + n_core_rxns] = -Vmax
            delG_cons_matrix[rxn_index + 1, rxn_index + n_core_rxns] = -Vmax
            # delG - k*zi <= k
            delG_cons_matrix[rxn_index + 2, 2 * n_core_rxns + rxn_index] = 1
            delG_cons_matrix[rxn_index + 3, 2 * n_core_rxns + rxn_index] = 1
            delG_cons_matrix[rxn_index + 2, n_core_rxns + rxn_index] = -K
            delG_cons_matrix[rxn_index + 3, n_core_rxns + rxn_index] = -K
            # delG - s.T RT ln(x) - S.T delGf = 0
            delG_cons_matrix[rxn_index + 4, 2 * n_core_rxns + rxn_index] = 1
            delG_cons_matrix[rxn_index + 5, 2 * n_core_rxns + rxn_index] = 1
            delG_cons_matrix[
                rxn_index + 4, 3 * n_core_rxns : 3 * n_core_rxns + n_mets
            ] = -S_rxn.T
            delG_cons_matrix[
                rxn_index + 5, 3 * n_core_rxns : 3 * n_core_rxns + n_mets
            ] = S_rxn.T
            delG_cons_matrix[
                rxn_index + 4, 3 * n_core_rxns + n_mets : 3 * n_core_rxns + 2 * n_mets
            ] = -S_rxn.T
            delG_cons_matrix[
                rxn_index + 5, 3 * n_core_rxns + n_mets : 3 * n_core_rxns + 2 * n_mets
            ] = S_rxn.T

            indicators.extend(
                [reaction.indicator_forward.name, reaction.indicator_reverse.name]
            )
            delGr.extend([reaction.delG_forward.name, reaction.delG_reverse.name])
            rhs_delG.extend([0, K, 0])
            sense.extend(["L", "L", "E"])

        lb_conc, lb_formation, ub_conc, ub_formation = ([], [], [], [])
        for met in self.metabolites:
            concentration.append(met.concentration_variable.name)
            formation.append(met.compound_variable.name)
            lb_conc.append(met.concentration_variable.lb)
            ub_conc.append(met.concentration_variable.ub)
            lb_formation.append(met.compound_variable.lb)
            ub_formation.append(met.compound_variable.ub)

        var_names = rxn_var + indicators + delGr + concentration + formation

        lhs = np.concatenate((mass_balance, delG_cons_matrix), axis=0)
        rhs = rhs_mass_bal + rhs_delG
        cons_sense = sense_mass_bal + sense
        lb = (
            [-1000] * rxn_tot
            + [0] * n_core_rxns
            + [-1e5] * n_core_rxns
            + lb_conc
            + lb_formation
        )
        ub = (
            [1000] * rxn_tot
            + [1] * n_core_rxns
            + [1e5] * n_core_rxns
            + ub_conc
            + ub_formation
        )

        return (lhs, rhs, var_names, np.array(lb), np.array(ub), cons_sense)
