from .reaction import thermo_reaction
import numpy as np
from ..util.thermo_constants import *
from ..util.constraints import (
    directionality,
    delG_indicator,
    concentration_exp,
    formation_exp,
    bounds_ellipsoid,
    quad_constraint,
    delG_constraint_expression,
)
from copy import deepcopy, copy
from ..util.dGf_calculation import calculate_dGf, cholesky_decomposition
from ..util.posdef import isPD, nearestPD
from ..util.util_func import Exclude_quadratic, correlated_pairs, quadratic_matrices
from .compound import Thermo_met
from warnings import warn
from six import iteritems
from cobra import Model
from .solution import get_solution, get_legacy_solution
from cobra.core.dictlist import DictList
import optlang
import itertools
import os
from random import choices
import string
from scipy import stats, linalg
from equilibrator_api import ComponentContribution, Q_

# present_dir = os.path.normpath(os.path.dirname(os.path.abspath(__file__)))
cwd = os.getcwd()

import logging

logs_dir = cwd + os.sep + "logs"
if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)

logging.basicConfig(
    filename=logs_dir + os.sep + "thermo_model.log", level=logging.DEBUG
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
            new_met = Thermo_met(metabolite=metabolite, updated_model=self,)
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
            new_reaction = thermo_reaction(cobra_rxn=reaction, updated_model=self,)
            self.reactions.append(new_reaction)

        try:
            self._solver = deepcopy(model.solver)
            self._miqc_solver = deepcopy(model.solver)
            # Cplex has an issue with deep copies
        except Exception:  # pragma: no cover
            self._solver = copy(model.solver)  # pragma: no cover
            self._miqc_solver = copy(
                model.solver
            )  # creating this for adding different delG constraint to copy later to gurobi and cplex interfaces

        self.Exclude_list = Exclude_list
        self.solver.configuration.tolerances.integrality = tolerance_integral

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
        """ Cplex interface to support QC

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
    def covariance_dG(self):
        """Covariance matrix of the reactants and groups from component contribution method. Borrowed from equilibrator-api. This covariance matrix is used for sampling from multi-dimensional space. 
        Note: This matrix is rank-deficient matrix

        Returns
        -------
        np.ndarray
            covariance matrix of the components
        """
        return covariance

    @property
    def problem_metabolites(self):
        """ Metabolites for which we can't calculate the Gibbs free energy of formation using component contribution method. If the metabolite is not covered by reactant/group contribution method, then we have to write the metabolite and corresponding reactions from tMFA analysis or have to lump the reactions.
        
        Returns:
            List
                List of metabolite not covered by component contribution
        """
        problematic_metabolites = []
        for met in self.metabolites:
            if met.is_proton:
                continue
            if not met.compound_vector:
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
            self._Exclude_reactions = self.Exclude_list + self.problematic_rxns
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
    def sphere_variables(self):
        """List of independent variables that represent the whole solution space. The covariance matrix is rank-deficient. Dependent variables can be written in linear combination of independent variables using cholesky matrix.

        Returns
        -------
        List
            List of optlang interface variables
        """
        try:
            return self._sphere_variables
        except AttributeError:
            self._sphere_variables = [
                var
                for var in self._miqc_solver.variables
                if var.name.startswith("Sphere_")
            ]
            return self._sphere_variables

    @property
    def correlated_metabolites(self):
        """Metabolites that are highly correlated thermodynamically, for example they share many common groups or co-occur all the time in training data.

        Returns
        -------
        Dict
            Dictionary of correlated metabolites
        """
        try:
            return self._correlated_metabolites
        except AttributeError:
            self._correlated_metabolites = correlated_pairs(self)
            return self._correlated_metabolites

    def update_thermo_variables(self):
        """ Generates reaction and metabolite variables required for thermodynamic analysis and adds to the model. We use two different methods to solve the tMFA problem. Traditional 'box' method employs MILP problem where components are allowed to vary between some s.d from mean. The other method uses MIQCP structure to use covariance matrix to capture covariance. Two methods share some common variables, where as MIQCP method requires independent variables to sample from original solution space.

        Common Variables: delG_reaction, indicator_reaction, concentration_metabolite
        MILP variables: metabolite error 
        MIQCP variables: independent component variables

        """
        # Add metabolite concentration variable
        conc_varibales = []
        for metabolite in self.metabolites:
            conc_variable = self.problem.Variable(
                "lnc_{}".format(metabolite.id),
                lb=np.log(metabolite.concentration_min),
                ub=np.log(metabolite.concentration_max),
            )
            conc_varibales.append(conc_variable)
        self.add_cons_vars(conc_varibales)
        self._miqc_solver.add(conc_varibales)

        # Adding the thermo variables for reactions, delG_reaction and indicator (binary)
        for rxn in self.reactions:
            if rxn.id in self.Exclude_reactions:
                continue
            delG_forward = self.problem.Variable(
                "dG_{}".format(rxn.forward_variable.name), lb=-1e5, ub=1e5
            )

            delG_reverse = self.problem.Variable(
                "dG_{}".format(rxn.reverse_variable.name), lb=-1e5, ub=1e5
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
            self.add_cons_vars(
                [delG_forward, delG_reverse, indicator_forward, indicator_reverse]
            )
            self._miqc_solver.add(
                [delG_forward, delG_reverse, indicator_forward, indicator_reverse]
            )  # for the MIQC constraints

        # Add group variables for BOX method (including dependent & independent)
        component_variables = np.array(
            [
                self.problem.Variable(
                    "component_{}".format(i),
                    lb=-np.sqrt(covariance[i, i]),
                    ub=np.sqrt(covariance[i, i]),
                )
                for i in range(len(covariance))
            ]
        )
        self.add_cons_vars(component_variables.tolist())

        # Add variables for the independent components, number will be equal to rank of component covariance matrix
        something = 669  # replace with rank of matrix
        sphere_vars = np.array(
            [
                self.problem.Variable("Sphere_{}".format(i), lb=-1, ub=1)
                for i in range(something)
            ]
        )
        self._miqc_solver.add(
            sphere_vars.tolist()
        )  # These variable are not used for box method, only for the MIQCP

    def _generate_constraints(self):
        """ Generates thermodynamic constraints for the model. See util/constraints.py for detailed explanation of constraints

        Vi - Vmax * Zi <= 0
        delGr - K + K * Zi <= 0
        delGr - RT * S.T * ln(x) - S.T @ delGf - delGtransport = 0

        
        Returns:
            List -- List of themrodynamic constraints
        """
        # First check if thermovariables are added to the model, if not update. Just being conservative on when to update variables
        if len(self.variables) <= 2.5 * len(self.reactions):
            self.update_thermo_variables()

        # Get the group variables for box method
        component_variables = [
            var for var in self.variables if var.name.startswith("component_")
        ]

        # Now we create common constraints between MIQC method and MIP methods and later add the delG constraint separately

        rxn_common_constraints = []
        delG_constraint_box = []
        delG_constraint_QC = []
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

            rxn_common_constraints.extend([dir_f, dir_r, ind_f, ind_r])

            # Create two different constraints for box method and MIQC method

            # delG constraint for box
            concentration_term = concentration_exp(rxn)
            met_term = formation_exp(rxn, component_variables=component_variables)

            lhs_forward = rxn.delG_forward - RT * concentration_term - met_term
            lhs_reverse = rxn.delG_reverse + RT * concentration_term + met_term
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
            delG_constraint_box.extend([delG_f, delG_r])

            # delG constraint for QC problem
            forward_qc, reverse_qc = delG_constraint_expression(rxn)
            delG_f_qc = self.problem.Constraint(
                forward_qc,
                lb=rhs,
                ub=rhs,
                name="delG_{}".format(rxn.forward_variable.name),
            )
            delG_r_qc = self.problem.Constraint(
                reverse_qc,
                lb=-rhs,
                ub=-rhs,
                name="delG_{}".format(rxn.reverse_variable.name),
            )
            delG_constraint_QC.extend([delG_f_qc, delG_r_qc])

        return (
            rxn_common_constraints + delG_constraint_box,
            rxn_common_constraints + delG_constraint_QC,
        )

    def update(self):
        """ Adds the generated thermo constaints to  model. Checks for duplication 
        """
        box_constraints, QC_constraints = self._generate_constraints()
        for cons in box_constraints:
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

        for cons in QC_constraints:
            if cons.name not in self._miqc_solver.constraints:
                self._miqc_solver.add([cons])
                logging.debug(
                    "Constraint {} added to the MIQC problem".format(cons.name)
                )
            else:
                logging.warning(
                    "Constraint {} already in the MIQC problem, removing previous entry".format(
                        cons.name
                    )
                )
                self._miqc_solver.remove(cons.name)
                self._miqc_solver.add([cons])

    def optimize(self, solve_method="MIQC", raise_error=False):
        """ solves the model with given constraints. By default, we try to solve the model with quadratic constraints. Note: Quadratic constraints are supported by Gurobi/Cplex currently. if either of two solvers are not found, one can solve 'box' type MILP problem.

        :param solve_method: Method to solve the problem, defaults to "MIQC", any other string input leades to solving with box MILP  
        :type solve_method: str, optional
        :param raise_error: , defaults to False
        :type raise_error: bool, optional
        :return: returns solution object
        :rtype: solution object (refer to Solution class)
        """

        if solve_method == "MIQC":
            if not (
                optlang.available_solvers["GUROBI"]
                or optlang.available_solvers["CPLEX"]
            ):
                logging.warning(
                    "GUROBI/CPLEX not available, Quadratic constraints are not supported by current solver"
                )
                return

            if self.solver.__class__.__module__ == "optlang.gurobi_interface":
                self.gurobi_interface.optimize()
                solution = get_legacy_solution(self, solver="gurobi")

                return solution

            elif self.solver.__class__.__module__ == "optlang.cplex_interface":
                solution = self.cplex_interface.solve()
                solution = get_legacy_solution(self, solver="cplex")

                return solution

        else:
            self.slim_optimize()
            solution = get_solution(self, raise_error=raise_error)

            return solution

    def Quadratic_constraint(self):
        """ Adds Quadratic constraint to the model's Gurobi/Cplex Interface. 
        (x-mu).T @ inv(cov) @ (x-mu) <= chi-square
        Note: This one creates one ellipsoidal constraint for all the metabolites that has non zero or non 'nan' formation energy, irrespective of the magnitude of variance. if the model is infeasible after adding this constraint, refer to util_func.py, find_correlated metabolites to add different ellipsoidal constraints to high variance and normal compounds to avoid possible numerical issues.
        
        Unable to retrieve quadratic constraints in Gurobi model, can see the QC when printed.

        :raises NotImplementedError: Implemented only for Gurobi/Cplex interfaces.
        :return: [description]
        :rtype: [type]
        """

        if self.solver.__class__.__module__ == "optlang.gurobi_interface":

            solver_interface = self._miqc_solver.copy()

            sphere_variables = [
                var
                for var in solver_interface.getVars()
                if var.VarName.startswith("Sphere_")
            ]

            solver_interface.addQConstr(
                np.sum(np.square(np.array(sphere_variables))) <= 1, name="unit_normal"
            )

            solver_interface.update()

            return solver_interface

        elif self.solver.__class__.__module__ == "optlang.cplex_interface":

            from cplex import Cplex, SparseTriple

            tmp_dir = present_dir + os.sep + os.pardir + os.sep + "tmp"

            if not os.path.exists(tmp_dir):
                os.makedirs(tmp_dir)

            rand_str = "".join(choices(string.ascii_lowercase + string.digits, k=6))
            # write cplex model to mps file and re read
            self.solver.problem.write(tmp_dir + os.sep + rand_str + ".mps")
            # Instantiate Cplex model
            cplex_model = Cplex()
            cplex_model.read(tmp_dir + os.sep + rand_str + ".mps")

            Sphere_var_names = [
                var_name
                for var_name in cplex_model.variables.get_names()
                if var_name.startswith("Sphere_")
            ]

            # Add the Sphere constraint
            cplex_model.quadratic_constraints.add(
                quad_expr=SparseTriple(
                    ind1=Sphere_var_names,
                    ind2=Sphere_var_names,
                    val=len(Sphere_var_names) * [1],
                ),
                rhs=1,
                name="unit_normal",
            )

            return cplex_model

        else:
            raise NotImplementedError("Current solver doesn't support QC")
            logging.error("Current solver doesnt support problesm of type MIQC")

    def calculate_S_matrix(self):
        """ Calculates the stoichiometric matrix (metabolites * Reactions)

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
        """ Function to add metabolite concentration ratio constraints to the model. E.g. ratio of redox pairs
        
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
        """ Reactions which can't be included in thermodynamic analysis
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
        """ Creates matrices structure of the MILP problem. Quadratic constraint is not exported.

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

