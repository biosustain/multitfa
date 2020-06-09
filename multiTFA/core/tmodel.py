from .reaction import thermo_reaction
import numpy as np
from ..util.thermo_constants import Vmax, K, RT
from ..util.constraints import (
    directionality,
    delG_indicator,
    concentration_exp,
    metabolite_variables,
    reaction_variables,
    met_exp_qp,
    bounds_ellipsoid,
    quad_constraint,
)
from copy import deepcopy, copy
from ..util.dGf_calculation import calculate_dGf, cholesky_decomposition
from .compound import Thermo_met
from numpy import array, dot, sqrt, diag, isnan
from warnings import warn
from six import iteritems
from cobra import Model
from .solution import get_solution
from cobra.core.dictlist import DictList
import optlang


class tmodel(Model):
    def __init__(
        self,
        model,
        Kegg_map={},
        Exclude_list=[],
        pH_I_T_dict={},
        concentration_dict={"min": {}, "max": {}},
        tolerance_integral=1e-9,
        del_psi_dict={},
        solve_method="box",
        debug=False,
    ):

        """ Class representation of tMFA model, dependeds on cobra model class.
        
        Arguments:
            model {Cobra model object} -- Cobra model 
        
        Keyword Arguments:
            Kegg_map {dict} -- Dictionary of metabolite id to Kegg/seed identifiers of all metabolites in the model (default: {{}})
            
            Exclude_list {list} -- Reaction ids that needs to be excluded from thermodynamic analysis (e.g: exchange/sinks)  (default: {[]})
            
            pH_I_T_dict {dict} -- Dictionary of pH, ionic strength of different compartments, parameters for all the compartments needs to be specified  (default: {{}})
            
            concentration_dict {dict} -- Dictionary of min/max concentrations of metabolites where available (default: {{'min':{},'max':{}}})
            
            tolerance_integral {float} -- Integral tolerance of the solver (We recemmond 1e-9 for tmfa problems) (default: {1e-9})
            
            del_psi_dict {dict} -- membrane potential dictionary for all the compartments in the model (default: {{}})
            
            debug {bool} -- Debugging flag (default: {False})
        """

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
                Kegg_map=Kegg_map,
                concentration_dict=concentration_dict,
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
                Kegg_map=Kegg_map,
                pH_I_T_dict=pH_I_T_dict,
                del_psi_dict=del_psi_dict,
            )
            self.reactions.append(new_reaction)

        try:
            self._solver = deepcopy(model.solver)
            # Cplex has an issue with deep copies
        except Exception:  # pragma: no cover
            self._solver = copy(model.solver)  # pragma: no cover

        self.Kegg_map = Kegg_map
        self.Exclude_list = Exclude_list
        self.pH_I_T_dict = pH_I_T_dict
        self.concentration_dict = concentration_dict
        self.del_psi_dict = del_psi_dict
        self.covariance_matrix()
        self.solver.configuration.tolerances.integrality = tolerance_integral
        self.Exclude_reactions = list(set(Exclude_list + self.problematic_rxns))
        self.update_thermo_variables()
        self.update()

        if not (
            optlang.available_solvers["GUROBI"] or optlang.available_solvers["CPLEX"]
        ):
            warn(
                'Quadratic constraints are only supported with GUROBI or CPLEX, Please use "box" or "sampling method"'
            )
            self.solve_method = "box"
        else:
            if self.solver.__class__.__module__ == "optlang.gurobi_interface":
                self.solve_method = "qc"
                self.solver.problem.update()
                self.gurobi_interface = self.solver.problem
                self.MIQP()
                self.gurobi_interface.update()
                self.gurobi_interface.write("gurobi_problem.lp")
            elif self.solver.__class__.__module__ == "optlang.cplex_interface":
                self.solve_method = "qc"
                self.cplex_interface = self.solver.problem.copy()
            else:
                self.solve_method = "box"
                warn(
                    "GUROBI/CPLEX not found. Quadratic constraints not supported for this model"
                )

    @property
    def cholskey_matrix(self):

        """Calculates cholesky matrix (square root of covariance matrix of compounds)
        
        Returns:
            np.ndarray 
        """

        std_dg, cov_dg = calculate_dGf(self.metabolites, self.Kegg_map)
        chol_matrix = cholesky_decomposition(std_dg, cov_dg)
        return chol_matrix

    @property
    def problem_metabolites(self):

        """ Metabolites for which we can't calculate the Gibbs free energy of formation using component contribution method

        Metabolites are considered problematic metabolites if whole row in cholesky matrix is zero
        
        Returns:
            List -- List of problematic metabolites
        """

        problematic_metabolites = []

        for met in self.metabolites:
            met_index = self.metabolites.index(met)
            if met.Kegg_id in ["C00080", "cpd00067"]:
                continue
            if np.count_nonzero(self.cholskey_matrix[:, met_index]) == 0 or isnan(
                met.delG_f
            ):  # or sqrt(diag(self.cov_dG)[self.metabolites.index(met)]) > 100:
                problematic_metabolites.append(met)

        return problematic_metabolites

    @property
    def problematic_rxns(self):

        """ Reactions which can't be included in thermodynamic analysis
            reactions involving problematic metabolites
        
        Returns:
            List -- List of reactions excluded, this combined with 'Exclude_list' gives us 'Exclude_reactions'
        """

        problematic_rxns = []
        for met in self.problem_metabolites:
            problematic_rxns.append(met.reactions)

        if len(problematic_rxns) > 0:
            problematic_rxns = frozenset.union(*problematic_rxns)
            problems = [i.id for i in problematic_rxns]
            return problems
        else:
            return []

    def covariance_matrix(self):
        """calculates the covariance matrix of Gibbs free energy of the metabolites
        
        Returns:
            std_dg [np.ndarray] -- std delg of all metabolites
            cov_dg [np.ndarray] -- covariance matrix
        """

        self.std_dG, self.cov_dG = calculate_dGf(self.metabolites, self.Kegg_map)

    def update_thermo_variables(self):

        """ Generates reaction and metabolite variables required for thermodynamic analysis and adds to the model

        Variables generated --

        Metabolite concentration variable
        Metabolite confidence variable
        Reaction delG variable
        Reaction flux indicator variable
        """

        # metabolite concentration, Ci and metabolite variables
        conc_variables = []
        met_variables = []

        for metabolite in self.metabolites:
            conc_var, met_variable = metabolite_variables(metabolite)
            conc_variables.append(conc_var)
            met_variables.append(met_variable)
        self.add_cons_vars(conc_variables + met_variables)

        # Now add reaction variables and generate remaining constraints
        for rxn in self.reactions:
            if rxn.id in self.Exclude_reactions:
                continue
            reaction_vars = reaction_variables(rxn)
            self.add_cons_vars(reaction_vars)

    def calculate_S_matrix(self):

        """ Calculates the stoichiometric matrix (metabolites * 2*Reactions)
        Each reaction has two column entries to represent forward and reverse.
        
        Note: S matrix is calculated only for the reactions involved in thermodynamic analysis.

        Returns:
            Tuple  -- Tuple of reaction order, np.ndarray of stoichiometric matrix
        """

        core_rxn = [
            rxn.id for rxn in self.reactions if rxn.id not in self.Exclude_reactions
        ]

        n_reactions = len(core_rxn)
        n_metabolites = len(self.metabolites)
        S_matrix = np.zeros((2 * n_reactions, n_metabolites))

        reaction_index = 0
        rxn_order = []
        for rxn in core_rxn:
            reaction = self.reactions.get_by_id(rxn)
            rxn_order.extend(
                [reaction.forward_variable.name, reaction.reverse_variable.name]
            )
            for metabolite, stoic in iteritems(reaction.metabolites):
                S_matrix[reaction_index, self.metabolites.index(metabolite)] = stoic
                S_matrix[
                    reaction_index + 1, self.metabolites.index(metabolite)
                ] = -stoic
            reaction_index = reaction_index + 2

        S = np.transpose(S_matrix)

        return rxn_order, S

    def _generate_constraints(self):

        """ Generates thermodynamic constraints for the model. See util/constraints.py for detailed explanation of constraints

        Vi - Vmax * Zi <= 0
        delGr - K + K * Zi <= 0
        delGr - RT * S.T * ln(x) - S.T @ cholesky @ Zf - S.T @ delGf - delGtransport = 0

        
        Returns:
            List -- List of themrodynamic constraints
        """

        rxn_constraints = []
        # Now add reaction variables and generate remaining constraints
        for rxn in self.reactions:
            if rxn.id in self.Exclude_reactions:
                continue

            # Direactionality constraint
            dir_f, dir_r = directionality(rxn)
            ind_f, ind_r = delG_indicator(rxn)

            # delG constraint
            concentration_term = concentration_exp(rxn)

            met_term = met_exp_qp(rxn)

            # temporarily removing ci expression and adding met_exp_qp for adding qcp. later we can check which is default and adjust accordingly

            lhs_forward = (
                rxn.delG_forward - RT * concentration_term - met_term
            )  # - ci_term
            lhs_reverse = (
                rxn.delG_reverse + RT * concentration_term + met_term
            )  # ci_term
            rhs = rxn.transform + rxn.transport_delG

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
            rxn_constraints.extend([dir_f, dir_r, ind_f, ind_r, delG_f, delG_r])

        return rxn_constraints

    def update(self):

        """ Adds the generated thermo constaints to  model. Checks for duplication 
        """
        constraints = self._generate_constraints()
        for cons in constraints:
            if cons.name not in self.constraints:
                self.add_cons_vars([cons])
            else:
                warn("Constraint {} already exists in the model".format(cons.name))
                self.solver.remove(cons.name)
                self.add_cons_vars([cons])

    def optimize(self):

        if self.solve_method == "box":
            solution = self.solver.optimize()
        elif self.solve_method == "qc":
            if self.solver.__class__.__module__ == "optlang.gurobi_interface":
                solution = self.gurobi_interface.optimize()
            elif self.solver.__class__.__module__ == "optlang.cplex_interface":
                solution = self.cplex_interface.solve()

        return solution

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

    def MIQP(self):

        if self.solver.__class__.__module__ == "optlang.gurobi_interface":

            # solver_interface = self.gurobi_interface

            # Get metabolite variable from gurobi interface
            metid_vars_dict = {}
            for var in self.gurobi_interface.getVars():
                if var.VarName.startswith("met_"):
                    metid_vars_dict[var.VarName[4:]] = var

        elif self.solver.__class__.__module__ == "optlang.cplex_interface":
            pass

        else:
            raise NotImplementedError(
                "Current solver does not support quadratic constraints, please use Gurobi or Cplex"
            )

        # Problem metabolites, if met.delGf == 0 or cholesky row is zeros then delete them
        delete_met, cov_mets, cov_met_inds = [], [], []

        for met in self.metabolites:
            if met.delG_f == 0 or isnan(met.delG_f):
                delete_met.append(met)
            else:
                cov_met_inds.append(self.metabolites.index(met))
                cov_mets.append(met)

        cov_dg = self.cov_dG
        # Pick indices of non zero non nan metabolites
        cov_dG = cov_dg[:, cov_met_inds]
        cov_dg = cov_dG[cov_met_inds, :]

        lhs, rhs = quad_constraint(
            cov_dg, cov_mets, metid_vars_dict
        )  # Calculate lhs, rhs for quadratic constraints

        # Calculate ellipsoid box bounds and set to variables
        bounds = bounds_ellipsoid(cov_dg)  # Check for posdef cov_dg
        for met in self.metabolites:
            if met in delete_met:
                continue
            metid_vars_dict[met.id].LB = met.delG_f - bounds[cov_mets.index(met)]
            metid_vars_dict[met.id].UB = met.delG_f + bounds[cov_mets.index(met)]

        # self.gurobi_interface.addConstr(lhs <= rhs, "qp_constraint")
        self.gurobi_interface.update()

        # self.gurobi_interface.write("QC_problem.lp")

        # return solver_interface

    def lp_matrices_matlab(self):
        """[Variable order:
            flux,excluded reactions, binary, delG, concentration, significance]
        """
        _, S = self.calculate_S_matrix()
        core_rxn = [
            rxn for rxn in self.reactions if rxn.id not in self.Exclude_reactions
        ]
        m, n = np.shape(S)

        S_Ex = np.zeros((len(self.Exclude_reactions), len(self.metabolites)))
        i = 0
        for rxn_id in self.Exclude_reactions:
            reaction = self.reactions.get_by_id(rxn_id)
            for metabolite, stoic in iteritems(reaction.metabolites):
                S_Ex[i, self.metabolites.index(metabolite.id)] = stoic
            i = i + 1
        S_Exclude = S_Ex.T

        _, y = np.shape(S_Exclude)
        mass_matrix = np.concatenate(
            (S, S_Exclude, np.zeros((m, 2 * m + 2 * n))), axis=1
        )
        flux_indicator = np.concatenate(
            (np.eye(n), np.zeros((n, y)), -1000 * np.eye(n), np.zeros((n, 2 * m + n))),
            axis=1,
        )
        delG_indicator_matrix = np.concatenate(
            (
                np.zeros((n, n)),
                np.zeros((n, y)),
                K * np.eye(n),
                np.eye(n),
                np.zeros((n, 2 * m)),
            ),
            axis=1,
        )
        delG_matrix = np.concatenate(
            (
                np.zeros((n, n)),
                np.zeros((n, y)),
                np.zeros((n, n)),
                np.eye(n),
                -RT * S.T,
                -S.T @ self.cholskey_matrix,
            ),
            axis=1,
        )

        lhs_matrix = np.concatenate(
            (mass_matrix, flux_indicator, delG_indicator_matrix, delG_matrix), axis=0
        )

        delG_rhs = []
        for rxn in core_rxn:
            delG_rhs.extend([rxn.calculate_rxn_delG(), -1 * rxn.calculate_rxn_delG()])

        rhs_matrix = []

        rhs_matrix.extend([0] * (m + n))
        rhs_matrix.extend([K] * n)
        rhs_matrix.extend(delG_rhs)
        rhs_matrix = np.array(rhs_matrix)

        core_rxn_ids = [rxn.id for rxn in core_rxn]

        binaries = []
        delGs = []
        flux_var = []
        for rxn_id in core_rxn_ids:
            flux_var.append(rxn_id)
            binaries.append(rxn_id + "_bi")
            delGs.append("delG_" + rxn_id)
            flux_var.append(rxn_id + "_rev")
            binaries.append(rxn_id + "_rev_bi")
            delGs.append("delG_" + rxn_id + "_rev")

        log_var = []
        sig_var = []
        for metid in self.metabolites:
            log_var.append("log_" + metid.id)
            sig_var.append("z_f_" + metid.id)
        var_list = (
            flux_var + self.Exclude_reactions + binaries + delGs + log_var + sig_var
        )

        core_flux_bounds_lb = []
        core_flux_bounds_ub = []
        for rxn in core_rxn:
            core_flux_bounds_lb.extend([rxn.lower_bound, -rxn.upper_bound])
            core_flux_bounds_ub.extend([rxn.upper_bound, -rxn.lower_bound])

        exclude_lb = []
        exclude_ub = []
        for rxn_id in self.Exclude_reactions:
            rxn = self.reactions.get_by_id(rxn_id)
            exclude_lb.append(rxn.lower_bound)
            exclude_ub.append(rxn.upper_bound)

        lower_bounds = (
            core_flux_bounds_lb
            + exclude_lb
            + [0] * n
            + [-1000] * n
            + [1e-5] * m
            + [-1.96] * m
        )
        upper_bounds = (
            core_flux_bounds_ub
            + exclude_ub
            + [1] * n
            + [1000] * n
            + [0.01] * m
            + [1.96] * m
        )

        return (
            lhs_matrix,
            rhs_matrix,
            var_list,
            np.array(lower_bounds),
            np.array(upper_bounds),
        )
