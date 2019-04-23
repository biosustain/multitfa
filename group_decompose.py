"""module for decomposing molecules into chemical groups."""
# The MIT License (MIT)
#
# Copyright (c) 2013 The Weizmann Institute of Science.
# Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
# Copyright (c) 2018 Institute for Molecular Systems Biology,
# ETH Zurich, Switzerland.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
import itertools
import logging
from typing import Dict, Iterator, List, Set, Tuple

import pandas as pd

from groups import (
    DEFAULT_GROUPS_DATA, Group, GroupsData, GroupVector,
    MalformedGroupDefinitionError)
from molecule import Molecule


class GroupDecomposition(object):
    """Class representing the group decomposition of a molecule.

    This class is similar to GroupVector, but contains more specific
    information about the decomposed molecule. For example, rather than just
    counting the number of instances of each group, GroupDecomposition also
    stored the list of sets of atoms in each of these instances.
    """

    def __init__(self,
                 groups_data: GroupsData,
                 mol: Molecule,
                 groups: List[Tuple[Group, List[Set[int]]]],
                 unassigned_nodes: Set[int]):
        """Create a GroupDecomposition object."""
        self.groups_data = groups_data
        self.mol = mol
        self.groups = groups
        self.unassigned_nodes = unassigned_nodes

    def ToDataFrames(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Return the decomposition as a tabular string.

        :return: a Pandas DataFrame
        """
        gr_data = []
        for group, node_sets in self.groups:
            for n_set in node_sets:
                gr_data.append((
                    group.name,
                    group.hydrogens or 0,
                    group.charge or 0,
                    group.nMg or 0,
                    ','.join([str(i) for i in n_set])
                ))
        gr_df = pd.DataFrame(
            data=gr_data,
            columns=['name', 'nH', 'charge', 'nMg', 'atoms']
        )

        unassigned_data = []
        if self.unassigned_nodes:
            all_atoms = self.mol.GetAtoms()
            for i in self.unassigned_nodes:
                a = all_atoms[i]
                unassigned_data.append((
                    i,
                    a.GetAtomicNum(),
                    Molecule.GetSymbol(a.GetAtomicNum()),
                    a.GetHvyValence(),
                    a.GetFormalCharge()
                ))

        unassigned_df = pd.DataFrame(
            data=unassigned_data,
            columns=['index', 'an', 'el', 'valence', 'charge']
        )

        return gr_df, unassigned_df

    def __str__(self) -> str:
        """Convert the groups to a string."""
        group_strs = []
        for gr, node_sets in self.NonEmptyGroups():
            if gr.hydrogens is None and gr.charge is None and gr.nMg is None:
                group_strs.append('%s x %d' % (gr.name, len(node_sets)))
            else:
                group_strs.append('%s [H%d %d %d] x %d' %
                                  (gr.name, gr.hydrogens, gr.charge,
                                   gr.nMg, len(node_sets)))
        return " | ".join(group_strs)

    def __len__(self) -> int:
        """Get the total number of groups."""
        counter = 0
        for _group, node_sets in self.NonEmptyGroups():
            counter += len(node_sets)
        return counter

    def AsVector(self) -> GroupVector:
        """Return the group in vector format.

        Note: self.groups contains an entry for *all possible* groups, which is
        why this function returns consistent values for all compounds.
        """
        return GroupVector(
            self.groups_data,
            [len(node_sets) for _, node_sets in self.groups] + [1]
        )

    def NonEmptyGroups(self) -> Iterator[Tuple[Group, list]]:
        """Generate non-empty groups."""
        for group, node_sets in self.groups:
            if node_sets:
                yield group, node_sets

    def UnassignedAtoms(self):
        """Generate unassigned atoms."""
        for i in self.unassigned_nodes:
            yield self.mol.GetAtoms()[i], i

    def NetCharge(self) -> float:
        """Return the net charge."""
        return self.AsVector().NetCharge()

    def Hydrogens(self) -> float:
        """Return the number of hydrogens."""
        return self.AsVector().Hydrogens()

    def Magnesiums(self) -> float:
        """Return the number of Mg2+ ions."""
        return self.AsVector().Magnesiums()

    def CountGroups(self) -> int:
        """Return the total number of groups in the decomposition."""
        return sum([len(gdata[-1]) for gdata in self.groups])

    @staticmethod
    def distribute(total: int, num_slots: int) -> List[List[int]]:
        """Distribution 'total' balls in 'num_slots' slots.

        Example: distribute(3, 2) = [[0, 3], [1, 2], [2, 1], [3, 0]]

        :return: a list of all distinct options for distributing the balls
        """
        if num_slots == 1:
            return [[total]]

        if total == 0:
            return [[0] * num_slots]

        all_options = []
        for i in range(total + 1):
            for opt in GroupDecomposition.distribute(total - i, num_slots - 1):
                all_options.append([i] + opt)

        return all_options

    @staticmethod
    def multi_distribute(
        total_slots_pairs: List[Tuple[int, int]]
    ) -> List[List[int]]:
        """Distribution 'total' balls in 'num_slots' slots multiple times.

        similar to distribute, but with more constraints on the sub-totals
        in each group of slots. Every pair in the input list represents
        the subtotal of the number of balls and the number of available balls
        for them. The total of the numbers in these slots will be equal to
        the subtotal.

        Example:
            multi_distribute([(1, 2), (2, 2)]) =
            [[0, 1, 0, 2], [0, 1, 1, 1], [0, 1, 2, 0],
             [1, 0, 0, 2], [1, 0, 1, 1], [1, 0, 2, 0]]

        in words, the subtotal of the two first slots must be 1, and
        the subtotal of the two last slots must be 2.
        """
        multilist_of_options = []
        for (total, num_slots) in total_slots_pairs:
            multilist_of_options.append(
                GroupDecomposition.distribute(total, num_slots)
            )

        return [sum(x) for x in itertools.product(*multilist_of_options)]

    def PseudoisomerVectors(self) -> List[GroupVector]:
        """Return a list of group vectors, one per pseudo-isomer."""
        if not self.CountGroups():
            logging.debug('No groups in this decomposition, not calculating '
                          'pseudoisomers.')
            return []

        # A map from each group name to its indices in the group vector.
        # Note that some groups appear more than once (since they can have
        # multiple protonation levels).
        group_name_to_index = {}

        # 'group_name_to_count' is a map from each group name to its number
        # of appearances in 'mol'
        group_name_to_count = {}
        for i, gdata in enumerate(self.groups):
            group, node_sets = gdata
            group_name_to_index.setdefault(group.name, []).append(i)
            group_name_to_count[group.name] = \
                group_name_to_count.get(group.name, 0) + len(node_sets)

        index_vector = []  # maps the new indices to the original ones that
        # are used in groupvec

        # A list of per-group pairs (count, # possible protonation levels).
        total_slots_pairs = []

        for group_name, groupvec_indices in group_name_to_index.items():
            index_vector += groupvec_indices
            total_slots_pairs.append((group_name_to_count[group_name],
                                      len(groupvec_indices)))

        # generate all possible assignments of protonations. Each group can
        # appear several times, and we can assign a different protonation
        # level to each of the instances.
        groupvec_list = []
        for assignment in self.multi_distribute(total_slots_pairs):
            v = [0] * len(index_vector)
            for i in range(len(v)):
                v[index_vector[i]] = assignment[i]
            v += [1]  # add 1 for the 'origin' group
            groupvec_list.append(GroupVector(self.groups_data, v))
        return groupvec_list

    # Various properties
    nonempty_groups = property(NonEmptyGroups)
    unassigned_atoms = property(UnassignedAtoms)
    hydrogens = property(Hydrogens)
    net_charge = property(NetCharge)
    magnesiums = property(Magnesiums)
    group_count = property(CountGroups)


class GroupDecompositionError(Exception):
    """General decomposition error."""

    def __init__(self, msg: str, decomposition: GroupDecomposition = None):
        """Create a GroupDecompositionError."""
        Exception.__init__(self, msg)
        self.decomposition = decomposition

    def __str__(self):
        """Return the exception error message."""
        return Exception.__str__(self)

    def GetDebugTable(self) -> str:
        """Get the decomposition info for debugging purposes."""
        if self.decomposition is not None:
            gr_df, unassigned_df = self.decomposition.ToDataFrames()
            s = "Groups:\n" + str(gr_df) + "\nUnassigned nodes:\n" + str(
                unassigned_df)
            return s
        else:
            return ''


class GroupDecomposer(object):
    """Decomposes compounds into their constituent groups."""

    def __init__(self, groups_data: GroupsData = None):
        """Construct a GroupDecomposer.

        Args:
            groups_data: a GroupsData object.
        """
        if groups_data is None:
            self.groups_data = DEFAULT_GROUPS_DATA
        else:
            self.groups_data = groups_data

    @staticmethod
    def FromGroupsFile(filename: str):
        """Create a GroupDecomposer from a CSV file."""
        assert filename
        gd = GroupsData.FromGroupsFile(filename)
        return GroupDecomposer(gd)

    @staticmethod
    def _RingedPChainSmarts(length: int) -> str:
        return ''.join(['[C,S][O;R1]',
                        '[P;R1](=O)([OH,O-])[O;R1]' * length, '[C,S]'])

    @staticmethod
    def _InternalPChainSmarts(length: int) -> str:
        return ''.join(['[C,S][O;R0]',
                        '[P;R0](=O)([OH,O-])[O;R0]' * length, '[C,S]'])

    @staticmethod
    def _TerminalPChainSmarts(length: int) -> str:
        return ''.join(['[OH,O-]',
                        'P(=O)([OH,O-])O' * length, '[C,S]'])

    @staticmethod
    def AttachMgToPhosphateChain(mol: Molecule, chain_map: Dict[Group, list],
                                 assigned_mgs: Set[int]) -> Set[int]:
        """Attach Mg2+ ions the appropriate groups in the chain.

        For each Mg2+ we see, we attribute it to a phosphate group if
        possible. We prefer to assign it to a terminal phosphate,
        but otherwise we assign it to a 'middle' group when there are 2 of
        them.

        :param mol: the molecule.
        :param chain_map: the groups in the chain.
        :param assigned_mgs: the set of Mg2+ ions that are already assigned.
        :return: The updated list of assigned Mg2+ ions.
        """
        def AddMg(p_group: Group, pmg_group: Group, mg: List[int]):
            node_set = chain_map[p_group].pop(0)
            mg_index = mg[0]
            node_set.add(mg_index)
            assigned_mgs.add(mg_index)
            chain_map[pmg_group].append(node_set)

        all_pmg_groups = (GroupsData.FINAL_PHOSPHATES_TO_MGS +
                          GroupsData.MIDDLE_PHOSPHATES_TO_MGS +
                          GroupsData.RING_PHOSPHATES_TO_MGS)
        for _mg in mol.FindSmarts('[Mg+2]'):
            if _mg[0] in assigned_mgs:
                continue

            for _p_group, _pmg_group in all_pmg_groups:
                if chain_map[_p_group]:
                    AddMg(_p_group, _pmg_group, _mg)
                    break

        return assigned_mgs

    @staticmethod
    def UpdateGroupMapFromChain(
        group_map: Dict[Group, List[Set[int]]],
        chain_map: Dict[Group, List[Set[int]]]
    ) -> Dict[Group, List[Set[int]]]:
        """Update the group_map by adding the chain."""
        for group, node_sets in chain_map.items():
            group_map.get(group, []).extend(node_sets)
        return group_map

    @staticmethod
    def FindPhosphateChains(
        mol: Molecule, max_length: int = 4,
        ignore_protonations: bool = False
    ) -> List[Tuple[Group, List[Set[int]]]]:
        """Find all phosphate chains.

        Chain end should be 'OC' for chains that do not really end, but link
        to carbons. Chain end should be '[O-1,OH]' for chains that end in an
        hydroxyl.

        :param mol: the molecule to decompose.
        :param max_length: the maximum length of a phosphate chain to consider.
        :param ignore_protonations: whether or not to ignore protonation values.

        :return: A list of 2-tuples (phosphate group, list of occurrences).
        """
        group_map = dict((pg, []) for pg in GroupsData.PHOSPHATE_GROUPS)
        v_charge = [a.GetFormalCharge() for a in mol.GetAtoms()]
        assigned_mgs = set()

        def pop_phosphate(
            pchain: List[int],
            p_size: int
        ) -> Tuple[Set[int], int]:
            if len(pchain) < p_size:
                raise Exception('trying to pop more atoms than are left in '
                                'the pchain')
            phosphate = pchain[0:p_size]
            charge = sum(v_charge[i] for i in phosphate)
            del pchain[0:p_size]
            return set(phosphate), charge

        def add_group(
            chain_map: Dict[Group, list],
            group_name: str,
            charge: int,
            atoms: Set[int]
        ) -> None:
            default = GroupsData.DEFAULTS[group_name]

            if ignore_protonations:
                chain_map[default].append(atoms)
            else:
                # NOTE(flamholz): We rely on the default number of
                # magnesiums being 0 (which it is).
                hydrogens = default.hydrogens + charge - default.charge
                group = Group(group_name, hydrogens, charge, default.nMg)
                if group not in chain_map:
                    # logging.warning('This protonation (%d) level is not
                    # allowed for terminal phosphate groups.' % hydrogens)
                    # logging.warning('Using the default protonation level (
                    # %d) for this name ("%s").' % (default.hydrogens,
                    # default.name))
                    raise GroupDecompositionError(
                        f"The group {group_name} cannot have nH = {hydrogens}"
                    )
                    # chain_map[default].append(atoms)
                else:
                    chain_map[group].append(atoms)

        # For each allowed length
        for length in range(1, max_length + 1):
            # Find internal phosphate chains (ones in the middle of the
            # molecule).
            smarts_str = GroupDecomposer._RingedPChainSmarts(length)
            chain_map = dict((k, []) for (k, _) in group_map.items())
            for pchain in mol.FindSmarts(smarts_str):
                working_pchain = list(pchain)
                working_pchain.pop()  # Lose the last carbon
                working_pchain.pop(0)  # Lose the first carbon

                if length % 2:
                    atoms, charge = pop_phosphate(working_pchain, 5)
                    add_group(chain_map, 'ring -OPO3-', charge, atoms)
                else:
                    atoms, charge = pop_phosphate(working_pchain, 9)
                    add_group(chain_map, 'ring -OPO3-OPO2-', charge, atoms)

                while working_pchain:
                    atoms, charge = pop_phosphate(working_pchain, 8)
                    add_group(chain_map, 'ring -OPO2-OPO2-', charge, atoms)

            assigned_mgs = GroupDecomposer.AttachMgToPhosphateChain(
                mol, chain_map, assigned_mgs)
            GroupDecomposer.UpdateGroupMapFromChain(group_map, chain_map)

            # Find internal phosphate chains (ones in the middle of the
            # molecule).
            smarts_str = GroupDecomposer._InternalPChainSmarts(length)
            chain_map = dict((k, []) for (k, _) in group_map.items())
            for pchain in mol.FindSmarts(smarts_str):
                working_pchain = list(pchain)
                working_pchain.pop()  # Lose the last carbon
                working_pchain.pop(0)  # Lose the first carbon

                if length % 2:
                    atoms, charge = pop_phosphate(working_pchain, 5)
                    add_group(chain_map, '-OPO3-', charge, atoms)
                else:
                    atoms, charge = pop_phosphate(working_pchain, 9)
                    add_group(chain_map, '-OPO3-OPO2-', charge, atoms)

                while working_pchain:
                    atoms, charge = pop_phosphate(working_pchain, 8)
                    add_group(chain_map, '-OPO2-OPO2-', charge, atoms)

            assigned_mgs = GroupDecomposer.AttachMgToPhosphateChain(
                mol, chain_map, assigned_mgs)
            GroupDecomposer.UpdateGroupMapFromChain(group_map, chain_map)

            # Find terminal phosphate chains.
            smarts_str = GroupDecomposer._TerminalPChainSmarts(length)
            chain_map = dict((k, []) for (k, _) in group_map.items())
            for pchain in mol.FindSmarts(smarts_str):
                working_pchain = list(pchain)
                working_pchain.pop()  # Lose the carbon

                atoms, charge = pop_phosphate(working_pchain, 5)
                add_group(chain_map, '-OPO3', charge, atoms)

                if not length % 2:
                    atoms, charge = pop_phosphate(working_pchain, 4)
                    add_group(chain_map, '-OPO2-', charge, atoms)

                while working_pchain:
                    atoms, charge = pop_phosphate(working_pchain, 8)
                    add_group(chain_map, '-OPO2-OPO2-', charge, atoms)

            assigned_mgs = GroupDecomposer.AttachMgToPhosphateChain(
                mol, chain_map, assigned_mgs)
            GroupDecomposer.UpdateGroupMapFromChain(group_map, chain_map)

        return [(pg, group_map[pg]) for pg in GroupsData.PHOSPHATE_GROUPS]

    def CreateEmptyGroupDecomposition(self) -> GroupDecomposition:
        """Create and empty group decomposition object."""
        emptymol = Molecule.FromSmiles("")
        decomposition = self.Decompose(emptymol,
                                       ignore_protonations=True,
                                       raise_exception=False)
        for i, (group, _) in enumerate(decomposition.groups):
            decomposition.groups[i] = (group, [])
        return decomposition

    def Decompose(
        self,
        mol: Molecule,
        ignore_protonations: bool = False,
        raise_exception: bool = False
    ) -> GroupDecomposition:
        """Decompose a molecule into groups.

        The flag 'ignore_protonations' should be used when decomposing a
        compound with lacing protonation representation (for example,
        the KEGG database doesn't posses this information). If this flag is
        set to True, it overrides the '(C)harge sensitive' flag in the
        groups file (i.e. - *PC)

        :param mol: the molecule to decompose.
        :param ignore_protonations: whether to ignore protonation levels.
        :param raise_exception: whether to assert that there are no unassigned
        atoms.

        :return: A GroupDecomposition object containing the decomposition.
        """
        unassigned_nodes = set(range(len(mol)))
        groups: List[Tuple[Group, List[Set[int]]]] = []

        def _AddCorrection(group, count):
            """Add empty sets for each 'correction' group found."""
            list_of_sets = [set() for _ in range(count)]
            groups.append((group, list_of_sets))

        for group in self.groups_data.groups:
            # Phosphate chains require a special treatment
            if group.IsPhosphate():
                pchain_groups = None
                if group.IgnoreCharges() or ignore_protonations:
                    pchain_groups = self.FindPhosphateChains(
                        mol, ignore_protonations=True)
                elif group.ChargeSensitive():
                    pchain_groups = self.FindPhosphateChains(
                        mol, ignore_protonations=False)
                else:
                    raise MalformedGroupDefinitionError(
                        'Unrecognized phosphate wildcard: %s' % group.name)

                for phosphate_group, group_nodesets in pchain_groups:
                    current_groups = []

                    for focal_set in group_nodesets:
                        if focal_set.issubset(unassigned_nodes):
                            # Check that the focal-set doesn't override an
                            # assigned node
                            current_groups.append(focal_set)
                            unassigned_nodes = unassigned_nodes - focal_set
                    groups.append((phosphate_group, current_groups))
            elif group.IsCodedCorrection():
                _AddCorrection(group, group.GetCorrection(mol))
            # Not a phosphate group or expanded correction.
            else:
                # TODO: if the 'ignore_protonation' flag is True,
                #  this should always use the pseudogroup with the lowest nH
                #  in each category regardless of the hydrogens in the given
                #  Mol.
                current_groups = []
                for nodes in mol.FindSmarts(group.smarts):
                    try:
                        focal_nodes = set(group.FilterFocalSet(nodes))
                    except IndexError as e:
                        logging.error(
                            'Focal set for group %s is out of range: %s' % (
                                str(group), str(group.focal_atoms)))
                        raise e

                    # check that the focal-set doesn't override an assigned
                    # node
                    if focal_nodes.issubset(unassigned_nodes):
                        current_groups.append(focal_nodes)
                        unassigned_nodes = unassigned_nodes - focal_nodes
                groups.append((group, current_groups))

        # Ignore the hydrogen atoms when checking which atom is unassigned
        for nodes in mol.FindSmarts('[H]'):
            unassigned_nodes = unassigned_nodes - set(nodes)

        decomposition = GroupDecomposition(self.groups_data, mol,
                                           groups, unassigned_nodes)

        if raise_exception and decomposition.unassigned_nodes:
            raise GroupDecompositionError(
                f"Unable to decompose {mol} into groups",
                decomposition)

        return decomposition
