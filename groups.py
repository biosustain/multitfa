"""module for managing groups."""
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
import json
import logging
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd
import quilt

from component_contribution import DEFAULT_QUILT_PKG

from molecule import Molecule


class GroupsDataError(Exception):
    """A general group data error."""

    pass


class MalformedGroupDefinitionError(GroupsDataError):
    """An error for malform group definitions."""

    pass


class FocalSet(set):
    """Class of the set of focal atoms."""

    def __init__(self, focal_atoms_str: str):
        """Create a FocalSet object."""
        if not focal_atoms_str:
            raise ValueError(
                'You must supply a non-empty focal atom string. '
                'You may use "None" or "All" in the obvious fashion.')

        self.focal_atoms_str = focal_atoms_str
        prepped_str = focal_atoms_str.strip().lower()

        if prepped_str == 'all':
            super(FocalSet, self).__init__()
            self.all_atoms = True
        elif prepped_str == 'none':
            super(FocalSet, self).__init__()
            self.all_atoms = False
        else:
            super(FocalSet, self).__init__(map(int, focal_atoms_str.split('|')))
            self.all_atoms = False

    def __str__(self) -> str:
        """Return as string."""
        return self.focal_atoms_str

    def __contains__(self, elt: int) -> bool:
        """Check if the set contains an element."""
        if self.all_atoms:
            return True
        else:
            return super(FocalSet, self).__contains__(elt)


class Group(object):
    """Representation of a single group."""

    def __init__(
        self,
        name: str,
        hydrogens: int,
        charge: int,
        nMg: int,
        smarts: str = None,
        focal_set: FocalSet = None
    ):
        """Create a Group object."""
        self.name = name
        self.hydrogens = hydrogens
        self.charge = charge
        self.nMg = nMg
        self.smarts = smarts
        self.focal_set = focal_set

    def _IsHydrocarbonGroup(self) -> bool:
        """Check if the group is a hydrocarbon."""
        return self.name.startswith('*Hc')

    def _IsSugarGroup(self) -> bool:
        """Check if the group is a sugar."""
        return self.name.startswith('*Su')

    def _IsAromaticRingGroup(self) -> bool:
        """Check if the group is aromatic."""
        return self.name.startswith('*Ar')

    def _IsHeteroaromaticRingGroup(self) -> bool:
        """Check if the group is a heteroaromatic ring."""
        return self.name.startswith('*Har')

    def IsPhosphate(self) -> bool:
        """Check if the group is a phosphate group."""
        return self.name.startswith('*P')

    def IgnoreCharges(self) -> bool:
        """Check if this group is set to ignore charges."""
        # (I)gnore charges
        return self.name[2] == 'I'

    def ChargeSensitive(self) -> bool:
        """Check if this group is charge sensitive."""
        # (C)harge sensitive
        return self.name[2] == 'C'

    def IsCodedCorrection(self) -> bool:
        """Return True for corrections with hand-written code."""
        return (self._IsHydrocarbonGroup() or
                self._IsAromaticRingGroup() or
                self._IsHeteroaromaticRingGroup())

    @staticmethod
    def _IsHydrocarbon(mol) -> int:
        """Test if a molecule is a simple hydrocarbon."""
        if mol.FindSmarts('[!C;!c]'):
            # If we find anything other than a carbon (w/ hydrogens)
            # then it's not a hydrocarbon.
            return 0
        return 1

    @staticmethod
    def _CountAromaticRings(mol) -> int:
        """Count the number of aromatic rings."""
        expressions = ['c1cccc1', 'c1ccccc1']
        count = 0
        for smarts_str in expressions:
            count += len(mol.FindSmarts(smarts_str))
        return count

    @staticmethod
    def _CountHeteroaromaticRings(mol) -> int:
        expressions = ['a1aaaa1', 'a1aaaaa1']
        count = 0
        all_atoms = mol.GetAtoms()
        for smarts_str in expressions:
            for match in mol.FindSmarts(smarts_str):
                atoms = set([all_atoms[i].atomicnum for i in match])
                atoms.discard(6)  # Ditch carbons
                if atoms:
                    count += 1
        return count

    def GetCorrection(self, mol) -> int:
        """Get the value of the correction for this molecule."""
        if self._IsHydrocarbonGroup():
            return self._IsHydrocarbon(mol)
        elif self._IsAromaticRingGroup():
            return self._CountAromaticRings(mol)
        elif self._IsHeteroaromaticRingGroup():
            return self._CountHeteroaromaticRings(mol)

        raise TypeError('This group is not a correction.')

    def FilterFocalSet(self, nodes: Iterable[int]) -> Iterable[int]:
        """Get the subset of Focal atoms in a list of nodes.

        :param nodes: the nodes matching this group.
        :return: The subset of focal atoms.
        """
        for i, node in enumerate(nodes):
            if i in self.focal_set:
                yield node

    def __str__(self) -> str:
        """Return a string representation of the group."""
        if self.hydrogens is None:
            return '%s' % self.name
        if self.charge is None:
            return '%s' % self.name
        if self.nMg is None:
            return '%s' % self.name
        return '%s [H%d Z%d Mg%d]' % \
            (self.name, self.hydrogens or 0, self.charge or 0, self.nMg or 0)

    def __eq__(self, other) -> bool:
        """Enable == checking.

        Only checks name, protons, charge, and nMg.
        """
        return (str(self.name) == str(other.name) and
                self.hydrogens == other.hydrogens and
                self.charge == other.charge and
                self.nMg == other.nMg)

    def __hash__(self):
        """We are HASHABLE.

        Note that the hash depends on the same attributes that are checked
        for equality.
        """
        return hash((self.name, self.hydrogens, self.charge, self.nMg))


class GroupsData(object):
    """Contains data about all groups."""

    ORIGIN = Group('Origin', hydrogens=0, charge=0, nMg=0)

    # Phosphate groups need special treatment, so they are defined in code...
    # TODO(flamholz): Define them in the groups file.

    # each tuple contains: (name, is_default, Group object)

    phosphate_groups = [
        ('initial H0', True,
         Group(name='-OPO3-', hydrogens=0, charge=-1, nMg=0)),
        ('initial H1', False,
         Group('-OPO3-', hydrogens=1, charge=0, nMg=0)),
        ('middle H0', True,
         Group('-OPO2-', hydrogens=0, charge=-1, nMg=0)),
        ('middle H1', False,
         Group('-OPO2-', hydrogens=1, charge=0, nMg=0)),
        ('final H0', True,
         Group('-OPO3', hydrogens=0, charge=-2, nMg=0)),
        ('final H1', False,
         Group('-OPO3', hydrogens=1, charge=-1, nMg=0)),
        ('final H2', False,
         Group('-OPO3', hydrogens=2, charge=0, nMg=0)),
        ('initial chain H0', True,
         Group('-OPO3-OPO2-', hydrogens=0, charge=-2, nMg=0)),
        ('initial chain H1', False,
         Group('-OPO3-OPO2-', hydrogens=1, charge=-1, nMg=0)),
        ('initial chain H2', False,
         Group('-OPO3-OPO2-', hydrogens=2, charge=0, nMg=0)),
        ('initial chain Mg1', False,
         Group('-OPO3-OPO2-', hydrogens=0, charge=0, nMg=1)),
        ('middle chain H0', True,
         Group('-OPO2-OPO2-', hydrogens=0, charge=-2, nMg=0)),
        ('middle chain H1', False,
         Group('-OPO2-OPO2-', hydrogens=1, charge=-1, nMg=0)),
        ('middle chain H2', False,
         Group('-OPO2-OPO2-', hydrogens=2, charge=0, nMg=0)),
        ('middle chain Mg1', False,
         Group('-OPO2-OPO2-', hydrogens=0, charge=0, nMg=1)),
        ('ring initial H0', True,
         Group('ring -OPO3-', hydrogens=0, charge=-1, nMg=0)),
        ('ring initial H1', False,
         Group('ring -OPO3-', hydrogens=1, charge=0, nMg=0)),
        ('ring initial chain H0', True,
         Group('ring -OPO3-OPO2-', hydrogens=0, charge=-2, nMg=0)),
        ('ring initial chain H1', False,
         Group('ring -OPO3-OPO2-', hydrogens=1, charge=-1, nMg=0)),
        ('ring initial chain H2', False,
         Group('ring -OPO3-OPO2-', hydrogens=2, charge=0, nMg=0)),
        ('ring middle chain H0', True,
         Group('ring -OPO2-OPO2-', hydrogens=0, charge=-2, nMg=0)),
        ('ring middle chain H1', False,
         Group('ring -OPO2-OPO2-', hydrogens=1, charge=-1, nMg=0)),
        ('ring middle chain H2', False,
         Group('ring -OPO2-OPO2-', hydrogens=2, charge=0, nMg=0)),
        ('ring initial chain Mg1', False,
         Group('ring -OPO2-OPO2-', hydrogens=0, charge=0, nMg=1))
    ]

    PHOSPHATE_GROUPS: List[Group] = []
    PHOSPHATE_DICT: Dict[str, Group] = {}
    DEFAULTS = {}
    for full_name, is_default, group in phosphate_groups:
        PHOSPHATE_GROUPS.append(group)
        PHOSPHATE_DICT[full_name] = group
        if is_default:
            DEFAULTS[group.name] = group

    RING_PHOSPHATES_TO_MGS = ((PHOSPHATE_DICT['ring initial chain H0'],
                               PHOSPHATE_DICT['ring initial chain Mg1']),)
    MIDDLE_PHOSPHATES_TO_MGS = ((PHOSPHATE_DICT['initial chain H0'],
                                 PHOSPHATE_DICT['initial chain Mg1']),)
    FINAL_PHOSPHATES_TO_MGS = ((PHOSPHATE_DICT['middle chain H0'],
                                PHOSPHATE_DICT['middle chain Mg1']),)

    def __init__(self, groups: List[Group]):
        """Construct GroupsData.

        :param groups: a list of Group objects.
        chemical groups.
        """
        self.groups = groups
        self.all_groups = self._GetAllGroups(self.groups)
        self.all_group_names = [str(g) for g in self.all_groups]
        self.all_group_hydrogens = np.array(
            [g.hydrogens or 0 for g in self.all_groups])
        self.all_group_charges = np.array(
            [g.charge or 0 for g in self.all_groups])
        self.all_group_mgs = np.array(
            [g.nMg or 0 for g in self.all_groups])

    def ToDataFrame(self) -> pd.DataFrame:
        """Convert to a DataFrame."""
        data = list(zip(self.all_group_charges,
                        self.all_group_hydrogens,
                        self.all_group_mgs,
                        self.all_group_names))
        df = pd.DataFrame(data=data, columns=['charge', 'num_protons',
                                              'num_mgs', 'full_name'])
        df.index = df.full_name.str.findall(r"^(.+) \[H\d Z-?\d Mg\d\]").str[0]
        df.index.name = 'group_name'
        return df

    def Count(self) -> int:
        """Count the number of groups."""
        return len(self.all_groups)

    count = property(Count)

    @staticmethod
    def _GetAllGroups(groups: Iterable[Group]) -> List[Group]:
        all_groups = []

        for group in groups:
            # Expand phosphate groups.
            if group.IsPhosphate():
                all_groups.extend(GroupsData.PHOSPHATE_GROUPS)
            else:
                all_groups.append(group)

        # Add the origin.
        all_groups.append(GroupsData.ORIGIN)
        return all_groups

    @staticmethod
    def FromDataFrame(gr_def_df: pd.DataFrame) -> object:
        """Initialize a GroupData from a DataFrame."""
        list_of_groups = []
        for row in gr_def_df.itertuples(index=True):
            logging.debug(f"Reading group definition for {row.NAME}")

            # Check that the smarts are good.
            if not Molecule.VerifySmarts(row.SMARTS):
                raise GroupsDataError('Cannot parse SMARTS expression: %s' %
                                      row.SMARTS)

            group = Group(row.NAME,
                          hydrogens=row.PROTONS,
                          charge=row.CHARGE,
                          nMg=row.MAGNESIUMS,
                          smarts=row.SMARTS,
                          focal_set=FocalSet(row.FOCAL_ATOMS))
            list_of_groups.append(group)

        logging.debug('Done reading groups data.')

        return GroupsData(list_of_groups)

    @staticmethod
    def FromQuilt() -> object:
        """Create a GroupsData object from quilt."""
        quilt.install(DEFAULT_QUILT_PKG, force=True)
        pkg = quilt.load(DEFAULT_QUILT_PKG)
        return GroupsData.FromDataFrame(pkg.train.group_definitions())

    @staticmethod
    def FromGroupsFile(file):
        """Initialize a GroupData from a CSV file."""
        assert file
        gr_def_df = pd.read_csv(file, index_col=None, header=0)
        return GroupsData.FromDataFrame(gr_def_df)

    def Index(self, gr):
        """Get the index of a group."""
        try:
            return self.all_groups.index(gr)
        except ValueError:
            raise ValueError('group %s is not defined' % str(gr))

    def GetGroupNames(self):
        """Get all group names."""
        return self.all_group_names


class GroupVector(list):
    """A vector of groups."""

    def __init__(
        self,
        groups_data: GroupsData,
        iterable: Iterable[float] = None
    ):
        """Construct a vector.

        Args:
            groups_data: a GroupsData object, contains data about all groups.
            iterable: data to load into the vector.
        """
        self.groups_data = groups_data

        if iterable is not None:
            super(GroupVector, self).__init__(iterable)
        else:
            super(GroupVector, self).__init__([0] * len(self))

    def __str__(self):
        """Return a sparse string representation of this group vector."""
        return " | ".join(
            [f"{name} x {count}" for
             name, count in zip(self.groups_data.GetGroupNames(), self)
             if count != 0]
        )

    def __len__(self) -> int:
        """Get the number of groups."""
        return len(self.groups_data.GetGroupNames())

    def __iadd__(self, other) -> object:
        """Add another group vector."""
        for i in range(len(self)):
            self[i] += other[i]
        return self

    def __isub__(self, other) -> object:
        """Subtract another group vector."""
        for i in range(len(self)):
            self[i] -= other[i]
        return self

    def __add__(self, other) -> object:
        """Add two group vectors."""
        result = GroupVector(self.groups_data)
        for i in range(len(self)):
            result[i] = self[i] + other[i]
        return result

    def __sub__(self, other) -> object:
        """Subtract two group vectors."""
        result = GroupVector(self.groups_data)
        for i in range(len(self)):
            result[i] = self[i] - other[i]
        return result

    def __eq__(self, other) -> bool:
        """Compare two group vectors."""
        for i in range(len(self)):
            if self[i] != other[i]:
                return False
        return True

    def __nonzero__(self) -> bool:
        """Check if the group vector is zero."""
        for i in range(len(self)):
            if self[i] != 0:
                return True
        return False

    def __mul__(self, other: float) -> object:
        """Multiply by a scalar."""
        try:
            c = float(other)
            return GroupVector(self.groups_data, [x*c for x in self])
        except ValueError:
            raise ValueError("A GroupVector can only be multiplied by a scalar"
                             ", given " + str(other))

    def NetCharge(self) -> float:
        """Return the net charge."""
        return float(np.dot(self, self.groups_data.all_group_charges))

    def Hydrogens(self) -> float:
        """Return the number of protons."""
        return float(np.dot(self, self.groups_data.all_group_hydrogens))

    def Magnesiums(self) -> float:
        """Return the number of Mg2+ ions."""
        return float(np.dot(self, self.groups_data.all_group_mgs))

    def RemoveEpsilonValues(self, epsilon=1e-10):
        """Set epsilon values to zero."""
        for i in range(len(self)):
            if abs(self[i]) < epsilon:
                self[i] = 0

    def ToJSONString(self):
        """Return a JSON representation of the group vector."""
        return json.dumps(dict([(i, x)
                                for (i, x) in enumerate(self) if x != 0]))

    @staticmethod
    def FromJSONString(groups_data, s):
        """Read from a JSON string."""
        v = [0] * groups_data.Count()
        for i, x in json.loads(s).items():
            v[int(i)] = x
        return GroupVector(groups_data, v)

    def __iter__(self) -> Iterable[float]:
        """Return a flat version of the group vector."""
        return super(GroupVector, self).__iter__()

    def as_array(self):
        """Return a NumPy array with the group counts."""
        return np.array(list(self.__iter__()))

    def as_dict(self) -> Dict[str, int]:
        """Convert the group vector into a dictionary."""
        return {name: count for name, count in
                zip(self.groups_data.GetGroupNames(), self.__iter__())
                if count != 0}


DEFAULT_GROUPS_DATA = GroupsData.FromQuilt()
