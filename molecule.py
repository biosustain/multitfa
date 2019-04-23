"""module for dealing with molecules."""
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
from typing import Iterable, List

import openbabel as ob


class OpenBabelError(Exception):
    """An error in Open Babel."""

    pass


class Molecule(object):
    """A sort of wrapper class for OpenBabel."""

    _obElements = ob.OBElementTable()
    _obSmarts = ob.OBSmartsPattern()

    def __init__(self):
        """Create an empty Molecule object."""
        self.obmol = ob.OBMol()

    def __len__(self) -> int:
        """Get the number of atoms."""
        return self.obmol.NumAtoms()

    def __str__(self) -> str:
        """Return a string representation of the Molecule."""
        return str(self.obmol)

    def GetAtoms(self) -> List[ob.OBAtom]:
        """Get the list of atoms."""
        return [self.obmol.GetAtom(i+1) for i in range(self.obmol.NumAtoms())]

    def FindSmarts(self, smarts: str) -> Iterable[List[int]]:
        """Use SMARTS to find substructures.

        Corrects the pyBel version of Smarts.findall() which returns results
        as tuples, with 1-based indices even though Molecule.atoms is 0-based.

        :param smarts: the SMARTS query to search for.
        :return: The re-mapped list of SMARTS matches.
        """
        Molecule._obSmarts.Init(smarts)
        if Molecule._obSmarts.Match(self.obmol):
            match_list = Molecule._obSmarts.GetMapList()
            return map(lambda m: [(n - 1) for n in m], match_list)
        else:
            return []

    @staticmethod
    def GetSymbol(atomic_num: int) -> str:
        """Get the symbol associated with an atomic number."""
        return Molecule._obElements.GetSymbol(atomic_num)

    @staticmethod
    def VerifySmarts(smarts):
        """Verify that the SMARTS string is valid."""
        return Molecule._obSmarts.Init(smarts)

    @staticmethod
    def FromSmiles(smiles: str):
        """Create a Molecule object from a SMILES string."""
        m = Molecule()
        obConversion = ob.OBConversion()
        obConversion.AddOption("w", obConversion.OUTOPTIONS)
        obConversion.SetInFormat("smiles")
        if not obConversion.ReadString(m.obmol, smiles):
            raise OpenBabelError("Cannot read the SMILES string: " + smiles)
        return m

    @staticmethod
    def FromInChI(inchi: str) -> object:
        """Create a Molecule object from an InChI string."""
        m = Molecule()
        m.inchi = inchi
        obConversion = ob.OBConversion()
        obConversion.AddOption("w", obConversion.OUTOPTIONS)
        obConversion.SetInFormat("inchi")
        obConversion.ReadString(m.obmol, m.inchi)
        return m
