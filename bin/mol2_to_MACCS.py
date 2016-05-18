#!/usr/bin/env python2

# output the MACCS bitstring of each molecule found in a MOL2 file

import os
import rdkit.Chem
import sys

from rdkit.Chem import MACCSkeys

def RetrieveMol2Block(fileLikeObject, delimiter="@<TRIPOS>MOLECULE"):
    """generator which retrieves one mol2 block at a time
    """
    mol2 = []
    for line in fileLikeObject:
        if line.startswith(delimiter) and mol2:
            yield "".join(mol2)
            mol2 = []
        mol2.append(line)
    if mol2:
        yield "".join(mol2)

input_fn = sys.argv[1]

with open(input_fn) as in_file:
    problems_fn = os.path.splitext(input_fn)[0] + ".errors.mol2"
    problem_mols = open(problems_fn, 'w')
    for mol2 in RetrieveMol2Block(in_file):
        mol = rdkit.Chem.MolFromMol2Block(mol2)
        mol_lines = str(mol2).split("\n")
        mol_name = mol_lines[1]
        bitstring = mol_name + " "
        try:
            maccs = MACCSkeys.GenMACCSKeys(mol)
            for bit in maccs:
                if bit:
                    bitstring += "1"
                else:
                    bitstring += "0"
            bitstring += "\n"
        except:
            bitstring = '0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000\n'
            problem_mols.write(mol2)
        sys.stdout.write(bitstring)
