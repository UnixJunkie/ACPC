#!/usr/bin/env python2

# output the MACCS bitstring of each molecule found in a MOL2 file

import rdkit.Chem
import sys

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

if __name__ == "__main__":
    import sys
    from rdkit.Chem import MACCSkeys
    with open(sys.argv[1]) as in_file:
        for mol2 in RetrieveMol2Block(in_file):
            mol = rdkit.Chem.MolFromMol2Block(mol2)
            try:
                maccs = MACCSkeys.GenMACCSKeys(mol)
                for bit in maccs:
                    if bit:
                        sys.stdout.write('1')
                    else:
                        sys.stdout.write('0')
                sys.stdout.write('\n')
            except:
                sys.stdout.write('0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000\n')
