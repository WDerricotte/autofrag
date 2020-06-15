#!/usr/bin/env python
"""
autofrag.py
A command line tool for parsing interacting fragments in noncovalent molecular complexes using graph theory

Handles the primary functions
"""
import sys
from molecule import *
from fragment import *
import numpy as np
    

def main():
    geom_file = sys.argv[1]

    mol = Molecule(geom_file)
    frag = Fragments()
    struct = mol.read_molecule()
    mol.store_molecule(struct)

    
    


    adj_matrix = frag.build_adj_matrix(mol.atoms, mol.geometry, mol.units)
    frags = frag.get_fragments(adj_matrix)
    print(frags)

    frag_output = open('frag_%s' %geom_file, "w+")
    frag_output.write('%d\n\n' %mol.natoms)
    for frag in frags:
        for atom_num in frag:
            atom_key = mol.atoms[atom_num-1]
            sym_only = "".join(sym for sym in atom_key if not sym.isdigit())
            frag_output.write(f"{sym_only}   ")
            coordinate = "  ".join(str(coord) for coord in mol.geometry[atom_key])
            frag_output.write("%s \n" %coordinate)
        index = frags.index(frag)
        if(index!=(len(frags)-1)):
            frag_output.write("--\n")
    frag_output.close()

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    main()
