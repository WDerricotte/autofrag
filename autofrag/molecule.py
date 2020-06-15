import sys
import os
import copy
import numpy as np

class Molecule:
    def __init__(self, filename):
        self.filename = filename
        self.natoms = 0
        self.atoms = []
        self.geometry = dict()

    def get_natoms(self):
        """
        Based on the molecule file this function returns the number of atoms
        
        Paramaters
        ----------
        None

        Returns
        -------
        natoms : int
            Number of atoms in the structure
        """
        xyz_file = open(self.filename, 'r')
        natoms = int(xyz_file.readline())
        return natoms

    def units(self):
        """
        Sets the units for the input structure, is ANGSTROM by default
        """
        self.units = 'ANGSTROM'

    def read_molecule(self):
        """
        Function to read Molecule from an xyz coordinate file

        Parameters
        ----------
        filename : str
            The name of the file containing the coordinates for your molecule

        Returns
        -------
        struct : str
            String containing only the atoms and their respective xyz coordinates
        """
        xyz_file = open(self.filename, 'r')
        self.natoms = self.get_natoms()
        struct_list = xyz_file.readlines()[1:]
        struct = ""
        struct = ' \n'.join([str(atom) for atom in struct_list])
        return struct 

    def get_symbols(self):
        """
        Function that returns a list of the unique atom labels with atom symbol and number
        for example the 3rd atom that is a carbon is listed as C3

        Parameters
        ----------
        None

        Returns
        -------
        geometry_keys : list
            A list of unique atom labels
        """
        return self.geometry.keys()

    def get_geometry(self):
        geom = self.geometry
        return copy.deepcopy(geom)

    def write_molecule(self, outfile, comment=""):
        """
        Writes the molecule to an output xyz file.

        Parameters
        ----------
        outfile : str
            The name of the file the molecule will be written to
        comment : str
            An optional comment that will appear on the second line of the file
            this is empty by default, writing a blank line

        Returns
        -------
        None
        """
        mol_file = open(outfile, "w+")
        atoms = self.atoms
        natoms = len(atoms)
        mol_file.write(f"{natoms}\n")
        mol_file.write(f"{comment}\n")
        for atom in atoms:
            x_coord = self.geometry[atom][0]
            y_coord = self.geometry[atom][1]
            z_coord = self.geometry[atom][2]
            mol_file.write(''.join(sym for sym in atom if not sym.isdigit()))
            mol_file.write(f"   {x_coord}   {y_coord}   {z_coord}\n")
        mol_file.close()


    def store_molecule(self, struct):
        """
        Takes the molecule that was read in from a file, and stores it in a dictionary.
        """
        struct = struct.split()  # Unparsed List of all symbols and coordinates
        for i in range(self.natoms):
            atom_key = "%s%d" %(str(struct[i*4]),(i+1))
            self.atoms.append(atom_key)
            atom_coord = [float(struct[(k+1) + i*4]) for k in range(3)]
            self.geometry[atom_key] = atom_coord

    def angstroms_to_bohr(self):
        bohr_geometry = {}
        for atom in self.atoms:
            bohr_geometry[atom] = [coord*1.88973 for coord in self.geometry[atom]]
        self.units = 'BOHR'
        self.geometry = bohr_geometry 

    def bohr_to_angstroms(self):
        atoms = self.get_symbols()
        angstrom_geometry = {}
        for atom in atoms:
            angstrom_geometry[atom] = [coord*0.529177 for coord in self.geometry[atom]]
        self.units = 'ANGSTROM'
        self.geometry = angstrom_geometry
             