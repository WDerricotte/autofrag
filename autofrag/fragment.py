import sys
import os
import numpy as np

periodic_table = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
    "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]
    
#Covalent radii taken from DOI: 10.1039/b801115j
#Everything beyond Cm was set to 1.80
covalent_radii = [0.00,0.32,0.28,1.28,0.96,0.84,0.76,0.71,0.66,0.57,0.58,1.66,1.41,1.21,1.11,1.07,1.05,1.02,1.06,2.03,1.76,1.70,1.60,1.53,1.39,1.61,1.52,1.50,1.24,1.32,1.22,1.22,1.20,1.19,1.20,1.20,1.16,2.20,1.95,1.90,1.75,
    1.64,1.54,1.47,1.46,1.42,1.39,1.45,1.44,1.42,1.39,1.39,1.38,1.39,1.40,2.44,2.15,2.07,2.04,2.03,2.01,1.99,1.98,1.98,1.96,1.94,1.92,1.92,1.89,1.90,1.87,1.87,1.75,1.70,1.62,1.51,1.44,1.41,1.36,1.36,1.32,1.45,
    1.46,1.48,1.40,1.50,1.50,2.60,2.21,2.15,2.06,2.00,1.96,1.90,1.87,180,169,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80,1.80]

vdw_radii = [0.00,1.10,1.40,1.81,1.53,1.92,1.70,1.55,1.52,1.47,1.54,2.27,1.73,1.84,2.10,1.80,1.75,1.88,2.75,2.31,1.87,2.11,1.85,1.90,1.83,2.02,
    3.03,2.49,1.93,2.17,2.06,2.06,1.98,2.16,3.43,2.68,1.96,2.02,2.07,1.97,2.02,2.20,3.48,2.83]

class Fragments:
    def __init__(self):
        pass

    @staticmethod
    def calculate_distance(atom1, atom2):
        distance = 0.00
        distance = np.sqrt((atom1[0]-atom2[0])**2+(atom1[1]-atom2[1])**2+(atom1[2]-atom2[2])**2)
        return distance	

    def check_bond(self,distance, atomtype1, atomtype2, covalent_radii, units):
        bond = True
        '''find bond threshold based on covalent radii
        taken from 10.1039/b801115j
        For everything beyond Cm 1.80 A was used.
        get radii from saved numpy array
        use atomtypes as numbers
        use factor to be sure to catch bonds. Bigger factor means more bonds are detected
        '''
        atomtype1 = int(periodic_table.index(atomtype1))
        atomtype2 = int(periodic_table.index(atomtype2))
        bond_threshold = 1.25 * (covalent_radii[atomtype1]+covalent_radii[atomtype2])
        if(units=="BOHR"):
            distance = distance*0.529177
        else:
            pass
        if distance > bond_threshold:
            bond = False
        return bond

    @staticmethod
    def build_dist_matrix(atoms, geometry):
        natoms = len(atoms)
        distance_matrix = np.zeros((natoms,natoms)) # Initialize distance matrix
        for i in range(natoms):
            atom1 = atoms[i]
            for j in range(natoms):
                atom2 = atoms[j]
                distance_matrix[i][j] = Fragments.calculate_distance(geometry[atom1],geometry[atom2])
        return distance_matrix


    def build_adj_matrix(self, atoms, geometry, units="ANGSTROM"):
        natoms = len(atoms)
        distance_matrix = Fragments.build_dist_matrix(atoms,geometry)
        adjacency_matrix = np.zeros((natoms,natoms))
        for i in range(natoms):
            atom1 = ''.join(sym for sym in atoms[i] if not sym.isdigit())
            for j in range(natoms):
                atom2 = ''.join(sym for sym in atoms[j] if not sym.isdigit())
                if self.check_bond(distance_matrix[i][j], atom1, atom2, covalent_radii, units):
                    adjacency_matrix[i][j] = 1
                else:
                    adjacency_matrix[i][j] = 0
                if(i==j):
                    adjacency_matrix[i][j] = 0
        return adjacency_matrix

    def build_deg_matrix(self, atoms, adjacency_matrix):
        natoms = len(atoms)
        degrees = np.sum(adjacency_matrix, axis=0) # Sums each column of Adjacency Matrix
        degree_matrix = np.diag(degrees)
        return degree_matrix

    def build_laplacian_matrix(self, atoms, adjacency_matrix):
        degree_matrix = self.build_deg_matrix(atoms, adjacency_matrix)
        laplacian = np.subtract(degree_matrix, adjacency_matrix)
        return laplacian
        
    def get_fragments(self, adjacency_matrix):
        fragments = []
        for i in range(len(adjacency_matrix)):
            bonded_atoms = [i+1]
            for j in range(len(adjacency_matrix)):
                element = adjacency_matrix[i][j]
                if(element==1):
                    bonded_atoms.append(j+1)
            fragments.append(bonded_atoms)
        # Join Fragments on the same molecule
        fragments_new = []
        for i in range(len(fragments)):
            for j in range(len(fragments)):
                if(set(fragments[i]) & set(fragments[j])): #See if the two lists share any values
                    combine = sorted(fragments[i] + list(set(fragments[j]) - set(fragments[i])))
                    if(len(combine)>0):
                        if(len(fragments_new) > 0):
                            # Does it share elements with any existing lists?
                            partner = True
                            for k in range(len(fragments_new)):
                                if(set(combine) & set(fragments_new[k])):
                                    # If So, append the new elements to that list
                                    fragments_new[k].extend(x for x in combine if x not in fragments_new[k])
                                    fragments_new[k] = sorted(fragments_new[k])
                                    partner = True
                                    break
                                else:
                                    partner=False
                            if(partner==False):
                                fragments_new.append(combine) 
                        if(len(fragments_new)==0):
                            # If Not, create a new fragment
                            fragments_new.append(combine)
        return fragments_new

