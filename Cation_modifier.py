from pymatgen.core import Structure
from pymatgen.core import Element
from glob import glob
import numpy as np
from scipy.spatial.transform import Rotation as R
from Functions import molecule_rotation, convert_site_index
import copy

struct = Structure.from_file("Test_structures/FPB_bulk_cubic.vasp")

#Need to do:
# Rotation of cations
#   1) Identify molecule
#       How? --> Bond information (maybe from VESTA)

# Read bond data
# Input data from user: what elements composing the molecule
elements_list=['C','H','N']



# When replacing inorganic cation to organic molecule,
# I need to operate with the center of mass
# command: Element.C.atomic_mass

n_atom_count_dict=dict()
label2site_index=dict()
for i in range(0,struct.num_sites):
    # Update label for each element
    if struct[i].specie in list(n_atom_count_dict.keys()):
        n_atom_count_dict.update({struct[i].specie:n_atom_count_dict[struct[i].specie]+1})
    else:
        n_atom_count_dict.update({struct[i].specie:1})
    #Example: BaTiO3 --> Ba1:0 Ti1:1 O1:2 O2:3 O3:4
    label2site_index.update({'{0}{1}'.format(struct.species[i], n_atom_count_dict[struct[i].specie]):i})
    
def Bond_length(pmg_struct, label2site_index, atom1_label,atom2_label):
    '''
    input: atom1_label, atom2_label -- str,  Ex) Fe5, O1, 1, 8 etc. 
            label2site_index -- dictionary defined in script
            pmg_struct -- pymatgen structure
    output: bondlength -- float, bond length of (atom1-atom2) will be calculated.
    '''
    atom1_pmg_index=convert_site_index(label2site_index, atom1_label)
    atom2_pmg_index=convert_site_index(label2site_index, atom2_label)
    bondlength=pmg_struct.get_distance(atom1_pmg_index,atom2_pmg_index)
    return bondlength

molecules_list=[]
molecule=['C1','N1','N2','H1','H2','H3','H4','H5']
molecules_list.append(molecule)
reference_point=struct.sites[label2site_index['C1']].frac_coords

# axis_vector and rotation angle
axis_vector=np.array([0,1,0])  ## To be modified
angle_degree = 45

struct=molecule_rotation(struct,molecule,label2site_index,
                         axis_vector,angle_degree,reference_point)

struct.to(fmt='poscar',filename='test.vasp')
    

    
    
    
    
    