from pymatgen.core import Structure
from pymatgen.core import Element
from glob import glob
import numpy as np
from scipy.spatial.transform import Rotation as R
from Functions import molecule_rotation, convert_site_index, Write_POSCAR, molecule_finder, create_df
import copy
import pandas as pd

struct = Structure.from_file("FPB_FA_arranged_multiple_rotation.vasp")
#struct = Structure.from_file("test.vasp")

# Read bond data
# Input data from user: what elements composing the molecule
elements_list=['C','H','N']



# When replacing inorganic cation to organic molecule,
# I need to operate with the center of mass
# command: Element.C.atomic_mass

# ############## USING PANDAS DATAFRAME ##########################
# n_atom_count_dict=dict()
# df=pd.DataFrame(columns=['site_index','atom_label','pmg_site','element'])
# for i in range(0,struct.num_sites):
#     # Update label for each element
#     if struct[i].specie in list(n_atom_count_dict.keys()):
#         n_atom_count_dict.update({struct[i].specie:n_atom_count_dict[struct[i].specie]+1})
#     else:
#         n_atom_count_dict.update({struct[i].specie:1})
        
#     label='{0}{1}'.format(struct.species[i], n_atom_count_dict[struct[i].specie])
#     # Append a site to the data frame
#     # If this part is costly, maybe using pd.concat would be faster (not sure yet)
#     df= df.append({'site_index':i, \
#                 'atom_label':'{0}{1}'.format(struct.species[i], n_atom_count_dict[struct[i].specie]), \
#                 'pmg_site':struct.sites[i],\
#                 'element':str((struct.sites[i]).specie)},ignore_index=True)
# ################################################################

def update_dataframe(struct):
    ############## USING PANDAS DATAFRAME ##########################
    n_atom_count_dict=dict()
    df=pd.DataFrame(columns=['site_index','atom_label','pmg_site','element'])
    for i in range(0,struct.num_sites):
        # Update label for each element
        if struct[i].specie in list(n_atom_count_dict.keys()):
            n_atom_count_dict.update({struct[i].specie:n_atom_count_dict[struct[i].specie]+1})
        else:
            n_atom_count_dict.update({struct[i].specie:1})
            
        label='{0}{1}'.format(struct.species[i], n_atom_count_dict[struct[i].specie])
        # Append a site to the data frame
        # If this part is costly, maybe using pd.concat would be faster (not sure yet)
        df= df.append({'site_index':i, \
                    'atom_label':'{0}{1}'.format(struct.species[i], n_atom_count_dict[struct[i].specie]), \
                    'pmg_site':struct.sites[i],\
                    'element':str((struct.sites[i]).specie)},ignore_index=True)
    return (df, n_atom_count_dict)
    ################################################################

df=create_df(struct)

molecules_list=[]
molecule=['C3','N5','N6','H11','H12','H13','H14','H15']
molecules_list.append(molecule)


### IDENTIFYING MOLECULE
## Start from any atom in a molecule, and spread the search
# Consider bond information
import pandas as pd
bond_csv=pd.read_csv('Bond_data_from_VESTA.csv')
molecule_formula='CN2H5'


# starting_atom=struct.sites[label2site_index['C1']]
# pd version
starting_atom=df[df['atom_label']=='C1']['pmg_site'].iloc[0]

# neighbors is a list of pmg_sites
neighbors=struct.get_neighbors(starting_atom,r=3)

### Let's make dictionaly with C-N, C-H entires.
### Make the data doulbe so that it can contain N-C, H-C entries with same bond length.
bond_dict=dict()
bond_data_file=open('Bond_data_from_VESTA.dat','r')
data= bond_data_file.readlines()
for i,line in enumerate(data):
    # line example
    # 'Zr     H    0.00000    2.29004  0  1  1  0  1'
    temp=line.split()
    bond_dict[temp[0]+'-'+temp[1]]=float(temp[3])
    bond_dict[temp[1]+'-'+temp[0]]=float(temp[3])

# list of atom labels composing a molecule (done for bond searching)
molecule=[] # ex) molecule=['C1','N1','N2','H1','H2']
# list of atom labels (neighbor searching needs to be done before added to molecule list)
atoms_to_search=[]  # ex) molecule=['C1','N1','N2','H1','H2']

# Add starting atom to the atoms_to_search list
# atoms_to_search.append(site_index2label[struct.index(starting_atom)])
# Pandas version
atoms_to_search.append(df[df['pmg_site']==starting_atom]['atom_label'].iloc[0])

# loop_flag will be 0 when there is not atoms to search in atoms_to_search
loop_flag=1
while loop_flag == 1:
    ## FINDING molecule by searching bond from a starting atom
    # In each loop, we search around A1 atom, called neighbors (A2)
    # If the discance between A1 and A2 is shorter than the bond length (defined in VESTA setting),
    # we identify A2 as part of molecule.
    # Then A2 will be the A1 in the next iteration    
    
    # A1_label is like 'C1'. 
    # A1 is pmg_site object
    # A1_element is string of element, like 'C'
    A1_label=atoms_to_search[0]
    A1=df[df['atom_label']==A1_label]['pmg_site'].iloc[0]
    A1_element=str(A1.specie)
    # A1=struct.sites[label2site_index[A1_label]]
    
    # neighbors is searched with radius of 3 ang. This value is conventionally enough
    #neighbors=struct.get_neighbors(struct.sites[label2site_index[A1_label]],r=3.0)
    # pd version
    neighbors=struct.get_neighbors(A1,r=3.0)
    
    for A2 in neighbors:
        # if the distance between two atoms is less then defined bond length,
        # then I add it to atoms_to_search
        
        # This line is required as A2 is often indicated to image outside of the cell
        # Otherwise, sometime this line rasis ValueError
        A2.to_unit_cell(in_place=True)
        A2_label = df[df['pmg_site']==A2]['atom_label'].iloc[0]
        
        # There are several conditions before adding to atoms_to_search
        # 1. the A1-A2 distance needs to be shorter than the bond information
        # 1-1. 'A1-A2' entry should be in bond_dict dictionary.
        #       If not, add the information in 'Bond_data_from_VESTA.dat' file.
        # 2. A2 should not be already in molecule or atoms_to_search list.
        # 
        if (A1_element+'-'+str(A2.specie)) in bond_dict.keys():
            if A1.distance(A2) < \
                bond_dict[A1_element+'-'+str(A2.specie)] and \
                A2_label not in molecule and \
                A2_label not in atoms_to_search:
                # Then, add to atoms_to_search
                atoms_to_search.append(A2_label)
                
    # After searching around A1 is done, add A1 to molecule member
    molecule.append(A1_label)
    # Then remove A1 from atoms_to_search list. The next item in the list will be A1 in the next iteration
    atoms_to_search.remove(A1_label)
    # Modify loop_flag to 0 if there is no item in atoms_to_search.
    if len(atoms_to_search) == 0:
        loop_flag=0

'''      
    Maybe I can further automate this process by starting with C atoms
    Molecules - list of followings:
        molecule - list of atomic labels
 Maybe I have to use pandas to connect label, site, and index?
 element as well
'''

'''
    Molecule data bases?
    LA, FA, MA, etc..?
    and possibly store their center of mass?:
        Some molecules used for RP phases, the center of mass would not be very good idea
    1. Center of mass from mass-weighted average
    2. Volumetric center of mass..? by using convex hull?
'''


# axis_vector and rotation angle
molecule=['C1', 'N2', 'H1', 'N1', 'H4', 'H5', 'H3', 'H2']
reference_point=(df[df['atom_label']=='C1']['pmg_site'].iloc[0]).frac_coords
axis_vector=np.array([1,1,0])  ## To be modified
angle_degree = 90 # degree
#struct=molecule_rotation(struct,molecule,axis_vector,angle_degree,reference_point)

reference_point=(df[df['atom_label']=='C3']['pmg_site'].iloc[0]).frac_coords
molecule=['C3','N5','N6','H11','H12','H13','H14','H15']
axis_vector=np.array([-1,1,0])  ## To be modified
angle_degree = 90 # degree
#struct=molecule_rotation(struct,molecule,axis_vector,angle_degree,reference_point)



struct.to(fmt='poscar',filename='test2.vasp')


# Translation
translation_vector=[0,0,0.2]
translation_vector=np.array(translation_vector)
translated_struct=copy.deepcopy(struct)
for atom in molecule:
    coords = copy.deepcopy((df[df['atom_label']==atom]['pmg_site'].iloc[0]).frac_coords)
    translated_coords = coords + translation_vector
    translated_struct.sites[df[df['atom_label']==atom]['site_index'].iloc[0]].frac_coords=translated_coords
    
#translated_struct.to(fmt='poscar',filename='test2.vasp')


## Write_POSCAR
# Write_POSCAR('write_test.vasp',struct)

# copy molecule
"""
Fucntionality:
    0. Remove a cation
    struct.remove_sites([0])
    struct.remove_sites([df[df['atom_label']==atom]['site_index'].iloc[0]])
    1. copy to specific site (copy and translate)
        
        add to n_atom_cound_dict
        struct.sites.append(A1)
    2. 

"""

temp=['C1','C2','C3','C4']
#for carbon in temp:
#    molecules_list.append(molecule_finder(struct,carbon,df))
    
    
# What if organic molecule has multiple carbon atoms?
carbon_list=df[df['element']=='C']['atom_label']
carbon_dict=dict()
# dictionary of bool. True: added to molecule. False: not added to molecule yet
for carbon in carbon_list:
    carbon_dict[carbon]=False
for carbon in carbon_list:
    if carbon_dict[carbon] == False:    
        molecule=molecule_finder(struct,carbon,df)
    else:
        continue
    
    for atom in molecule:
        if atom in carbon_dict.keys():
            carbon_dict[atom] = True
    molecules_list.append(molecule)
    
def list_all_molecules(pmg_struct):
    """

    Parameters
    ----------
    pmg_struct : pymatgen Structure
        VASP POSCAR structure with organic molecule.

    Returns
    -------
    molecules : list
        List of molecule. molecule is a list of atom labels.

    """
    df=create_df(struct)
    # Make dictionary of carbon atoms
    carbon_list=df[df['element']=='C']['atom_label']
    # Each item is True or False. True if it is already captured by molecule finder
    carbon_dict=dict()
    for carbon in carbon_list:
        carbon_dict[carbon]=False
    
    
    molecules=[]
    for carbon in carbon_list:
        if carbon_dict[carbon] == False:    
            # Calls molecule finder based on carbon atom
            molecule=molecule_finder(struct,carbon,df)
        else:
            continue
        
        for atom in molecule:
            if atom in carbon_dict.keys():
                carbon_dict[atom] = True
        molecules.append(molecule)
    return molecules

