from pymatgen.core import Structure
from pymatgen.core import Element
from glob import glob
import numpy as np
from scipy.spatial.transform import Rotation as R
from Functions import molecule_rotation, convert_site_index
import copy
import pandas as pd

struct = Structure.from_file("Test_structures/FPB_bulk_cubic.vasp")
#struct = Structure.from_file("test.vasp")

# Read bond data
# Input data from user: what elements composing the molecule
elements_list=['C','H','N']



# When replacing inorganic cation to organic molecule,
# I need to operate with the center of mass
# command: Element.C.atomic_mass

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
################################################################


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
reference_point=(df[df['atom_label']=='C1']['pmg_site'].iloc[0]).frac_coords
axis_vector=np.array([0,1,0])  ## To be modified
angle_degree = 45 # degree
struct=molecule_rotation(struct,molecule,df,
                          axis_vector,angle_degree,reference_point)
struct.to(fmt='poscar',filename='test2.vasp')


# Translation
translation_vector=[0,0,0.2]
translation_vector=np.array(translation_vector)
translated_struct=copy.deepcopy(struct)
for atom in molecule:
    coords = copy.deepcopy((df[df['atom_label']==atom]['pmg_site'].iloc[0]).frac_coords)
    translated_coords = coords + translation_vector
    translated_struct.sites[df[df['atom_label']==atom]['site_index'].iloc[0]].frac_coords=translated_coords
    
translated_struct.to(fmt='poscar',filename='test2.vasp')



######### FUNCTIONS ###########
def Write_POSCAR(filename,pmg_struct,element_sequence=None):
    """
    Parameters
    ----------
    filename : str
        name of POSCAR file. PatternLength_PatternNumber_TransVector.vasp
    CN_str : str
        Sequence of numbers indicating coordination number of B-cations
    lattice : numpy array
        3x3 array describing lattice vectors
    A_array : numpy array
        n x 3 array; atomic position of A-caitons.
    B_array : numpy array
        n x 3 array; atomic position of B-caitons.
    O_array : numpy array
        list of lists. Sub-lists have 3 numbers representing position of oxygen atoms
    A : str
        Name of A-cation if not defined, Sr is used.
    B : str
        Name of B-cation if not defined, Fe is used.
    Returns
    -------
    None.
    """
    
    out_file=open(filename,'w')
    out_file.write("Generated POSCAR file\n") #first comment line
    out_file.write("1.0\n") # scale 
    for i in range(np.shape(lattice)[0]):
        out_file.write("{0:20.10f} {1:20.10f} {2:20.10f}\n".format(lattice[i,0],lattice[i,1],lattice[i,2]))
    out_file.write("   {0}   {1}   O\n".format(A,B)) # Atom labels
    out_file.write("   {0}    {1}   {2}\n".format(np.shape(A_array)[0],
                                              np.shape(B_array)[0],
                                              np.shape(O_array)[0]))
    out_file.write("   Direct\n")
    for i in range(np.shape(A_array)[0]):
        out_file.write("{0:16.10f} {1:19.10f} {2:19.10f}\n".format(A_array[i,0],A_array[i,1],A_array[i,2]))
    # for i,B_site in enumerate(B_list):
    #     out_file.write("    "+"        ".join('%.10f' % entry for entry in B_site))
    #     out_file.write("\n")
    #     out_file.close()
    for i in range(np.shape(B_array)[0]):
        out_file.write("{0:16.10f} {1:19.10f} {2:19.10f}\n".format(B_array[i,0],B_array[i,1],B_array[i,2]))
    for i in range(np.shape(O_array)[0]):
        out_file.write("{0:16.10f} {1:19.10f} {2:19.10f}\n".format(O_array[i,0],O_array[i,1],O_array[i,2]))
    # Printing oxygen sites is a bit different as it is with list format
    # for i,O_site in enumerate(O_list):
    #     out_file.write("    "+"        ".join('%.10f' % entry for entry in O_site))
    #     out_file.write("\n")
    #     out_file.close()




    ## Definition of lattice vectors
# lattice = np.array([[a*np.sqrt(2), 0, 0],
#                     [0, a*np.sqrt(6), 0],
#                     [a*np.sqrt(2)*((t*3)%1),a*np.sqrt(6)*t,a/np.sqrt(3)]])
lattice = struct.lattice.matrix
out_file=open('test.vasp','w')
out_file.write("Generated POSCAR file\n") #first comment line
out_file.write("1.0\n") # scale 
for i in range(np.shape(lattice)[0]):
    out_file.write("{0:20.10f} {1:20.10f} {2:20.10f}\n".format(lattice[i,0],lattice[i,1],lattice[i,2]))
species_list=struct.species
reduced_species=[str(i) for n, i in enumerate(species_list) if i not in species_list[:n]]
out_file.write("  "+"  ".join('%3s' % entry for entry in reduced_species))
out_file.write("\n")
num_each_element=[species_list.count(Element(i)) for i in reduced_species]
out_file.write("  "+"  ".join('%3d' % entry for entry in num_each_element))
out_file.write("\n")
out_file.write("Direct\n")
for element in reduced_species:
    site_list=df[df['element']=='H']['pmg_site']
    for site in site_list:
        out_file.write("  "+"        ".join('%.10f' % entry for entry in site.frac_coords)+'\n')
out_file.close()
## Write_POSCAR
#Write_POSCAR('write_test.vasp',struct)
    
    
        
    