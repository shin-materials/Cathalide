# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 23:57:09 2021

@author: yongjin
"""

import numpy as np
from pymatgen.core import Structure
import copy
from scipy.spatial.transform import Rotation as R
from pymatgen.core import Element
import pandas as pd

def convert_site_index(label2site_index,str_atom):
    """
    input: label2site_index -- dictionary defined in script
        str_atom -- str, atom label like Fe1, or could be user index.
    output: pmg_site_index -- int, pymatgen index, starting from 0
    User index starts from 1, while python site index starts from 0.
    Based on the format, this function returns the pymatgen index
    """
    if str_atom.isnumeric():
        # if is numeric, python index is -1 from user index.
        # user index starts from 1, python index from 0
        pmg_site_index=int(str_atom)-1
    else:
        pmg_site_index=label2site_index[str_atom]
    return pmg_site_index





def molecule_rotation(pmg_struct,molecule,label_df,axis_vector,angle,reference_point):
	"""
	input:
	pmg_struct: pymatgen structure. Structure containing organic molecule to rotate
	molecule: list of strings, containing atomic labels composing molecule to rotate
				ex) molecule=['C1','N1','N2','H1','H2','H3','H4','H5']
	label_df: Pandas DataFrame with columns=['site_index','atom_label','pmg_site','element']
	axis_vector: numpy array with size (3,). This defines the direction of rotation axis
				This function will normalize the vector.
	angle - float: degree angles to rotate the molecule with respect to the rotation axis.
				The rotation will be done anti-clockwise
	reference_point - numpy array with size (3,). This is a fractional coordinate of point 
				that the rotation axis goes through
	
	output:
	rotated_struct: pymatgen structure 
	"""
	rotated_struct=copy.deepcopy(pmg_struct)
	
	# convert reference point to cartesian
	a=pmg_struct.lattice.a
	b=pmg_struct.lattice.b
	c=pmg_struct.lattice.c
	reference_point = reference_point * np.array([a,b,c])

	## Transformation (2)
	## ROTATION VECTOR
	## Vector from np.array (should be normalized),
	## and multiply the angle of rotation (in radian)
	## Ex) r = R.from_rotvec(np.pi/2 * np.array([0, 0, 1]))
	# convert to real space axis_vector by introducing lattice vector
	
	axis_vector = axis_vector * np.array([a,b,c])
	#axis_vector=axis_vector * np.array([a,b,c])
	# normalize
	axis_vector = axis_vector/np.linalg.norm(axis_vector)
	# make degree to radian
	angle_radian = angle / 180 * np.pi
	# create rotation operator
	rotation=R.from_rotvec(angle_radian * axis_vector)

	for atom in molecule:
		# Get index from label of atom,
		# call sites from pymatgen struct
		# and parse the fractional coordinates
		#coords=copy.deepcopy((label_df[label_df['atom_label']==atom]['pmg_site'].iloc[0]).frac_coords)
		coords=copy.deepcopy((label_df[label_df['atom_label']==atom]['pmg_site'].iloc[0]).coords)
		
		## Transformation (1)
		## subtract the coordinate of rotation-axis-position  = reference point
		coords_diff = coords-reference_point
		# adjust for when the operation includes the periodic boundary
		adjust_array= -np.array([a,b,c])*(coords_diff > np.array([a*0.5,b*0.5,c*0.5])).astype(int) + \
					np.array([a,b,c])*(coords_diff < np.array([-a*0.5,-b*0.5,-c*0.5])).astype(int)
		coords_temp = coords_diff + adjust_array
		shift_vector = coords_temp - coords

		# apply rotation
		coords_temp=rotation.apply(coords_temp)
		# Move to reference point again
		coords = coords_temp - shift_vector
		rotated_struct.sites[label_df[label_df['atom_label']==atom]['site_index'].iloc[0]].coords=coords

	return rotated_struct



def molecule_translation(pmg_struct,molecule,label_df,translation_vector):
	"""
	input:
	pmg_struct: pymatgen structure. 
		Structure containing organic molecule to rotate
	molecule: list of strings
		 containing atomic labels composing molecule to rotate
				ex) molecule=['C1','N1','N2','H1','H2','H3','H4','H5']
	label2site_index: dictionary defined in script
	translation_vector: List of size 3, or numpy array with size (3,)
		This defines the vector of translation of the molecule
	angle - float: degree angles to rotate the molecule with respect to the rotation axis.
				The rotation will be done anti-clockwise
	reference_point - numpy array with size (3,). This is a fractional coordinate of point 
				that the rotation axis goes through
	
	output:
	translated_struct: pymatgen structure 

	"""
	# convert to numpy array
	translation_vector=np.array(translation_vector)

	translated_struct=copy.deepcopy(pmg_struct)
	for atom in molecule:
	    coords = copy.deepcopy((label_df[label_df['atom_label']==atom]['pmg_site'].iloc[0]).frac_coords)
	    translated_coords = coords + translation_vector
	    translated_struct.sites[label_df[label_df['atom_label']==atom]['site_index'].iloc[0]].frac_coords=translated_coords
	    
	rotated_struct=copy.deepcopy(pmg_struct)

	return translated_struct

######### FUNCTIONS ###########
def Write_POSCAR(filename,struct, label_df, element_sequence=None):
    """
    Parameters
    ----------
    filename : str
        name of POSCAR file. PatternLength_PatternNumber_TransVector.vasp
    struct : pymatgen Structure object
        structure object to write as POSCAR
    DF : Pandas dataframe
        Pandas DataFrame with columns=['site_index','atom_label','pmg_site','element']
    element_sequence : str
        Sequence of elements to be written in POSCAR
        ex) 'CNHPbI' or 'C N H Pb I' will affect the 7th line in the POSCAR file
    Returns
    -------
    None.
    """
    # Idenfity lattice vectors
    lattice = struct.lattice.matrix
    # Get element of each atoms
    species_list=struct.species
    # Remove repeated entry in the element list
    reduced_species=[str(i) for n, i in enumerate(species_list) if i not in species_list[:n]]
    
    out_file=open(filename,'w')
    out_file.write("Generated POSCAR file\n") #first comment line
    out_file.write("1.0\n") # scale 
    # Print lattice part
    for i in range(np.shape(lattice)[0]):
        out_file.write("{0:20.10f} {1:20.10f} {2:20.10f}\n".format(lattice[i,0],lattice[i,1],lattice[i,2]))
    # Print elements
    out_file.write("  "+"  ".join('%3s' % entry for entry in reduced_species))
    out_file.write("\n")
    # Print the number of atoms for each element
    num_each_element=[species_list.count(Element(i)) for i in reduced_species]
    out_file.write("  "+"  ".join('%3d' % entry for entry in num_each_element))
    out_file.write("\n")
    
    out_file.write("Direct\n")
    for element in reduced_species:
        site_list=label_df[label_df['element']=='H']['pmg_site']
        for site in site_list:
            out_file.write("  "+"        ".join('%.10f' % entry for entry in site.frac_coords)+'\n')
    out_file.close()
    
    return


def molecule_rotation_bak(pmg_struct,molecule,label2site_index,axis_vector,angle,reference_point):
	"""
	input:
	pmg_struct: pymatgen structure. Structure containing organic molecule to rotate
	molecule: list of strings, containing atomic labels composing molecule to rotate
				ex) molecule=['C1','N1','N2','H1','H2','H3','H4','H5']
	label2site_index -- dictionary defined in script
	axis_vector: numpy array with size (3,). This defines the direction of rotation axis
				This function will normalize the vector.
	angle - float: degree angles to rotate the molecule with respect to the rotation axis.
				The rotation will be done anti-clockwise
	reference_point - numpy array with size (3,). This is a fractional coordinate of point 
				that the rotation axis goes through
	
	output:
	rotated_struct: pymatgen structure 
	"""
	rotated_struct=copy.deepcopy(pmg_struct)

	for atom in molecule:
	    # Get index from label of atom,
	    # call sites from pymatgen struct
	    # and parse the fractional coordinates
	    coords=copy.deepcopy(pmg_struct.sites[label2site_index[atom]].frac_coords)
	    
	    ## Transformation (1)
	    ## subtract the coordinate of rotation-axis-position  = reference point
	    coords_temp=coords-reference_point
	    
	    ## Transformation (2)
	    ## ROTATION VECTOR
	    ## Vector from np.array (should be normalized),
	    ## and multiply the angle of rotation (in radian)
	    ## Ex) r = R.from_rotvec(np.pi/2 * np.array([0, 0, 1]))

	    # normalize
	    axis_vector=axis_vector/np.linalg.norm(axis_vector)
	    # make degree to radian
	    angle_radian = angle / 180 * np.pi
	    rotation=R.from_rotvec(angle_radian * axis_vector)
	    # apply rotation
	    coords_temp=rotation.apply(coords_temp)
	    # Move to reference point again
	    coords=coords_temp+reference_point
	    rotated_struct.sites[label2site_index[atom]].frac_coords=coords

	return rotated_struct

def molecule_finder(pmg_struct,starting_atom_label,label_df):
	### IDENTIFYING MOLECULE
	## Start from any atom in a molecule, and spread the search
	# Consider bond information
	
	bond_csv=pd.read_csv('Bond_data_from_VESTA.csv')
	# molecule_formula='CN2H5'


	# starting_atom=pmg_struct.sites[label2site_index['C1']]
	# pd version
	#starting_atom=label_df[label_df['atom_label']==starting_atom_label]['pmg_site'].iloc[0]

	# neighbors is a list of pmg_sites
	#neighbors=pmg_struct.get_neighbors(starting_atom,r=3)

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
	# atoms_to_search.append(site_index2label[pmg_struct.index(starting_atom)])
	#atoms_to_search.append(label_df[label_df['pmg_site']==starting_atom]['atom_label'].iloc[0])
	atoms_to_search.append(starting_atom_label)
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
	    A1=label_df[label_df['atom_label']==A1_label]['pmg_site'].iloc[0]
	    A1_element=str(A1.specie)
	    # A1=pmg_struct.sites[label2site_index[A1_label]]
	    
	    # neighbors is searched with radius of 3 ang. This value is conventionally enough
	    #neighbors=pmg_struct.get_neighbors(pmg_struct.sites[label2site_index[A1_label]],r=3.0)
	    # pd version
	    neighbors=pmg_struct.get_neighbors(A1,r=3.0)
	    
	    for A2 in neighbors:
	        # if the distance between two atoms is less then defined bond length,
	        # then I add it to atoms_to_search
	        
	        # This line is required as A2 is often indicated to image outside of the cell
	        # Otherwise, sometime this line rasis ValueError
	        A2.to_unit_cell(in_place=True)
	        A2_label = label_df[label_df['pmg_site']==A2]['atom_label'].iloc[0]
	        
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
	return molecule
	# Pandas version