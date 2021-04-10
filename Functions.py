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
