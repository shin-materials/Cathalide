# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 23:57:09 2021

@author: yongjin
"""

import numpy as np
from pymatgen.core import Structure
import copy
from scipy.spatial.transform import Rotation as R


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


def molecule_rotation(pmg_struct,molecule,label2site_index,axis_vector,angle,reference_point):
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
