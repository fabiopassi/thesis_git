""" This module contains the functions to evaluate the centroid and the average retain
    probability of an amminoacid """

# Importing modules
import numpy as np
from numba import jit



# Dictionary that takes in input the element letter and returns the approximated mass in u.m.a.

elem_to_mass = {
    'C' : 12,
    'N' : 14,
    'O' : 16,
    'S' : 32
}



# Functions

def eval_centroid(amminoacid) :
    """ This function takes in input a list of atoms (i.e. a list of lists) and evaluates their
        mass-weighted average positions """
    
    avg_position = np.zeros(3)
    mass_amminoacid=0

    for atom in amminoacid :   
        mass = elem_to_mass[atom[2][0]]         # Assign the correct mass to a certain atom
        mass_amminoacid += mass
        avg_position[0] += atom[4] * mass       # x
        avg_position[1] += atom[5] * mass       # y
        avg_position[2] += atom[6] * mass       # z

    return avg_position/(mass_amminoacid)



def eval_avg_probability(amminoacid):
    """ This function takes in input a list of atoms (i.e. a list of lists) and evaluates their
        average retain probability """
    
    avg_probability = 0

    for atom in amminoacid :
        avg_probability += atom[-1]

    return avg_probability/len(amminoacid)



def eval_var_probability(amminoacid):
    """ This function takes in input a list of atoms (i.e. a list of lists) and evaluates the
        variance of their retain probability """
    
    avg_probability = eval_avg_probability(amminoacid)
    var_probability = 0

    for atom in amminoacid :
        var_probability += np.power(atom[-1] - avg_probability ,2)

    return var_probability/len(amminoacid)



@jit(nopython=True)
def eval_shell_numba(at_residues, array_residues, radius) :
    """ Function that takes in input two 2D-np arrays and a float (radius).
        This function returns the total number of amminoacids in array_residues whose distance from one
        component of at_residues is < radius, excluding other atomistic residues """
    
    num_mg_residues = 0

    if at_residues.ndim == 1 :
        for j in range(array_residues.shape[0]) :
            if np.rint(array_residues[j][-1]) == -1 :
                dx = array_residues[j][1] - at_residues[1]
                dy = array_residues[j][2] - at_residues[2]
                dz = array_residues[j][3] - at_residues[3]
                dr = np.sqrt( np.square(dx) + np.square(dy) + np.square(dz) )
                if dr < radius : 
                    num_mg_residues += 1
                    array_residues[j][-1] = 0

    else :
        for i in range(at_residues.shape[0]) :
            for j in range(array_residues.shape[0]) :
                if np.rint(array_residues[j][-1]) == -1 :
                    dx = array_residues[j][1] - at_residues[i][1]
                    dy = array_residues[j][2] - at_residues[i][2]
                    dz = array_residues[j][3] - at_residues[i][3]
                    dr = np.sqrt( np.square(dx) + np.square(dy) + np.square(dz) )
                    if dr < radius : 
                        num_mg_residues += 1
                        array_residues[j][-1] = 0
        

    return num_mg_residues, array_residues



@jit(nopython = True)
def eval_num_at_model(arr_resid, at_res, radius_canvas) :
    """ Function to evaluate the number of atoms in the model """
    
    # Reset all MG and CG atoms to CG
    for i in range(arr_resid.shape[0]) :
        if np.rint(arr_resid[i][-1]) != 1 : arr_resid[i][-1] = -1

    # Calculate the number of atoms from AA, MG and CG residues and sum them up
    if at_res.ndim <= 1 :
        num_retained_atoms = np.rint(at_res[4])
    else :
        num_retained_atoms = np.rint(at_res[:,4]).sum()
    
    num_mg_residues, arr_resid = eval_shell_numba(at_res, arr_resid, radius_canvas)
    num_cg_residues = arr_resid.shape[0] - num_mg_residues - at_res.shape[0]
    num_at_model = num_retained_atoms + 4 * num_mg_residues + num_cg_residues

    return num_at_model, arr_resid