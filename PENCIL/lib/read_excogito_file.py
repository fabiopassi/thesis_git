""" This module contains the function to read the retained atoms from the output 
    files of excogito """

# Importing modules
import numpy as np



# Functions

def read_optimized_mapping (filename) :
    """ This function takes in input the name of a file and extracts the ideal mapping 
        produced by excogito. The function returns ...  """
    
    with open(filename, "r") as map_file:
        
        check_string = map_file.readline()

        while check_string.strip() != "last_mapping" :
            check_string = map_file.readline()

        binary_map = map_file.readline()
        retained_atoms = map_file.readline()
    
        # Removing extra words
        binary_map = binary_map.replace("conv", "")
        binary_map = binary_map.replace("mapping\n", "")
        # Converting the string of char into array of int
        binary_map = binary_map.strip().split(" ")
        binary_map = np.array([ int(bit) for bit in binary_map ])

        return binary_map, retained_atoms

