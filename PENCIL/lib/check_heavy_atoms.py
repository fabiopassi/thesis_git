""" This module contains a function to check that the extraction of heavy atoms from the .gro file is 
    compatible with MDTraj library results. """

# Importing modules
import pandas as pd
import mdtraj

# Function

def check_output_heavy_atoms(heavy_atoms, N_heavy_atoms, file_gro_path) :

    # Printing output of excogito_2_canvas.py in pandas dataframe format (more readable)
    print("Output of excogito_2_canvas:\n")
    print("The number of excogito_2_canvas heavy atoms is ", N_heavy_atoms)
    print("Heavy atoms :\n", pd.DataFrame(heavy_atoms))
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    print("Output of mdtraj:\n")

    # Printing the same procedure but obtained using MDTraj
    mdtraj_gro = mdtraj.load(file_gro_path)
    mdtraj_gro_topology = mdtraj_gro.topology
    no_h = mdtraj_gro_topology.select('type != H')      # Removing hydrogens
    print("The number of mdtraj heavy atoms is ", len(no_h))
    mdtraj_gro_heavy = mdtraj.load(file_gro_path, atom_indices = list(no_h))    # Load heavy atoms
    table, bonds = mdtraj_gro_heavy.topology.to_dataframe()         # Writing in dataframe style to compare easily the outputs
    print(table.head())
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    return 0
