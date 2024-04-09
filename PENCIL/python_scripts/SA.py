#**********************************************************************************************************************************************#
#***                                                                                                                                        ***#
#***              Script that suggests how to create a CANVAS model compatible with the output of propre and excogito analysis              ***#
#***                                                                                                                                        ***#
#**********************************************************************************************************************************************#


# Importing standard modules
import sys
import os
import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functools import partial
from multiprocessing import Pool, cpu_count


# Looking for excogito_2_canvas directory ...
python_scripts_dir = os.path.dirname(os.path.realpath(__file__)).split("/")
del python_scripts_dir[0]
del python_scripts_dir[-1]
excogito_2_canvas_dir = ""
for x in python_scripts_dir : excogito_2_canvas_dir = excogito_2_canvas_dir + "/" + x

# ... and importing my modules from "lib" directory
sys.path.append(excogito_2_canvas_dir)
from lib.read_excogito_file import *
from lib.check_heavy_atoms import *
from lib.input_control import *
from lib.write_output_files import *
from lib.simulated_annealing import *
from lib.manipulate_amminoacid import *
from lib.read_gromacs import *

# record start time
start_time = time.perf_counter()



# STEP 0: Reading input files from command line 

#         Necessary input files:
#           1) .gro file of the protein under study
#           2) Path to the directory with optimized excogito output files
#           3) Code with which the output files will be indicated (e.g. the file for coloring vmd residues)

#         Optional input files:
#           4) Path to directory with inside read_gromacs.py module (default: "/home/fabio/thesis/canvas/lib")
#           5) Number of atoms suggested by propre for this protein 

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False) 

input_args=parser.add_argument_group("Required Arguments") 

input_args.add_argument('-c', '--code', dest='prot_code', action='store', metavar = 'str',help = argparse.SUPPRESS)                 # Mandatory
input_args.add_argument('-d', '--directory', dest='path', action='store', metavar = 'str',help = argparse.SUPPRESS)                 # Mandatory
input_args.add_argument('-g', '--gro', dest='file_gro_path', action='store', metavar = 'FILE', help = argparse.SUPPRESS)            # Mandatory

input_args.add_argument('-i', '--initial-T', dest='T_0', action="store", metavar = 'float', help = argparse.SUPPRESS)               # Optional
input_args.add_argument('-n', '--num-propre', dest='num_propre', action="store", metavar = 'int', help = argparse.SUPPRESS)         # Optional
input_args.add_argument('-p', '--parallel', dest='parallel_bool', action="store", metavar = 'int', help = argparse.SUPPRESS)        # Optional
input_args.add_argument('-r', '--radius', dest='radius_canvas', action="store", metavar = 'int', help = argparse.SUPPRESS)          # Optional
input_args.add_argument('-s', '--steps', dest='SA_steps', action="store", metavar = 'int', help = argparse.SUPPRESS)                # Optional
input_args.add_argument('-t', '--decay-time', dest='tau', action="store", metavar = 'int', help = argparse.SUPPRESS)                # Optional

input_args.add_argument('-h', '--help', action='help', help = argparse.SUPPRESS)                                                    # Optional


# ------------------------------------------------------------------------------------------------------------------------------------------

# Printing help message if the script does not have any arguments  
if len(sys.argv)==1:
    print_help()
    sys.exit(1)

# Check if help option is used
if(sys.argv[1].strip() == "--help" or sys.argv[1].strip() == "-h"):
    print_help()
    quit()

# Printing help message if the script does not present valid arguments (e.g. "--")
check_argv_errors()

# Parsing arguments, printing error and help message if the script presents more than one task, and checking if mandatory files are present  
try:
    args  = parser.parse_args()
except SystemExit:
    print("\nOnly one task is allowed. Check that each flag (e.g. -g) is followed by its specific argument")
    print("Look below for more information.\n")
    print_help() 
    sys.exit()  

file_gro_path   = args.file_gro_path                            # Mandatory
path            = args.path                                     # Mandatory 
prot_code       = args.prot_code                                # Mandatory

T_0             = args.T_0                                      # Optional
num_propre      = args.num_propre                               # Optional
parallel        = args.parallel_bool                            # Optional
radius_canvas   = args.radius_canvas                            # Optional
SA_steps        = args.SA_steps                                 # Optional
tau             = args.tau                                      # Optional

# Check if mandatory files are present and exist in the system
mandatory_files_present(file_gro_path, path, prot_code)

# If T_0 is present, check that it is a positive float
T_0 = check_float_T_0_opt_arg(T_0)

# If num_propre is present, check that it is a positive integer
num_propre = check_int_num_propre_opt_arg(num_propre, file=0)

# If parallel is present, check that it is 0 or 1
parallel = check_parallel_bool_opt_arg(parallel)

# If radius_canvas is present, check that it is a positive float
radius_canvas = check_float_radius_canvas_opt_arg(radius_canvas, file=0)

# If SA steps is present, check that it is a positive integer
SA_steps = check_int_SA_steps_opt_arg(SA_steps)

# If tau is present, check that it is a positive float
tau = check_float_tau_opt_arg(tau)

# Create a directory to store output files with the path: _path-to-github-dir_/output_files/prot_code
output_dir = excogito_2_canvas_dir + "/" + "output_files" + "/" + prot_code
if not os.path.isdir(output_dir) :
    os.system(f"mkdir {output_dir}")



# STEP 1: Read GRO of the protein,save positions of atoms, atom types and atom numbers

# Import module to read .gro file from the directory specified in input or from my CANVAS directory

# Reading .gro file
with open(file_gro_path, "r") as myfile :
    file_gro = readgro_all(myfile)

# Removing hydrogen atoms from the list of atoms
heavy_atoms = []                    # List of heavy atoms
for atom in file_gro[0] :
    if atom[2][0] != "H"  :
        heavy_atoms.append(atom)
N_heavy_atoms = len(heavy_atoms)                   # Number of heavy atoms (no hydrogen)

# Checking output via mdtraj
#check_output_heavy_atoms(heavy_atoms, N_heavy_atoms, file_gro_path)
# NOTE: uncomment if you want to check that the .gro file is read correctly and hydrogens are all removed

# Reindexing heavy_atoms such that atom numbers go from 0 to 810
for (i,atom) in enumerate(heavy_atoms) :
    atom[3] = i



# STEP 2: Taking the files produced by excogito in the "path" directory

files_in_path = os.listdir(path)
excogito_files = []
for file in files_in_path :
    if "SA" in file and ".dat" in file:
        excogito_files.append(file)
# Check if the files have been found correctly
num_files = len(excogito_files)
if num_files == 0 : 
    print("\nERROR. There are no files in the format *SA*.dat produced by excogito in the directory given in input. Please check the directory content.\n")
    quit()

print("\nThe number of excogito files used for the calculation of probabilities is", num_files,"\n")



# STEP 3: For each file in the list excogito_files I extract the optimized mapping

# I define the vector with the probability of retaining an atom during the mapping
retain_probability = []      

for file in excogito_files :
    # NOTE: Maybe retained atoms is also not important
    binary_map, retained_atoms = read_optimized_mapping(path + "/" + file)
    retain_probability.append(binary_map)

retain_probability = np.array(retain_probability)
retain_probability = np.mean(retain_probability, axis=0)



# STEP 4: I add the retain probability at each atom in the list heavy_atoms as last entry
[atom.append(retain_probability[i]) for (i, atom) in enumerate(heavy_atoms)]



# STEP 5: I group together atoms belonging to the same residue

# Define list of residues
list_residue_with_atoms = []
# I put all the atoms belonging to the same residue in the same sublist cycling on residue number
for i in range(heavy_atoms[0][0], heavy_atoms[-1][0]+1):
    amminoacid = filter(lambda atom : atom[0] == i, heavy_atoms)
    list_residue_with_atoms.append(list(amminoacid))

# Small recap on data structure:
# list_residue_with_atoms: list that contains amminoacids (lists)
# |
# |
# --> amminoacid: list that contains atom (lists)
#           |
#           |
#           --> atom: list that contains info on atoms read from .gro file



# STEP 6: Create a list of residues. Each element in the list has the following elements:
#   1) Number of residue
#   2) Name of the residue
#   3) np.array of size 3 for the CM position
#   4) Number of heavy atoms in the residue
#   5) Average retain probability
#   6) Variance of the retain probability
#   7) Maximum probability of keeping an atom of the residue  

list_residues = []

for amminoacid in list_residue_with_atoms :
    # I build the list with the info on the residue, extracted from the consitutent atoms
    residue = []
    residue.append(amminoacid[0][0])                            # Number of residue
    residue.append(amminoacid[0][1])                            # Name of residue
    residue.append(eval_centroid(amminoacid))                   # Centroid
    residue.append(len(amminoacid))                             # Number of heavy atoms in residue
    residue.append(eval_avg_probability(amminoacid))            # Average retain probability
    residue.append(eval_var_probability(amminoacid))            # Variance of retain probability
    residue.append(max( [ atom[-1] for atom in amminoacid ] ))  # Maximum probability of keeping an atom belonging to the reisdue
    list_residues.append(residue)

print("First residue of the structure : ", list_residues[0], "\n")

# Plot of the maximum probability of keeping an atom belonging to a residue
avg_probability = []
for amminoacid in list_residues :
    avg_probability.append(amminoacid[-1]) 

plt.figure()
plt.hist(np.asarray(avg_probability), color="red", ec="orange", lw=1 )
plt.xlabel("Probability value of the most conserved atom")
plt.ylabel("Number of residues")
plt.show()

# Write files to draw the protein in VMD according to the probability
write_vmd_coloring_values(list_residues, output_dir)
write_vmd_coloring_script(output_dir)



# STEP 7: Sort the residues by probability. The sorting has two phases:
#
#        1) Sort by maximum probability: residues are sorted by the highest probability 
#           of keeping one atom of the residue
#
#        2) Sort by average probability: residues with the same value of maximum probability, are sorted by the 
#           highest average retain probability of the whole residue 
#
# *NOTE* : One could also introduce a third sorting layer based on variance or a sort of signal-to-noise ratio, but it is not necessary in my opinion

# Phase 1: Sorting by maximum single atom probability
list_residues_sorted_max = sorted(list_residues, key=lambda resid: resid[-1], reverse=True) 

# Phase 2: Subsorting by maximum average probability
# Group together in sublists the residues with the same value of maximum probability
tmp = list_residues_sorted_max[0]
list_to_be_sorted = []          # List that contains sublists of residues with equal maximum probability
sublist = []

for residue in list_residues_sorted_max :
    if residue[-1] != tmp[-1] :
        list_to_be_sorted.append(sublist)
        sublist = []
    sublist.append(residue)
    tmp = residue 

# Sort sublists
list_to_be_sorted = [ sorted(sublist, key=lambda resid: resid[-3], reverse=True) for sublist in list_to_be_sorted ]

# Restore just one list with the atoms sorted correctly
list_residues_sorted = []
for sublist in list_to_be_sorted :
    [ list_residues_sorted.append(x) for x in sublist ]

# Printing residues sorted by probability
print("Residues sorted by probability:\n")
print(pd.DataFrame(list_residues_sorted)[:25])
print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# INTERMEDIATE STEP: Converting the list list_residues_sorted in a numpy array. This is necessary in order to use numba framework to speed up calculations

# The 2D array is called "array_residues_sorted", it has dimensions (num_residues, ...) and each row has the following composition:

#   -> Column 0: Residue number
#   -> Column 1: x coordinate of the centroid position
#   -> Column 2: y coordinate of the centroid position
#   -> Column 3: z coordinate of the centroid position
#   -> Column 4: Number of heavy atoms in the residue
#   -> Column 5: Average retain probability
#   -> Column 6: Maximum probability of keeping an atom of the residue
#   -> Column 7: State of the atom (1 = AA, 0 = MG, -1 = CG) (this column is useful for the simulated annealing procedure)

# NOTE: The column 7 is initialized to -1 for all amminoacids

array_residues_sorted = np.empty( (len(list_residues_sorted), 8) )

for i in range(len(list_residues_sorted)) :
    array_residues_sorted[i][0] = list_residues_sorted[i][0]
    array_residues_sorted[i][1] = list_residues_sorted[i][2][0]
    array_residues_sorted[i][2] = list_residues_sorted[i][2][1]
    array_residues_sorted[i][3] = list_residues_sorted[i][2][2]
    array_residues_sorted[i][4] = list_residues_sorted[i][3]
    array_residues_sorted[i][5] = list_residues_sorted[i][4]
    array_residues_sorted[i][6] = list_residues_sorted[i][6]
    array_residues_sorted[i][7] = -1


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



###############################################################################################
#######################         SIMULATED ANNEALING PART            ###########################
###############################################################################################


# NOTE: This step is executed only if -n flag is supplied!

# STEP 8: SYSTEM INITIALIZATION

#         We initialize the system for the simulated annealing by taking the first atoms from the sorted list and considering them atomistic
#         (with their MG shells) until we have a number of atoms higher than num_propre
#         (this must be optained by PROPRE and given in input to the script).
#         Workflow:
#           0) Put all atoms in state -1 (i.e. CG)
#           1) Take the first residue in the sorted list and put it in state 1 (i.e. atomistic)
#           2) Put the residues closer than radius_canvas to it in state 0, if they are not in state 1 (i.e. MG)  
#           3) Evaluate the number of atoms in the hypotetical model using the formula:
#           num_at_model = #atoms_in_kept_residue + 4 * num_atoms_in_state_0 + 1 * num_atoms_in_state_-1
#           4) check if num_at_model exceeds the number of atoms given by propre
#           5) Repeat steps 1 to 4 adding one more atom every time until num_at_model > num_propre

if (num_propre is None) :

    print("\nThe simulated annealing part for the CANVAS model construction is not calculated.")
    print("To do this calculation, please supply PROPRE output via -n flag.\n")
    file_residues_sorted = "file_resid_sort.txt"
    print(f"Writing to file {output_dir}/{file_residues_sorted} the list of residues sorted by probability required to run the C routine ...")
    write_residue_sorted_file(output_dir, file_residues_sorted, array_residues_sorted)
    print("Done.\n")

    # Record end time
    end_time = time.perf_counter()
    # Calculate elapsed time and print
    elapsed_time = round((end_time - start_time), 3)
    print(f"\nElapsed time: {elapsed_time} s\n")
    quit() 

file_residues_sorted = "file_resid_sort.txt"
print(f"Writing to file {output_dir}/{file_residues_sorted} the list of residues sorted by probability, required for postprocessing...")
write_residue_sorted_file(output_dir, file_residues_sorted, array_residues_sorted)
print("Done.\n")

print("Initializing the system for the simulated annealing calculation:\n")

num_at_model = 0                    # Number of atoms in the CANVAS model
num_at_residues = 0                 # Number of residues modeled atomistic
dr = 0.1                            # Offset to take into account that here we work with centroids, CANVAS works with single-atom distances, hence the radius is underestimate here
at_residues = None                  # Atomistic residues (state 1)

print("The calculation is done considering the radius of the canvas model equal to ", radius_canvas, "nm\n")
radius_canvas += dr                 # Radius correction

# Initialize the model

while num_at_model < num_propre :

    # Increase the number of atomistic residues
    num_at_residues += 1
    
    # Reset atom count
    num_at_model = 0
    
    # I make the residue atomistic and add it to the list
    array_residues_sorted[num_at_residues-1][-1] = 1

    if at_residues is None :
        at_residues = array_residues_sorted[num_at_residues-1].copy()
    elif at_residues.ndim == 1 :
        at_residues = ( np.append(at_residues, array_residues_sorted[num_at_residues-1]) )
        at_residues = np.reshape(at_residues, ( 2, array_residues_sorted.shape[1] ) )
    else :
        at_residues = np.append(at_residues, [array_residues_sorted[num_at_residues-1]], axis=0)
    
    # I evaluate the total number of heavy atoms in the model
    num_at_model, array_residues_sorted = eval_num_at_model(array_residues_sorted, at_residues, radius_canvas)

    # Status message
    print("Number of atomistic residues : ", num_at_residues, "\tNumber of atoms in CANVAS model : ", num_at_model)

print("\nAtomistic residues in the initial configuration:\n")
print(pd.DataFrame(list_residues_sorted)[:num_at_residues])

# Evaluating average probability of the most probable atom for the cost function
avg_retain_probability = np.asarray(avg_probability).mean()



# STEP 9: SA OPTIMIZATION (probability, first neighbours and constraint on atom number)

#   Now I apply simulated annealing to the system, to find a mapping that minimizes the cost function. 
#   Workflow:
#           0) Propose a move (remove an AA residue, add a residue actually in state 0 or -1)
#           1) Recompute the MG or CG atoms based on new AA residues
#           2) Evaluate the new cost function  
#           3) Accept the move with probability exp(-dE/T), where T decreases exponentially with SA duration
#           4) If the number of sites of the model becomes lower than num_propre, add a new AA residue

# Cost function parameters
k_1 = 1/avg_retain_probability                          # Probability term
k_2 = 0.1*num_at_residues/0.04                          # Number of atoms term (compared to propre number)

# Evaluating initial cost function
E = cost_function_SA(array_residues_sorted, at_residues, num_at_model, num_propre, k_1, k_2)

# Parameters
if tau is None : tau = - SA_steps / np.log(0.001)                                                               # Initialize decay time if it is not
num_trials = 1000                                                                                               # Number of proposed moves for T_0 estimate
print("\n" + "~" * 60)
print("\nEvaluating starting temperature with", num_trials, "proposed MC moves...\n")

# Evaluating starting temperature (if not set by the user)0.1*at_res_global.shape[0]/0.04
if T_0 is None :
    T_0 = eval_start_T_numba(array_residues_sorted, at_residues, E, num_propre, radius_canvas, num_trials, k_1, k_2)

# Printing status message
print("~" * 60)
print("\nStarting simulated annealing procedure with the following parameters\n")
print(f"\t-> MC_steps =", SA_steps)
print(f"\t-> T_0 = {T_0:.3f}")
print(f"\t-> tau = {tau:.3f}\n")
print(f"\t-> k_1 = {k_1:.3f}")
print(f"\t-> k_2 = {k_2:.3f}\n")


# Performing SA optimization in serial or parallel mode, depending on parallel flag

if parallel :

    ncpus = cpu_count()
    print("Starting parallel execution of SA optimization...")
    print(f"The number of optimizations performed is {ncpus}")
    n = range(ncpus)

    pool = Pool(ncpus)

    temp = partial(SA_optimization, output_dir, SA_steps, T_0, tau, array_residues_sorted, at_residues, radius_canvas, num_propre, list_residues_sorted, E, num_at_model, k_1, k_2)
    check = pool.map(temp, iterable=n)

    pool.close()
    pool.join()

else :
    n = 0
    print("Starting serial execution of SA optimization...")
    check = SA_optimization(output_dir, SA_steps, T_0, tau, array_residues_sorted, at_residues, radius_canvas, num_propre, list_residues_sorted, E, num_at_model, k_1, k_2, n)


# Record end time
end_time = time.perf_counter()
 
# calculate elapsed time and print
elapsed_time = round((end_time - start_time), 3)
print(f"\nElapsed time: {elapsed_time} s\n")
