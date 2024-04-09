#**********************************************************************************************************************************************#
#***                                                                                                                                        ***#
#***                               Script that suggests one single mapping out of all the possible optimizaitons                            ***#
#***                                                                                                                                        ***#
#**********************************************************************************************************************************************#


# Importing standard modules
import sys
import os
import argparse
import time
import numpy as np


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
#           1) Code with directory in output_files where data produced by NAME_SCRIPT.py are stored
#           2) Number of atoms given by propre

#         Optional input files:
#           3) Alternative path for file with sorted residues

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False) 

input_args=parser.add_argument_group("Required Arguments") 

input_args.add_argument('-d', '--directory', dest='SA_dir', action='store', metavar = 'str',help = argparse.SUPPRESS)               # Mandatory
input_args.add_argument('-n', '--num-propre', dest='num_propre', action="store", metavar = 'int', help = argparse.SUPPRESS)         # Mandatory

input_args.add_argument('-r', '--radius', dest='radius_canvas', action="store", metavar = 'float', help = argparse.SUPPRESS)        # Optional
input_args.add_argument('-h', '--help', action='help', help = argparse.SUPPRESS)                                                    # Optional
input_args.add_argument('-f', '--txt_file', dest='file_resid', action="store", metavar = 'str', help = argparse.SUPPRESS)           # Optional

# ------------------------------------------------------------------------------------------------------------------------------------------

# Printing help message if the script does not have any arguments  
if len(sys.argv)==1:
    print_help_postprocessing()
    sys.exit(1)

# Check if help option is used
if(sys.argv[1].strip() == "--help" or sys.argv[1].strip() == "-h"):
    print_help_postprocessing()
    quit()

# Printing help message if the script does not present valid arguments (e.g. "--")
check_argv_errors()

# Parsing arguments, printing error and help message if the script presents more than one task, and checking if mandatory files are present  
try:
    args  = parser.parse_args()
except SystemExit:
    print("\nOnly one task is allowed. Check that each flag (e.g. -g) is followed by its specific argument")
    print("Look below for more information.\n")
    print_help_postprocessing() 
    sys.exit()  

SA_dir          = args.SA_dir                                   # Mandatory
num_propre      = args.num_propre                               # Mandatory

radius_canvas   = args.radius_canvas                            # Optional
file_resid      = args.file_resid                               # Optional

# If num_propre is present, check that it is a positive integer
num_propre = check_int_num_propre_opt_arg(num_propre, file=1)

# If radius_canvas is present, check that it is a positive float
radius_canvas = check_float_radius_canvas_opt_arg(radius_canvas, file=1)

# Check if mandatory directory is present and exist in the system
SA_dir = excogito_2_canvas_dir + "/" + "output_files" + "/" + SA_dir

if (SA_dir is None):
    print("\n####################################################################################")
    print("ERROR. The mandatory directory is missing.")
    print("       Look below for further help.")
    print("####################################################################################\n\n")
    print_help_postprocessing()
    quit() 

if not os.path.isdir(SA_dir) :
    print("\n####################################################################################")
    print("ERROR. The mandatory directory name is wrong. Directory not found.")
    print("       Please check the name and try again.")
    print("####################################################################################\n\n")
    quit()

# If file_resid is present, then that path is used. Otherwise the default file "file_resid_sort.txt" will be assumed to be in SA_dir
if file_resid is None :
    
    file_resid = SA_dir + "/" + "file_resid_sort.txt"

# Check that the suggested path exists and that the file is not empty

if not os.path.isfile(file_resid) :
    print("\n####################################################################################")
    print("ERROR. The residue file path is wrong. File not found.")
    print("       Please check the path and try again, or leave the default option.")
    print(f"       If you are usign the default option, make sure to have a file named file_resid_sort.txt in directory {SA_dir}")
    print("####################################################################################\n\n")
    quit()

if(os.path.getsize(file_resid) == 0):
    print("\n####################################################################################")
    print(f"ERROR. Your residue file is empty. Please, fill out {file_resid} with significant data or use another file.")
    print("####################################################################################\n\n")
    quit()



# STEP 1 : Reading the sorted residues from the file and storing them in a list
with open(file_resid, "r") as res_file :
    
    first_line = res_file.readline()
    num_prot_residues = int(first_line.replace("\n", ""))
    list_residues = []

    for line in res_file :
        line = line.split("\t")
        del line[-1]
        residue = []
        residue.append(int(line[0]))                        # Residue number
        residue.append(float(line[1]))                      # x CM
        residue.append(float(line[2]))                      # y CM
        residue.append(float(line[3]))                      # z CM
        residue.append(int(line[4]))                        # Number of atoms in the residue
        residue.append(float(line[5]))                      # Average probability
        residue.append(float(line[6]))                      # Maximum single atom probability

        list_residues.append(residue)



# STEP 2 : Reading all the optimized structures suggested from NAME_PREVIOUS_PROGRAM and save them into a list of lists (optimal structures)

files_in_path = os.listdir(SA_dir)
data_files = []

for file in files_in_path :
    if ".log" in file and "SA" in file :
        data_files.append(file)

num_files = len(data_files)

if num_files == 0 : 
    print("\nERROR. There are no files in the format SA*.log produced by PREVIOUS_SCRIPT in the directory given in input. Please check the directory content.\n")
    quit()

print(f"\nThe number of data files found in {SA_dir} is {num_files}\n")

optimal_structures = []

for file in data_files:

    with open(SA_dir + "/" + file, "r") as log_file :
        check = "START"
        while "Atomistic residues after simulated annealing" not in check :
            check = log_file.readline()
        log_file.readline()
        log_file.readline()
        
        retained_res = []

        
        for line in log_file :
            # Reading python output file
            if "\t" not in line :
                line = line.split(" ")
                while '' in line : line.remove('')
                if line[0] == '\n' : break
                retained_res.append( int(line[1]) )
            # Reading C output file
            else :
                line = line.split("\t")
                if line[0] == '\n' : break
                retained_res.append( int(line[1]) )

    optimal_structures.append(retained_res)



# STEP 3 : Find all the residues contained in AT LEAST one maping and evaluate the frequency with which they appear in all mappings

# Firstly I make all the "optimal_structures" list entries homogeneous
max_len = max( [ len(x) for x in optimal_structures ] )
for structure in optimal_structures :
    while len(structure) < max_len :
        structure.append(-1)

# Count occurrencies
map_residues, residue_frequency = np.unique(np.array(optimal_structures), return_counts=True)



# STEP 4 : Keep just as atomistic in the final mapping the residues with appearence >= 75 % and sort them, so that it is easier to fill up holes in the sequence
at_res_final_map = sorted( [ x for (i,x) in enumerate(map_residues) if residue_frequency[i]/num_files >= 0.75 and x != -1 ] )

if not at_res_final_map :
    print("\nNo atom was found with retain probability >= 75 %")
    print("The suggestion is to examine all the output mapping from PREVIOUS_PROGRAM yourself, because the variability is high.\n")
    print("We encourage you to choose one among the output mappings with one of the following criteria:")
    print("\t-> Choose the mapping with the lowest value of mapping entropy")
    print("\t-> Choose the mapping which contains regions biologically relevant (if you have this information a priori) or the mapping which contains the residue with highest probability\n")    
    quit()

print("Residues with appearence >= 75% from optimized mappings :\n")
print(f"Residues :\t{at_res_final_map}")



# STEP 5 : If in at_res_final_map there is a hole of 1 or 2 residues in a chain of sequential AA residues, I fill it up, so that all that sequence is modelled AA

holes_2 = [ at_res_final_map[i] + 2 for i in range(len(at_res_final_map) - 1) if at_res_final_map[i+1] == at_res_final_map[i] + 3 ]

at_res_map_holes_filled = sorted(at_res_final_map + holes_2)

holes_1 = [ at_res_map_holes_filled[i] + 1 for i in range(len(at_res_map_holes_filled) - 1) if at_res_map_holes_filled[i+1] == at_res_map_holes_filled[i] + 2 ]

at_res_map_holes_filled = sorted(at_res_map_holes_filled + holes_1)

print("\n\nAtomistic residues in the model after filling the holes :")
print(at_res_map_holes_filled, "\n")



# STEP 6 : If the number of atoms in this way is lower than propre number, then do the following:

# Take all the sequential chains and look at the two MG atoms before and after each chain: take the one with the highest probability and add it to the model 

# Example :

#       Residues with probability > 75% : 44 45 46 101 102 103
#       If num_at_model < num_propre, then look at the probabilities of residues 43 47 100 104: pick the one with highest probability and add it to the model
#       Continue until num_at_model > num_propres


# Building list with AA residues and assign state for num_at_model calculation

num_at_model = 0                    # Number of atoms in the CANVAS model
dr = 0.1                            # Offset to take into account that here we work with centroids, CANVAS works with single-atom distances, hence the radius is underestimate here
radius_canvas += dr
list_at_res_filename = "list_at_res.dat"                        # Name for the output file
block_py_out_file = SA_dir + "/" + list_at_res_filename         # Path to the output file
at_residues = []

for num_res_aa in at_res_map_holes_filled :
    for residue in list_residues :
        
        if residue[0] == num_res_aa : 
            residue.append(1)
            at_residues.append(residue)

for residue in list_residues :
    if residue[-1] != 1 : residue.append(-1)

# Convert list_residues and at_residues in arrays, so functions from libraries can be used
array_residues_sorted = np.array(list_residues)
at_residues = np.array(at_residues)

# Calculate the number of atoms in the model considering just AA residues found at step 4 (with holes filled) and quit if it is greater than num_propre
num_at_model, array_residues_sorted = eval_num_at_model(array_residues_sorted, at_residues, radius_canvas)

if np.rint(num_at_model) > num_propre :

    print("\nThe suggested choice of AA residues based on the optimization results is :\n")
    for i in range(at_residues.shape[0]) : print(int(at_residues[i][0]))
    print(f"\nNumber of atoms in the model : {int(num_at_model)} (propre number of atoms : {num_propre}).\n")
    write_canvas_list_at_res(at_residues, block_py_out_file)
    
else :
    # If it is not, then add AA residues with the criterion explained above until num_at_model > num_propre
    print(f"\nActual number of atoms in the model : {int(num_at_model)} ( propre number of atoms : {num_propre} )")
    print(f"\nAdding residues adjacent to AA sequences until reaching the number of atoms indicated by propre ...\n\n")

    while num_at_model < num_propre :

        # Find the number of the residues at the beginning and end of each AA sequence
        # (in the above example, we look for 
        #       -> 44 and 101 , stored in list "begin"
        #       -> 46 and 103 , stored in list "end" 
        # )
        begin = [ at_residues[0] ]
        end = [ at_residues[-1] ]

        for i in range(at_residues.shape[0]-1) :
            if np.rint(at_residues[i+1][0]) != np.rint(at_residues[i][0])+1 :
                end.append(at_residues[i])
                begin.append(at_residues[i+1])

        # Now look in the list of all the residues for the MG residues placed immediately before or after the begin or end residues
        # (in the above example, we look for 
        #       -> 43 and 100
        #       -> 47 and 104 
        # )

        candidates = []

        for extremum in begin :
            idx = np.where( np.rint(array_residues_sorted[:,0]) == np.rint(extremum[0]) - 1 )
            candidates.append(array_residues_sorted[idx[0][0]])

        for extremum in end :
            idx = np.where( np.rint(array_residues_sorted[:,0]) == np.rint(extremum[0]) + 1 )
            candidates.append(array_residues_sorted[idx[0][0]])

        # Status message
        print("Finding candidate residues...")
        num_res_candidates = [ int(x[0]) for x in candidates ]
        print(f"Candidate residues found : {num_res_candidates}")

        # Take the one with the highest probability
        best_res = np.argmax(np.array(candidates)[:,6])

        # Add it to AA residues and compute the number of atoms in the model
        at_residues = np.append(at_residues, [candidates[best_res]], axis = 0)
        at_residues = np.array( sorted(at_residues, key=lambda resid: resid[0]) )
        new_aa_res = np.where( array_residues_sorted[:,0] == np.rint(candidates[best_res][0]) )
        array_residues_sorted[new_aa_res[0][0]][-1] = 1
        print(f"Adding residue : {int(array_residues_sorted[new_aa_res[0][0]][0])}")

        # Recalculate the number of atoms in the model
        num_at_model, array_residues_sorted = eval_num_at_model(array_residues_sorted, at_residues, radius_canvas)

        # Status message
        print(f"Number of atoms in the model : {int(num_at_model)} ( propre number of atoms : {num_propre} )\n\n")

    # Print final result and output file
    print("\nThe suggested choice of AA residues based on the optimization results is :\n")
    for i in range(at_residues.shape[0]) : print(int(at_residues[i][0]))
    print(f"\nNumber of atoms in the model : {int(num_at_model)} (propre number of atoms : {num_propre}).\n")
    write_canvas_list_at_res(at_residues, block_py_out_file)


# Record end time
end_time = time.perf_counter()
 
# calculate elapsed time and print
elapsed_time = round((end_time - start_time), 3)
print(f"\nElapsed time: {elapsed_time} s\n")