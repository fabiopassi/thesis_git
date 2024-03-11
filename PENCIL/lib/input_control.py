# Importing modules
import os
import sys



####################################        Checking flags        ####################################

def check_argv_errors():
    """ Function to check if non accepted flags are present """

    for i in range(1, len(sys.argv)):

        if(i%2 + 1 == 0):

            if(sys.argv[i][0] == '-' and len(sys.argv[i]) == 1):
                print("Error. '-' is not accepted as flag. Use, for example, '-g' instead of '-'. Look below for further help.\n");
                print_help()
                quit()

            if(sys.argv[i][0] == '-' and sys.argv[i][1] != '-' and  len(sys.argv[i]) > 2):
                print("Error. Each flag must contain '-' plus ONLY ONE letter. Example: -g. Look below for further help.\n");
                print_help()
                quit()

            if(sys.argv[i][0] == '-' and sys.argv[i][1] == '-' and  len(sys.argv[i]) == 2):
                print("Error. '--' is not an allowed flag, it must be followed by a word. Example: --gro")
                print_help()
                quit()



####################################        Checking mandatory arguments        ####################################

def mandatory_files_present(file_gro_path, path, prot_code): 
    """ Function to check if all mandatory arguments are present and if files exist"""

    # Check if a file path was passed
    if (file_gro_path is None):
        print("\n####################################################################################")
        print("ERROR. The Coordinate file is missing.")
        print("       Look below for further help.")
        print("####################################################################################\n\n")
        print_help()
        quit() 

    # Check for file existence in the machine
    if not os.path.isfile(file_gro_path) :
        print("\n####################################################################################")
        print("ERROR. The Coordinate file path is wrong. File not found.")
        print("       Please check the path and try again.")
        print("####################################################################################\n\n")
        quit()

    if(os.path.getsize(file_gro_path) == 0):
        print("\n####################################################################################")
        print(f"ERROR. Your coordinate file is empty. Please, fill out {file_gro_path} with significant data or use another file.")
        print("####################################################################################\n\n")
        quit()

    # Check if a directory path was passed
    if (path is None):
        print("\n####################################################################################")
        print("ERROR. The path to the optimized excogito files is missing.")
        print("       Look below for further help.")
        print("####################################################################################\n\n")
        print_help()
        quit() 

    # Check for directory existence in the machine
    if not os.path.isdir(path) :
        print("\n####################################################################################")
        print("ERROR. The excogito directory path is wrong. Directory not found.")
        print("       Please check the path and try again.")
        print("####################################################################################\n\n")
        quit() 

    if (prot_code is None) :
        print("\n####################################################################################")
        print("ERROR. The code of the protein with which output will be named is missing.")
        print("       Look below for further help.")
        print("####################################################################################\n\n")
        print_help()
        quit()



####################################        Utilities        ####################################

def isIntFloat(x):
    """ It checks if x is an integer or a float """
    try:
        float(x)
        return True
    except ValueError:
        return False
        

def check_Int_Float_Str(x):
    """ It checks if x if an integer, a float or a string """

    if(isIntFloat(x)):
        try:
            x = int(x)
        except:
            x = float(x)
    return x



####################################        Checking optional arguments        ####################################


def check_float_T_0_opt_arg(T_0) :
    """ Function to check that the T_0 optional argument is a float """
    
    if(T_0 is not None):   
        T_0 = check_Int_Float_Str(T_0)
              
        if(not isinstance(T_0, float)):
            print("\n####################################################################################")
            print("ERROR. '-i/--initial-T {}' set, but not recognized. Only a positive float number number is allowed. 'strings', or negative numbers are not permitted.".format(T_0))
            print("       Please, insert a positive float number (-i/--initial-T <FLOAT>) or ignore this flag.")
            print("       Look below for further help.") 
            print("####################################################################################\n\n")
            print_help()
            quit()
        else:
            if(T_0 <= 0):
                print("\n####################################################################################")
                print("ERROR. '-i/--initial-T {}' set, but not recognized. Only a float number higher than 0 is allowed".format(T_0))  
                print("       A negative number or zero are meaningless.")
                print("       Please, insert a positive float number (-i/--initial-T <FLOAT>) or ignore this flag.") 
                print("       Look below for further help.")
                print("####################################################################################\n\n")
                print_help()
                quit()

    return T_0



def check_int_num_propre_opt_arg(num_propre, file) :
    """ Function to check that the num_propre optional argument is an integer
        The variable file depends on which script calls the function:
            0 : SCRIPT SENZA NOME
            1 : post_processing.py
     """
    
    if(num_propre is not None):   
        num_propre = check_Int_Float_Str(num_propre)
              
        if(not isinstance(num_propre, int)):
            print("\n####################################################################################")
            print("ERROR. '-n/--num_propre {}' set, but not recognized. Only a positive integer number is allowed. 'float', 'strings', or negative numbers are not permitted.".format(num_propre))
            print("       Please, insert a positive integer number (-n/--num_propre <INT>) or ignore this flag.")
            print("       Look below for further help.") 
            print("####################################################################################\n\n")
            if file == 0 : print_help()
            if file == 1 : print_help_postprocessing()
            quit()
        else:
            if(num_propre <= 0):
                print("\n####################################################################################")
                print("ERROR. '-n/--num_propre {}' set, but not recognized. Only a positive integer number higher than 0 is allowed".format(num_propre))  
                print("       A negative number or zero are meaningless.")
                print("       Please, insert a positive integer number (-n/--num_propre <INT>) or ignore this flag.") 
                print("       Look below for further help.")
                print("####################################################################################\n\n")
                if file == 0 : print_help()
                if file == 1 : print_help_postprocessing()
                quit() 
    
    return num_propre



def check_parallel_bool_opt_arg(parallel_bool):
    """ Function to check that parallel_bool optional argument is only 0 or 1 """

    if(parallel_bool is not None):
        parallel_bool = check_Int_Float_Str(parallel_bool)

        if(not isinstance(parallel_bool, int)):  # if flag_mdtraj is NOT integer (float or string)
            print("\n####################################################################################")
            print("ERROR. '-p/--parallel {}' set, but not recognized. Only 0 or 1 are allowed. 'float', 'strings', or negative numbers are not permitted.".format(parallel_bool))
            print("       Please, insert 0 or 1 (-p/--parallel <INT>) or ignore this flag leaving the default value 1 (-p/--parallel 1).")
            print("       Look below for further help.") 
            print("####################################################################################\n\n")
            print_help()
            quit()
        else:
            if(parallel_bool != 0 and parallel_bool != 1) :
                print("\n####################################################################################")
                print("ERROR. '-p/--parallel {}' set, but not recognized. Only 0 or 1 are allowed".format(parallel_bool))  
                print("       Different integers are meaningless.")
                print("       Please, insert 0 or 1 (-p/--parallel <INT>) or ignore this flag leaving the default value of 1 (-p/--parallel 1).") 
                print("       Look below for further help.")
                print("####################################################################################\n\n")
                print_help()
                quit()

    if(parallel_bool is None):
        parallel_bool = 1
    
    return parallel_bool



def check_float_radius_canvas_opt_arg(radius_canvas, file) :
    """ Function to check that the radius_canvas optional argument is an float
        The variable file depends on which script calls the function:
            0 : SCRIPT SENZA NOME
            1 : post_processing.py
    """
    
    if(radius_canvas is not None):   
        radius_canvas = check_Int_Float_Str(radius_canvas)
              
        if(not isinstance(radius_canvas, float)):
            print("\n####################################################################################")
            print("ERROR. '-r/--radius {}' set, but not recognized. Only a positive float number number is allowed. 'strings', or negative numbers are not permitted.".format(radius_canvas))
            print("       Please, insert a positive float number (-r/--radius <FLOAT>) or ignore this flag.")
            print("       Look below for further help.") 
            print("####################################################################################\n\n")
            if file == 0 : print_help()
            if file == 1 : print_help_postprocessing()
            quit()
        else:
            if(radius_canvas <= 0):
                print("\n####################################################################################")
                print("ERROR. '-r/--radius {}' set, but not recognized. Only a float number higher than 0 is allowed".format(radius_canvas))  
                print("       A negative number or zero are meaningless.")
                print("       Please, insert a positive float number (-r/--radius <FLOAT>) or ignore this flag.") 
                print("       Look below for further help.")
                print("####################################################################################\n\n")
                if file == 0 : print_help()
                if file == 1 : print_help_postprocessing()
                quit()

    if(radius_canvas is None) :
        radius_canvas = 0.7

    return radius_canvas



def check_int_SA_steps_opt_arg(SA_steps) :
    """ Function to check that the SA_steps optional argument is an integer """
    
    if(SA_steps is not None):   
        SA_steps = check_Int_Float_Str(SA_steps)
              
        if(not isinstance(SA_steps, int)):
            print("\n####################################################################################")
            print("ERROR. '-s/--steps {}' set, but not recognized. Only a positive integer number is allowed. 'float', 'strings', or negative numbers are not permitted.".format(SA_steps))
            print("       Please, insert a positive integer number (-s/--steps <INT>) or ignore this flag.")
            print("       Look below for further help.") 
            print("####################################################################################\n\n")
            print_help()
            quit()
        else:
            if(SA_steps <= 0):
                print("\n####################################################################################")
                print("ERROR. '-s/--steps {}' set, but not recognized. Only a positive integer number higher than 0 is allowed".format(SA_steps))  
                print("       A negative number or zero are meaningless.")
                print("       Please, insert a positive integer number (-s/--steps <INT>) or ignore this flag.") 
                print("       Look below for further help.")
                print("####################################################################################\n\n")
                print_help()
                quit()

    if(SA_steps is None) :
        SA_steps = 40000

    return SA_steps



def check_float_tau_opt_arg(tau) :
    """ Function to check that the tau optional argument is an float """
    
    if(tau is not None):   
        tau = check_Int_Float_Str(tau)

        if(not isinstance(tau, float)):
            print("\n####################################################################################")
            print("ERROR. '-t/--decay-time {}' set, but not recognized. Only a positive float number number is allowed. 'strings', or negative numbers are not permitted.".format(tau))
            print("       Please, insert a positive float number (-t/--decay-time <FLOAT>) or ignore this flag.")
            print("       Look below for further help.") 
            print("####################################################################################\n\n")
            print_help()
            quit()
        else:
            if(tau <= 0):
                print("\n####################################################################################")
                print("ERROR. '-t/--decay-time {}' set, but not recognized. Only a float number higher than 0 is allowed".format(tau))  
                print("       A negative number or zero are meaningless.")
                print("       Please, insert a float positive number (-t/--decay-time <FLOAT>) or ignore this flag.") 
                print("       Look below for further help.")
                print("####################################################################################\n\n")
                print_help()
                quit()

    return tau



####################################        Help functions        ####################################

def print_help(): 
    """ help function """

    print("\nUsage: $~ python3 {} -c <prot_code> -d <excogito_file_directory> -g <.gro file> [-i <T_0>] [-n <num_propre>] [-p <parallel_bool>] [-r <radius_canvas>] [-s num_steps] [-t decay_time]".format(sys.argv[0]))
   
    print("\n-----------------------------------------------------------------------------------------------------\n")
    print("Mandatory arguments :\n")
    print("   prot_code                     MANDATORY          Code with which the output files are labelled")
    print("   excogito_file_directory       MANDATORY          Directory with excogito 'optimize' command output")
    print("   .gro file                     MANDATORY          File of atom Coordinates in .gro format\n")

    print("Optional arguments :\n")
    print("   T_0                           OPTIONAL           Initial temperature of SA calculation (default: chosen in order to accept 75% of initial moves)")
    print("   num_propre                    OPTIONAL           Number of atoms to keep in CG model, output of propre")
    print("   parallel_bool                 OPTIONAL           Decide if you want more SA minimization in parallel {1 : parallel | 0 : serial} (default: 1)")
    print("   radius_canvas                 OPTIONAL           Value for the radius of CANVAS model (choice 2) (default : 0.7 nm)")
    print("   num_steps                     OPTIONAL           Number of MC steps for the simulated annealing minimization (default: 40 000)")
    print("   decay_time                    OPTIONAL           Decay time for temperature in the simulated annealing minimization");
    print("                                                    (default: calculated such that the final temperature is 0.001 * initial temperature).\n")
    
    print("-----------------------------------------------------------------------------------------------------\n")

    print("Hereafter the list of flags:\n")

    print("   -c  --code                    string             Protein code to label output")
    print("   -d  --directory               string             Directory with output of excogito optimize")
    print("   -g  --gro                     FILE               Coordinate FILE (gro format)\n")

    print("  [-i]  [--initial-T]            float              Initial temperature for SA calculation.")
    print("  [-n]  [--num-propre]           int                Number of atoms obtained by propre.")
    print("  [-p]  [--parallel]             int                Serial(0) or parallel(1) execution. Just 0 and 1 are allowed.")
    print("  [-r]  [--radius]               float              Radius of CANVAS model.")
    print("  [-s]  [--steps]                int                Number of MC steps.")
    print("  [-t]  [--decay-time]           float              Decay time.\n");
    
    print("  [-h]  [--help]                                    Give this help list\n")



def print_help_postprocessing(): 
    """ help function """

    print("\nUsage: $~ python3 {} -d <SA_dir> -n <num_propre> [-r <radius_canvas>] [-f <file_resid.txt>]".format(sys.argv[0]))
   
    print("\n-----------------------------------------------------------------------------------------------------\n")
    print("Mandatory arguments :\n");  
    print("   SA_dir                        MANDATORY          Name of the directory in output_files which contains the data produced by the NAME_OF_SCRIPT")
    print("   num_propre                    MANDATORY          Number of atoms to keep in CG model, output of propre\n")
    
    print("Optional arguments :\n")
    
    print("   radius_canvas                 OPTIONAL           Value for the radius of CANVAS model (choice 2) (default : 0.7 nm)")
    print("   file_resid.txt                OPTIONAL           Path to the files which contains all the residues sorted by probability (produced by NAME_OF_SCRIPT)")
    print('                                                    Default : a file named "file_resid_sort.txt" will be searched in SA_dir directory.\n')
    
    print("-----------------------------------------------------------------------------------------------------\n");
    print("Hereafter the list of flags:\n");

    print("   -d  --directory               string             Directory with the optimized mappings for CANVAS by NAME_OF_SCRIPT")
    print("   -n  --num-propre              int                Number of atoms obtained by propre.\n")
    
    print("  [-r]  [--radius]               float              Radius of CANVAS model.")
    print("  [-f]  [--txt-file]             string             Path to the file with residues sorted.\n")
    
    print("  [-h]  [--help]                                    Give this help list\n")

