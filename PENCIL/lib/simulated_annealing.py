# Importing modules
import numpy as np
import time
import sys
import pandas as pd
from lib.manipulate_amminoacid import eval_num_at_model
from numba import jit



@jit(nopython = True)
def replace_random_atom(initial_array, at_residues):
    """ Function that takes an AA residue and it exchanges it with a CG or MG residue (MC move proposal).
        The initial status of the system is described from arrays initial_array and at_residues (2D-arrays).
        The outputs arr_new and at_residues_new contain the updated list of residues where one has been removed from AA list and one has been added """

    # Generating random number for the residue to remove
    rand_remove = np.random.randint(at_residues.shape[0])
    num_resid_to_remove = at_residues[rand_remove][0]
    
    # Generating random number for the residue to add until I get a residue not AA
    flag = 1
    while flag :
        rand_add = np.random.randint(initial_array.shape[0])
        if np.rint(initial_array[rand_add][-1]) != 1 : flag = 0

    # Making the old AA residue CG in the global list
    for i in range(initial_array.shape[0]) :
        if initial_array[i][0] == num_resid_to_remove : initial_array[i][-1] = -1

    # Adding an AA residue
    initial_array[rand_add][-1] = 1
    at_residues[rand_remove] = initial_array[rand_add].copy()

    return initial_array, at_residues



# def eval_start_T_numba(initial_array, at_residues, E, num_propre, radius_canvas, num_trials, k_1, k_2) :
#     """ Function to evaluate the starting temperature: I propose num_trials moves from the initial configuration and I fix the temperature so that the acceptance
#         probability is 0.75 """
    
#     T = np.empty(num_trials)                            # Array with temperature values
#     n = 0            

#     for n in range(num_trials) :

#         # Copy the initial system
#         arr_new = initial_array.copy()                  # Atomistic residues in the proposed move
#         at_residues_new = at_residues.copy()            # Array of all residues in the proposed move

#         # Generating move
#         arr_new, at_residues_new = replace_random_atom(arr_new, at_residues_new)               

#         # Reset all MG and CG residues to CG and calculating MG and CG atoms with the new configuration
#         num_at_model_new, arr_new = eval_num_at_model(arr_new, at_residues_new, radius_canvas)

#         # Evaluating the new cost function
#         E_new = cost_function_SA(arr_new, at_residues_new, num_at_model_new, num_propre, k_1, k_2)

#         # Accepting the move or rejecting 
#         if E_new < E :
#             T[n] = 0
#         else :
#             T[n] = (E - E_new) / np.log(0.75)

#     # Now I sort my T vector and I take the element at 0.75*num_trials(i.e. a temperature that make me prbably accept 75% of the moves)
#     T_sorted = np.sort(T)
#     T_0 = T_sorted[ int(0.75*num_trials) ]

#     return T_0


def eval_start_T_numba(initial_array, at_residues, E, num_propre, radius_canvas, num_trials, k_1, k_2) :
    """ Function to evaluate the starting temperature: I propose num_trials moves from the initial configuration and I fix the temperature so that the acceptance
        probability is 0.75 """
    
    delta = np.empty(num_trials)                            # Array with deltas
    initial_acceptance_ratio = 0.75
    n = 0            

    for n in range(num_trials) :

        # Copy the initial system
        arr_new = initial_array.copy()                  # Atomistic residues in the proposed move
        at_residues_new = at_residues.copy()            # Array of all residues in the proposed move

        # Generating move
        arr_new, at_residues_new = replace_random_atom(arr_new, at_residues_new)               

        # Reset all MG and CG residues to CG and calculating MG and CG atoms with the new configuration
        num_at_model_new, arr_new = eval_num_at_model(arr_new, at_residues_new, radius_canvas)

        # Evaluating the new cost function
        E_new = cost_function_SA(arr_new, at_residues_new, num_at_model_new, num_propre, k_1, k_2)

        # Calculating cost function difference
        delta[n] = np.abs(E - E_new)

    # Temperature calculated imposing average acceptance 0.75
    T_0 = - np.mean(delta) / np.log(initial_acceptance_ratio)
    
    return T_0



@jit(nopython=True)
def cost_function_SA(array_residues, at_residues, num_at_model, num_propre, k_1, k_2) :
    """ SA cost function:
        -> Probability term
        -> First neighbours
        -> Spring for atom number """
    

    E = 0

    for i in range(at_residues.shape[0]) :
        # Probability term
        E -= k_1 * at_residues[i][6]
        for j in range(array_residues.shape[0]) :
            # First neighbours
            if np.rint(array_residues[j][0]) == np.rint(at_residues[i][0])+1 : 
                E -= np.rint(array_residues[j][-1])
                break

    if num_at_model > num_propre : E += k_2 * np.power((num_at_model - num_propre) / num_propre, 2)

    return E



def cost_function_MC(list_residues, at_residues) :
    """ MC cost function:
        -> First neighbours
        -> Second neighbours """

    E = 0

    for amminoacid in at_residues :
        for residue in list_residues :
            # First neighbours
            if residue[0] == amminoacid[0]+1 : 
                E -= residue[-1]
                break
           
    return E



def progress_bar(prog, actual_prog) :
    """ Function to print progress bar on terminal """
    
    idx = 1
    max_ast = np.rint(100*idx).astype(int)             # The maximum number of asterisks
    num_ast = np.floor (prog*idx).astype(int)          # Actual number of asterisks
    if prog != actual_prog :
        sys.stdout.write ("\r" + "Calculating" + "\t" + "[ " + "#" * num_ast + "." * (max_ast -num_ast) + " ]" + "  " + prog.astype(str) + "%" )
        if prog % 10 != 0 : actual_prog = prog
    
    return actual_prog




def SA_optimization(output_dir, SA_steps, T_0, tau, arr_resid_global, at_res_global, radius_canvas, num_propre, list_residues_sorted, E, num_at_model, k_1, k_2, n) :
    """ This function performs the simulated annealing optimization procedure """

    # Measuring time taken to execute the function
    start_time = time.perf_counter()

    # Fixing spacing on terminal before progress bar
    if n == 0 : print()
    
    with open(f"{output_dir}/SA_data_{n}.txt", "w") as data_file :

        data_file.write(f"# File with data from the Simulated Annealing optimization number {n}\n\n")
        data_file.write("#" + "~" * 60 + "\n\n")
        data_file.write("# Starting simulated annealing procedure with the following parameters:\n")
        data_file.write(f"#\t-> MC_steps = {SA_steps}\n")
        data_file.write(f"#\t-> T_0 = {T_0:.3f}\n")
        data_file.write(f"#\t-> tau = {tau:.3f}\n")
        data_file.write(f"#\t-> k_1 = {k_1:.3f}\n")
        data_file.write(f"#\t-> k_2 = {k_2:.3f}\n\n")
        data_file.write("#" + "~" * 60 + "\n\n")
        data_file.write("# The data written here are the following:\n")
        data_file.write("# \tcolumn 0 : Number of the move\n")
        data_file.write("# \tcolumn 1 : Acceptance bit (1 accepted, 0 rejected)\n")
        data_file.write("# \tcolumn 2 : Temperature value\n")
        data_file.write("# \tcolumn 3 : Value of the cost function after the move\n")
        data_file.write("# \tcolumn 4 : Number of atomistic residues\n\n")

        with open(f"{output_dir}/SA_{n}.log", "w") as log_file:

            log_file.write(f"# Log file for the Simulated Annealing optimization number {n}\n\n")
            log_file.write("#" + "~" * 60 + "\n\n")
            log_file.write("# Starting simulated annealing procedure with the following parameters:\n")
            log_file.write(f"#\t-> MC_steps = {SA_steps}\n")
            log_file.write(f"#\t-> T_0 = {T_0:.3f}\n")
            log_file.write(f"#\t-> tau = {tau:.3f}\n")
            log_file.write(f"#\t-> k_1 = {k_1:.3f}\n")
            log_file.write(f"#\t-> k_2 = {k_2:.3f}\n\n")
            log_file.write("#" + "~" * 60 + "\n\n")

            actual_prog = -1        # For progress status

            # Initializing the arrays for proposed moves equal to the initial ones
            array_residues_sorted = arr_resid_global.copy()     # Local variable for the list of sorted residues
            at_residues = at_res_global.copy()                  # Local variable for the list of atomistic residues
            arr_new = array_residues_sorted.copy()              # Atomistic residues in the proposed move
            at_residues_new = at_residues.copy()                # Array of all residues in the proposed move

            for mc_step in range(SA_steps) :

                # Propose a move
                arr_new, at_residues_new = replace_random_atom(arr_new, at_residues_new)

                # Reset all MG and CG residues to CG and calculating MG and CG atoms with the new configuration
                num_at_model_new, arr_new = eval_num_at_model(arr_new, at_residues_new, radius_canvas)

                # Evaluating the new cost function
                E_new = cost_function_SA(arr_new, at_residues_new, num_at_model_new, num_propre, k_1, k_2)

                # Accept or reject the move (and upgrade all data if accepted)
                acc_flag = 0
                if E_new < E :
                    E = E_new
                    at_residues = at_residues_new.copy()
                    array_residues_sorted = arr_new.copy()
                    num_at_model = num_at_model_new
                    acc_flag = 1
                else :
                    dE = E_new - E
                    T = T_0 * np.exp(-mc_step/tau)
                    acc_prob = np.exp(-dE / T)
                    if np.random.rand() < acc_prob :
                        E = E_new
                        at_residues = at_residues_new.copy()
                        array_residues_sorted = arr_new.copy()
                        num_at_model = num_at_model_new
                        acc_flag = 1
                    else :
                        # If move not accepted, I copy initial values
                        at_residues_new = at_residues.copy()
                        arr_new = array_residues_sorted.copy()

                # Progress status (on .log file) and progress bar (on terminal, printd only by process with n = 0)
                prog = np.floor(100*mc_step/(SA_steps-1)).astype(int)

                if n == 0 : actual_prog = progress_bar(prog, actual_prog)

                if prog != actual_prog and prog % 10 == 0 :
                    intermediate_time = time.perf_counter()
                    elapsed_time = round((intermediate_time - start_time), 3)
                    log_file.write(f"\nElapsed time: {elapsed_time} s\n")
                    log_file.write(f"Calculation progress : {prog} % ...\n\n")
                    actual_prog = prog

                
                # Check if num_at_model < num_propre. If yes, add a random AA residue
                if num_at_model < num_propre :

                    # Pick one MG or CG residue at random
                    log_file.write(f"\n\nNumber of atoms in the model : {num_at_model}\n")
                    flag = 1
                    while flag :
                        rand_add = np.random.randint(array_residues_sorted.shape[0])
                        if np.rint(array_residues_sorted[rand_add][-1]) != 1 : flag = 0

                    # Adding an AA residue
                    log_file.write(f"Added residue : {list_residues_sorted[rand_add]}\n")
                    array_residues_sorted[rand_add][-1] = 1
                    #num_at_residues += 1
                    at_residues = np.append(at_residues, [array_residues_sorted[rand_add]], axis=0)

                    # Recalculating num_at_model and correct states after adding the atom
                    num_at_model, array_residues_sorted = eval_num_at_model(array_residues_sorted, at_residues, radius_canvas)

                    # Updating cost function
                    E = cost_function_SA(array_residues_sorted, at_residues, num_at_model, num_propre, k_1, k_2)

                    # Copyinig the new initial configuration
                    arr_new = array_residues_sorted.copy()
                    at_residues_new = at_residues.copy()

                    # Status message
                    log_file.write(f"Number of atomistic residues: {at_residues.shape[0]}\n")
                    log_file.write("Number of atoms in the model : {num_at_model}\n\n")

                # Writing data to file
                data_file.write(f"{mc_step}\t{acc_flag}\t{T_0 * np.exp(-mc_step/tau):.3f}\t{E:.3f}\t{len(at_residues)}\n")

            # Fixing spacing on terminal after progress bar
            if n == 0 : print()

            # Now we have to print the output in a readable and appealing form: to do so, we retrieve convert out at_residues array in a list with the same format of the
            # lists used in the first part of the code
            at_residues_list_output = [ residue for residue in list_residues_sorted if residue[0] in np.rint(at_residues[:,0]) ]


            # Sort final atomistic residues and print them
            at_residues_list_output = sorted(at_residues_list_output, key=lambda resid: resid[0])
            log_file.write("\nAtomistic residues after simulated annealing:\n\n")
            log_file.write(f"{pd.DataFrame(at_residues_list_output)}")
            log_file.write(f"\n\nNumber of sites in the model : {num_at_model} ( propre atom number : {num_propre} )")

    return 0

