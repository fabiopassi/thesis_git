# Functions

def write_vmd_coloring_values(list_residues, output_dir) :
    """ Function to write on .txt file the average retain probability and the maximum retain probability for vmd representation """

    with open(f"""{output_dir}/vmd_probability_values.txt""", "w") as vmd_plot_values :
        vmd_plot_values.write("; Residue_number   Residue_type    Max_prob    Avg_prob\n")
        for residue in list_residues :
            vmd_plot_values.write(f"""{residue[0]}\t{residue[1]}\t{residue[-1]}\t{residue[-3]}\n""")



def write_vmd_coloring_script(output_dir) :
    """ Function to write on .tcl file the instructions for the VMD tk console to create a representation for the probability of residues """

    with open(f"""{output_dir}/vmd_probability_col.tcl""", "w") as vmd_tcl_script :
       vmd_tcl_script.write(f"""set fp [open "{output_dir}/vmd_probability_values.txt" r]\n\n""")
       vmd_tcl_script.write('# Count the number of lines in a text file\n\n')
       vmd_tcl_script.write('set number 0\n')
       vmd_tcl_script.write('while { [gets $fp line] >= 0 } {\n')
       vmd_tcl_script.write('\tincr number\n')
       vmd_tcl_script.write('}\n\n')
       vmd_tcl_script.write('close $fp\n\n')
       vmd_tcl_script.write('############################################\n\n')
       vmd_tcl_script.write(f"""set fp [open "{output_dir}/vmd_probability_values.txt" r]\n\n""")
       vmd_tcl_script.write('set file_data [read $fp]\n\n')
       vmd_tcl_script.write('close $fp\n\n')
       vmd_tcl_script.write('set data [split $file_data "\\n"]\n\n')
       vmd_tcl_script.write('for {set i 1} {$i < $number} {incr i} {\n\n')
       vmd_tcl_script.write('\tset a [lindex $data $i 0];\n')
       vmd_tcl_script.write('\tset b [lindex $data $i 1];\n')
       vmd_tcl_script.write('\tset c [lindex $data $i 2];\n')
       vmd_tcl_script.write('\tset d [lindex $data $i 3];\n')
       vmd_tcl_script.write('\tset e [atomselect top "residue [expr $i - 1]"];\n\n')
       vmd_tcl_script.write('\t$e set beta $c\n')
       vmd_tcl_script.write('\t$e set radius $d\n')
       vmd_tcl_script.write('}\n')



def write_residue_sorted_file(output_dir, filename, array_residues_sorted) :
    """ Function that creates the file with all the residues sorted by probability, required to run the C program.
        The first line contains the total number of residues, to make the memory allocation in C easier.
        NOTE: I print variables respecting the types I want in C."""
    
    with open(f"{output_dir}/{filename}", "w") as res_file :
        res_file.write(f"{array_residues_sorted.shape[0]}\n")
        for i in range(array_residues_sorted.shape[0]) :
            res_file.write(f"{array_residues_sorted[i][0].astype(int)}\t")
            res_file.write(f"{array_residues_sorted[i][1]}\t")
            res_file.write(f"{array_residues_sorted[i][2]}\t")
            res_file.write(f"{array_residues_sorted[i][3]}\t")
            res_file.write(f"{array_residues_sorted[i][4].astype(int)}\t")
            res_file.write(f"{array_residues_sorted[i][5]}\t")
            res_file.write(f"{array_residues_sorted[i][6]}\t")
            res_file.write(f"{array_residues_sorted[i][7].astype(int)}\n")
       

def write_canvas_list_at_res(at_res, file) :
    """ Function to write the file with the list of AA residues for the creation of CANVAS model """

    with open(file, "w") as list_at_res_file :

        for i in range(at_res.shape[0]) : list_at_res_file.write(f"{int(at_res[i][0])}\n")

    print(f"\nThe file {file} containing the list of all atom residues has been written correctly.")
    print("This file must be given to CANVAS block.py script via the flag -l.")
