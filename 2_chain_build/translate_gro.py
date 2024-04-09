#!/bin/python3

with open("misfolded_chain_1.gro", "r") as initial_file :
    with open("misfolded_chain_2.gro", "w") as final_file :

        # First lines
        line_1 = initial_file.readline()
        final_file.write(line_1)

        line_2 = initial_file.readline()
        final_file.write(line_2)
        num_atoms = int(line_2.split("\n")[0][1:])

        check = "start"

        while check != "" and check is not None :

            # Split fields and remove extra char
            check = initial_file.readline()

            check = check.split(" ")
            while "" in check :
                check.remove("")

            # If the first field has a ".", then ignore it
            if "." in check[0] :
                quit()

            # Split the first two fields (residue type and residue number)
            res_type = check[0][-3] + check[0][-2] + check[0][-1]
            res_num = ""
            for i in range( len(check[0]) - 3 ) : res_num += check[0][i]

            # Translation along x coordinate
            dx = 0                                      # Translation along x in nm
            dy = 0                                      # Translation along y in nm
            dz = 4                                      # Translation along z in nm
            check[3] = float(check[3]) + dx             # x positions
            check[4] = float(check[4]) + dy				# y positions
            check[5] = float(check[5]) + dz				# z positions

            

            # Writing output
            final_file.write(f"{res_num:>5}{res_type:<5}{check[1]:>5}{check[2]:>5}{check[3]:8.3f}{check[4]:8.3f}{check[5]:8.3f}\n")
