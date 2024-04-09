with open("radius_charge-allatom.dat", "r") as canvas_file :
    
    charges = []
    radii = []

    canvas_file.readline()

    for line in canvas_file :
        split_line = line.split()
        while "" in split_line :
            split_line.remove("")

        radii.append(float(split_line[2]))
        charges.append(float(split_line[3]))



with open("initial-frame.pdb", "r") as pdbfile :
    with open("initial-frame.pqr", "w") as pqrfile :

        for _ in range(4) :
            pdbfile.readline()

        i = 0
        for line in pdbfile :

            line = line.split()
            while "" in line :
                line.remove()

            if len(line) > 7 :
                line[8] = charges[i]
                line[9] = radii[i]

                pqrfile.write(f"{line[0]:<6}{line[1]:>5}{line[2]:>5} {line[3]:<4}{line[4]:>6}{line[5]:>12}{line[6]:>8}{line[7]:>8}{line[8]:>8.4f}{line[9]:>7.4f}\n")

                i += 1
    