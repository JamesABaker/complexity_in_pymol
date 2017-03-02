import sys
from pymol import cmd, stored


def complexity(selection='all'):
    # DESCRIPTION
    '''
    In pymol type `run Path/To/complexity_in_pymol/complexity.py` and hit enter. Then type `complexity` in the pymol terminal. This function will extract and recolour transmembrane helices according to the complexity score of the helix.
    '''

    write_input_fasta_file = open("inputseq.fasta", "w")
    write_segment_positions_file = open("TMsegments.txt", "w")

    cmd.hide("all")
    cmd.show('cartoon', "all")
    cmd.color('gray', "all")

    # writing the fasta sequence to the input
    sequence = cmd.get_fastastr('all')
    write_input_fasta_file.writelines([sequence])

    # Do sequence predition. For now, let's use the predefined ones for 1MT5,
    # which is probably a SP and a TMH.
    tmh_list = [[410, 438]]

    # It is much less prone to error if only 1 TMH is fed at a time.
    for a_tmh in tmh_list:
        start = int(a_tmh[0])
        stop = int(a_tmh[1])
        print("Checking tmh at positions %s-%s" % (start, stop))
        output_locations = str(str(start) + " , " + str(stop))
        write_segment_positions_file.writelines([output_locations])

        # run tmsoc
        tmh_complexity = "NONE"
        commandstring = 'bash perl.sh'
        os.system(commandstring)

        # reading the TMSOC result

        # This static link should be made dynamic so that it works after pymol
        # imports the python.
        tmsoc_result = open("TMSOCoutput.txt", 'r')

        for line in tmsoc_result:
            if ";simple" in line:
                complexity = "Simple"
            elif ";complex" in line:
                complexity = "Complex"
            elif ";twilight" in line:
                complexity = "Twilight"
            else:
                print(line)
                pass

        cmd.select("1MT5_B", "resi %s-%s" % (start, stop))
        cmd.show('cartoon', "1MT5_B")
        if tmh_complexity == "Complex":
            cmd.color('red', "resi %s-%s" % (start, stop))
        elif tmh_complexity == "Simple":
            cmd.color('blue', "resi %s-%s" % (start, stop))
        elif tmh_complexity == "Twilight":
            cmd.color('purple', "resi %s-%s" % (start, stop))
        elif tmh_complexity == "NONE":
            print("Complexity calculation failed.")
        print("%s-%s complexity=" % (start, stop), tmh_complexity)


cmd.extend("complexity", complexity)
