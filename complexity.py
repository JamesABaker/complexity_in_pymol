import sys
import subprocess
from subprocess import call
from pymol import cmd, stored

def complexity(selection='all'):
    '''
    In pymol type `run Path/To/complexity_in_pymol/complexity.py` and hit enter. Then type `complexity` in the pymol terminal. This function will extract and recolour transmembrane helices according to the complexity score of the helix.
    '''

    write_input_fasta_file = open("sequence.fasta", "w")
    write_segment_positions_file = open("TMsegments.txt", "w")

    cmd.hide("all")
    cmd.show('cartoon', "all")
    cmd.color('gray', "all")

    # writing the fasta sequence to the input
    sequence = cmd.get_fastastr('all')
    print(sequence)
    write_input_fasta_file.writelines([sequence])

    # Do sequence predition. For now, let's use the predefined ones for 1MT5,
    tmh_list = [[5,25]]

    # It is much less prone to error if only 1 TMH is fed at a time.
    for a_tmh in tmh_list:
        start = int(a_tmh[0])
        stop = int(a_tmh[1])
        print(("Checking tmh at positions %s-%s" % (start, stop)))
        output_locations = str(str(start) + "," + str(stop))
        write_segment_positions_file.writelines([output_locations])

        #Let's set the TMH complexity to a null value incase of an error.
        tmh_complexity = "NONE"

        # run tmsoc
        commandstring = 'perl TMSOC.pl ./sequence.fasta ./TMsegments.txt >> TMSOCoutput.txt'
        subprocess.call([commandstring], shell=True)


        # reading the TMSOC result

        # This static link should be made dynamic so that it works after pymol
        # imports the python.
        tmsoc_result = open("TMSOCoutput.txt", 'r')

        for line in tmsoc_result:
            print(line)
            if ";simple" in line:
                tmh_complexity = "Simple"
            elif ";complex" in line:
                tmh_complexity = "Complex"
            elif ";twilight" in line:
                tmh_complexity = "Twilight"
            else:
                pass

        #Remove the output to stop it intefering with other iterations.
        subprocess.call("rm sequence.fasta", shell=True)
        subprocess.call("rm TMsegments.txt", shell=True)
        subprocess.call("rm TMSOCoutput.txt", shell=True)


        cmd.select("new_selection", "resi %s-%s" % (start, stop))
        if tmh_complexity == "Complex":
            cmd.color('red', "resi %s-%s" % (start, stop))
        elif tmh_complexity == "Simple":
            cmd.color('blue', "resi %s-%s" % (start, stop))
        elif tmh_complexity == "Twilight":
            cmd.color('purple', "resi %s-%s" % (start, stop))
        elif tmh_complexity == "NONE":
            cmd.color('orange', "resi %s-%s" % (start, stop))
            print("Complexity calculation failed.")

        print(("%s-%s complexity=" % (start, stop), tmh_complexity))
    print("Red is complex, blue is simple, and purple is twighlight. Orange helices indicate an error.")

cmd.extend("complexity", complexity)
