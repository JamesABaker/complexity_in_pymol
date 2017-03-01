import sys
from Bio import SeqIO
from Bio.PDB import *
from pymol import cmd, stored

write_input_fasta_file = open("inputseq.fasta", "w")
write_segment_positions_file = open("TMsegments.txt", "w")
cmd.color('white', "all")

def complexity(selection='all'):

    # writing the fasta sequence to the input
    sequence = cmd.get_fastastr('all')
    #I need to get this to read actual ID
    pdb_id = "1MT5_B"
    write_input_fasta_file.writelines([str(">" + pdb_id), sequence])

    # Do sequence predition. For now, let's use the predefined ones for 1MT5,
    # which is probably a SP and a TMH.
    tmh_list = [[410 , 438]]

    #It is much less prone to error if only 1 TMH is fed at a time.
    for a_tmh in tmh_list:
        start=int(a_tmh[0])
        stop=int(a_tmh[1])
        output_locations = str(str(start) + " , " + str(stop))
        print(output_locations)
        write_segment_positions_file.writelines([output_locations])

        # run tmsoc
        complexity = "NONE"
        commandString = 'bash perl.sh'
        os.system(commandString)

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
                #print(line)
                pass


        cmd.select("1MT5_B", "resi a_tmh[0]-a_tmh[1]")
        cmd.show('sphere', "1MT5_B")
        if complexity=="Complex":
            cmd.color('red', "resi a_tmh[0]-a_tmh[1]")
        elif complexity=="Simple":
            cmd.color('blue', "resi a_tmh[0]-a_tmh[1]")
        elif complexity=="Twilight":
            cmd.color('purple', "resi a_tmh[0]-a_tmh[1]")
        elif complexity=="NONE":
            print("Complexity calculation failed.")


    # DESCRIPTION
    # In pymol type `run Path/To/complexity_in_pymol/complexity.py` and hit
    # enter. Then type `complexity` in the pymol terminal. This function will
    # extract and recolour transmembrane helices according to the complexity
    # score of the helix.

    #
    # Your code goes here
    #
    return (sequence)

cmd.extend("complexity", complexity)
