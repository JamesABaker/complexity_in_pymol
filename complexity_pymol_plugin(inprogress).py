from __future__ import absolute_import
import sys
import subprocess
from subprocess import call
from subprocess import check_output
from pymol import cmd, stored
from io import open


def complexity(selection=u'all'):
    '''
    In pymol type `run Path/To/complexity_in_pymol/complexity.py` and hit enter. Then type `complexity` in the pymol terminal.
    This function will extract and recolour transmembrane helices according to the complexity score of the helix.
    '''
    print sys.version
    # might be useful for dynamic file locaitons
    #pwd = check_output("pwd", shell=True)
    #pwd = str(pwd.rstrip()).encode('utf-8')
    #print("File locations:", pwd)

    write_input_fasta_file = open("sequence.fasta", "w", encoding="utf-8")
    write_segment_positions_file = open(
        "TMsegments.txt", "w", encoding="utf-8")

    cmd.hide("all")
    cmd.show('cartoon', "all")
    cmd.color('gray', "all")

    # writing the fasta sequence to the input
    sequence = cmd.get_fastastr('all')

    print "Writing sequence.fasta\n"
    write_input_fasta_file.write(sequence)

    # Do sequence predition. For now, let's use the predefined ones for 1MT5,
    tmh_list = [[5, 25]]

    # It is much less prone to error if only 1 TMH is fed at a time.
    for a_tmh in tmh_list:
        start = int(a_tmh[0])
        stop = int(a_tmh[1])
        print ("Checking tmh at positions %s-%s" % (start, stop))
        output_locations = str(start) + "," + str(stop)
        print "Writing TMsegments.txt\n"
        write_segment_positions_file.write(unicode([output_locations]))

        # Checks input files are in place.
        file_list = check_output("ls", shell=True)
        fasta_file = False
        segment_file = False
        for line in file_list.split('\n'):
            if "sequence.fasta" in line:
                fasta_file = True
            if "TMsegments.txt" in line:
                segment_file = True

        if fasta_file == False or segment_file == False:
            print 'Input files for TMSOC are not present, but, they should be.'

        elif fasta_file == True and segment_file == True:
            print 'Input files for TMSOC found.'

            print u"Selecting TM region.\n"
            cmd.select("new_selection", "resi %s-%s" % (start, stop))

            # Let's set the TMH complexity to a null value incase of an error.
            tmh_complexity = "NONE"

            # run tmsoc
            #'perl TMSOC.pl sequence.fasta TMsegments.txt '

            command = "perl TMSOC.pl sequence.fasta TMsegments.txt"
            print "Running", command
            results = check_output(
                [command], shell=True)

            print "Result:"
            for line in results.split('\n'):
                print line
                if ";simple" in line:
                    tmh_complexity = u"Simple"
                elif ";complex" in line:
                    tmh_complexity = u"Complex"
                elif ";twilight" in line:
                    tmh_complexity = u"Twilight"
                else:
                    pass

            if tmh_complexity == "Complex":
                cmd.color('red', "resi %s-%s" % (start, stop))
            elif tmh_complexity == "Simple":
                cmd.color('blue', "resi %s-%s" % (start, stop))
            elif tmh_complexity == "Twilight":
                cmd.color('purple', "resi %s-%s" % (start, stop))
            elif tmh_complexity == "NONE":
                cmd.color('orange', "resi %s-%s" % (start, stop))
                print "Complexity calculation failed.\n"
            else:
                print "Catastrophic failure. The complexity variable behaved in a completely unexpected way."

            print "%s-%s complexity=" % (start, stop), tmh_complexity, "\n"

            # Remove the output to stop it intefering with other iterations.
            #print("Cleaning input and output files.")
            #subprocess.call("rm sequence.fasta", shell=True)
            #subprocess.call("rm TMsegments.txt", shell=True)
            #subprocess.call("rm TMSOCoutput.txt", shell=True)

        print "Red is complex, blue is simple, and purple is twighlight. Orange indicates an error reported by TMSOC."

cmd.extend("complexity", complexity)
