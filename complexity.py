#!/usr/bin/python

import sys
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
from Bio import SeqIO
import wget
import os
import subprocess
from subprocess import check_output

p = PDBParser(PERMISSIVE=1)


if len(sys.argv) != 2:
    print("------------------------------------------")
    print("Usage : python3 ./complexity.py <YOUR PDB>")
    print("------------------------------------------")
    sys.exit()
else:
    pdb_id_filename = sys.argv[1]
    pdb_id = os.path.splitext(pdb_id_filename)[0]

    # Check if a local file already exists.
    print("Checking local files for %s" % pdb_id_filename)
    try:
        with open(pdb_id_filename, 'r') as file:
            print("Local file, %s, found." % pdb_id_filename)
            pass

    # Exceptions MUST NOT be left naked.
    except(FileNotFoundError):
        print("Local file not found.")
        print("Downloading %s from the PDB..." % pdb_id)
        # Download the file
        pdburl = ("https://files.rcsb.org/download/%s" % pdb_id_filename)
        file_download = wget.download(pdburl)

    # Get the sequence
    print('\nConverting %s to fasta file...\n' % pdb_id)
    pdb_fasta_filename = str(pdb_id + ".fasta")
    count = SeqIO.convert(pdb_id_filename, "pdb-seqres",
                          pdb_fasta_filename, "fasta")
    # arguments to get_structure(title,filename)
    with open(pdb_fasta_filename, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            print("id:", record.id)
            print("sequence:", record.seq)

            ### Writing TMsegments.txt ###

            # Tries running the TMHMM program. If it is found the coordinates
            # will automatically be entered. If the program fails, the user
            # manually enters the TMH boundaries and these are passed to TMSOC.
            tmh_results = subprocess.getoutput(["tmhmm", pdb_fasta_filename])

            if "tmhmm: not found" in tmh_results:
                print(
                    "\nTMHMM not installed properly. You must enter the TMH boundaries manually.\n")
                user_tmh_boundaries = input(
                    "Please enter the start and stop locations of the TMH.\nSeparate start and stop by a comma, and a new TMH by a space.\nExample: 2,22 48,68 1023,1046\n")
                # tmh_boundaries_for_file_export=user_tmh_boundaries.rsplit(sep="
                # ", maxsplit=-1)
                with open('TMsegments.txt', 'w') as f:
                    write_data = f.write(user_tmh_boundaries)
            else:
                print(tmh_results)

            ### Writing FASTA for the 1 chain only ###

            with open('sequence.fasta', 'w') as f:
                write_fasta = f.write(
                    "> " + record.id + ":\n" + str(record.seq))

            # Let's set the TMH complexity to a null value incase of an
            # error.
            tmh_complexity = "NONE"

            # run tmsoc
            #'perl TMSOC.pl sequence.fasta TMsegments.txt '

            command = "perl TMSOC.pl sequence.fasta TMsegments.txt"
            # print "Running", command
            print("Running TMSOC..."
            results=str(check_output([command], shell=True))

            print("\nTMSOC results:\n" + results)

            print("\nTMH TMSOC outcome:\n")
            # This will only work for single pass. needs to take into account
            # multiple results else it will only report the last result.
            if ";simple" in results:
                tmh_complexity="Simple"
            elif ";complex" in results:
                tmh_complexity="Complex"
            elif ";twilight" in results:
                tmh_complexity="Twilight"
            else:
                pass
            print(tmh_complexity)
