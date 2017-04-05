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

with open("tmh_list.txt") as f:
    content = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end
# of each line
content = [x.strip() for x in content]

for pdb_id in content:
    pdb_id_filename = pdb_id.split('_')[0] + ".pdb"

    # Check if a local file already exists.
    print("Checking local files for %s" % pdb_id_filename)
    try:
        with open("original_pdbs/%s" % pdb_id_filename, 'r') as file:
            print("Local file, %s, found." % pdb_id_filename)
            #os.rename("original_pdb/%s" % pdb_id_filename, pdb_id_filename)
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
    pdb_fasta_filename = str(pdb_id.split('_')[0] + ".fasta")
    count = SeqIO.convert(pdb_id_filename, "pdb-seqres",
                          pdb_fasta_filename, "fasta")

    # arguments to get_structure(title,filename)
    with open(pdb_fasta_filename, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            print("id:", record.id)
            print("sequence:", record.seq)

            ### Writing FASTA for the 1 chain only ###

            with open('sequence.fasta', 'w') as f:
                write_fasta = f.write(
                    ">" + record.id + ":\n" + str(record.seq))

            ### Writing TMsegments.txt ###

            ### Attempts to find the ID in the PDBTM XML file


                # Tries running the phobius program. If it is found the coordinates
                # will automatically be entered. If the program fails, the user
                # manually enters the TMH boundaries and these are passed to TMSOC.
                tmh_results = subprocess.getoutput(
                    ["perl phobius.pl sequence.fasta"])

                if "No such file or directory" in tmh_results:
                    print(tmh_results)
                    print("\nPhobius not installed properly.\n")
                    user_tmh_boundaries = input(
                        "Please enter the start and stop locations of the TMH.\nSeparate start and stop by a comma, and a new TMH by a space.\nExample: 2,22 48,68 1023,1046\n")
                    with open('TMsegments.txt', 'w') as f:
                        write_data = f.write(user_tmh_boundaries)
                else:
                    # print(tmh_results)

                    protein_coordinates = []
                    for line in tmh_results.splitlines():
                        if "TRANSMEM" in line:
                            line_coordinates = []
                            for item in line.split(" "):
                                item.replace(" ", "")
                                # Phobius output contains something like: TRANSMEM detected: FT   TRANSMEM     30     56
                                # The two numbers can be extracted since they are
                                # the only integers.
                                try:
                                    line_coordinates.append(int(item))
                                except(ValueError):
                                    pass
                            # Now we check that there is only a start and a stop.
                            if len(line_coordinates) == 2:
                                protein_coordinates.append(line_coordinates)
                            else:
                                print(
                                    "ERROR in", pdb_id, ". The transmem line does not contain only start and stop coordinates.")
                                print(line)
                    with open('TMsegments.txt', 'w') as f:
                        for coordinates in protein_coordinates:
                            for number, vector in enumerate(coordinates):
                                write_data = f.write(str(vector))
                                if number == 0:
                                    write_data = f.write(",")
                            write_data = f.write(" ")



            # Let's set the TMH complexity to a null value incase of an
            # error.
            tmh_complexity = "NONE"

            # run tmsoc
            #'perl TMSOC.pl sequence.fasta TMsegments.txt '

            command = "perl TMSOC.pl sequence.fasta TMsegments.txt"
            # print "Running", command
            print("Running TMSOC...")
            results = str(check_output([command], shell=True))

            print("\nTMSOC results:\n" + results + "\n")

            #Move files to download folders.

            #os.rename(pdb_fasta_filename, "original_fasta/%s" % pdb_fasta_filename)
            #os.rename(pdb_id_filename, "original_pdb/%s" % pdb_id_filename)
