#!/usr/bin/python

import sys
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
from Bio import SeqIO
import wget
import os

p = PDBParser(PERMISSIVE=1)



if len(sys.argv) != 2:
    print("------------------------------------------------------")
    print("Usage : ./complexity.py <YOUR PDB> <TMH start> <TMH stop>")
    print("------------------------------------------------------")
    sys.exit()
else:
    pdb_id_filename = sys.argv[1]
    pdb_id = os.path.splitext(pdb_id_filename)[0]

    #Check the old file isn't there!!!!
    print("Checking local files for %s" % pdb_id_filename)
    try:
        with open(pdb_id_filename, 'r') as file:
            pass
    except():
        print("Local file not found.")
        print("Downloading %s from the PDB..." % pdb_id)
        # Download the file
        pdburl = ("https://files.rcsb.org/download/%s" % pdb_id_filename)
        file_download = wget.download(pdburl)

    #Get the sequence
    print('\nConverting %s to fasta file...' % pdb_id)
    pdb_fasta_filename=str(pdb_id+".fasta")
    count = SeqIO.convert(pdb_id_filename, "pdb-seqres", pdb_fasta_filename, "fasta")
    # arguments to get_structure(title,filename)
    with open(pdb_fasta_filename, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            print(record.seq)
