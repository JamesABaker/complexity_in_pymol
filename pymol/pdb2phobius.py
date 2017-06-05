import sys
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
from Bio import SeqIO
import os
import subprocess
from subprocess import check_output

pdb_fasta_filename = sys.argv[1]

print("Opening %s..." % pdb_fasta_filename)

# with open(pdb_fasta_filename, "rU") as input_handle:
# with open("sequence.fasta", "w") as output_handle:
#    sequences = SeqIO.parse(input_handle, "pdb-seqres")
#    count = SeqIO.write(sequences, output_handle, "fasta")
count = SeqIO.convert(pdb_fasta_filename, "pdb-seqres",
                      "sequence.fasta", "fasta")
print("Converted %i records" % count)

### Writing TMsegments.txt ###

# Tries running the phobius program. If it is found the coordinates
# will automatically be entered. If the program fails, the user
# manually enters the TMH boundaries and these are passed to
# TMSOC.
print("Running Phobius on sequence to estimate TMH regions...")
with open('phobius_output.txt', 'wb', 0) as output_file:
    subprocess.Popen(
        ["perl", "phobius.pl", "sequence.fasta"], stdout=output_file)

with open("phobius_output.txt", "rt") as in_file:

    text = in_file.read()
    print("Phobius complete.")
# print(text)
