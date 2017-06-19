#    This file is part of 3D TMH Complexity.
#
#    3D TMH Complexity is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    3D TMH Complexity is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with 3D TMH Complexity.  If not, see <http://www.gnu.org/licenses/>.

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


one_letter = {'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q',
              'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y',
              'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A',
              'GLY': 'G', 'PRO': 'P', 'CYS': 'C'}

parser = PDBParser()
structure = parser.get_structure('pdb_structure', pdb_fasta_filename)
for model in structure:
    for chain in model:
        print(chain.id, "\n")
        sequence = []
        residue_numbers = []
        new_chain = True
        discontinuation_events = []
        for residue in chain:
            if residue.get_resname() in one_letter:
                print(one_letter[residue.get_resname()], residue.id[1])
                sequence.append(one_letter[residue.get_resname()])
                residue_numbers.append(residue.id[1])
                if new_chain == True:
                    continuous_sequence_number = residue.id[1]
                new_chain = False
                if residue.id[1] == continuous_sequence_number:
                    continuous_sequence_number = continuous_sequence_number + 1
                else:
                    discontinuation_event = int(residue.id[1])
                    discontinuation_events.append(discontinuation_event)
                    continuous_sequence_number = residue.id[1] + 1
                    print("Discontinuation detected. Skips and missing regions are usually because because parts of the structure could not be confidently resolved.")
                if residue.id[1] < 0:
                    print(
                        "Position below 0 detected. This could be the result of imperfect tag cleavage.")
            else:
                print("Skipping", residue.get_resname(),
                      "at position", residue.id[1])
        with open("sequence.fasta", "a") as fasta_file:
            header = ">" + ";" + pdb_fasta_filename.split('.')[0] + ";" + chain.id + ";" + str(min(
                residue_numbers)) + ";" + str(max(residue_numbers)) + ";" + str(discontinuation_events) + ";"
            fasta_file.write(header)
            fasta_file.write("\n")
            fasta_file.write("".join(sequence))
            fasta_file.write("\n")


# count = SeqIO.convert(pdb_fasta_filename, "pdb-seqres",
#                      "sequence.fasta", "fasta")
#print("Converted %i records" % count)

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
