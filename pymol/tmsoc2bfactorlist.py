from Bio.PDB import *
from Bio import SeqIO

import sys

pdb_id = sys.argv[1]


def b_factor_complexity(complexity):
    if str('complex') in str(complexity):
        return 1
    if str('twilight') in str(complexity):
        return 2
    if  str('simple') in str(complexity):
        return 3
    else:
        return 4
        print("ERROR in reading TMSOC output for TMH.")


with open("tmsoc_output.txt", 'r') as tmsoc_output:
    content = tmsoc_output.readlines()

    for line_number, line in enumerate(content):
        if "2. Masked FASTA sequence:" in line:
            max_line = line_number

    p = PDBParser()
    structure = p.get_structure('X', pdb_id)

    b_factors = []


    #Change this to read from fasta
    for record in SeqIO.parse("sequence.fasta", "fasta"):
        for n in range(len(record.seq)):
            b_factors.append(0)

    for new_line_number, line in enumerate(content):

        if new_line_number > 0 and new_line_number < max_line:
            line = line.replace(',', ";")
            line = line.split(';')
            tmh_sequence = line[0]
            start_position = int(line[1])
            end_position = int(line[2])
            score1 = line[3]
            score2 = line[4]
            score3 = line[5]
            complexity = line[6]

            p = PDBParser()
            structure = p.get_structure('X', pdb_id)
            for model in structure:
                for chain in model:
                    for residue_number, residue in enumerate(chain):
                        if residue_number <= end_position and residue_number >= start_position:
                            b_factors[residue_number] = b_factor_complexity(str(complexity))
                        else:
                            pass

with open("bfactors.txt", 'w') as bfactors_output_file:
    for i in b_factors:
        bfactors_output_file.write("%s\n" % i)
