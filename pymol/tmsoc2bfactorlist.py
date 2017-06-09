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


from Bio.PDB import PDBParser
from Bio import SeqIO

import sys


def b_factor_complexity(complexity):
    if str('complex') in str(complexity):
        return 1
    elif str('twilight') in str(complexity):
        return 2
    elif str('simple') in str(complexity):
        return 3
    else:
        return 4
        print("ERROR in reading TMSOC output for TMH.")


new_tmsoc_result = "1. TM segment(s) summary:"
all_results = []
current_result = []

for line in open('tmsoc_output.txt'):
    if line.startswith(new_tmsoc_result) and current_result:
        # if line starts with thew new TMSOC result token and the current chunk
        # is not empty
        all_results.append(current_result[:])  # add not empty chunk to chunks
        current_result = []  # make current chunk blank
    # just append a line to the current chunk on each iteration
    current_result.append(line)

all_results.append(current_result)  # append the last chunk outside the loop

for current_result in all_results:
    b_factors = []
    for record in SeqIO.parse("sequence.fasta", "fasta"):
        print(record.seq)
        for n in range(len(record.seq)):
            b_factors.append(0)

    for item_position, line in enumerate(current_result):
        line = line.replace(',', ";")
        line = line.split(';')

        if str("1. TM segment(s) summary:") not in str(line) and str('2. Masked FASTA sequence:') not in str(line) and len(line) == 7:
            tmh_sequence = line[0]
            start_position = int(line[1])
            end_position = int(line[2])
            score1 = line[3]
            score2 = line[4]
            score3 = line[5]
            complexity = line[6]

            for record in SeqIO.parse("sequence.fasta", "fasta"):

                for residue_number, residue in enumerate(record.seq):
                    if residue_number <= end_position and residue_number >= start_position and end_position - start_position > 1:
                        b_factors[residue_number] = b_factor_complexity(
                            str(complexity))
                    else:
                        pass
        elif str(">") in str(line):
            subchain_id = line[0].replace('\n', '')
            subchain_id = subchain_id.replace('>', '_')
            subchain_id = subchain_id.replace(':', '_')

    with open("bfactors%s.txt" % subchain_id, 'w') as bfactors_output_file:
        for i in b_factors:
            bfactors_output_file.write("%s\n" % i)
