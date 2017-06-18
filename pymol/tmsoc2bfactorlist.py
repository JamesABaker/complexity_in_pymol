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
    print("Processing TMSOC result chunk:", current_result)
    for tmsoc_result_line in current_result:
        if ">" in tmsoc_result_line:
            tmsoc_result_line=tmsoc_result_line.split(";")
            #Here we ignore the first entry, which is simply ">"
            tmsoc_result_line_pdb_code=str(tmsoc_result_line[1])
            tmsoc_result_line_chain_id=str(tmsoc_result_line[2])
            tmsoc_result_line_start_location=int(tmsoc_result_line[3])
            tmsoc_result_line_end_location=int(tmsoc_result_line[4])
            tmsoc_result_line_breaks=list(tmsoc_result_line[5])
            tmsoc_result_line_masked_sequence=str(tmsoc_result_line[6])


    #Change this to proper length of the segment and add compensation factors for + AND - start positions
    #for record in SeqIO.parse("sequence.fasta", "fasta"):
    #    print(record.seq)
    #    for n in range(len(record.seq)):
    #        b_factors.append(0)

    with open('sequence.fasta', "r") as sequence_file:
        lines=sequence_file.read().splitlines()
        #[] needs a ; at the end
        records_from_lines=''.join(lines).split('>;')
        print("Looking for chunk match in FASTA file...")
        for record in records_from_lines:
            print("Checking chunk:", record)
            record=record.split(";")
            if len(record) == 6:
                record_pdb_code=str(record[0])
                record_chain_id=str(record[1])
                record_start_location=int(record[2])
                record_end_location=int(record[3])
                record_breaks=list(record[4])
                record_sequence=str(record[5])


                #Integrity checking.
                if len(record_sequence)==abs(record_start_location-record_end_location):
                    print("WARNING! Length mismatch in", record_pdb_code, "chain ", record_chain_id)
                if abs(record_start_location-record_end_location) == len(record_breaks):
                    print("Mismatch distance matches number of breaks")

                # Checks if the record matches the current TMSOC chunk ind terms of PDB code and chain ID.
                if record_pdb_code == tmsoc_result_line_pdb_code and record_chain_id == tmsoc_result_line_chain_id:
                    # Sets a new variable to hold the current result chunk if it finds a match. This needs to be checked and reset.
                    current_result_record=record
                    current_result_record_chain_id=str(current_result_record[1])
                    current_result_record_start_location=int(current_result_record[2])
                    current_result_record_end_location=int(current_result_record[3])
                    current_result_record_breaks=list(current_result_record[4])
                    current_result_record_sequence=str(current_result_record[5])

                    #Needs testing to accomodate premature/- start positions
                    for n in range(len(record_sequence)+record_start_location):
                        b_factors.append(0)
            else:
                print("Record contained the wrong number of items, skipping to next record.")
                pass

    for item_position, line in enumerate(current_result):
        line = line.replace(',', ";")
        line = line.split(';')

        if str("1. TM segment(s) summary:") not in str(line) and str('2. Masked FASTA sequence:') not in str(line) and str('>') not in str(line) and len(line) == 7:
            print(line)
            tmh_sequence = line[0]
            #These start and end positions are NOT the real start positions, but rather distances from the beining of the first amino acid in the fasta file sequence.
            start_position = int(line[1])
            end_position = int(line[2])
            complexity_score = line[3]
            hydrophobicity_score = line[4]
            z_score = line[5]
            complexity = line[6]

            for residue_number, residue in enumerate(current_result_record_sequence):
                if residue_number <= end_position and residue_number >= start_position and end_position - start_position > 1:
                    #need to add option for complexity score
                    b_factors[residue_number+current_result_record_start_location] = b_factor_complexity(str(complexity))
                else:
                    pass

        #elif str(">") in str(line):
            #subchain_id = line[0].replace('\n', '')
            #subchain_id = subchain_id.replace('>', '_')
            #subchain_id = subchain_id.replace(':', '_')

    with open("bfactors_%s_%s.txt" % (record_pdb_code, current_result_record_chain_id), 'w') as bfactors_output_file:
        for i in b_factors:
            bfactors_output_file.write("%s\n" % i)
