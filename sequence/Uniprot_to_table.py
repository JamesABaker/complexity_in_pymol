from __future__ import division
from Bio import SeqIO
import numpy as np
import os
import subprocess
import re
import sys
from Bio.PDB import PDBParser
from Bio import SeqIO
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator
import scipy
from scipy import stats
import numpy as np
import pylab
import sys
from decimal import *


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


#### Complexity code ####

def complexity_extraction(sequence, tmh_locations):
    '''
    Fills the empty lists with values of complexities of TMHs.
    '''

    # Generates the TMSOC input files from the source file.

    # Runs TMSOC.

    # Processes the TMSOC input files.
    new_tmsoc_result = "1. TM segment(s) summary:"
    all_results = []
    current_result = []

    for line in open('tmsoc_output.txt'):
        if line.startswith(new_tmsoc_result) and current_result:
            # if line starts with thew new TMSOC result token and the current chunk
            # is not empty
            # add not empty chunk to chunks
            all_results.append(current_result[:])
            current_result = []  # make current chunk blank
        # just append a line to the current chunk on each iteration
        current_result.append(line)

    # append the last chunk outside the loop
    all_results.append(current_result)

    for current_result in all_results:
        b_factors = []
        print("Processing TMSOC result chunk:", current_result)
        for tmsoc_result_line in current_result:
            if ">" in tmsoc_result_line:
                tmsoc_result_line = tmsoc_result_line.split(";")
                # Here we ignore the first entry, which is simply ">"
                tmsoc_result_line_pdb_code = str(tmsoc_result_line[1])
                tmsoc_result_line_chain_id = str(tmsoc_result_line[2])
                tmsoc_result_line_start_location = int(tmsoc_result_line[3])
                tmsoc_result_line_end_location = int(tmsoc_result_line[4])
                tmsoc_result_line_breaks = list(tmsoc_result_line[5])
                # masked fast doesn't always exist
                # tmsoc_result_line_masked_sequence=str(tmsoc_result_line[6])

        # Change this to proper length of the segment and add compensation factors for + AND - start positions
        # for record in SeqIO.parse("sequence.fasta", "fasta"):
        #    print(record.seq)
        #    for n in range(len(record.seq)):
        #        b_factors.append(0)

        with open('sequence.fasta', "r") as sequence_file:
            lines = sequence_file.read().splitlines()
            #[] needs a ; at the end
            records_from_lines = ''.join(lines).split('>;')
            print("Looking for chunk match in FASTA file...")
            for record in records_from_lines:
                print("Checking chunk:", record)
                record = record.split(";")
                if len(record) == 6:
                    record_pdb_code = str(record[0])
                    record_chain_id = str(record[1])
                    record_start_location = int(record[2])
                    record_end_location = int(record[3])
                    record_breaks = list(record[4])
                    record_sequence = str(record[5])

                    # Integrity checking.
                    if len(record_sequence) == abs(record_start_location - record_end_location):
                        print("WARNING! Length mismatch in",
                              record_pdb_code, "chain ", record_chain_id)
                    if abs(record_start_location - record_end_location) == len(record_breaks):
                        print("Mismatch distance matches number of breaks")

                    # Checks if the record matches the current TMSOC chunk ind
                    # terms of PDB code and chain ID.
                    if record_pdb_code == tmsoc_result_line_pdb_code and record_chain_id == tmsoc_result_line_chain_id:
                        # Sets a new variable to hold the current result chunk if
                        # it finds a match. This needs to be checked and reset.
                        current_result_record = record
                        current_result_record_chain_id = str(
                            current_result_record[1])
                        current_result_record_start_location = int(
                            current_result_record[2])
                        current_result_record_end_location = int(
                            current_result_record[3])
                        current_result_record_breaks = list(
                            current_result_record[4])
                        current_result_record_sequence = str(
                            current_result_record[5])

                        # Needs testing to accomodate premature/- start
                        # positions
                        for n in range(len(record_sequence) + record_start_location):
                            b_factors.append(0)
                else:
                    print(
                        "Record contained the wrong number of items, skipping to next record.")
                    pass

        # This is the list that will be printed as a graph.
        complexity_scores_for_plot = []

        for item_position, line in enumerate(current_result):
            line = line.replace(',', ";")
            line = line.split(';')

            if str("1. TM segment(s) summary:") not in str(line) and str('2. Masked FASTA sequence:') not in str(line) and str('>') not in str(line) and len(line) == 7:
                print(line)
                tmh_sequence = line[0]
                # These start and end positions are NOT the real start positions,
                # but rather distances from the beining of the first amino acid in
                # the fasta file sequence.
                start_position = int(line[1])
                end_position = int(line[2])
                complexity_score = float(line[3])
                hydrophobicity_score = float(line[4])
                z_score = line[5]
                complexity = line[6]
                if start_position == -1 and end_position == -1:
                    print("Skipping placeholder entry")
                    pass
                else:
                    complexity_scores_for_plot.append(complexity_score)

                    for residue_number, residue in enumerate(current_result_record_sequence):
                        if residue_number <= end_position and residue_number >= start_position and end_position - start_position > 1:
                            # need to add option for complexity interpretation
                            b_factors[
                                residue_number + current_result_record_start_location] = str(complexity_score)
                            # b_factors[residue_number+current_result_record_start_location] = b_factor_complexity(str(complexity))
                        else:
                            pass


# Input files should be obtained in text format downloaded from Uniprot
# and moved to the same directory as this script.
input_filenames = [
    "GPCR_UniRef50.txt"
]

# Parameters for tmh allowances
minimum_tmd_length = 16
maximum_tmd_length = 38
length_excluded_tmds = []
feature_type = "TRANSMEM"
alternative_feature = "INTRAMEM"

# Before we conduct the analysis we must determine how many empty TMH
# lists to make. For example if the protein with the most TMHs in the
# dataset has 13 TMHs, then there will be 13 empty lists made in a list
# representing the 1-13 possible TMH locations.
for input_file in input_filenames:
    # These are the parameters used by the biopython Seq.IO module
    filename = input_file
    input_format = "swiss"

    # We need to check against nearby features to prevent overlapping
    # flanking regions. Note here we want to avoid clashing with INTRAMEM
    # regions, since their flanking regions may be similar, however INTRAMEM regions
    # should not be included in the logged transmembrane features since
    # they may have extremely short transmembrane sequences which would
    # disrupt very much so the alignment of the flanking regions.

    # We iterate through each record, parsed by biopython to find the most
    # number of TMDs.
    tmd_count = 0
    for record in SeqIO.parse(filename, input_format):
        this_record_tmd_count = 0
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                this_record_tmd_count = this_record_tmd_count + 1
        if this_record_tmd_count > tmd_count:
            tmd_count = this_record_tmd_count

    print("Maximum tmh count in", input_file, "is", tmd_count)
    list_of_complexity_scores_in_tmh = []
    for n in range(tmd_count):
        list_of_complexity_scores_in_tmh.append([])

    # Check lists are empty.
    for complexity_list in list_of_complexity_scores_in_tmh:
        if not complexity_list:
            print("Lists not empty.")
    # Checks the max_tmh_count is the same length as the number of lists
    if max_tmh_count != len(list_of_complexity_scores_in_tmh):
        print("Lists do not match the max number of TMHs in this dataset")

    # Now we can iterate through the records inserting complexity scores into
    # the list.
    for record in SeqIO.parse(filename, input_format):
        this_record_tmd_count = 0
        #Sequence fasta file
        header=str(">",record.id)
        sequence=record.seq
        
        #TMH positions file
        for i, f in enumerate(record.features):
            if f.type == feature_type:
