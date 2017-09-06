from __future__ import division
from decimal import *
from Bio import SeqIO
import subprocess
import scipy
from scipy import stats
import numpy as np


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


# To do:
# Separate hydrophobicity calculation - more efficient and flexible
# Length restrictions and statistics for final results
# Control datasets
# Generate graphics

# Input files should be obtained in text format downloaded from Uniprot
# and moved to the same directory as this script.
input_filenames = [
    "GPCR_UniRef50.txt"
]

# Parameters for tmh allowances
minimum_tmd_length = 16
maximum_tmd_length = 38
feature_type = "TRANSMEM"
alternative_feature = "INTRAMEM"

#### Complexity calculation code ####


def tmsoc_calculation(sequence, tmh_locations):
    '''
    Returns a list with values of complexities of TMHs.
    '''

    list_of_tmh_features = []
    # Generates the TMSOC input files from the source file.
    with open("TMsegments.txt", 'w') as tm_segments_file:
        tm_segments_file.write(tmh_locations)
    with open("sequence.fasta", 'w') as fasta_file:
        fasta_file.write(">placeholderheader\n")
        fasta_file.write(str(sequence))
    # Runs TMSOC.
    perl_script_output = subprocess.check_output(
        ["perl", "TMSOC.pl", "sequence.fasta", "TMsegments.txt"])
    # Processes the TMSOC input files. Normally the TMSOC output would be a
    # multi-line file, not a string.
    list_of_output_lines = str(perl_script_output.decode('UTF-8'))
    list_of_output_lines = list_of_output_lines.replace("\n", ":")
    list_of_output_lines = list_of_output_lines.split(":")
    for line in list_of_output_lines:
        line = str(line.replace(',', ";"))
        line = line.split(";")
        # First we can ignore the "junk" lines
        if len(line) == 7:
            tmh_sequence = line[0]
            start_position = int(line[1])
            end_position = int(line[2])
            complexity_score = float(line[3])
            hydrophobicity_score = float(line[4])
            z_score = line[5]
            complexity = line[6]
            # Lenght restrictions-these should obey parameters given at
            # begining of script.
            if abs(start_position - end_position) < minimum_tmd_length or abs(start_position - end_position) > maximum_tmd_length:
                complexity_score = str("null")

            list_of_tmh_features.append(complexity_score)

    return(list_of_tmh_features)

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
    # number of TMDs. This avoids list index exceeded errors later on.
    tmd_count = 0
    for record in SeqIO.parse(filename, input_format):
        this_record_tmd_count = 0
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                this_record_tmd_count = this_record_tmd_count + 1
        if this_record_tmd_count > tmd_count:
            tmd_count = this_record_tmd_count
            if this_record_tmd_count > 7:
                print(record.id, " contained more than 7 TMHs.")

    # Generate a list of empty lists, one for each tmh set.
    print("Maximum tmh count in", input_file, "is", tmd_count)

    list_of_complexity_scores_in_tmh = []
    list_of_hydrophobicity_scores_in_tmh = []

    for n in range(tmd_count):
        list_of_complexity_scores_in_tmh.append([])
        list_of_hydrophobicity_scores_in_tmh.append([])

    # Now we can iterate through the records inserting complexity scores into
    # the empty lists.
    for record in SeqIO.parse(filename, input_format):

        # Sequence fasta file
        sequence = record.seq

        # TMH positions file
        tmh_positions = str("")
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                tmh_positions = tmh_positions + \
                    (str(f.location.start) + "," + str(f.location.end) + " ")

        # Adds the complexity score of a helix at (for example the 4th helix)
        # to the complexity list of lists (for example in the 4th position)
        for n, i in enumerate(tmsoc_calculation(sequence, tmh_positions)):
            # Null entries are added for TMHs that are not within length
            # restrictions.
            if i == "null":
                pass
            else:
                list_of_complexity_scores_in_tmh[n].append(int(i))

        # Adds to the hydrophobicity list of lists
        # for n, i in enumerate(hydrophobicity_calculation(sequence, tmh_positions)):
        #    list_of_hydrophobicity_scores_in_tmh[n].append(int(i))

    # Reporting the results
    print(filename, "results\n")

    # stats for complexity
    for n, i in enumerate(list_of_complexity_scores_in_tmh):

        # n is the index, so for human readable numbers we need to add 1. i.e
        # the first helix is n=0, so we report it as n+1.
        print("TMH ", n + 1)
        print("Mean complexity:", np.mean(i), ", N:", len(i))

        if n + 1 < len(list_of_complexity_scores_in_tmh):
            print("KS of TMH ", n + 1, " to ", n + 2, ":", scipy.stats.ks_2samp(
                list_of_complexity_scores_in_tmh[n], list_of_complexity_scores_in_tmh[n + 1]))

    print("\n")

    '''
    # stats for hydrophobicity
    for n, i in enumerate(list_of_hydrophobicity_scores_in_tmh):
        # n is the index, so for human readable numbers we need to add 1. i.e
        # the first helix is n=0, so we report it as n+1.

        print("TMH ", n + 1)
        print("Mean Hydrophobicity:", np.mean(i), ", N:", len(i))

        if n + 1 < len(list_of_hydrophobicity_scores_in_tmh):
            print("KS of TMH ", n + 1, " to ", n + 2, ":", scipy.stats.kruskal(
                list_of_hydrophobicity_scores_in_tmh[n], list_of_hydrophobicity_scores_in_tmh[n + 1]))
    '''
