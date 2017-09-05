from __future__ import division
from Bio import SeqIO
import subprocess
import re
import sys
from Bio.PDB import PDBParser
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator
import scipy
from scipy import stats
import numpy as np
import pylab
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

def complexity_extraction(sequence, tmh_locations):
    '''
    Returns a list with values of complexities of TMHs.
    '''
    list_of_complexities=[]
    # Generates the TMSOC input files from the source file.
    with open("TMsegments.txt", 'w') as tm_segments_file:
        tm_segments_file.write(tmh_locations)

    with open("sequence.fasta", 'w') as fasta_file:
        fasta_file.write(">placeholderheader\n")
        fasta_file.write(str(sequence))

    # Runs TMSOC.
    perl_script_output=subprocess.check_output(["perl", "TMSOC.pl", "sequence.fasta", "TMsegments.txt"])

    # Processes the TMSOC input files. Normally the TMSOC output would be a multi-line file, not a string.
    list_of_output_lines=str(perl_script_output.decode('UTF-8'))
    list_of_output_lines=list_of_output_lines.replace("\n", ":")
    list_of_output_lines=list_of_output_lines.split(":")
    for line in list_of_output_lines:

        line = str(line.replace(',', ";"))
        line = line.split(";")
        #First we can ignore the "junk" lines
        if len(line)==7:
            tmh_sequence = line[0]
            start_position = int(line[1])
            end_position = int(line[2])
            complexity_score = float(line[3])
            hydrophobicity_score = float(line[4])
            z_score = line[5]
            complexity = line[6]
            #Lenght restrictions
            if start_position == -1 and end_position == -1:
                print("Skipping placeholder entry")
                pass
            else:
                list_of_complexities.append(complexity_score)
    return(list_of_complexities)

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

    # Generate a list of empty lists, one for each tmh set. This avoids issues
    # of exceeding list indices in advance.
    print("Maximum tmh count in", input_file, "is", tmd_count)
    list_of_complexity_scores_in_tmh = []
    for n in range(tmd_count):
        list_of_complexity_scores_in_tmh.append([])

    # Now we can iterate through the records inserting complexity scores into
    # the empty lists.
    for record in SeqIO.parse(filename, input_format):
        this_record_tmd_count = 0
        # Sequence fasta file
        sequence = record.seq
        # TMH positions file
        tmh_positions = str("")
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                tmh_positions=tmh_positions+(str(f.location.start)+ ","+ str(f.location.end) + " ")

        for n, i in enumerate(complexity_extraction(sequence, tmh_positions)):
            list_of_complexity_scores_in_tmh[n].append(i)
