from __future__ import division
from Bio import SeqIO
import numpy as np
import os
import subprocess
import re
import sys

from __future__ import division


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



# Input files should be obtained in text format downloaded from Uniprot
# and moved to the same directory as this script.
input_filenames = [
    "UniArch.txt",
    "UniBacilli.txt",
    "UniCress.txt",
    "UniEcoli.txt",
    "UniER.txt",
    "UniFungi.txt",
    "UniGolgi.txt",
    "UniHuman.txt",
    "UniPM.txt"]


# Parameters for tables
flank_lengths = [5, 10, 20]
minimum_tmd_length = 16
maximum_tmd_length = 38
length_excluded_tmds = []
flank_clash_amendment = [True, False]

# For each file, a table is generated for each of the flank lengths set.
for input_file in input_filenames:
    for flank_clash_amendment_status in flank_clash_amendment:
        for flank_length in flank_lengths:

            # Single-pass and mult-pass are treated differently because at the end
            # we report on the numbers so we avoid re-iterating through the larger
            # datasets, and throughout the process they are handled differently for
            # example when looking for clashes of flanking regions.
            number_of_records = 0
            number_of_records_correct_length = 0
            number_of_records_single = 0
            number_of_records_correct_length_single = 0
            number_of_records_multi = 0
            number_of_records_correct_length_multi = 0

            # File output name
            output_filename = input_file.replace(
                ".txt", "_%s_flanklength_flankclash%s.csv" % (flank_length, str(flank_clash_amendment_status)))
            half_flanks_output_filename = input_file.replace(
                ".txt", "_%s_flanklength_flankclash%s_only_half_flanks.csv" % (flank_length, str(flank_clash_amendment_status)))
            full_flanks_output_filename = input_file.replace(
                ".txt", "_%s_flanklength_flankclash%s_only_full_flanks.csv" % (flank_length, str(flank_clash_amendment_status)))

            # The header row in the file.
            with open(output_filename, 'w') as my_file:
                my_file.write("Name and description, ID, N terminal inside/outside, tmh start location, tmh end location, full protein sequence, tmh sequence, N flank sequence, C flank sequence, transmembrane helix sequential number, number of transmembrane helices in protein\n")
            my_file.closed
            with open(half_flanks_output_filename, 'w') as my_file:
                my_file.write("Name and description, ID, N terminal inside/outside, tmh start location, tmh end location, full protein sequence, tmh sequence, N flank sequence, C flank sequence, transmembrane helix sequential number, number of transmembrane helices in protein\n")
            my_file.closed
            with open(full_flanks_output_filename, 'w') as my_file:
                my_file.write("Name and description, ID, N terminal inside/outside, tmh start location, tmh end location, full protein sequence, tmh sequence, N flank sequence, C flank sequence, transmembrane helix sequential number, number of transmembrane helices in protein\n")
            my_file.closed

            # These are the parameters used by the biopython Seq.IO module
            filename = input_file
            input_format = "swiss"
            feature_type = "TRANSMEM"
            subcellular_location = "TOPO_DOM"
            # We need to check against nearby features to prevent overlapping
            # flanking regions. Note here we want to avoid clashing with INTRAMEM
            # regions, since their flanking regions may be similar, however INTRAMEM regions
            # should not be included in the logged transmembrane features since
            # they may have extremely short transmembrane sequences which would
            # disrupt very much so the alignment of the flanking regions.
            avoid_features = ["TRANSMEM", "INTRAMEM"]
            unknown = 0
            list_of_multipass_total_tmh_counts = []
            reliable_flank_length = flank_length

            # We iterate through each record, parsed by biopython.
            for record in SeqIO.parse(filename, input_format):
                new_record = True
                tmd_count = 0
                for i, f in enumerate(record.features):
                    if f.type == feature_type:

                        n_terminal_start = "None"
                        tmh_record = []
                        tmd_count = tmd_count + 1

                        # Calls the parsing functions as strings for later use.
                        id_of_record = record.id
                        name_of_record = str(
                            record.description).replace(",", "")

                        # The first separation of single or multi-pass is based on
                        # annotation itself, however this is prone to error, so is
                        # checked later.
                        if "; Single-pass" in str(record):
                            single_or_multi_topology = "Single-pass"
                        if "; Multi-pass" in str(record):
                            single_or_multi_topology = "Multi-pass"

                        # Some transmembrane annotations have unknown sequence
                        # positions where it is ambiguous. These features are
                        # discounted.
                        if "UnknownPosition" in str(f.location):
                            pass
                            unknown = unknown + 1
                            print id_of_record, "had an unknown position TMH."
                        else:
                            full_sequence = str(record.seq)
                            tmh_start = int(f.location.start)
                            tmh_stop = int(f.location.end)
                            tmh_sequence = str(
                                record.seq[(f.location.start):(f.location.end)])

                            # If a clash is detected in either flank, the
                            # script will compensate for the clash and stop
                            # checking to see if there are additional
                            # clashes at that flank. This is reset per TMH.
                            C_terminal_flank_clash = False
                            N_terminal_flank_clash = False

                            # We iterate through each feature in the record, to
                            # assertain if there are any clashes.
                            for each_features in record.features:
                                # Another check for unknown positions to rule
                                # out.
                                if "UnknownPosition" in str(each_features.location):
                                    pass
                                else:
                                    #  We identify if any of these features are within the flanking residue distance.
                                    # If they are clashes, then the flank length is half the distance between the features.
                                    # Note now, we only need to worry about the flank length of the logged feature, since the other
                                    # flank length will be dealt with when
                                    # iterating through the features looking for TM
                                    # regions.

                                    # Let's first check if the N terminus of this
                                    # feature comes after the C terminus of the
                                    # feature we're logging.
                                    if C_terminal_flank_clash is not True:
                                        if each_features.location.start > f.location.end:
                                            # Is the feature and the N terminal flank starting within the flank
                                            # of the C terminus of the feature that
                                            # we're logging?
                                            if (each_features.location.start - reliable_flank_length) < (f.location.end + reliable_flank_length)and flank_clash_amendment_status == True:
                                                if str(each_features.type) in avoid_features:
                                                    # There will be a clash/overlapping
                                                    # of transmembrane flank lengths at
                                                    # the C flank.
                                                    max_flank_size = each_features.location.start - f.location.end
                                                    flank2_length = max_flank_size / 2
                                                    C_terminal_flank_clash = True
                                                else:
                                                    # There is a feature, but it's not
                                                    # a transmembrane region, so
                                                    # flanking regions won't be double
                                                    # counted
                                                    flank2_length = reliable_flank_length
                                            else:
                                                flank2_length = reliable_flank_length
                                        else:
                                            flank2_length = reliable_flank_length

                                    # Now we will check the same as above, however
                                    # for features that the C terminas falls before
                                    # the feature being logged N-terminus
                                    if N_terminal_flank_clash is not True:
                                        if each_features.location.end < f.location.start:
                                            if (each_features.location.end + reliable_flank_length) > (f.location.start - reliable_flank_length) and flank_clash_amendment_status == True:
                                                if str(each_features.type) in avoid_features:
                                                    max_flank_size = f.location.start - each_features.location.end
                                                    flank1_length = max_flank_size / 2
                                                    N_terminal_flank_clash = True
                                                else:
                                                    flank1_length = reliable_flank_length
                                            else:
                                                flank1_length = reliable_flank_length
                                        else:
                                            flank1_length = reliable_flank_length

                                    # Now we will amend the flank length incase it
                                    # exceeds the maximum protein length, or the
                                    # 0th sequence item (the start of the
                                    # sequence)
                                    if (f.location.start - int(flank1_length)) >= 0:
                                        N_terminal_flank = str(
                                            record.seq[(f.location.start - int(flank1_length)):(f.location.start)])
                                    elif (f.location.start - int(flank1_length)) < 0:
                                        N_terminal_flank = str(
                                            record.seq[0:(f.location.start)])
                                    if (f.location.end + int(flank2_length)) > len(record.seq):
                                        C_terminal_flank = str(
                                            record.seq[(f.location.end):(int(len(record.seq)))])
                                    else:
                                        C_terminal_flank = str(
                                            record.seq[(f.location.end):(f.location.end + int(flank2_length))])

                            # Now, the orientation is determined for all features
                            # according to the cytoplasmic annotation. The total
                            # number of cytoplasmic annotations, and non
                            # cytoplasmic annotations are counted, and depending on
                            # the starting locations falling "inside" or "outside" the cytoplasm,
                            # the annotation is determined.

                            list_of_cyto_starts = []
                            list_of_non_cyto_starts = []

                            for feature_number, other_features in enumerate(record.features):

                                # Some features contain unknown locations. These are
                                # discarded.
                                if "UnknownPosition" in str(other_features.location):
                                    pass
                                else:
                                    # If the feature matches the cytoplasm...
                                    if other_features.type == subcellular_location:
                                        if "Cytoplasmic" in str(other_features.qualifiers):

                                            # The feature start locations are
                                            # recorded as being cytoplasmic or not.
                                            # This, combined with knowing the start
                                            # location, allows us to go through the
                                            # TMHs and figure out which way each is
                                            # oriented quickly.
                                            list_of_cyto_starts.append(
                                                int(other_features.location.start))

                                            total_tmd_count = 0
                                            for number_of_features, a_feature in enumerate(record.features):

                                                if a_feature.type == feature_type:
                                                    total_tmd_count = total_tmd_count + 1
                                            if total_tmd_count > 1 and new_record == True:
                                                list_of_multipass_total_tmh_counts.append(
                                                    total_tmd_count)
                                            new_record = False

                                        else:
                                            list_of_non_cyto_starts.append(
                                                int(other_features.location.start))
                                    else:
                                        pass

                            # We will now check if the residues preceding the
                            # TMD are intra, or extra cytoplasmic. Because this
                            # calls locations below the feature location start,
                            # it must be called only when above 1. If this
                            # cannot be done, the topology can be infered in
                            # the next script.
                            if n_terminal_start == "none" and f.location.start > 1:
                                previous_feautre_location = f.location.start - 1
                                for index, a_features in enumerate(record.features):
                                    if a_features.type == subcellular_location and a_features.location.start < previous_feautre_location and a_features.location.end > previous_feautre_location:
                                        if "Cytoplasmic" in str(a_features.qualifiers):
                                            n_terminal_start == "Inside"
                                        else:
                                            n_terminal_start = "Outside"

                            # If the previous method did not work to identify topology, we can still
                            # imply the topology.
                            if n_terminal_start == "None":
                                # This checks that the subcellular starting
                                # locations exist.
                                if len(list_of_cyto_starts) > 0 and len(list_of_non_cyto_starts) > 0:

                                    # We need to furthermore check whether the
                                    # starting location was inside or outside.
                                    if min(list_of_non_cyto_starts) < min(list_of_cyto_starts):

                                        # The modulo is used to check if the
                                        # current tmd number is even or odd. From
                                        # this and the first TMD we know whether
                                        # this TMD starts inside or outside.
                                        if tmd_count % 2 == 0:
                                            n_terminal_start = "Inside"
                                        elif tmd_count % 2 == 1:
                                            n_terminal_start = "Outside"

                                    elif min(list_of_cyto_starts) < min(list_of_non_cyto_starts):
                                        if tmd_count % 2 == 0:
                                            n_terminal_start = "Outside"
                                        elif tmd_count % 2 == 1:
                                            n_terminal_start = "Inside"
                                    else:
                                        pass

                            # Checks if the single_or_multi_topology variable
                            # exists (it would have been created earlier)
                            if 'single_or_multi_topology' not in locals() or 'single_or_multi_topology' not in globals():
                                if total_tmd_count == 1:
                                    single_or_multi_topology = "Single-pass"
                                if total_tmd_count > 1:
                                    single_or_multi_topology = "Multi-pass"

                            if "Inside" in n_terminal_start or "Outside" in n_terminal_start:
                                #+1s are used since slices originally call how many steps to iterate rather than the sequence postion. This matches the Uniprot sequence numbering
                                tmh_record = [name_of_record, id_of_record, n_terminal_start, tmh_start + 1, tmh_stop,
                                              full_sequence, tmh_sequence, N_terminal_flank, C_terminal_flank, tmd_count, total_tmd_count]

                                number_of_records = number_of_records + 1

                                if total_tmd_count == 1:
                                    number_of_records_single = number_of_records_single + 1
                                if total_tmd_count > 1:
                                    number_of_records_multi = number_of_records_multi + 1

                                if len(tmh_sequence) >= minimum_tmd_length and len(tmh_sequence) <= maximum_tmd_length:
                                    with open(output_filename, 'a') as my_file:
                                        for i in tmh_record:
                                            my_file.write(str(i))
                                            my_file.write(",")
                                        my_file.write("\n")
                                    number_of_records_correct_length = number_of_records_correct_length + 1

                                    # Now we see if the flanks are either half
                                    # length, or full length.
                                    if len(C_terminal_flank) == max_flank_size and len(N_terminal_flank) == max_flank_size:
                                        with open(full_flanks_output_filename, 'a') as my_file:
                                            for i in tmh_record:
                                                my_file.write(str(i))
                                                my_file.write(",")
                                            my_file.write("\n")

                                    if len(C_terminal_flank) >= max_flank_size / 2 and len(N_terminal_flank) >= max_flank_size / 2:
                                        with open(half_flanks_output_filename, 'a') as my_file:
                                            for i in tmh_record:
                                                my_file.write(str(i))
                                                my_file.write(",")
                                            my_file.write("\n")

                                    if total_tmd_count == 1:
                                        number_of_records_correct_length_single = number_of_records_correct_length_single + 1
                                    if total_tmd_count > 1:
                                        number_of_records_correct_length_multi = number_of_records_correct_length_multi + 1
                                else:
                                    length_exclusion_info = str(
                                        id_of_record) + "_" + str(tmd_count)
                                    length_excluded_tmds.append(
                                        length_exclusion_info)

                            # No records should be here with none. This is for
                            # debugging only.
                            elif "None" in n_terminal_start:
                                pass
                            else:
                                pass
        # General information regarding the output file useful as a log.
        print input_file, ", flank_clash_amendment_status:", flank_clash_amendment_status
        print "Number of TMHs in dataset:", number_of_records
        print "Number of TMHs after dumping incorrect lengths:", number_of_records_correct_length
        print "Mean number of tmhs in multipass proteins:", np.mean(list_of_multipass_total_tmh_counts)
        print "S.D of tmhs in multipass proteins:", np.std(list_of_multipass_total_tmh_counts)
        print "Records"
        print number_of_records
        print number_of_records_correct_length
        print "Single-pass"
        print "Total:", number_of_records_single
        print "After length exclusion:", number_of_records_correct_length_single
        print "Multi-pass"
        print "Total:", number_of_records_multi
        print "After length exclusion:", number_of_records_correct_length_multi

        exclusion_ids_output_filename = input_file.replace(
            ".txt", "_%s_flanklength_flankclash%s_logged_lengthexclusionIDs.txt" % (flank_length, str(flank_clash_amendment_status)))
        with open(exclusion_ids_output_filename, 'a') as my_file:
            my_file.write("# ID_HelixNumber\n")
        my_file.closed
        for item in length_excluded_tmds:
            with open(exclusion_ids_output_filename, 'a') as my_file:
                my_file.write(item)
                my_file.write("\n")
            my_file.closed



#### Complexity code ####

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

                    # Needs testing to accomodate premature/- start positions
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
                        #b_factors[residue_number+current_result_record_start_location] = b_factor_complexity(str(complexity))
                    else:
                        pass

        # elif str(">") in str(line):
            #subchain_id = line[0].replace('\n', '')
            #subchain_id = subchain_id.replace('>', '_')
            #subchain_id = subchain_id.replace(':', '_')

    # Plot the complexity histogram.
    if len(complexity_scores_for_plot) == 0:
        print("Skipping chain %s since no complexity scores were found. This is most likely because there are no TMHs in this chain." %
              current_result_record_chain_id)
        pass
    else:
        print("Complexity scores")
        start_count = 0
        complexity_count = []
        for i in complexity_scores_for_plot:
            start_count = start_count + 1
            print(start_count, i)
            complexity_count.append(start_count)
        print("Complexity scores for vector of %s\n" % record_pdb_code)
        print(','.join(str(i) for i in complexity_scores_for_plot))

        plt.plot(complexity_count, complexity_scores_for_plot)
        plt.xlabel('Helix Number for %s_%s' %
                   (record_pdb_code, current_result_record_chain_id))
        plt.ylabel('Complexity')
        filename = record_pdb_code.replace(
            '.pdb', '') + "_" + current_result_record_chain_id + ".pdf"
        plt.gcf().subplots_adjust(bottom=0.2)
        plt.savefig(filename)
        # plt.show()

        plt.clf()
        plt.cla()

    with open("bfactors_%s_%s.txt" % (record_pdb_code, current_result_record_chain_id), 'w') as bfactors_output_file:
        for i in b_factors:
            bfactors_output_file.write("%s\n" % i)
