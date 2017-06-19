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
from io import open

print("Processing Phobius output...")
with open("phobius_output.txt", 'r') as phobius_output:
    content = phobius_output.readlines()
    for line in content:

        content_for_id_written = False
        # This needs tidying up.
        line = line.replace('  ', ' ')
        line = line.replace('  ', ' ')
        line = line.replace('  ', ' ')
        line = line.replace('  ', ' ')

        line = line.split(' ')
        # print(line)
        row_type = line[0]

        if row_type == "FT":
            structure = line[1]
            start_position = line[2]
            end_position = line[3]
            if structure == "TRANSMEM":
                # Writes the segment location file
                print("Checking tmh at positions %s-%s..." %
                      (start_position, end_position))
                output_locations = str(start_position) + \
                    "," + str(end_position)
                print("Writing segments.txt file...")

                # The shell script removes this file between runs, so append should
                # not be an issue.
                with open("TMsegments.txt", 'a') as sequence_output:
                    sequence_output.write(output_locations + " ")
                content_for_id_written = True
            elif structure == "SIGNAL":
                # Writes the segment location file
                print("Checking tmh of signal sequence at positions %s-%s..." %
                      (start_position, end_position))
                output_locations = str(start_position) + \
                    "," + str(end_position)
                print("Writing segments.txt file...")

                # The shell script removes this file between runs, so append should
                # not be an issue.
                with open("TMsegments.txt", 'a') as sequence_output:
                    sequence_output.write(output_locations + " ")
                content_for_id_written = True
            else:
                pass

        elif row_type == "ID":
            sequence_id = line[1]
            structure = "NULL"
            start_position = "NULL"
            end_position = "NULL"

        elif "//" in row_type:

            if content_for_id_written is False:
                # This is a work around. TMSOC cannot interpret NONE type
                # input, so a placeholder is needed.
                with open("TMsegments.txt", 'a') as sequence_output:
                    sequence_output.write("0")
            else:
                pass

            with open("TMsegments.txt", 'a') as sequence_output:
                sequence_output.write("\n")

        else:
            pass
