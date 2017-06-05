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
