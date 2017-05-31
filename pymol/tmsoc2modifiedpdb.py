with open("tmsoc_output.txt", 'r') as tmsoc_output:
    content=tmsoc_output.readlines()
    for line_number, line in enumerate(content):
        if "2. Masked FASTA sequence:" in line:
            max_line=line_number
    for new_line_number, line in enumerate(content):

        if new_line_number>0 and new_line_number<max_line:
            line=line.replace(',',";")
            line=line.split(';')
            tmh_sequence=line[0]
            start_position=line[1]
            end_position=line[2]
            score1=line[3]
            score2=line[4]
            score3=line[5]
            complexity=line[6]
            print(complexity)
