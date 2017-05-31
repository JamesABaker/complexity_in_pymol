# Clearing files from previous run.
rm TMsegments.txt
rm sequence.fasta
rm temp.fasta
rm tmsoc_output.txt
rm phobius_output.txt


python3 pdb2phobius.py
sleep 0.5 #The pdb2phobius script takes a split second to run on some files, but for some reason, wait doesn't actually allow it to finish.
python3 phobius2tmsoc.py
perl TMSOC.pl sequence.fasta TMsegments.txt > tmsoc_output.txt
python3 tmsoc2modifiedpdb.py
