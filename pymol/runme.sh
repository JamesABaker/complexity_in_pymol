echo -n 'What is the file name of the PDB? Example answer: 1a91.pdb'
echo
read pdb_id

# Clearing files from previous run.
echo 'Removing old log files.'
rm TMsegments.txt
rm sequence.fasta
rm tmsoc_output.txt
rm phobius_output.txt

# Regenerating empty log files. Some python setups need a file to exist in order to write to it.
echo 'Creating new log files.'
touch TMsegments.txt
touch sequence.fasta
touch tmsoc_output.txt
touch phobius_output.txt


# Scripts that process the inputs and outputs in order. If you only have python installed, python3 can be renamed as python.
echo 'Parsing' pdb_id 'to Phobius.'
python3 pdb2phobius.py "$pdb_id"
sleep 0.5 #The pdb2phobius script takes a split second to run on some files, but for some reason, wait doesn't actually allow it to finish.
echo 'Parsing Phobius output to TMSOC.'
python3 phobius2tmsoc.py
echo 'Running TMSOC.'
perl TMSOC.pl sequence.fasta TMsegments.txt > tmsoc_output.txt
echo 'Transforming TMSOC output to bfactor list.'
python3 tmsoc2bfactorlist.py
