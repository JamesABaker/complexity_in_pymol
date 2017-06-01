echo -n 'What is the file name of the PDB? Example answer: 1a91.pdb'
echo
read pdb_id

# Clearing files from previous run.
rm TMsegments.txt
rm sequence.fasta
rm tmsoc_output.txt
rm phobius_output.txt
rm bfactors.txt

# Scripts in order. If you only have python installed, python3 can be renamed as python.
python3 pdb2phobius.py "$pdb_id"
sleep 0.5 #The pdb2phobius script takes a split second to run on some files, but for some reason, wait doesn't actually allow it to finish.
python3 phobius2tmsoc.py
perl TMSOC.pl sequence.fasta TMsegments.txt > tmsoc_output.txt
python3 tmsoc2bfactorlist.py "$pdb_id"
