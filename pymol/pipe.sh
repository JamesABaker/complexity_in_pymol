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


echo "Complexity For Pymol  Copyright (C) 2017 James Alexander Baker"
echo "This program comes with ABSOLUTELY NO WARRANTY."
echo "This is free software, and you are welcome to redistribute it"
echo "under certain conditions."

# The current parameter should be just the pdb code.
pdb_id=$1

# Clearing files from previous run.
echo 'Removing old log files...'
rm TMsegments.txt
rm sequence.fasta
rm tmsoc_output.txt
rm phobius_output.txt

# Regenerating empty log files. Some python setups need a file to exist in order to write to it.
echo 'Creating new log files...'
touch TMsegments.txt
touch sequence.fasta
touch tmsoc_output.txt
touch phobius_output.txt


# Scripts that process the inputs and outputs in order. If you only have python installed, python3 can be renamed as python.
echo 'Parsing' pdb_id 'to Phobius.'
python3 pdb2phobius.py "$pdb_id"
sleep 1 #The pdb2phobius script takes a split second to run on some files, but for some reason, wait doesn't actually allow it to finish.
echo 'Parsing Phobius output to TMSOC.'
python3 phobius2tmsoc.py
echo 'Running TMSOC.'
perl TMSOC.pl sequence.fasta TMsegments.txt > tmsoc_output.txt
echo 'Transforming TMSOC output to bfactor list.'
python3 tmsoc2bfactorlist.py
