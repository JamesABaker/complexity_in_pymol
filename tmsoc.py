# Import PyMOL's stored module.  This will allow us with a
# way to pull out the PyMOL data and modify it in our script.
# See below.
import sys
from Bio.PDB import *
from pymol import stored

#Notes for bfactor editing
#get_bfactor(self)
# set_bfactor(self, bfactor)

#Convert pdb to fasta

fasta = cmd.get_fastastr('all')

# Writing input files for TMSOC
with open("inputseq.fasta", 'w') as temp_tmh_fasta:
    temp_tmh_fasta.write(fasta)

with open("TMsegments.txt", 'w') as temp_tmh_seq:
    pass

#Needs to be rewritten to include pdb, not swissprot format.
#The
for i, f in enumerate(record.features):
    if f.type == feature_type:
        try:
            # Coordinates from the original TMSOC project look like
            # 37,63 74,96
            #These need to be rewritten to be compatible with pdb files
            with open("TMsegments.txt", 'a') as temp_tmh_seq:
                temp_tmh_seq.write(str(str(f.location.start + 1) + "," + str((f.location.end + 1)) + " "))
        except TypeError:
            pass

# Running TMSOC
complexity = "NONE"
commandString = 'bash perl.sh'
os.system(commandString)
tmsoc_result = open("TMSOCoutput.txt", 'r')

for line in tmsoc_result:
    if ";simple" in line:
        complexity = "Simple"
        complexity_list.append(complexity)
    elif ";complex" in line:
        complexity = "Complex"
        complexity_list.append(complexity)
    elif ";twilight" in line:
        complexity = "Twilight"
        complexity_list.append(complexity)

#Each residue needs a colour that describes it
#white = NONE 1
#red = complex 2
#blue = simple 3
#purple = twilight 4
#yellow = non-tmh helix
# Residues will be replaced by the corresponding numbers 1 by 1 which will be placed into a file

# load the protein
cmd.load("protA.pdb")

# open the file of new values (just 1 column of numbers, one for each alpha carbon)
inFile = open("newBFactors", 'r')

# create the global, stored array
stored.newB = []

# read the new B factors from file
for line in inFile.readlines(): stored.newB.append( float(line) )

# close the input file
inFile.close()

# clear out the old B Factors
alter protA, b=0.0

# update the B Factors with new properties
alter protA and n. CA, b=stored.newB.pop(0)

# color the protein based on the new B Factors of the alpha carbons
cmd.spectrum("b", "protA and n. CA")
