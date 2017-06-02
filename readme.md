# Status

[![Code Issues](https://www.quantifiedcode.com/api/v1/project/8a4ca942e31146de8448bb69a75c384f/badge.svg)](https://www.quantifiedcode.com/app/project/8a4ca942e31146de8448bb69a75c384f)

- [x] TM prediction integration.
- [x] Manual TMH location definitions.
- [x] TMSOC integration.
- [x] PyMol compatible scipt.
- [ ] PDBTM integration.
- [x] b-factor PDB generator indicating complexity.
- [ ] Database builder.

These scripts are in development and are provided as is. Use at your own risk.

# About

This is a development project. The aim is to roll out [TMSOC](http://tmsoc.bii.a-star.edu.sg/) like function into 3D protein structures. The project has two components. First is a PyMol script for visualising sequence complexity in the 3D structures. Second is a structural database that describes 3D arrangements of membrane protein intramembrane spaces.

Currently, **the software is available for Linux distributions only** since it is dependent on the `decodeanhmm` file underlying [Phobius](http://software.sbc.su.se/cgi-bin/request.cgi?project=phobius).

## PyMol

### About

The script uses Phobius to estimate TMH boundaries. TMSOC is then used to assess complexity. The TMH regions are coloured accordingly:

- non-TMH loops = dark blue
- complex TMH = light blue
- twilight TMH = green
- simple TMH = orange

If a region is in red, Phobius identified a TMH, however, TMSOC experienced an error.

### Installation

This software requires linux, python3, biopython, and Pymol.

### Running the script
#### Overview
Move the pdb file into the complexity_in_pymol/pymol/ directory and then run `bash runme.sh` in a terminal from the complexity_in_pymol/pymol/ folder. In Pymol, navigate to the complexity_in_pymol/pymol/ and after loading the structure enter `run bfactors.py` and then `complexb YOURSTRUCTURENAME`.

#### Step by step instructions
 1. Move the pdb file you wish to work on into the pymol directory. For exampl"maximum=4"ds/complexity_in_pymol/pymol/1a91.pdb`
 2. In a terminal, navigate to the pymol folder. For example:
 `cd Downloads/complexity_in_pymol/pymol`
 3. In a terminal, run the `runme.sh`. This runs a series of commands in order and manages and logs the various input and output files. For example:
 `bash runme.sh`
 4. Enter the pdb filename when prompted. If you experience an error during this process, please log it by creating a [new issue](https://github.com/jbkr/complexity_in_pymol/issues/new) on github.
 5. Load up the Pymol application.
 6. Navigate in the application Pymol to the complexity_in_pymol/pymol directory. For example:
 `cd Downloads/complexity_in_pymol/pymol`
 7. In Pymol load the structure. For example:
 `load 1a91.pdb`
 8. In Pymol, load the bfactors script. Use `run bfactors.py`
 9. In Pymol run the bfactor function on your desired molecule. For example `complexb 1a91`


## Database

The script takes a list of PDB IDs. It opens or downloads the PDB file. It then checks if the PDBTM database contains TMD boundary information. Phobius is also used as a backup. If Phobius fails, and there is no PDBTM entry, then manual co-ordinates can be entered. After that, the TMHs run through TMSOC. The PDB file is copied and modified so that the B-factor column reflects the TMH complexity. All the information and logs are then added to a database.

### PDBTM dataset preparation

The PDBTM dataset is not entirely ready to use from the [download page](http://pdbtm.enzim.hu/?_=/download/files). In order to make it compatible with these scripts, as well as other populat XML modules, you must:

1. remove all cases of `<?xml version="1.0"?>`, except the first instance.
2. run `xmlint --format PDBTM_DOWNLOAD.xml > NEW_XML_FILE.xml`
