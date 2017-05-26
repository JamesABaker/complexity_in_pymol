# About

This is a development project. The aim is to roll out [TMSOC](http://tmsoc.bii.a-star.edu.sg/) like function into 3D protein structures. The project has two components. First is a PyMol script for visualising sequence complexity in the 3D structures. Second is a structural database that describes 3D arrangements of membrane protein intra membrane spaces.

## PyMol

### Usage
- Open pymol.
- Navigate to the folder containing `complexity_pymol_plugin.py` in pymol.
- Enter `run complexity_pymol_plugin.py` and enter in pymol.
- load the 3D structure
- Run `complexity` in pymol.

### About
The script uses phobius to estimate TMH boundaries. TMSOC is then used to assess complexity. The TMH regions are coloured accordingly.

## Database
The script takes a list of PDB IDs. It opens or downloads the PDB file. It then checks if the PDBTM database contains TMD boundary information. Phobius is also used as a backup. If phobius fails, and there is no PDBTM entry, then manual co-ordinates can be entered. After that, the TMHs are ran through TMSOC. The PDB file is copied and modified so that the B-factor column reflects the TMH complexity. All the information and logs are then added to a database.

# To do

- [x] TM prediction integration.
- [x] Manual TMH location definitions.
- [x] TMSOC integration.
- [ ] PyMol compatible scipt.
- [ ] PDBTM integration.
- [ ] b-factor PDB generator indicating complexity.
- [ ] Database builder.

[![Code Issues](https://www.quantifiedcode.com/api/v1/project/8a4ca942e31146de8448bb69a75c384f/badge.svg)](https://www.quantifiedcode.com/app/project/8a4ca942e31146de8448bb69a75c384f)

# PDBTM dataset preparation

The PDBTM dataset is not entirely ready to use from the [download page](http://pdbtm.enzim.hu/?_=/download/files). In order to make it compatible with these scripts, as well as other populat XML modules, you must:

1. remove all cases of `<?xml version="1.0"?>`, except the first instance.
2. run `xmlint --format PDBTM_DOWNLOAD.xml > NEW_XML_FILE.xml`
