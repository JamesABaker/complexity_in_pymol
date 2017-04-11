# About

This is a development project. The aim is to roll out [TMSOC](http://tmsoc.bii.a-star.edu.sg/) like function into 3D protein structures.

# To do

- [x] TM prediction integration.
- [x] Manual TMH location definitions.
- [x] TMSOC integration.
- [ ] PDBTM integration.
- [ ] b-factor PDB generator indicating complexity.
- [ ] Database builder.

[![Code Issues](https://www.quantifiedcode.com/api/v1/project/8a4ca942e31146de8448bb69a75c384f/badge.svg)](https://www.quantifiedcode.com/app/project/8a4ca942e31146de8448bb69a75c384f)

# PDBTM dataset preparation

The PDBTM dataset is not entirely ready to use from the [download page](http://pdbtm.enzim.hu/?_=/download/files). In order to make it compatible with these scripts, as well as other populat XML modules, you must:

1. remove all cases of `<?xml version="1.0"?>`, except the first instance.
2. run `xmlint --format PDBTM_DOWNLOAD.xml > NEW_XML_FILE.xml`
