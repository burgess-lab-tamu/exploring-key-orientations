#!/bin/tcsh

# should I backup PDBinfo.txt and PPcomplexes first?

echo retrieving list of all PDB ids in RCSB

wget https://data.rcsb.org/rest/v1/holdings/current/entry_ids -O PDBids.json

python identify_missing_pdbs.py PDBids.json PDBinfo.txt >! missing_pdbs.txt

echo retrieving metadata attributes for missing PDBs

python get_pdb_info.py -f missing_pdbs.txt >! PDBinfo.json

echo appending info on new PDBs to PDBinfo.txt

python process_pdb_info.py PDBinfo.json >> PDBinfo.txt

echo reselecting PPcomplexes.txt

python select_complexes.py PDBinfo.txt >! PPcomplexes.txt

