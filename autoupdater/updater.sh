#!/bin/bash

make interface_sizes
mv PDBinfo oldPDBInfo
mv PDBids oldPDBIds
wget 'http://www.rcsb.org/pdb/rest/getEntityInfo' -O PDBinfo
wget 'http://www.rcsb.org/pdb/rest/getCurrent' -O PDBids
python extractDiffPDBInfo.py oldPDBIds PDBinfo updatedPdbs.txt > tempPDBInfo  # temp File - tempPDBInfo
python select_PDBs.py tempPDBInfo | tee >> PPcomplexes tempPPComplexes           # temp File - tempPPComplexes
python download.py tempPPComplexes /apps/burgess/pdb
python interfaces.py tempPPComplexes hetero_fast | tee >> interface_pairs.hetero tempInterfacePairsHetero
./interface_sizes tempInterfacePairsHetero | tee >> interface_pairs.hetero.sizes  tempInterfacePairsHeteroSizes
python symmetry_reduction.py tempInterfacePairsHeteroSizes | tee >> interface_pairs.hetero.sym.sizes  tempInterfacePairsHeteroSymSizes
rm oldPDBInfo oldPDBIds tempPDBInfo tempPPComplexes tempInterfacePairsHetero tempInterfacePairsHeteroSizes tempInterfacePairsHeteroSymSizes
# Create a database file.
