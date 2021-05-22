rm PDBinfo
rm PDBids
cd  pdb 
rm *.pdb
cd ..
rm PPcomplexes interface_pairs.hetero interface_pairs.hetero.sizes interface_pairs.hetero.sym.sizes updatedPdbs.txt 
cp tempUpdatedPdbs.txt updatedPdbs.txt
cp /home/ioerger/burgess/PDBinfo-2-3-17 PDBinfo
cp /home/ioerger/burgess/PDBids-2-3-17 PDBids
cp /home/ioerger/burgess/PPcomplexes-2-3-17 PPcomplexes
cp /home/ioerger/burgess/interface_pairs.hetero-2-3-17 interface_pairs.hetero
cp /home/ioerger/burgess/interface_pairs.hetero.sizes-2-3-17 interface_pairs.hetero.sizes
cp /home/ioerger/burgess/interface_pairs.hetero.sym.sizes-2-3-17 interface_pairs.hetero.sym.sizes
