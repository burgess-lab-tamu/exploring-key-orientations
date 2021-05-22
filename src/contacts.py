import sys,os,math
from pdb import *
from KDTree import *

pdb = sys.argv[1]
ch1,ch2 = sys.argv[2],sys.argv[3]

chains = read_pdb(pdb)
hashchains = {}
for chain in chains: 
  if len(chain)>0:
    hashchains[chain[0].chn] = chain # based on chn of first Atom
A,B = hashchains[ch1],hashchains[ch2] # error if not found

tree = KDTree([b.co for b in B])

for a in A:
  d,co = tree.nearest(a.co)
  b = tree.index_of(co)
  cutoff = 4.0
  if d<cutoff*cutoff:
    b = B[b]
    print a.rnam,a.rnum,a.chn,a.anam,b.rnam,b.rnum,b.chn,b.anam
