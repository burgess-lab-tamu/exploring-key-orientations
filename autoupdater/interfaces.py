import sys,os
sys.path.append(os.path.dirname(os.path.realpath(__file__))+"/../..")
from PPI.src.utils import *
from PPI.src.pdb import *

# this does not handle alternative conformations
# return angle (or -1 if not defined)

PDBdir = "/home/ioerger/burgess/PP/"
#PDBdir = "/ocean/oahu/pdb/all/"

def side_chain_angle_to_contact(atom1,prot1,atom2):
  CA,CB = None,None
  for a in prot1:  
    if a.rnum==atom1.rnum and a.chn==atom1.chn:
      if a.anam==" CA ": CA = a
      if a.anam==" CB ": CB = a
  if CA==None or CB==None: return -1
  ang = angle3(CB.co,CA.co,atom2.co)
  return ang

def contact_residues(A,B,cutoff=4.0,pointing=45.0):
  tree = KDTree([b.co for b in B])
  atoms = []  
  for a in A:
    #if a.anam in [" C  "," N  "," O  "," CA "]: continue
    d,co = tree.nearest(a.co)
    if d<cutoff*cutoff:
      b = tree.index_of(co)
      ang = side_chain_angle_to_contact(a,A,B[b])
      if ang!=-1 and ang<pointing: atoms.append(a)
  return atoms

# previous version returned rnums...
  residues = [x.rnum for x in atoms]
  residues = list(set(residues)) # remove duplicates
  residues.sort()
  return residues

def show_chains(complexes_file):
  cnt = 0
  for line in open(complexes_file):
    #if cnt==6: break 
    cnt += 1
    w = line.rstrip().split()
    fname = PDBdir+"pdb%s.ent" % (w[0].lower())
    prot = read_pdb(fname)
    allatoms = []
    for ch in prot: allatoms += ch
    chaintypes = w[5:]
    print '#',line,

    for i in range(len(chaintypes)):
      ch = chaintypes[i]
      for c in ch[5:]:
        print w[0],c,ch[:4],len(filter(lambda arec: arec.chn==c,allatoms))

def unique_rnums(atoms):
  residues = [x.rnum for x in atoms]
  residues = list(set(residues)) # remove duplicates
  residues.sort()
  return residues

def id(atom):
  return "%s%s%s" % (atom.rnam,atom.rnum,atom.chn)

def unique_ids(atoms):
  residues = [id(x) for x in atoms]
  residues = list(set(residues)) # remove duplicates
  residues.sort()
  return residues

# evaluate interface of each pair of chains

def hetero_complexes(complexes_files):
 cnt = 0
 for line in open(complexes_file):
  #if cnt==6: break 
  cnt += 1
  w = line.rstrip().split()
  fname = PDBdir+"pdb%s.ent" % (w[0].lower())
  prot = read_pdb(fname)
  allatoms = []
  for ch in prot: allatoms += ch
  chaintypes = w[5:]
  print "#",line,

  for i in range(len(chaintypes)):
    for j in range(len(chaintypes)):
      if i<>j and chaintypes[i][:4]=="prot" and chaintypes[j][:4]=="prot":
        chains1 = chaintypes[i][5:]
        chains2 = chaintypes[j][5:]
        print w[0],chains1[0],chains2[0] ; continue #########
        
        for a in chains1:
          for b in chains2:
            atoms1 = filter(lambda arec: arec.chn==a,allatoms)
            atoms2 = filter(lambda arec: arec.chn==b,allatoms)
            if len(atoms2)==0: continue
            c = contact_residues(atoms1,atoms2,pointing=90)
            c = unique_rnums(c)
            #print a,len(atoms1),b,len(atoms2),len(c)
            print w[0],a,b,
            for rnum in c: print rnum,
            print

# don't actually read PDBs and extract interface residues

def hetero_complexes_fast(complexes_files):
 cnt = 0
 for line in open(complexes_file):
  #if cnt==6: break 
  cnt += 1
  w = line.rstrip().split()
  fname = PDBdir+"pdb%s.ent" % (w[0].lower())
  chaintypes = w[5:]
  print "#",line,

  if len(chaintypes)>20: print "# %s: skip - too many chains" % w[0]; continue
  for i in range(len(chaintypes)):
    for j in range(len(chaintypes)):
      if i<>j and chaintypes[i][:4]=="prot" and chaintypes[j][:4]=="prot":
        #chains1 = chaintypes[i][5:]
        #chains2 = chaintypes[j][5:]
        #print w[0],chains1[0],chains2[0]
        for chains1 in chaintypes[i][5:]:
          for chains2 in chaintypes[j][5:]:
            print w[0],chains1,chains2

def NA_complexes(complexes_file):
 cnt = 0
 for line in open(complexes_file):
  #if cnt==6: break 
  cnt += 1
  w = line.rstrip().split()
  fname = PDBdir+"pdb%s.ent" % (w[0].lower())
  prot = read_pdb(fname)
  allatoms = []
  for ch in prot: allatoms += ch
  chaintypes = w[5:]
  print "#",line,

  NAatoms = []
  for i in range(len(chaintypes)):
    if chaintypes[i][:4]=="nuc_":
      for a in chaintypes[i][5:]: NAatoms += filter(lambda arec: arec.chn==a,allatoms)
  if len(NAatoms)==0: continue

  for i in range(len(chaintypes)):
    if chaintypes[i][:4]=="prot":
      for a in chaintypes[i][5:]:
        PRatoms = filter(lambda arec: arec.chn==a,allatoms)
        c = contact_residues(PRatoms,NAatoms,pointing=90)
        c = unique_rnums(c)
        print w[0],a,
        for rnum in c: print rnum,
        print

def NA_complexes_fast(complexes_file):
 cnt = 0
 for line in open(complexes_file):
  #if cnt==6: break 
  cnt += 1
  w = line.rstrip().split()
  fname = PDBdir+"pdb%s.ent" % (w[0].lower())
  #prot = read_pdb(fname)
  #allatoms = []
  #for ch in prot: allatoms += ch
  chaintypes = w[5:]
  print "#",line,

  if len(chaintypes)>20: print "# %s: skip - too many chains" % w[0]; continue
  for i in range(len(chaintypes)):
    if chaintypes[i][:4]=="prot":
      for a in chaintypes[i][5:]:
        print w[0],a,"NA"

if __name__=="__main__":
  complexes_file = sys.argv[1] # PPcomplexes or PNAcomplexes
  if sys.argv[2]=='show': show_chains(complexes_file)
  elif sys.argv[2]=='hetero': hetero_complexes(complexes_file)
  elif sys.argv[2]=='hetero_fast': hetero_complexes_fast(complexes_file)
  elif sys.argv[2]=='NA': NA_complexes(complexes_file)
  elif sys.argv[2]=='NA_fast': NA_complexes_fast(complexes_file)
