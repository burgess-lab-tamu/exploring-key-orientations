from pdb import *
from KDTree import *
from utils import *
from superpose_svd import *

# return mean abs diff of dist between each pair

def pairwise_rmsd(X,Y):
  if len(X)!=len(Y): error("unequal size lists in pairwise_rmsd()")
  s,n = 0,0
  for i in range(len(X)):
    for j in range(len(Y)):
      if i<j: s += abs(distsq(X[i],X[j])-distsq(Y[i],Y[j])); n += 1
  return math.sqrt(s/float(n))

def rmsd(X,Y):
  s = 0
  for i in range(len(X)):
    s += distsq(X[i],Y[i])
  return math.sqrt(s/float(len(X)))

def get_chain(chains,chainid):
  for chain in chains:
    if len(chain)>=1 and chain[0].chn==chainid: return chain
  return None

AAs = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO ASN ARG SER THR VAL TRP".split()
NAs = "DA DG DC DT A G C T U".split()

def get_nucleic_acids(chains):
  atoms = []
  for chain in chains:
    for atom in chain:
      if atom.rnam in NAs: atoms.append(atom)
  return atoms

def select_atoms(pdbatoms,ids):
  dict = {}
  for a in pdbatoms: dict[str(a.anum)] = a; dict[a.anam.strip()] = a # nums and names should be distinct
  return [dict[i] for i in ids]

# this does not handle alternative conformations
# return angle (or -1 if not defined)

def side_chain_angle_to_contact(atom1,prot1,atom2):
  CA,CB = None,None
  for a in prot1:  
    if a.rnum==atom1.rnum and a.chn==atom1.chn:
      if a.anam==" CA ": CA = a
      if a.anam==" CB ": CB = a
  if CA==None or CB==None: return -1
  ang = angle3(CB.co,CA.co,atom2.co)
  return ang

# given 2 lists of pdbatoms
#   find those in first that are within 4 A of an atom in second
# return list of residues

# filter out those residues whose CA-CB vector does not point within
#   45 degrees of contact point

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
  residues = [x.rnum for x in atoms]
  residues = list(set(residues)) # remove duplicates
  residues.sort()
  return residues

# in degrees

def angle4(A,B,C,D):
  E = diff(B,A) # E=B-A
  F = diff(D,C) # F=D-C
  return angle(E,F) # defined in bioinfo/utils.py; in degrees

# each list contains 4 coords
# return weighted combination of:
#   diff between CA dists
#   diff between CA-CB angles

def compare_2vecs(coordsA,coordsB,dist_weight,angle_weight):
  a = dist(coordsA[0],coordsA[2])
  b = dist(coordsB[0],coordsB[2])
  c = angle4(coordsA[0],coordsA[1],coordsA[2],coordsA[3])
  d = angle4(coordsB[0],coordsB[1],coordsB[2],coordsB[3])
  return (dist_weight*abs(a-b)+angle_weight*abs(c-d))

def center(coords):
  center1 = []
  for i in range(3): center1.append(sum([x[i] for x in coords])/float(len(coords)))
  return center1

def compare_3vecs_orig(coordsA,coordsB,dist_weight,angle_weight):
  CAs1 = [coordsA[x] for x in [0,2,4]]
  CAs2 = [coordsB[x] for x in [0,2,4]]
  dst = pairwise_rmsd(CAs1,CAs2)
  angs = 0
  for i in range(3):
    for j in range(3):
      if i<j: # 0,1 ; 0,2 ; 1,2
        ang1 = angle4(coordsA[2*i],coordsA[2*i+1],coordsA[2*j],coordsA[2*j+1])
        ang2 = angle4(coordsB[2*i],coordsB[2*i+1],coordsB[2*j],coordsB[2*j+1])
        angs += abs(ang1-ang2)
  ang = angs/3.0
  return (dist_weight*dst+angle_weight*ang)

# angle is in degrees

def compare_3vecs_no_transform(coordsA,transformed,dist_weight,angle_weight):
  #dst = pairwise_rmsd(coordsA,transformed) # this compares internal pairs
  dst = rmsd(coordsA,transformed) # standard 6-atom RMSD
  ang = 0 # in degrees
  for i in range(3):
    ang += angle4(coordsA[2*i],coordsA[2*i+1],transformed[2*i],transformed[2*i+1])
  ang /= 3.0
  score = dist_weight*dst+angle_weight*ang
  return (score,dst,ang)

def compare_3vecs(coordsA,coordsB,dist_weight,angle_weight):
  transformed = superpose_svd(coordsA,coordsB)
  return compare_3vecs_no_transform(coordsA,transformed,dist_weight,angle_weight)

# angle is in degrees

def compare_2vecs_no_transform(coordsA,transformed,dist_weight,angle_weight):
  dst = rmsd(coordsA,transformed) # standard 4-atom RMSD
  ang = 0 # in degrees
  for i in range(2): ang += angle4(coordsA[2*i],coordsA[2*i+1],transformed[2*i],transformed[2*i+1])
  ang /= 2.0
  score = dist_weight*dst+angle_weight*ang
  return (score,dst,ang)

def compare_2vecs(coordsA,coordsB,dist_weight,angle_weight):
  transformed = superpose_svd(coordsA,coordsB)
  return compare_2vecs_no_transform(coordsA,transformed,dist_weight,angle_weight)

# returns a list of (score,i,j,k,transformed lig coords) sorted by score
#   where i,j,k and resnums of matching residues in prot

def search_3ab(prot,res,lig,patterncoords,dist_weight,ang_weight,cutoff,measure=False):
  # make dictionaries to look up CA/CB atoms by resnum
  CAs,CBs = {},{}
  for a in prot:
    if a.anam==" CA ": CAs[a.rnum] = a
    if a.anam==" CB ": CBs[a.rnum] = a

  a,b = 0,0
  hits = []
  if measure==True:
    i,j,k = res[0],res[1],res[2]
    protcoords = [CAs[i].co,CBs[i].co,CAs[j].co,CBs[j].co,CAs[k].co,CBs[k].co]
    (score,dst,ang) = compare_3vecs_no_transform(protcoords,patterncoords,dist_weight,ang_weight)
    hits.append((score,dst,ang,i,j,k))

  else:
    for i in res:
      for j in res:
        if i==j: continue
        d = abs(dist(CAs[i].co,CAs[j].co)-dist(patterncoords[0],patterncoords[2]))
        if d<cutoff:
          for k in res:
            if k==i or k==j: continue
            b += 1
            protCAs = [CAs[i].co,CAs[j].co,CAs[k].co]
            ligCAs = [patterncoords[0],patterncoords[2],patterncoords[4]]
            d = pairwise_rmsd(protCAs,ligCAs)
            if d<cutoff:
              protcoords = [CAs[i].co,CBs[i].co,CAs[j].co,CBs[j].co,CAs[k].co,CBs[k].co]
              #(score,dst,ang) = compare_3vecs_orig(protcoords,patterncoords,dist_weight,ang_weight),0,0
              (score,dst,ang) = compare_3vecs(protcoords,patterncoords,dist_weight,ang_weight)
              hits.append((score,dst,ang,i,j,k))

  hits.sort()
  results = []
  for (score,dst,ang,i,j,k) in hits[:3]: # top 3?
    coords = [CAs[i].co,CBs[i].co,CAs[j].co,CBs[j].co,CAs[k].co,CBs[k].co]
    (U,c1,c2) = get_transform(coords,patterncoords)  
    ligcoords = [x.co for x in lig]
    ligcoords = apply_transform(ligcoords,U,c1,c2)
    results.append((score,dst,ang,i,j,k,ligcoords))
  n = len(res)
  if measure==False: sys.stderr.write("REMARK residue combinations searched: %d (filtered) %d (total)\n" % (b,(n*(n-1)*(n-2))))
  return results

def search_2ab(prot,res,lig,patterncoords,dist_weight,ang_weight,cutoff,measure=False):
  # make dictionaries to look up CA/CB atoms by resnum
  CAs,CBs = {},{}
  for a in prot:
    if a.anam==" CA ": CAs[a.rnum] = a
    if a.anam==" CB ": CBs[a.rnum] = a

  a,b = 0,0
  hits = []
  if measure==True:
    i,j = res[0],res[1]
    protcoords = [CAs[i].co,CBs[i].co,CAs[j].co,CBs[j].co]
    (score,dst,ang) = compare_2vecs_no_transform(protcoords,patterncoords,dist_weight,ang_weight)
    hits.append((score,dst,ang,i,j))

  else:
    for i in res:
      for j in res:
        if i==j: continue
        protCAs = [CAs[i].co,CAs[j].co]
        ligCAs = [patterncoords[0],patterncoords[2]]
        d = pairwise_rmsd(protCAs,ligCAs)
        if d<cutoff:
          protcoords = [CAs[i].co,CBs[i].co,CAs[j].co,CBs[j].co]
          (score,dst,ang) = compare_2vecs(protcoords,patterncoords,dist_weight,ang_weight)
          hits.append((score,dst,ang,i,j))

  hits.sort()
  results = []
  for (score,dst,ang,i,j) in hits[:3]:
    coords = [CAs[i].co,CBs[i].co,CAs[j].co,CBs[j].co]
    (U,c1,c2) = get_transform(coords,patterncoords)  
    ligcoords = [x.co for x in lig]
    ligcoords = apply_transform(ligcoords,U,c1,c2)
    results.append((score,dst,ang,i,j,ligcoords))
  n = len(res)
  if measure==False: sys.stderr.write("REMARK residue combinations searched: %d (filtered) %d (total)\n" % (b,(n*(n-1))))
  return results

def read_pdb_models(fname):
  models,remarks = [],[]
  lines = open(fname).readlines()
  lines = [x.strip() for x in lines]
  ends = []
  for i in range(len(lines)):
    if lines[i]=="ENDMDL": ends.append(i)
  if lines[-1]!="ENDMDL": ends.append(len(lines))
  for i in range(len(ends)):
    a,b = 0,ends[i]
    if i!=0: a = ends[i-1]
    rems,atoms = [],[]
    for j in range(a,b):
      if lines[j][:6]=="ATOM  ": atoms.append(atom_rec(lines[j]))
      if lines[j][:6]=="REMARK": rems.append(lines[j])
    models.append(atoms)
    remarks.append(rems)
  return models,remarks

