# returns a list of chains, which are a list of ATOM records

def read_pdb(filename,hetatms=False,hydrogens=False):
  file = open(filename,"r")
  lines = file.readlines()
  file.close()
  return parse_pdb(lines,hetatms)

def parse_pdb(lines,hetatms=False,hydrogens=False):
  chains = []
  atoms = []
  for line in lines:
    if line.find("ATOM")==0:
      if line[13]!='H' or hydrogens==True:
        if line[16]==' ' or line[16]=='A': # only first atom if multiple altLoc
          atoms.append(atom_rec(line))
    if line.find("HETATM")==0 and hetatms==True: atoms.append(atom_rec(line))
    if line.find("TER")==0 or line.find("END")==0: 
      chains.append(atoms)
      atoms = []
  if len(atoms)>0: chains.append(atoms)
  return chains

# parse ATOM records into (anum,anam,rnum,rnam,chn,x,y,z,occ,bf)

class atom_rec:
  def __init__(self,line):
    self.line = line
    self.anum = int(line[6:11])
    self.anam = line[12:16]
    self.rnam = line[17:20].strip()
    self.chn = line[21:22]
    self.rnum = int(line[22:26])
    self.x = float(line[30:38])
    self.y = float(line[38:46])
    self.z = float(line[46:54])
    self.co = [self.x,self.y,self.z]
    self.occ = float(line[54:60])
    self.bf = float(line[60:66])
  def toString(self):
    return "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" % \
      (self.anum,self.anam,self.rnam,self.chn,self.rnum%10000, \
       self.co[0],self.co[1],self.co[2],self.occ,self.bf)
  def copy(self): return atom_rec(self.toString())

if __name__=="__main__":
  test = "4HHB.pdb"
  pdb = read_pdb(test)
  print "num_chains("+test+")="+str(len(pdb))
