import sys,os

dir = sys.argv[2] # target director

i = 1
for line in open(sys.argv[1]): # input file with pdb ids at front of each line
  w = line.split()
  pdb = w[0]
  fname = "pdb%s.ent.gz" % pdb.lower()
  outfile = "%s.pdb" % pdb
  # if not os.path.exists("%s/%s" (dir,outfile):
  print "downloading %s (%s)" % (outfile,i)
  os.system("wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/%s" % fname)
  os.system("gunzip -f %s" % fname)
  os.system("mv pdb%s.ent %s/%s" % (pdb.lower(),dir,outfile))
  i += 1
  
