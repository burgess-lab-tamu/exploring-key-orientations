import sys,os
pdb = sys.argv[1]

#fname = "pdb%s.ent.Z" % pdb.lower()
#os.system("wget ftp://ftp.rcsb.org/pub/pdb/data/structures/all/pdb/%s" % fname)
#os.system("zcat %s > %s.pdb" % (fname,pdb))
#os.system("rm %s" % fname)

fname = "pdb%s.ent.gz" % pdb.lower()
os.system("wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/%s" % fname)
os.system("gunzip -f %s" % fname)
os.system("mv pdb%s.ent %s.pdb" % (pdb.lower(),pdb))
