import sys

pdbs = []
pairs = {}

for line in open(sys.argv[1]):
  if "prot" in line:
    w = line.split()
    pdbs.append(w)
  if line[0]!='#':
    w = line.split()
    id = "%s-%s-%s" % (w[0],w[1],w[2])
    pairs[id] = int(w[3])

for w in pdbs:
  pdb = w[1]
  groups = filter(lambda x: x[:5]=="prot:",w)
  groups = [x[5:] for x in groups]
  print ' '.join(w)
  n = len(groups)
  for i in range(n):
    for j in range(n):
      if i!=j: # include both directions, A->B and B->A
        A,B = groups[i],groups[j]
        bestsize,bestA,bestB,bestid = -1,-1,-1,-1
        for ch1 in A:
          for ch2 in B:
            id = "%s-%s-%s" % (pdb,ch1,ch2)
            size = pairs.get(id,-1)
            print "# groups=%s,%s" % (i,j),id,size
            if size>bestsize: bestsize,bestA,bestB,bestid = size,ch1,ch2,id
        if bestsize!=-1: 
          print "# best(%s,%s):" % (bestA,bestB),bestid,bestsize
          print "%s %s %s %s" % (pdb,bestA,bestB,bestsize)

