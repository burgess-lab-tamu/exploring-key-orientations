import sys,json

known = {}
for line in open(sys.argv[2]): # PDBinfo.txt
  if line[0]=='#': continue
  w = line.rstrip().split('\t')
  id = w[1]
  known[id] = 1

allpdbs = json.load(open(sys.argv[1])) # PDBids.json
allpdbs = [a.encode('ascii') for a in allpdbs]

for id in allpdbs:
  if id not in known: print id
