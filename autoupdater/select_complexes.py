import sys

for line in open(sys.argv[1]):
  if line[0]=='#': continue
  w = line.rstrip().split('\t')
  id,meth = w[1],w[3]
  if meth!="X-ray": continue
  res = float(w[4])
  if res>3.0: continue
  nchains,ngroups,nprot = int(w[5]),int(w[6]),int(w[7])
  if nprot<2 or nchains>20: continue
  vals = [id,'xray']+w[5:]
  print ' '.join(vals)
