import sys,json,datetime

infile = sys.argv[1]

data = json.load(open(infile))

pdbs = []
for subset in data["subsets"]:
  pdbs += subset["data"]["entries"]

n,c = len(pdbs),0

for i in range(n):
  info = pdbs[i]
  id = info["entry"]["id"]
  meth = info["rcsb_entry_info"]["experimental_method"]
  res = info["rcsb_entry_info"]["resolution_combined"] # could be None [1.3] or [4.2,4.0]
  year = info["rcsb_accession_info"]["initial_release_date"][:4]
  if res==None: res = 0
  else: res = min(res)
  #if meth!="X-ray" or res==0 or res>3.0: continue
  date = str(datetime.datetime.now())
  date = date[:date.index('.')]
  vals = [date,id,year,meth,res]
  chains = info["polymer_entities"]
  #if chains==None: sys.stderr.write("skipping "+id+" because it has no chains\n"); continue
  if chains==None: chains = []
  nchains,nprot = 0,0
  groups = []
  translate = {"Protein":"prot","DNA":"nuc_","RNA":"nuc_"} # to match old format of PPcomplexes
  for chain in chains:
    polytype = chain["entity_poly"]["rcsb_entity_polymer_type"]
    if polytype not in translate: continue # NA-hybrid, Other
    chainids = chain["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"]
    groups.append("%s:%s" % (translate[polytype],''.join(chainids)))
    nchains += len(chainids)
    if polytype=="Protein": nprot += 1
  vals = vals+[nchains,len(groups),nprot]+groups
  print '\t'.join([str(x) for x in vals])
  c += 1
  sys.stderr.write("adding data for %s\n" % id)

sys.stderr.write("processed %d PDB structures\n" % c)



