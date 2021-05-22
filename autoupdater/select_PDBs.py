import sys
import xml.dom.minidom
#from xml.dom.minidom import parse, Node

#doc = xml.dom.minidom.parse("PDBinfo-12-20-12")
#doc = xml.dom.minidom.parse("temp.info")
doc = xml.dom.minidom.parse(sys.argv[1])

tot,xray = 0,0
for node in doc.getElementsByTagName("PDB"):
  tot += 1
  #fam = node1.getElementsByTagName("Family")[0].childNodes[0].data
  id = node.getAttribute("structureId")
  method = node.getElementsByTagName("Method")[0].getAttribute("name")
  entities = node.getElementsByTagName("Entity")
  subunits = []
  badchains = False
  for entity in entities:
    typ = str(entity.getAttribute('type'))
    chains = entity.getElementsByTagName("Chain")
    chainids = [str(chain.getAttribute('id')) for chain in chains]
    # exclude files with 2-letter chain ids
    for c in chainids:
      if len(c)!=1: 
	badchains = True
	break
    if badchains==True:
	break
    subunits.append((typ,chainids))
  if badchains==True:
    print "# skipping %s because some chains have ids with >1 letter" % id
    continue
  flag = ''
  if method=='xray': 
    xray += 1
    prot_subunits = 0
    for subunit in subunits:
      if subunit[0]=='protein': prot_subunits += 1
    if prot_subunits>=2: # for hetero complexes
      nchains = len(node.getElementsByTagName('Chain'))
      print id,method,nchains,len(entities),prot_subunits,
      for subunit in subunits:
        s = "nuc_:"
        if subunit[0]=="protein": s = "prot:"
        s += ''.join(subunit[1])
        print s,
      print

