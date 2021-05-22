import sys
sys.path.append("/pacific/home/ioerger/genomics")
from utils import *
import json
import urllib3
import urllib

http = urllib3.PoolManager()

search_URI = "https://search.rcsb.org/rcsbsearch/v1/query"
data_URI = "https://data.rcsb.org/graphql"

#    citation { rcsb_journal_abbrev year title rcsb_authors }

#    polymer_entities {
#      rcsb_entity_source_organism { 
#        common_name 
#        scientific_name 
#        rcsb_gene_name { value } } 

template = """{
  entries(entry_ids: %s) {
    entry { id }
    struct { title }
    rcsb_accession_info { initial_release_date }

    rcsb_entry_info {
      experimental_method
      assembly_count
      resolution_combined
    }

    polymer_entities {
      rcsb_polymer_entity_container_identifiers {
        entity_id
        auth_asym_ids
      }
      entity_poly { rcsb_entity_polymer_type rcsb_sample_sequence_length }
      rcsb_polymer_entity { pdbx_description }
    }

    nonpolymer_entities {
      pdbx_entity_nonpoly { name }
      rcsb_nonpolymer_entity { formula_weight }
      rcsb_nonpolymer_entity_container_identifiers { nonpolymer_comp_id } 
    }
  }
}"""

def get_pdb_info(pdbcodes):
  query = template % pdbcodes
  query = query.replace("'",'"')

  r = http.request('POST',data_URI,body=query,headers={"Content-Type": "application/graphql"})
  res_str = r.data.decode('utf-8')
  res_json = json.loads(res_str) # dictionary
  return res_json

###############

# 10000 pdbs takes ~1min and 20Mb

if __name__=="__main__":
  pdbcodes = []
  if "-f" in sys.argv:
    fname = sys.argv[sys.argv.index("-f")+1]
    for line in open(fname):
      w = line.split()
      pdb = w[0]
      pdbcodes.append(pdb)
  else: pdbcodes = sys.argv[1:]

  i,n = 0,len(pdbcodes)
  W = 10000
  print "{ \"subsets\": ["
  first = True
  while i<n:
    if first: first = False
    else: print ","
    sys.stderr.write("working on pdbs %s-%s\n" % (i,i+W))
    subset = pdbcodes[i:i+W]
    res_json = get_pdb_info(subset)
    print json.dumps(res_json,indent=2)
    i += W
  print "] }"
