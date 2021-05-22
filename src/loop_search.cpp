#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "triplet_lib.hpp"
#include "gsl_utils.h" // in bioinfo/

Atom* get_first_atom(vector<Atom*>& atoms,int rnum)
{
  for (int i=0 ; i<atoms.size() ; i++)
    if (rnum==atoms[i]->rnum) return atoms[i];
}

// caller should check that 4 atoms were found (no missing or duplicates)

vector<Atom*> get_bb_atoms(vector<Atom*>& chain,int rnum)
{
  vector<Atom*> atoms;
  // could make this faster with a hash table
  for (int i=0 ; i<chain.size() ; i++)
    if (chain[i]->rnum==rnum && strcmp(chain[i]->anam," N  ")==0) atoms.push_back(chain[i]);
  for (int i=0 ; i<chain.size() ; i++)
    if (chain[i]->rnum==rnum && strcmp(chain[i]->anam," CA ")==0) atoms.push_back(chain[i]); 
  for (int i=0 ; i<chain.size() ; i++)
    if (chain[i]->rnum==rnum && strcmp(chain[i]->anam," C  ")==0) atoms.push_back(chain[i]);
  for (int i=0 ; i<chain.size() ; i++)
    if (chain[i]->rnum==rnum && strcmp(chain[i]->anam," O  ")==0) atoms.push_back(chain[i]);
  return atoms; 
}

Atom* find_atom(vector<Atom*>* atoms,char* anam)
{
  char buf[5];
  sprintf(buf," %-3s",anam);
  for (int i=0 ; i<atoms->size() ; i++) 
    if (strcmp(atoms->at(i)->anam,buf)==0) return atoms->at(i);
  return NULL;
}

#define MAXATOMS 20

Match* compare_atom_sets(vector<Atom*>& a,vector<Atom*>& b)
{
  if (a.size()!=b.size()) { printf("error: atom sets must be same size\n"); exit(0); } // should also be between 3 and MAXATOMS
  int N=a.size();
  Match* match=new Match();
  double co1[MAXATOMS][3],co2[MAXATOMS][3];
  for (int i=0 ; i<N ; i++)
    for (int j=0 ; j<3 ; j++) { co1[i][j] = a[i]->co[j]; co2[i][j] = b[i]->co[j]; }

  for (int i=0 ; i<3 ; i++) match->cen1[i] = match->cen2[i] = 0;
  for (int i=0 ; i<N ; i++)
    for (int j=0 ; j<3 ; j++) { match->cen1[j] += co1[i][j]; match->cen2[j] += co2[i][j]; }
  for (int i=0 ; i<3 ; i++) { match->cen1[i] /= (float)N; match->cen2[i] /= (float)N; }
  for (int i=0 ; i<N ; i++)
    for (int j=0 ; j<3 ; j++) { co1[i][j] -= match->cen1[j]; co2[i][j] -= match->cen2[j]; }

  // assume coords are pre-centered, rotation of first onto second
  match->rms = superpose_rms((double*)co1,(double*)co2,N,match->rot);

/*
  // this score is the max distance
  double co3[MAXATOMS][3];
  for (int i=0 ; i<N ; i++) match->transform(a[i]->co,co3[i]);
  double maxdsq=0;
  for (int i=0 ; i<N ; i++) {
    double d=dist2(b[i]->co,co3[i]);
    if (d>maxdsq) maxdsq = d; }
  match->score = sqrt(maxdsq);
*/
  match->score = match->rms;

  return match;
}

//--------------------------------------------------

class LoopMatch 
{
 public:
  Match* match;
  int rnum,a,b;
  vector<Atom*>* conformer;
  LoopMatch(Match* m,int r,int wina,int winb,vector<Atom*>* conf) { match=m; rnum=r; a=wina; b=winb; conformer=conf; }
  LoopMatch() {} // empty constructor
};

int LMcomparator(LoopMatch& lm1,LoopMatch& lm2) { return lm1.match->score < lm2.match->score; }

//--------------------------------------------------

int main(int argc, char** argv)
{
  if (argc<8) { printf("usage: loop_search <protein.pdb> <chain1> <chain2> <conformers.pdb> <atomlist1> <atomlist2> <loopsize>\n  atomlist is a comma-separated list of atom names in the scaffold, like N5,C4,C3,O2 (where C4 is the C-alpha)\n"); exit(0); }

  char* pdb_fname=argv[1];
  char ch1=argv[2][0];
  char ch2=argv[3][0];
  char* conformers_fname=argv[4];
  istringstream sslignames1(argv[5]); // comma-separated list of 4 atom names, like N5,C4,C3,O2
  istringstream sslignames2(argv[6]); 
  int W=atoi(argv[7]); // loop size

  vector<vector<Atom*>*>* conformers=read_conformers(conformers_fname); // separated by ENDMDL
  printf("REMARK scaffold: models=%d, atoms=%d\n",(int)conformers->size(),(int)conformers->at(0)->size());
  vector<Atom*>* pdb=read_pdb(pdb_fname);

  vector<Atom*> chain1,chain2;
  for (int i=0 ; i<pdb->size() ; i++)
  {
    if (pdb->at(i)->chn==ch1) chain1.push_back(pdb->at(i));
    if (pdb->at(i)->chn==ch2) chain2.push_back(pdb->at(i));
  }

  printf("REMARK complex: %s, natoms(chain=%c)=%d, natoms(chain=%c)=%d\n",argv[2],ch1,(int)chain1.size(),ch2,(int)chain2.size());

  vector<int>* resnumsp=interface_residues(chain1,chain2);
  vector<int> resnums;
  for (int i=0 ; i<resnumsp->size() ; i++) resnums.push_back(resnumsp->at(i));
  printf("REMARK interface resnums (chain %c): ",ch1);
  for (int i=0 ; i<resnums.size() ; i++) printf("%d,",resnums[i]);
  printf("\n");
  printf("REMARK\n");

  vector<string> vlignames;
  string name;
  while (getline(sslignames1,name,',')) { vlignames.push_back(name); }
  while (getline(sslignames2,name,',')) { vlignames.push_back(name); }

  vector<vector<Atom*>* > ligatoms4allconfs;
  for (int k=0 ; k<conformers->size() ; k++) 
  {
    vector<Atom*>* conformer=conformers->at(k);
    vector<Atom*>* ligatoms=new vector<Atom*>(); // 8 atoms: 2 sets of <N,CA,C,O>
    for (int i=0 ; i<vlignames.size() ; i++) 
      ligatoms->push_back(find_atom(conformers->at(0),(char*)vlignames[i].c_str()));
    ligatoms4allconfs.push_back(ligatoms);
  }

  vector<LoopMatch> matches;

  for (int i=0 ; i<resnums.size() ; i++)
  { 
    int rnum=resnums[i];

    for (int j=0 ; j<W ; j++) // include ends of window
    {
      int a=rnum-W+j+1,b=rnum+j;
      vector<Atom*> res1=get_bb_atoms(chain1,a);
      vector<Atom*> res2=get_bb_atoms(chain1,b);
      if (res1.size()!=4) { printf("REMARK not enough backbone atoms, rnum=%d\n",a); continue; }
      if (res2.size()!=4) { printf("REMARK not enough backbone atoms, rnum=%d\n",b); continue; }
      for (int k=0 ; k<res2.size() ; k++) res1.push_back(res2[k]); // use append, add, or insert?

      float bestscore=-1;
      LoopMatch bestconf;
      for (int k=0 ; k<conformers->size() ; k++)
      {
        vector<Atom*>* conformer=conformers->at(k);
        vector<Atom*>* ligatoms=ligatoms4allconfs[k]; // just the 8 fiducials

        Match* m=compare_atom_sets(*ligatoms,res1); 
        if (bestscore==-1 || m->score < bestscore) { 
          bestscore = m->score; 
          bestconf = LoopMatch(m,rnum,a,b,conformer); }
      }
      matches.push_back(bestconf);
    }
  }

  sort(matches.begin(),matches.end(),LMcomparator);

  for (int i=0 ; i<matches.size() ; i++)
  { 
    LoopMatch lm=matches[i];
    Match* m=lm.match;
    int rnum=lm.rnum;
    int a=lm.a,b=lm.b;
    vector<Atom*>* conformer=lm.conformer;

    printf("REMARK best match for pdb %s, chain %c, %s%d, res window %d-%d, rms=%f\n",pdb_fname,ch1,get_first_atom(chain1,rnum)->rnam,rnum,a,b,m->score);
    for (int k=0 ; k<conformer->size() ; k++)
    {
      double temp[3]; // make temp copy of coords so we don't overwrite them
      for (int d=0 ; d<3 ; d++) temp[d] = conformer->at(k)->co[d]; 
      m->transform(conformer->at(k)->co,conformer->at(k)->co);
      conformer->at(k)->print();
      for (int d=0 ; d<3 ; d++) conformer->at(k)->co[d] = temp[d];
    }
    printf("ENDMDL\n");
  }
}
