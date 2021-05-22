#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <time.h>
#include "gsl_utils.h"
#include "getopt.h"
#include "triplet_lib.hpp"

using namespace std;

double pairwise_rmsd(vector<Atom*>& X,vector<Atom*>& Y)
{
  double r=0;
  int n=X.size(); // should be same length as Y
  for (int i=0 ; i<n ; i++)
    for (int j=0 ; j<n ; j++)
      if (i<j)
        r += fabs(dist2(X[i]->co,X[j]->co)-dist2(Y[i]->co,Y[j]->co));
  return sqrt(r/float(n));
}

//---------------------------------------------------------------

// atomlist looks like "6,20,10,25,16,32", should be alternating Ca,Cb order

void usage()
{
  printf("usage: triplet_search [flags] <conformers_pdb> <anum_list> <complex_pdb> <chain1> <chain2>\n");
  exit(-1);
}

float FILTER=4.0;
int show_best=0;

int main(int argc,char** argv)
{
  static struct option long_options[] = {
    {"best", 0, 0, 0},
    {"filter", 1, 0, 0},
    {"Wdist", 1, 0, 0},
    {"Wang", 1, 0, 0},
    {"contact", 1, 0, 0},
    {"pointing", 1, 0, 0},
    {NULL, 0, NULL, 0}
  };
  int option_index=0;
  int c;
  while ((c=getopt_long(argc,argv,"",long_options,&option_index))!=-1)
  {
    int this_option_optind = optind ? optind : 1;
    switch (c) 
    {
      case 0:
        if (strcmp("best",long_options[option_index].name)==0) show_best = 1;
        else if (strcmp("filter",long_options[option_index].name)==0) FILTER = atof(optarg);
        else if (strcmp("Wang",long_options[option_index].name)==0) Wang = atof(optarg);
        else if (strcmp("Wdist",long_options[option_index].name)==0) Wdist = atof(optarg);
        else if (strcmp("contact",long_options[option_index].name)==0) CONTACT = atof(optarg);
        else if (strcmp("pointing",long_options[option_index].name)==0) POINTING = atof(optarg);
        else { printf ("error: invalid flag: -%c\n", c); exit(-1); }
        break;
      default:
        printf ("error: invalid flag: --%s\n",long_options[option_index].name); exit(-1);
    }
  }

  if (argc-optind<5) usage();
  char* ligfile=argv[optind];
  char* anumlist=argv[optind+1];

  vector<vector<Atom*>*>* conformers=read_conformers(ligfile);
  vector<char*>* remarks=read_frames(ligfile);
  if (conformers->size()!=remarks->size()) { printf("error: unequal number of frames REMARKS as conformers\n"); exit(-1); }
  vector<int> frames;
  for (int i=0 ; i<remarks->size() ; i++)
  {
    vector<char*>* words=split(remarks->at(i),' ');
    frames.push_back(atoi(words->at(2)));
  }

  char* pdbfile=argv[optind+2];
  char ch1=argv[optind+3][0];
  char ch2=argv[optind+4][0];
  vector<Atom*>* pdb=read_pdb(pdbfile);

  time_t curtime;
  struct tm *loctime;
  curtime = time (NULL);
  loctime = localtime (&curtime);


  printf("REMARK # command: "); for (int i=0 ; i<argc ; i++) printf("%s ",argv[i]); printf("\n");
  //char* version=split(version_string,' ')->at(2);
  //printf("REMARK # version: %s\n",version);
  printf("REMARK # run date: "); fputs (asctime (loctime), stdout);
  printf("REMARK # params: filter=%0.1f, contact=%0.1f, pointing=%0.1f, Wang=%.1f, Wdist=%0.1f\n",FILTER,CONTACT,POINTING,Wang,Wdist);

  // get anums of 6 ligand Ca and Cb atoms
  vector<char*>* words=split(anumlist,',');
  vector<int> liganums;
  for (int i=0 ; i<words->size() ; i++)
  {
    int a;
    sscanf(words->at(i),"%d",&a);
    liganums.push_back(a);
  }

  vector<Atom*> chain1,chain2;
  for (int i=0 ; i<pdb->size() ; i++)
    if (pdb->at(i)->chn==ch1) chain1.push_back(pdb->at(i));
  for (int i=0 ; i<pdb->size() ; i++)
  {
    if (strcmp(argv[optind+4],"NA")==0) { if (nucleic(pdb->at(i))==1) chain2.push_back(pdb->at(i)); }
    else { if (pdb->at(i)->chn==ch2) chain2.push_back(pdb->at(i)); }
  }
  vector<int>* resnumsp=interface_residues(chain1,chain2);
  vector<int> resnums;
  for (int i=0 ; i<resnumsp->size() ; i++) resnums.push_back(resnumsp->at(i));
  printf("REMARK # interface resnums: ");
  for (int i=0 ; i<resnums.size() ; i++) printf("%d,",resnums[i]);
  printf("\n");
  printf("REMARK\n");


  unordered_map<int,Atom*> Calphas;
  unordered_map<int,Atom*> Cbetas;
  for (int i=0 ; i<pdb->size() ; i++)
  {
    Atom* atom=pdb->at(i);
    if (atom->chn==ch1)
    {
      if (strcmp(atom->anam," CA ")==0) Calphas[atom->rnum] = atom;
      if (strcmp(atom->anam," CB ")==0) Cbetas[atom->rnum] = atom;
    }
  }

  // begin search (for each conformer...)

  vector<Match*> hits;
  for (int m=0 ; m<conformers->size() ; m++)
  {
    vector<Atom*>* conformer=conformers->at(m);
    vector<Atom*> ligatoms; // 6 Ca and Cb atoms, alternating
    
    for (int j=0 ; j<6 ; j++)
      for (int i=0 ; i<conformer->size() ; i++)
        if (conformer->at(i)->anum==liganums[j])
          ligatoms.push_back(conformer->at(i));
    if (ligatoms.size()!=6) { printf("error: not all 6 lig atoms found\n"); exit(-1); }

    // loop through possible triplets
    Match* best=NULL;
    for (int i=0 ; i<resnums.size() ; i++)
    {
      if (Cbetas.find(resnums[i])==Cbetas.end()) continue;
      for (int j=0 ; j<resnums.size() ; j++)
      {
        if (Cbetas.find(resnums[j])==Cbetas.end()) continue;
        float d12=sqrt(dist2(Calphas[resnums[i]]->co,Calphas[resnums[j]]->co));
        float e12=sqrt(dist2(ligatoms[0]->co,ligatoms[2]->co));
        if (fabs(d12-e12)>FILTER) continue;

        for (int k=0 ; k<resnums.size() ; k++)
        {
         vector<Atom*> X;
         vector<Atom*> Y;
         X.push_back(Calphas[resnums[i]]); X.push_back(Calphas[resnums[j]]); X.push_back(Calphas[resnums[k]]);
         Y.push_back(ligatoms[0]); Y.push_back(ligatoms[2]); Y.push_back(ligatoms[4]);
         double r=pairwise_rmsd(X,Y);
         if (r<FILTER)
         {
          if (Cbetas.find(resnums[k])==Cbetas.end()) continue;
          if (i==j || i==k || j==k) continue;
          int a=resnums[i],b=resnums[j],c=resnums[k];
          vector<Atom*> triplet;
          triplet.push_back(Calphas[a]);
          triplet.push_back(Cbetas[a]);
          triplet.push_back(Calphas[b]);
          triplet.push_back(Cbetas[b]);
          triplet.push_back(Calphas[c]);
          triplet.push_back(Cbetas[c]);
          Match* match=compare_3vecs(ligatoms,triplet);
          match->i=i; match->j=j; match->k=k;
          if (best==NULL || match->score<best->score) best = match;
         }
        }
      }
    }
    if (best==NULL) continue;
    best->conf = m;
    hits.push_back(best);
  }

  if (hits.size()==0) { printf("no matches found\n"); exit(0); }
  if (show_best==1) sort(hits.begin(),hits.end(),match_compare);
  printf("num hits = %ld\n",hits.size());

  for (int j=0 ; j<hits.size() ; j++)
  {
    Match* best=hits[j];
    int m=best->conf; // index of conformer
    vector<Atom*>* conformer=conformers->at(m);

    printf("%s",remarks->at(m));
    //printf("REMARK conformer %d: best match rnums=%d,%d,%d, rms=%0.2f, angle=%0.2f, score=%0.2f\n",m,resnums[best->i],resnums[best->j],resnums[best->k],best->rms,best->angle,best->score);
    printf("REMARK match: chain=%c resids=(%d, %d, %d)\n",ch1,resnums[best->i],resnums[best->j],resnums[best->k]);
    printf("REMARK index=%d, score=%0.2f, rms=%0.2f, angle=%0.2f\n",m,best->score,best->rms,best->angle);
    printf("REMARK EKO pdb: %s chain: %c origframe: %d residues: %d %d %d rms: %f ang: %f score: %f\n",pdbfile,ch1,frames[m],resnums[best->i],resnums[best->j],resnums[best->k],best->rms,best->angle,best->score);
    for (int n=0 ; n<conformer->size() ; n++)
    {
      double co1[3],co2[3];
      Atom* atom=conformer->at(n);
      for (int i=0 ; i<3 ; i++) co1[i] = atom->co[i];
      best->transform(co1,co2);
      for (int i=0 ; i<3 ; i++) atom->co[i] = co2[i];
      atom->print();
      for (int i=0 ; i<3 ; i++) atom->co[i] = co1[i];
    }
    printf("END\n");
    if (show_best==1) break; 
  }
}

