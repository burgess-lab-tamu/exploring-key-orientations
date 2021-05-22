#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <set>
#include <unordered_map>
#include <time.h>
#include "gsl_utils.h"
#include "getopt.h"
#include "triplet_lib.hpp"

using namespace std;

float FILTER=2.0;
int show_best=0;

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

class Triplet
{
public:
  char* pdb;
  char chn;
  vector<char*> resid;
  vector<Atom*> atoms;
  float com[3];

  Triplet(char* buf)
  {
    vector<char*>* w=split(buf,' '); // free these to save memory?
    pdb = w->at(0);
    chn = w->at(1)[0];
    for (int i=0 ; i<3 ; i++) resid.push_back(w->at(3+i));
    for (int i=0 ; i<3 ; i++) com[i] = atof(w->at(6+i));
    for (int i=0 ; i<6 ; i++)
    {
      float x=atof(w->at(9+i*3)),y=atof(w->at(10+i*3)),z=atof(w->at(11+i*3));
      Atom* a=new Atom(x,y,z);
      atoms.push_back(a);
    } 
  }
  void print() { printf("%s %c %s %s %s\n",pdb,chn,resid[0],resid[1],resid[2]); }
};

// calculate pairwise distance among 3 Calphas for each triplet
// sort them so shortest is first and long is last
// return false (do not pass filter) if either the
//   shortest side of each triangle differs by more than 2 A (dist^2<4.0)
//   or the longest sides differ by more than 2 A
// pass in sorted list of triangle side lengths for ligatoms

bool filter(vector<Atom*> B,vector<float>& d)
{
  vector<float>e;
  e.push_back(sqrt(dist2(B[0]->co,B[2]->co)));
  e.push_back(sqrt(dist2(B[2]->co,B[4]->co)));
  e.push_back(sqrt(dist2(B[0]->co,B[4]->co)));
  sort(e.begin(),e.end());
  //printf("d=(%0.2f,%0.2f,%0.2f) e=(%0.2f,%0.2f,%0.2f)\n",d[0],d[1],d[2],e[0],e[1],e[2]);
  for (int i=0 ; i<3 ; i++) if (fabs(d[i]-e[i])>FILTER) return false;
  return true;
}

bool filter2(char* buf,vector<float>& d)
{
  double a[3];
  double b[3];
  double c[3];
  int n=0;
  char* p=buf;
  while (n++<9) { while (*(p++)!=' '); } n--;
  sscanf(p,"%lf %lf %lf",&a[0],&a[1],&a[2]);
  while (n++<15) { while ((*p++)!=' '); } n--;
  sscanf(p,"%lf %lf %lf",&b[0],&b[1],&b[2]);
  while (n++<21) { while ((*p++)!=' '); } n--;
  sscanf(p,"%lf %lf %lf",&c[0],&c[1],&c[2]);
  printf("<%0.2f,%0.2f,%0.2f> <%0.2f,%0.2f,%0.2f> <%0.2f,%0.2f,%0.2f>\n",a[0],a[1],a[2],b[0],b[1],b[2],c[0],c[1],c[2]);

  double r=sqrt(dist2(a,b));
  double s=sqrt(dist2(a,c));
  double t=sqrt(dist2(b,c));
  
  double mn = min(r,min(s,t));
  double mx = max(r,max(s,t));
  if (fabs(mn-d[0])>FILTER) return false;
  if (fabs(mx-d[2])>FILTER) return false;
  return true;
}

//---------------------------------------------------------------

// atomlist looks like "6,20,10,25,16,32", should be alternating Ca,Cb order

void usage()
{
  printf("usage: triplet_search [flags] <conformers_pdb> <anum_list> <triplet_db>\n");
  exit(-1);
}

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

  if (argc-optind<3) usage();
  char* ligfile=argv[optind];
  char* anumlist=argv[optind+1];

  //------------------------------------

  vector<vector<Atom*>*>* conformers=read_conformers(ligfile);
  vector<char*>* remarks=read_frames(ligfile);
  if (conformers->size()!=remarks->size()) { printf("error: unequal number of frames REMARKS as conformers\n"); exit(-1); }
  vector<int> frames;
  for (int i=0 ; i<remarks->size() ; i++)
  {
    vector<char*>* words=split(remarks->at(i),' ');
    frames.push_back(atoi(words->at(2)));
  }

  time_t curtime;
  struct tm *loctime;
  curtime = time (NULL);
  loctime = localtime (&curtime);

  printf("REMARK # command: "); for (int i=0 ; i<argc ; i++) printf("%s ",argv[i]); printf("\n");
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

  //------------------------------------

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
    vector<float> d; // triangle side lengths, for filtering
    d.push_back(sqrt(dist2(ligatoms[0]->co,ligatoms[2]->co)));
    d.push_back(sqrt(dist2(ligatoms[2]->co,ligatoms[4]->co)));
    d.push_back(sqrt(dist2(ligatoms[0]->co,ligatoms[4]->co)));
    sort(d.begin(),d.end());

    // loop through triplets database
    vector<Triplet> triplets;
    vector<Match*> matches;
    int a=0,b=0,c=0;
    char buf[10000];

    FILE* fil=fopen(argv[optind+2],"r");
    while (fgets(buf,10000,fil)!=NULL)
    {
      if (buf[0]=='#') continue;
      a++;
      if (a%100000==0) fprintf(stderr,"total triplets=%d, filtered out=%d, hits found=%d\n",a,c,b);
      //if (filter2(buf,d)==false) { c++; continue; }
      Triplet triplet(buf);
      //if (filter(triplet.atoms,d)==false) { c++; continue; }
      Match* match=compare_3vecs(ligatoms,triplet.atoms);
      if (match->rms<1.0)
      {
        triplets.push_back(triplet);
        matches.push_back(match);
        b++; 
//printf("rms=%0.2f\n",match->rms);
      }
    }
    fclose(fil);

    vector<pair<float,int> > scores;
    for (int i=0 ; i<triplets.size() ; i++)
      scores.push_back(pair<float,int>(matches[i]->rms,i));
    sort(scores.begin(),scores.end());

    for (int i=0 ; i<scores.size() ; i++)
    {
      float rms=scores[i].first;
      int j=scores[i].second;
      printf("rank=%d, rms=%0.3f, ",i,rms);
      triplets[j].print();
      printf("transform: ");
      for (int k=0 ; k<3 ; k++) printf("%0.3f ",matches[j]->cen1[k]);
      for (int k=0 ; k<9 ; k++) printf("%0.3f ",matches[j]->rot[k]);
      for (int k=0 ; k<3 ; k++) printf("%0.3f ",triplets[j].com[k]);
      printf("\n");
      // void Match::transform(double* A,double* B); // B=rot.(A-cen1)+cen2
    }

    break; // just search first conformer
  }

}

