#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <time.h>
#include "triplet_lib.hpp"
#include <sys/stat.h>

using namespace std;

extern "C" {

char* lowercase(char* r)
{
  int n=strlen(r);
  char* s=(char*)malloc(n+1);
  for (int i=0 ; i<n ; i++)
  {
    if (r[i]>='A' && r[i]<='Z') s[i] = r[i]-'A'+'a';
    else s[i] = r[i];
  }
  s[n] = '\0';
  return s;
}

// this is a bad idea ; creates a copy, and loses ptr so can't delete
vector<char*> split2(char* s1,char c) { return *split(s1,c); }

int main(int argc,char** argv)
{
  // make a first pass to ensure all files exist
  FILE* fil=fopen(argv[1],"r");
  char buf[1000];
  int flag=0;

  // input: interface_pairs.hetero or interface_pairs.NA
  fil = fopen(argv[1],"r");
  while (fgets(buf,1000,fil)!=NULL)
  {
    if (strstr(buf,"xray")!=NULL) printf("%s",buf);
    if (buf[0]=='#') continue;
    buf[strlen(buf)-1] = '\0'; // strip CRLF
    vector<char*> words = split2(buf,' ');
    //char* name=lowercase(words[0]);
    char* name=words[0];
    char fname[1000];
    //sprintf(fname,"/ocean/bali/ioerger/burgess/PP/pdb%s.ent",name);
    //sprintf(fname,"/ocean/oahu/pdb/all/pdb%s.ent",name);
    sprintf(fname,"/apps/burgess/pdb/%s.pdb",name);

    struct stat statbuf;
    if (stat(fname,&statbuf)==-1) { 
      printf("# error: %s not found\n",fname);
      continue; }

    char ch1=words[1][0],ch2=words[2][0];
    int NA=0;
    if (strcmp(words[2],"NA")==0) NA = 1;
     
    vector<Atom*>* pdb=read_pdb(fname); 
    vector<Atom*> A,B;
    for (int i=0 ; i<pdb->size() ; i++)
    {
      Atom* atom=pdb->at(i);
      if (atom->chn==ch1) A.push_back(atom);
      if (NA==1 && nucleic(atom)==1) B.push_back(atom);
      else if (NA==0 && atom->chn==ch2) B.push_back(atom); 
    }

    vector<int>* contacts=interface_residues(A,B);
    if (contacts->size()>=3) printf("%s %c %s %ld\n",words[0],ch1,words[2],contacts->size());
    else printf("# %s %c %s too few contacts\n",name,ch1,words[2]);
    fflush(stdout);
    for (int i=0 ; i<pdb->size() ; i++) delete pdb->at(i);
    delete pdb;
    delete contacts; 
  }
}


}
