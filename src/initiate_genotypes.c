
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"


/* =========================================================================================== */
int initiate_genotypes(t_topology* top,t_genotype** parent,t_genotype** child,
                       t_codonlimits* limits, t_parameters* pars) {

  int i,j;
  int n_codons;
  int gidx[9999];
  int maxn_codon;

  char dovdw[]="v";
  char dobon[]="b";
  char doang[]="a";
  char dodih[]="d";
  char doimp[]="i";

/***** determine genome size *****/

  for(i=0;i<9999;i++) gidx[i]=0;

  if(strstr(pars->opt,dovdw)!=NULL) {
    for(i=0;i<top->n_nb;i++) {
      gidx[top->egrpidx[i]]=1;
      gidx[top->sgrpidx[i]]=1;
    }
  }

  if(strstr(pars->opt,dobon)!=NULL) {
    for(i=0;i<top->n_bonds;i++) {
      gidx[top->bgrpidx[i]]=1;
    }
  }

  if(strstr(pars->opt,doang)!=NULL) {
    for(i=0;i<top->n_angles;i++) {
      gidx[top->agrpidx[i]]=1;
      gidx[top->cgrpidx[i]]=1;
    }
  }

  if(strstr(pars->opt,dodih)!=NULL) {
    for(i=0;i<top->n_dihedrals;i++) {
      gidx[top->dgrpidx[i]]=1;
    }
  }

  if(strstr(pars->opt,doimp)!=NULL) {
    for(i=0;i<top->n_impropers;i++) {
      gidx[top->igrpidx[i]]=1;
    }
  }

  n_codons=0;

  for(i=1;i<9999;i++) if(gidx[i]==1) n_codons++;

  maxn_codon= 2*top->n_nb + top->n_bonds + 2*(top->n_angles) + top->n_dihedrals + top->n_impropers;

/***** allocate memory for genomes *****/

  limits->lower = (double*)calloc((size_t)(n_codons+1),sizeof(double));
  limits->mean  = (double*)calloc((size_t)(n_codons+1),sizeof(double));
  limits->upper = (double*)calloc((size_t)(n_codons+1),sizeof(double));
  limits->codonidx = (int*)calloc((size_t)(maxn_codon+1),sizeof(int));
  limits->topolidx = (int*)calloc((size_t)(n_codons+1),sizeof(int));

  n_codons=0;
  for(i=1;i<9999;i++) {
    if(gidx[i]==1) {
      n_codons++;
      limits->topolidx[n_codons]=i;
      limits->codonidx[i]=n_codons;
    }
  }

  *parent = (t_genotype*)calloc(pars->n_individuals_max, sizeof(t_genotype));
  *child  = (t_genotype*)calloc(pars->n_individuals_max, sizeof(t_genotype));

  for(i=0;i<pars->n_individuals_max;i++) {
    (*parent)[i].codon = calloc((size_t)(n_codons+1),sizeof(double));
    (*parent)[i].index = i;
    (*parent)[i].length = n_codons;
    (*parent)[i].fitness=0.0;

    (*child)[i].codon = calloc((size_t)(n_codons+1),sizeof(double));
    (*child)[i].index = i;
    (*child)[i].length = n_codons;
    (*child)[i].fitness=0.0;
  }

/***** assign values to inital genome *****/

  if(strstr(pars->opt,dovdw)!=NULL) {

    for(i=0;i<top->n_nb;i++) {

      for(j=0;j<pars->n_individuals_max;j++) {
        (*parent)[j].codon[limits->codonidx[top->egrpidx[i]]] = top->epsilon[i];
        (*parent)[j].codon[limits->codonidx[top->sgrpidx[i]]] = top->sigma[i];
      }

      limits->mean[limits->codonidx[top->egrpidx[i]]]  = top->epsilon[i];
      limits->lower[limits->codonidx[top->egrpidx[i]]] = top->epsilon[i]*FRACDEPSMIN;
      limits->upper[limits->codonidx[top->egrpidx[i]]] = top->epsilon[i]*FRACDEPSMAX;

      limits->mean[limits->codonidx[top->sgrpidx[i]]]  = top->sigma[i];
      limits->lower[limits->codonidx[top->sgrpidx[i]]] = top->sigma[i]*FRACDSIGMIN;
      limits->upper[limits->codonidx[top->sgrpidx[i]]] = top->sigma[i]*FRACDSIGMAX;
    }
  }

  if(strstr(pars->opt,dobon)!=NULL) {

    for(i=0;i<top->n_bonds;i++) {

      for(j=0;j<pars->n_individuals_max;j++) {
        (*parent)[j].codon[limits->codonidx[top->bgrpidx[i]]] = top->bondk[i];
      }

      limits->mean[limits->codonidx[top->bgrpidx[i]]]  = top->bondk[i];
      limits->lower[limits->codonidx[top->bgrpidx[i]]] = BONDFCMIN;
      limits->upper[limits->codonidx[top->bgrpidx[i]]] = BONDFCMAX;
    }
  }


  if(strstr(pars->opt,doang)!=NULL) {

    for(i=0;i<top->n_angles;i++) {

      for(j=0;j<pars->n_individuals_max;j++) {
        (*parent)[j].codon[limits->codonidx[top->agrpidx[i]]] = top->angleeq[i];
        (*parent)[j].codon[limits->codonidx[top->cgrpidx[i]]] = top->anglefc[i];
      }

      limits->mean[limits->codonidx[top->agrpidx[i]]]  = top->angleeq[i];
      limits->lower[limits->codonidx[top->agrpidx[i]]] = top->angleeq[i]-DANGMAX;
      limits->upper[limits->codonidx[top->agrpidx[i]]] = top->angleeq[i]+DANGMAX;

      limits->mean[limits->codonidx[top->cgrpidx[i]]]  = top->anglefc[i];
      limits->lower[limits->codonidx[top->cgrpidx[i]]] = 0.0;
      limits->upper[limits->codonidx[top->cgrpidx[i]]] = ANGLEFCMAX;
    }
  }


  if(strstr(pars->opt,dodih)!=NULL) {

    for(i=0;i<top->n_dihedrals;i++) {

      for(j=0;j<pars->n_individuals_max;j++) {
        (*parent)[j].codon[limits->codonidx[top->dgrpidx[i]]] = top->dihfc[i];
      }

      limits->mean[limits->codonidx[top->dgrpidx[i]]]  = top->dihfc[i];
      limits->lower[limits->codonidx[top->dgrpidx[i]]] = DIHEDRALFCMIN;
      limits->upper[limits->codonidx[top->dgrpidx[i]]] = DIHEDRALFCMAX;
    }
  }

  if(strstr(pars->opt,doimp)!=NULL) {

    for(i=0;i<top->n_impropers;i++) {

      for(j=0;j<pars->n_individuals_max;j++) {
        (*parent)[j].codon[limits->codonidx[top->igrpidx[i]]] = top->impfc[i];
      }

      limits->mean[limits->codonidx[top->igrpidx[i]]]  = top->impfc[i];
      limits->lower[limits->codonidx[top->igrpidx[i]]] = 0.0;
      limits->upper[limits->codonidx[top->igrpidx[i]]] = IMPROPERFCMAX;
    }
  }

  printf("# We have %d codons in each of %d individual genomes\n", n_codons,pars->n_individuals_max);

  return(n_codons);
}
