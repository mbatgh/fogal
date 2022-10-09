#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"
#include "vec.h"


int write_mme(const char* filename, const t_parameters* pars, t_genotype* genome, t_topology* top, 
              t_codonlimits* limits, double** coords, double* qme, double* w_qme, double* mme, double* mme0) {


  int i,nread;
  double sume,sum,avge,f,f1;

  FILE* fef;

  nread = fitness(pars, genome, top, limits, coords, qme, w_qme, mme, mme0);

/***** determine average of energies *****/

  if(pars->n_qme>0) {
    sume = 0.0;
    for(i=0; i<pars->n_qme; i++) sume += mme[i];
    avge=sume/(double)pars->n_qme;
    for(i=0; i<pars->n_qme; i++) mme[i] = (mme[i]-avge)/pars->stdeve;
  }
  for(i=pars->n_qme; i<pars->n_qme+pars->n_qmf; i++) mme[i] = mme[i]/pars->stdevf;

/***** open output file *****/

  if((fef=fopen(filename,"w"))==NULL) error_msg("cannot open ef-init file");

  for(i=0;i<pars->n_qme;i++) {
    fprintf(fef,"%8d %14.6le %14.6le %14.6le %14.6le %14.6le\n", i, qme[i], mme0[i], mme[i], qme[i]*pars->stdeve, mme[i]*pars->stdeve);
  }
  for(i=pars->n_qme; i<pars->n_qme+pars->n_qmf;i++) {
    fprintf(fef,"%8d %14.6le %14.6le %14.6le %14.6le %14.6le\n", i, qme[i], mme0[i], mme[i], qme[i]*pars->stdevf, mme[i]*pars->stdevf);
  }

/***** calculate fitness *****/

  f = 0.0;
  sum = 0.0;
  for(i=0;i<pars->n_qme;i++) {
    f1=(qme[i]-mme[i]);
    f += f1*f1*w_qme[i];
    sum += w_qme[i];
  }
  for(i=pars->n_qme; i<pars->n_qme+pars->n_qmf;i++) {
    f1=(qme[i]-mme[i]);
    f += f1*f1*w_qme[i];
    sum += w_qme[i];
  }
  genome->fitness = sqrt(f/sum);

  printf("# Best fitness: %lf\n", genome->fitness);
  printf("# saved QM vs GMX energies in %s\n", filename);

  fclose(fef);

  return(0);
}
