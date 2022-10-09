
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"



/* =========================================================================================== */
int read_qme(t_parameters* pars, t_topology* top, double* qme, double* w_qme) {

  FILE *fqme;
  double sumqme, sumqmf, avgqme, avgqmf;
  int i;
  int n_qmeread;


/* open qme file */

  if(  (fqme=fopen(pars->qme_filename,"r"))==NULL) error_msg("Can not open qme file");
  else printf("# opened %s\n", pars->qme_filename);

/* read energies and forces */

  n_qmeread=0;
  rewind(fqme);
  while(fscanf(fqme,"%le %le", &(qme[n_qmeread]),&(w_qme[n_qmeread]))==2) n_qmeread++;
  fclose(fqme);

  if(n_qmeread == pars->n_qme + pars->n_qmf ) {
    printf("# Read %d pairs of values (energies/forces and weights) from %s\n", n_qmeread, pars->qme_filename);
    printf("# assuming these include %d energies and %d forces\n", pars->n_qme, pars->n_qmf);
  } else {
    error_msg("Unexpected number of items in qm file");
  }

/* convert units from Hartree/Bohr to kJ.mol and nm, and gradient to force (including sign) */

  for(i=0;           i<pars->n_qme;             i++) qme[i] *= 2625.5002;
  for(i=pars->n_qme; i<pars->n_qme+pars->n_qmf; i++) qme[i] *= -49614.759609591609;

/* determine average and stdev of energies and forces */

  if(pars->n_qme>0) {

    sumqme = 0.0;
    for(i=0; i<pars->n_qme; i++) sumqme+=qme[i];
    avgqme=sumqme/(double)pars->n_qme;
    sumqme = 0.0;
    for(i=0; i<pars->n_qme; i++) sumqme+=(qme[i]-avgqme)*(qme[i]-avgqme);
    pars->stdeve=sqrt(sumqme/(double)pars->n_qme);
    if(pars->stdeve==0.0) {
      printf("# Warning stdeve = 0. setting as 1.0");
      pars->stdeve=1.0;
    }
    for(i=0; i<pars->n_qme; i++) qme[i] = (qme[i]-avgqme)/pars->stdeve;
  }

  if(pars->n_qmf>0) {

    sumqmf= 0.0;
    for(i=pars->n_qme; i<pars->n_qme+pars->n_qmf; i++) sumqmf+=qme[i];
    avgqmf=sumqmf/(double)pars->n_qmf;
    sumqmf= 0.0;
    for(i=pars->n_qme; i<pars->n_qme+pars->n_qmf; i++) sumqmf+=(qme[i]-avgqmf)*(qme[i]-avgqmf);
    pars->stdevf=sqrt(sumqmf/(double)pars->n_qmf);
    if(pars->stdevf==0.0) {
      printf("# Warning stdevf = 0. setting as 1.0");
      pars->stdevf=1.0;
    }
    for(i=pars->n_qme; i<pars->n_qme+pars->n_qmf; i++) qme[i] = qme[i]/pars->stdevf;
  }

  printf("# Normalized QM input\n");
  if(pars->n_qme>0) printf("# <E>   = %le +/- %le\n", avgqme, pars->stdeve);
  else printf("# <|E|> = NA, no energies provided\n");
  if(pars->n_qmf>0) printf("# <|F|> = %le +/- %le\n", avgqmf, pars->stdevf);
  else printf("# <|F|> = NA, no forces provided\n");
  fflush(NULL);
  
  return(0);
}
