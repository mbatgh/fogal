
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"



/* =========================================================================================== */
int mutate(t_genotype* genotype, t_codonlimits* limits, double rate) {

/*
on average for 100*rate percent codons modify with a new value
chosen within limits using a triangular prob distribution centered on the input parameter value
*/

  int i;
  double myrand,a,b,c,F,U;
  int nmut;

  nmut=0;

  for(i=1;i<=genotype->length;i++) {

    myrand=(double)rand()/(double)RAND_MAX;

    if(myrand<=rate) {

      nmut++;

      a=limits->lower[i];
      b=limits->upper[i];
      c=limits->mean[i];

      U=(double)rand()/(double)RAND_MAX;

      F=(c-a)/(b-a);

      if(U<F) genotype->codon[i] = a+sqrt(U*(b-a)*(c-a));
      else    genotype->codon[i] = b-sqrt((1-U)*(b-a)*(b-c));

    }
  }

  return(nmut);
}
