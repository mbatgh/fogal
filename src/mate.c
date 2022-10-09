
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"


int mate(t_genotype* parent1, t_genotype* parent2, t_genotype* child, int n, double p_crossover) {

  int i, counted;
  double ra;

  counted=0;

  for(i=1;i<=parent1->length;i++) {

    ra=((double)rand()/RAND_MAX);
    if(ra<p_crossover) counted++;

/* if counted is incremented change to the other parent */
/* the larger p_crossover the more often this will occur */

    if (counted & 1 == 1) child->codon[i]=parent1->codon[i];
    else                  child->codon[i]=parent2->codon[i];

  }

  child->length=parent1->length;
  child->index=n;
  child->fitness=0.0;

  return(counted-1);
}
