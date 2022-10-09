
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"

int copy_gene(t_genotype* source, t_genotype* dest) {

  int i;

  dest->length=source->length;
  dest->fitness=source->fitness;
  dest->index=source->index;    
  for(i=0;i<source->length+1;i++) {
    dest->codon[i]=source->codon[i];
  }

  return(i-1);
}
