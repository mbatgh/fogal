#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"

int read_coords(t_parameters* pars,t_topology* top, double** coords) {

  char* buffer;
  size_t ncread;
  size_t bufsize;
  int nread;
  FILE* fcor;
  int i,j;
  char sdummy[8];
  int n_coords_read;

  bufsize = 500;
  buffer = (char *)malloc(bufsize * sizeof(char));

/***** open coord file *****/

  if( ( fcor=fopen(pars->xyz_filename,"r") )==NULL) error_msg("Can not open xyz file");
  else printf("# reading corrdinate data from %s\n", pars->xyz_filename);

  n_coords_read=0;

  for(i=0;i<pars->n_conf;i++) {

/***** read xyz file coordinates *****/

    for(j=0; j < top->n_atoms*pars->n_mol; j++) {
      ncread = getline(&buffer,&bufsize,fcor);
      nread = sscanf(buffer, " %s %lf %lf %lf", sdummy, &(coords[n_coords_read][XX]), &(coords[n_coords_read][YY]), &(coords[n_coords_read][ZZ]));
      if(nread!=4) error_msg("trouble reading coordinates 2");

/***** convert to nano meter *****/

      coords[n_coords_read][XX] *= 0.1;
      coords[n_coords_read][YY] *= 0.1;
      coords[n_coords_read][ZZ] *= 0.1;

      n_coords_read++;
    }
  }

  fclose(fcor);

  printf("# Allocated memory and read %d (%d*%d*%d) sets of 3D coordinates from %s\n", 
         n_coords_read, top->n_atoms,pars->n_mol,pars->n_conf,pars->xyz_filename);
  free(buffer);
  return(0);
}
