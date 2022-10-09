#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"


/* =========================================================================================== */
int read_topology(const t_parameters* pars, t_topology* top) {

  int section=0;
  char tmpchar1[256];
  char tmpchar2[256];
  double tmplf;
  double* ptmplf=&tmplf;
  int tmpint;
  int* ptmpint=&tmpint;
  char c1;
  char* pc1=&c1;
  int tmpint1;
  int* ptmpint1=&tmpint1;
  int tmpint2;
  int* ptmpint2=&tmpint2;

  int idum1, idum2;
  int nitems;

  int counted_types=0;
  int counted_atoms=0;
  int counted_pairs=0;
  int counted_bonds=0;
  int counted_angles=0;
  int counted_dihedrals=0;
  int counted_impropers=0;
  int counted_nb=0;

  char buf[256];
  FILE *ftop;

  if( ( ftop=fopen(pars->topology_filename,"r") )==NULL)
    error_msg("Can not open topology file");

  while(fgets(buf, 256, (FILE*)ftop)) {

    if(strstr(buf,"nonbond_params")!=NULL) {section=8; if(fgets(buf, 256, (FILE*)ftop)==NULL) error_msg("problem reading topology");}
    if(strstr(buf,"atomtypes")!=NULL)      {section=1; if(fgets(buf, 256, (FILE*)ftop)==NULL) error_msg("problem reading topology");}
    if(strstr(buf,"atoms")!=NULL)          {section=2; if(fgets(buf, 256, (FILE*)ftop)==NULL) error_msg("problem reading topology");}
    if(strstr(buf,"bonds")!=NULL)          {section=3; if(fgets(buf, 256, (FILE*)ftop)==NULL) error_msg("problem reading topology");}
    if(strstr(buf,"pairs")!=NULL)          {section=4; if(fgets(buf, 256, (FILE*)ftop)==NULL) error_msg("problem reading topology");}
    if(strstr(buf,"angles")!=NULL)         {section=5; if(fgets(buf, 256, (FILE*)ftop)==NULL) error_msg("problem reading topology");}
    if(strstr(buf,"dihedrals")!=NULL)      {
      if(section==6) section=7; else section=6;
      if(fgets(buf, 256, (FILE*)ftop)==NULL) error_msg("problem reading topology");
    }

    nitems=words(buf);

/***** types *****/

    if(section==1 && nitems>=7) {

      if(sscanf(buf,"%s %s %lf %lf %s %lf %lf %c %d %d", top->atomtype[counted_types], tmpchar1,
        ptmplf, ptmplf, tmpchar2, &(top->sigma0[counted_types]), &(top->epsilon0[counted_types]), pc1, &idum1, &idum2)>=7) counted_types++;
    }

/***** nonbonded LJ  *****/

    if(section==8 && nitems>=5) {

      if(sscanf(buf,"%s %s %d %lf %lf %c %d %d", top->nbat1[counted_nb], top->nbat2[counted_nb], 
            &idum1, &(top->sigma[counted_nb]), &(top->epsilon[counted_nb]), pc1, 
            &(top->sgrpidx[counted_nb]), &(top->egrpidx[counted_nb])) >= 5) counted_nb++;
    }

/***** atoms *****/

    if(section==2 && nitems>=8) {

      if(sscanf(buf,"%d %s %d %s %s %d %lf %lf",
        &(top->ndex[counted_atoms]),
        top->type[counted_atoms],
        ptmpint1,
        top->resname,
        top->atomname[counted_atoms],
        ptmpint2,
    	&(top->charge[counted_atoms]),
        &(top->mass[counted_atoms]))>=8) {
          if(top->type[counted_atoms][0]=='h')
            top->ish[top->ndex[counted_atoms]] = 1;
          else top->ish[top->ndex[counted_atoms]]=0;
          counted_atoms++;
      }

    }

/***** bonds *****/

    if(section==3 && nitems>=5) {

      if(sscanf(buf,"%d %d %d %lf %lf %c %d",
           &(top->bondi[counted_bonds]),
           &(top->bondj[counted_bonds]),
           ptmpint,
           &(top->bondr[counted_bonds]),
           &(top->bondk[counted_bonds]),
           pc1,
           &(top->bgrpidx[counted_bonds]))>=5 ) {
             if((top->ish[top->bondi[counted_bonds]]==1 || top->ish[top->bondj[counted_bonds]]==1) && pars->hcons==1) {
               top->isconstrained[counted_bonds]=1;
             } else { top->isconstrained[counted_bonds]=0; }
             counted_bonds++;
           }
    }

/***** pairs *****/

    if(section==4 && nitems>=3) {
      if(sscanf(buf,"%d %d %d",&(top->pairi[counted_pairs]),&(top->pairj[counted_pairs]),ptmpint)>=3) counted_pairs++;
    }

/***** angles *****/

    if(section==5 && nitems>=6) {

      if(sscanf(buf,"%d %d %d %d %lf %lf %c %d %d",
        &(top->anglei[counted_angles]),
        &(top->anglej[counted_angles]),
        &(top->anglek[counted_angles]),
        ptmpint,
        &(top->angleeq[counted_angles]),
        &(top->anglefc[counted_angles]),
        pc1,
        &(top->agrpidx[counted_angles]),
        &(top->cgrpidx[counted_angles]))>=6) counted_angles++;
    }

/***** dihedrals *****/

    if(section==6 && nitems>=8) {

      if(sscanf(buf,"%d %d %d %d %d %lf %lf %lf %c %d",
        &(top->dihi[counted_dihedrals]),
        &(top->dihj[counted_dihedrals]),
        &(top->dihk[counted_dihedrals]),
        &(top->dihl[counted_dihedrals]),
        ptmpint,
        &(top->dihphase[counted_dihedrals]),
        &(top->dihfc[counted_dihedrals]),
        &(top->dihpn[counted_dihedrals]),
        pc1,
        &(top->dgrpidx[counted_dihedrals]))>=8) counted_dihedrals++;
    }

/***** impropers *****/

    if(section==7 && nitems>=8) {

      if(sscanf(buf,"%d %d %d %d %d %lf %lf %lf %c %d",
        &(top->impi[counted_impropers]),
        &(top->impj[counted_impropers]),
        &(top->impk[counted_impropers]),
        &(top->impl[counted_impropers]),
        ptmpint,
        &(top->impphase[counted_impropers]),
        &(top->impfc[counted_impropers]),
        &(top->imppn[counted_impropers]),
        pc1,
        &(top->igrpidx[counted_impropers]))>=8) counted_impropers++;
    }
  }

  fclose(ftop);

  printf("# Read initial topology from %s\n", pars->topology_filename);
  printf("#   number of types:     %8d (%d)\n", top->n_types, counted_types);
  printf("#   number of atoms:     %8d (%d)\n", top->n_atoms, counted_atoms);
  printf("#   number of bonds:     %8d (%d)\n", top->n_bonds, counted_bonds);
  printf("#   number of pairs:     %8d (%d)\n", top->n_pairs, counted_pairs);
  printf("#   number of angles:    %8d (%d)\n", top->n_angles, counted_angles);
  printf("#   number of dihedrals: %8d (%d)\n", top->n_dihedrals, counted_dihedrals);
  printf("#   number of impropers: %8d (%d)\n", top->n_impropers, counted_impropers);
  printf("#   number of vdw pars:  %8d (%d)\n", top->n_nb, counted_nb);

  fflush(stdout);

  return(0);
}
