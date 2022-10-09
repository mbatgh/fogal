#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"



/* =========================================================================================== */

int allocate_memory(t_parameters* pars,t_topology* top, double*** coords,
                    double** qme, double** w_qme, double** mme, double** mme0) {

  int section=0;
  int n_qmeread;
  char buf[256];
  double ddummy1,ddummy2;
  FILE* fqme;
  FILE* ftop;


/* open topology file and determine size */

  if( ( ftop=fopen(pars->topology_filename,"r") )==NULL)
    error_msg("Can not open topology file");

  top->n_types=0;
  top->n_atoms=0;
  top->n_pairs=0;
  top->n_bonds=0;
  top->n_angles=0;
  top->n_dihedrals=0;
  top->n_impropers=0;
  top->n_nb=0;

  while(fgets(buf, 256, (FILE*)ftop)) {

    if(strstr(buf,"atomtypes")!=NULL) {section=1;if(fgets(buf, 256, (FILE*)ftop)==NULL) error_msg("problem reading");}
    if(strstr(buf,"atoms")!=NULL)     {section=2;if(fgets(buf, 256, (FILE*)ftop)==NULL) error_msg("problem reading");}
    if(strstr(buf,"bonds")!=NULL)     {section=3;if(fgets(buf, 256, (FILE*)ftop)==NULL) error_msg("problem reading");}
    if(strstr(buf,"pairs")!=NULL)     {section=4;if(fgets(buf, 256, (FILE*)ftop)==NULL) error_msg("problem reading");}
    if(strstr(buf,"angles")!=NULL)    {section=5;if(fgets(buf, 256, (FILE*)ftop)==NULL) error_msg("problem reading");}
    if(strstr(buf,"dihedrals")!=NULL) {
      if(section==6) section=7; else section=6;
      if(fgets(buf, 256, (FILE*)ftop)==NULL) error_msg("problem reading");
    }
    if(strstr(buf,"nonbond_params")!=NULL) {section=8; if(fgets(buf, 256, (FILE*)ftop)==NULL) error_msg("problem reading");}

    if(section==1 && words(buf)>=7) (top->n_types)++;
    if(section==2 && words(buf)>=8) (top->n_atoms)++;
    if(section==3 && words(buf)>=5) (top->n_bonds)++;
    if(section==4 && words(buf)>=3) (top->n_pairs)++;
    if(section==5 && words(buf)>=6) (top->n_angles)++;
    if(section==6 && words(buf)>=8) (top->n_dihedrals)++;
    if(section==7 && words(buf)>=8) (top->n_impropers)++;
    if(section==8 && words(buf)>=5) (top->n_nb)++;
  }

  top->atomtype  = c0mat(top->n_types,4);
  top->sigma0    = d0vec(top->n_types);
  top->epsilon0  = d0vec(top->n_types);

  top->nbat1     = c0mat(top->n_types*top->n_types,4);
  top->nbat2     = c0mat(top->n_types*top->n_types,4);
  top->sigma     = d0vec(top->n_nb);
  top->epsilon   = d0vec(top->n_nb);
  top->sgrpidx   = i0vec(top->n_nb);
  top->egrpidx   = i0vec(top->n_nb);
  top->lj12      = d0vec(top->n_atoms*pars->n_mol*top->n_atoms*pars->n_mol);
  top->lj6       = d0vec(top->n_atoms*pars->n_mol*top->n_atoms*pars->n_mol);

  top->type     = c0mat(top->n_atoms,4);
  top->atomname = c0mat(top->n_atoms,5);
  top->charge   = d0vec(top->n_atoms);
  top->mass     = d0vec(top->n_atoms);
  top->ndex     = i0vec(top->n_atoms);
  top->ish      = i0vec(top->n_atoms + 1);

  top->bondi   = i0vec(top->n_bonds);
  top->bondj   = i0vec(top->n_bonds);
  top->bondr   = d0vec(top->n_bonds);
  top->bondk   = d0vec(top->n_bonds);
  top->bgrpidx = i0vec(top->n_bonds);
  top->isconstrained = i0vec(top->n_bonds);

  top->pairi = i0vec(top->n_pairs);
  top->pairj = i0vec(top->n_pairs);

  top->anglei          = i0vec(top->n_angles);
  top->anglej          = i0vec(top->n_angles);
  top->anglek          = i0vec(top->n_angles);
  top->angleeq         = d0vec(top->n_angles);
  top->anglefc         = d0vec(top->n_angles);
  top->agrpidx         = i0vec(top->n_angles);
  top->cgrpidx         = i0vec(top->n_angles);

  top->dihi          = i0vec(top->n_dihedrals);
  top->dihj          = i0vec(top->n_dihedrals);
  top->dihk          = i0vec(top->n_dihedrals);
  top->dihl          = i0vec(top->n_dihedrals);
  top->dihphase      = d0vec(top->n_dihedrals);
  top->dihfc         = d0vec(top->n_dihedrals);
  top->dihpn         = d0vec(top->n_dihedrals);
  top->dgrpidx       = i0vec(top->n_dihedrals);

  top->impi          = i0vec(top->n_impropers);
  top->impj          = i0vec(top->n_impropers);
  top->impk          = i0vec(top->n_impropers);
  top->impl          = i0vec(top->n_impropers);
  top->impphase      = d0vec(top->n_impropers);
  top->impfc         = d0vec(top->n_impropers);
  top->imppn         = d0vec(top->n_impropers);
  top->igrpidx       = i0vec(top->n_impropers);

  top->neighbor1     = i0vec(pars->n_mol*pars->n_mol*top->n_atoms*top->n_atoms);
  top->neighbor2     = i0vec(pars->n_mol*pars->n_mol*top->n_atoms*top->n_atoms);

  top->nb_el_fac     = d0vec(pars->n_mol*top->n_atoms*pars->n_mol*top->n_atoms);
  top->nb_lj_fac     = d0vec(pars->n_mol*top->n_atoms*pars->n_mol*top->n_atoms);

  top->nbidx         = i0vec(top->n_atoms*pars->n_mol*top->n_atoms*pars->n_mol);

/* open qme file and count entries */  

  if(  (fqme=fopen(pars->qme_filename,"r"))==NULL) error_msg("Can not open qme file");
  else printf("# opened %s\n", pars->qme_filename);

  n_qmeread=0;
  while(fscanf(fqme,"%le %le", &ddummy1,&ddummy2)==2) n_qmeread++;
  if(n_qmeread != pars->n_qme + pars->n_qmf) {
    printf("expected: %d, found: %d", pars->n_qme+pars->n_qmf, n_qmeread);
    error_msg("Unexpected number of values in qme file");
  }

/* allocate memory for coordinates, QM and MM eneries and forces */

  *coords            = d0mat(top->n_atoms*pars->n_mol*pars->n_conf,3);

  *qme               = (double*)calloc(n_qmeread, sizeof(double));
  *w_qme             = (double*)calloc(n_qmeread, sizeof(double));
  *mme               = (double*)calloc(n_qmeread, sizeof(double));
  *mme0              = (double*)calloc(n_qmeread, sizeof(double));

  fclose(ftop);
  fclose(fqme);

  return(n_qmeread);

}
