
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"


/* =========================================================================================== */
int write_top(const char* filename, const t_parameters* pars, t_topology* top, t_genotype* genotype, t_codonlimits* limits) {

  int i;
  FILE* ofil;
  double eq0, fc0;
  double eps, sig;

  char dovdw[]="v";
  char dobon[]="b";
  char doang[]="a";
  char dodih[]="d";
  char doimp[]="i";

  if(access(filename, F_OK)==0) {
    fprintf(stderr,"in write_mod_top: tried writing to %s\n", filename);
    error_msg("topology file exists, abort ...");
  }

  if((ofil=fopen(filename,"w"))==NULL) error_msg("cannot open topology file for writing");

  fprintf(ofil,"[ defaults ]\n");
  fprintf(ofil,"1               2               yes             0.5     0.8333\n");

/*################# types ##################*/

  fprintf(ofil,"\n[ atomtypes ]\n");

  for(i=0;i<top->n_types;i++) {

    eps=top->epsilon0[i];
    sig=top->sigma0[i];

    fprintf(ofil," %-8s %-8s 0.00000  0.00000   A  %14.5e %14.5e\n",
      top->atomtype[i],top->atomtype[i],sig,eps);
  }

/*################# vdw parameters ##################*/

  fprintf(ofil,"\n[ nonbond_params ]\n");

  for(i=0;i<top->n_nb; i++) {

    eps=top->epsilon[i];
    sig=top->sigma[i];
    if(top->egrpidx[i]!=0 && strstr(pars->opt,dovdw)!=NULL) eps=genotype->codon[limits->codonidx[top->egrpidx[i]]];
    if(top->sgrpidx[i]!=0 && strstr(pars->opt,dovdw)!=NULL) sig=genotype->codon[limits->codonidx[top->sgrpidx[i]]];

    fprintf(ofil," %-8s %-8s  1  %14.6e %14.6e ; %d %d\n",
      top->nbat1[i],top->nbat2[i],sig,eps,top->sgrpidx[i],top->egrpidx[i]);
  }

/*################# molecules ##################*/

  fprintf(ofil,"\n[ moleculetype ]\n");
  fprintf(ofil," %s   3\n\n", top->resname);

  fprintf(ofil,"\n[ atoms ]\n");

  for(i=0;i<top->n_atoms;i++) {
    fprintf(ofil,"%8d  %6s  1  %3s  %4s %5d %12.6lf %12.5lf\n", 
      top->ndex[i],
      top->type[i],
      top->resname,
      top->atomname[i],
      top->ndex[i],
      top->charge[i],
      top->mass[i]);
  }


/*################# bonds ##################*/

  fprintf(ofil,"\n[ bonds ]\n");

  for(i=0;i<top->n_bonds;i++) {
      
    eq0=top->bondr[i];
    fc0=top->bondk[i];  

    if(top->bgrpidx[i]!=0 && strstr(pars->opt,dobon)!=NULL) fc0=genotype->codon[limits->codonidx[top->bgrpidx[i]]];
    
    fprintf(ofil,"%8d%8d   1 %12.4e %12.4e ; %d\n",
      top->bondi[i],top->bondj[i],eq0,fc0,top->bgrpidx[i]);
  }

/*################# pairs ##################*/

  fprintf(ofil,"\n[ pairs ]\n");

  for(i=0;i<top->n_pairs;i++) {
    fprintf(ofil,"%8d%8d       1\n", top->pairi[i], top->pairj[i]);
  }


/*################# angles ##################*/

  fprintf(ofil,"\n[ angles ]\n");

  for(i=0;i<top->n_angles;i++) {

    eq0=top->angleeq[i];
    fc0=top->anglefc[i];
    
    if(top->agrpidx[i]!=0 && strstr(pars->opt,doang)!=NULL) eq0=genotype->codon[limits->codonidx[top->agrpidx[i]]];
    if(top->cgrpidx[i]!=0 && strstr(pars->opt,doang)!=NULL) fc0=genotype->codon[limits->codonidx[top->cgrpidx[i]]];

    fprintf(ofil,"%8d%8d%8d      1%12.4e %12.4e ; %d %d\n",
      top->anglei[i],top->anglej[i],top->anglek[i],eq0,fc0,top->agrpidx[i],top->cgrpidx[i]);
  }


/*################# dihedrals ##################*/

  fprintf(ofil,"\n[ dihedrals ]\n");

  for(i=0;i<top->n_dihedrals;i++) {

    fc0=top->dihfc[i];
    if(top->dgrpidx[i]!=0 && strstr(pars->opt,dodih)!=NULL) fc0=genotype->codon[limits->codonidx[top->dgrpidx[i]]];

    fprintf(ofil,"%8d%8d%8d%8d     9  %12.2f %12.5f %5.0f ; %d \n",
      top->dihi[i],top->dihj[i],top->dihk[i],top->dihl[i],top->dihphase[i],fc0,top->dihpn[i], top->dgrpidx[i]);
  }


/*################# impropers ##################*/

  fprintf(ofil,"\n[ dihedrals ]\n");

  for(i=0;i<top->n_impropers;i++) {

    fc0=top->impfc[i];
    if(top->igrpidx[i]!=0 && strstr(pars->opt,doimp)!=NULL) fc0=genotype->codon[limits->codonidx[top->igrpidx[i]]];
      
    fprintf(ofil,"%8d%8d%8d%8d     4  %12.2f %12.5f %5.0f ; %d\n",
       top->impi[i],top->impj[i],top->impk[i],top->impl[i],top->impphase[i],fc0,top->imppn[i], top->igrpidx[i]);
  }


  fprintf(ofil,"\n[ system ]\n");
  fprintf(ofil," system\n");
  fprintf(ofil,"\n[ molecules ]\n");
  fprintf(ofil," %s %d\n", top->resname, pars->n_mol);
  
  fclose(ofil);
  
  return(0);
}
