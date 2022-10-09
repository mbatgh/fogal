
/*
available: b d j w y x

h   error_msg("")
z   PFILE

a   pars->initial_mutrate
c   pars->p_crossover
e   IQMEFILE
f   pars->f_elite
g   pars->n_generations_max
i   pars->n_individuals_max
k   pars->n_conf
l   pars->f_lucky
m   pars->codon_mutrate
n   pars->n_mol
o   IOPT
p   ITOPFILE
q   pars->n_qme
v   pars->n_qmf
r   IRUNID
s   pars->rseed
t   XYZFILE
u   pars->individual_mutrate
w   pars->vdw_cutoff
b   pars->hcons
 
initial_mutationrate     pars->initial_mutrate
p_crossover              pars->p_crossover
qmefile                  IQMEFILE
f_elite                  pars->f_elite
n_individuals            pars->n_individuals_max
f_lucky                  pars->f_lucky
mutationrate_genome      pars->codon_mutrate
n_mol                    pars->n_mol
n_generations            pars->n_generations_max
opt                      IOPT
topfile                  ITOPFILE
n_qme                    pars->n_qme
n_conformers             pars->n_conf
runid                    IRUNID
rseed                    pars->rseed
xyzfile                  XYZFILE
mutationrate_population  pars->individual_mutrate
n_qmf                    pars->n_qmf
hcons                    pars->hcons
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"


int get_parameters(int ARGC,char** ARGV, t_parameters* pars) {

/***** variables set in input file or on the command line *****/

  int i;
  int irf,irc;

  FILE* fp;
  char buf[256];
  char tmpstring[256];

  int IRUNID              = 0;
  int IOPT                = 0;
  int ITOPFILE            = 0;
  int IQMEFILE            = 0;
  int XYZFILE             = 0;
  int PFILE               = 0;

  time_t t;

  pars->vdw_cutoff         = -1.0;

  pars->n_generations_max  = -1;
  pars->n_individuals_max  = -1;
  pars->n_qme              = -1;
  pars->n_qmf              = -1;
  pars->n_mol              = -1;
  pars->n_conf             = -1;
  pars->hcons              = -1;

  pars->initial_mutrate    = -1.0;
  pars->codon_mutrate      = -1.0;
  pars->individual_mutrate = -1.0;
  pars->p_crossover        = -1.0;
  pars->f_elite            = -1.0;
  pars->f_lucky            = -1.0;

  pars->rseed              = -1;

  pars->topology_filename  = calloc(256, sizeof(char));
  pars->qme_filename       = calloc(256, sizeof(char));
  pars->xyz_filename       = calloc(256, sizeof(char));
  pars->runid              = calloc(256, sizeof(char));
  pars->opt                = calloc(256, sizeof(char));


/***** ... analyse the command line arguments ... *****/

  irc=0;


  for(i=1;i<ARGC;i++) {

    if(ARGV[i][0]!='-') error_msg("wrong format in command line");
    else switch(ARGV[i][1]) {

        case 'z': PFILE = ++i;
                  irc++;
                  break;
        case 'b': pars->hcons = atoi(ARGV[++i]);
                  irc++;
                  break;
        case 'o': IOPT = ++i;
                  irc++;
                  break;
        case 'g': pars->n_generations_max = atoi(ARGV[++i]);
                  irc++;
                  break;
        case 'i': pars->n_individuals_max = atoi(ARGV[++i]);
                  irc++;
                  break;
        case 'a': pars->initial_mutrate = atof(ARGV[++i]);
                  irc++;
                  break;
        case 'm': pars->codon_mutrate = atof(ARGV[++i]);
                  irc++;
                  break;
        case 'u': pars->individual_mutrate = atof(ARGV[++i]);
                  irc++;
                  break;
        case 'c': pars->p_crossover = atof(ARGV[++i]);
                  irc++;
                  break;
        case 'f': pars->f_elite = atof(ARGV[++i]);
                  irc++;
                  break;
        case 'l': pars->f_lucky = atof(ARGV[++i]);
                  irc++;
                  break;
        case 's': pars->rseed = atoi(ARGV[++i]);
                  irc++;
                  break;
        case 'n': pars->n_mol = atoi(ARGV[++i]);
                  irc++;
                  break;
        case 'q': pars->n_qme = atoi(ARGV[++i]);
                  irc++;
                  break;
        case 'v': pars->n_qmf = atoi(ARGV[++i]);
                  irc++;
                  break;
        case 'p': ITOPFILE = ++i;
                  irc++;
                  break;
        case 'e': IQMEFILE = ++i;
                  irc++;
                  break;
        case 't': XYZFILE = ++i;
                  irc++;
                  break;
        case 'r': IRUNID = ++i;
                  irc++;
                  break;
        case 'k': pars->n_conf = atoi(ARGV[++i]);
                  irc++;
                  break;
        case 'w': pars->vdw_cutoff = atof(ARGV[++i]);
                  irc++;
                  break;
        case 'h': error_msg("");
                  break;
        default:  error_msg("unrecognized command line parameter");
      }
  }


  if(ITOPFILE>0) strncpy(pars->topology_filename,ARGV[ITOPFILE],255);
  if(IQMEFILE>0) strncpy(pars->qme_filename,ARGV[IQMEFILE],255);
  if(XYZFILE>0)  strncpy(pars->xyz_filename,ARGV[XYZFILE],255);
  if(IRUNID>0)   strncpy(pars->runid,ARGV[IRUNID],255);
  if(IOPT>0)     strncpy(pars->opt,ARGV[IOPT],255);


  if(PFILE!=0) {

    if(  (fp=fopen(ARGV[PFILE],"r"))==NULL  ) error_msg("# cannot open parameter file");

    else {

      printf("# opened %s\n", ARGV[PFILE]);

      irf=0;

/***** remove comments *****/

      while(fgets(buf, 256, (FILE*)fp)) {

        i=0;
        do {
          if(buf[i]=='#' || buf[i]==';') buf[i]='\0';
          i++;
        } while (buf[i]!='\n' && buf[i-1]!='\0');


/***** how many tokens are left? *****/

        if(words(buf)==2) {

          if(strstr(buf,"n_conformers") !=NULL && pars->n_conf==-1)
            if(sscanf(buf,"%s %d",  tmpstring, &(pars->n_conf))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"n_generations") !=NULL && pars->n_generations_max==-1)
            if(sscanf(buf,"%s %d",  tmpstring, &(pars->n_generations_max))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"n_individuals")    !=NULL && pars->n_individuals_max==-1)
            if(sscanf(buf,"%s %d",  tmpstring, &(pars->n_individuals_max ))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"n_qme") !=NULL && pars->n_qme==-1)
             if(sscanf(buf,"%s %d",  tmpstring, &(pars->n_qme))==2) irf++;
             else error_msg("Problem reading parameter file");

          if(strstr(buf,"n_qmf") !=NULL && pars->n_qmf==-1)
             if(sscanf(buf,"%s %d",  tmpstring, &(pars->n_qmf))==2) irf++;
             else error_msg("Problem reading parameter file");

          if(strstr(buf,"n_mol") !=NULL && pars->n_mol==-1)
             if(sscanf(buf,"%s %d",  tmpstring, &(pars->n_mol))==2) irf++;
             else error_msg("Problem reading parameter file");

          if(strstr(buf,"initial_mutationrate") !=NULL && pars->initial_mutrate==-1.0)
            if(sscanf(buf,"%s %lf", tmpstring, &(pars->initial_mutrate))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"mutationrate_genome")  !=NULL && pars->codon_mutrate==-1.0) 
            if(sscanf(buf,"%s %lf", tmpstring, &(pars->codon_mutrate))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"mutationrate_population")  !=NULL && pars->individual_mutrate==-1.0) 
            if(sscanf(buf,"%s %lf", tmpstring, &(pars->individual_mutrate))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"p_crossover")  !=NULL && pars->p_crossover==-1.0) 
            if(sscanf(buf,"%s %lf", tmpstring, &(pars->p_crossover))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"f_elite")  !=NULL && pars->f_elite==-1.0) 
            if(sscanf(buf,"%s %lf", tmpstring, &(pars->f_elite))==2) irf++; 
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"f_lucky")  !=NULL && pars->f_lucky==-1.0)
            if(sscanf(buf,"%s %lf", tmpstring, &(pars->f_lucky))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"rseed")   !=NULL && pars->rseed==-1)
            if(sscanf(buf,"%s %d",  tmpstring, &(pars->rseed))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"vdw_cutoff") !=NULL && pars->vdw_cutoff==-1.0)
            if(sscanf(buf,"%s %lf",  tmpstring, &(pars->vdw_cutoff))==2) irf++;
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"topfile") !=NULL && ITOPFILE==0)
            if(sscanf(buf,"%s %s", tmpstring, pars->topology_filename)==2) {ITOPFILE=1; irf++;}
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"qmefile") !=NULL && IQMEFILE==0)
            if(sscanf(buf,"%s %s", tmpstring, pars->qme_filename)==2)      {IQMEFILE=1; irf++;}
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"xyzfile")!=NULL && XYZFILE==0)
            if(sscanf(buf,"%s %s", tmpstring, pars->xyz_filename)==2)     {XYZFILE=1; irf++;}
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"runid")   !=NULL && IRUNID==0)
            if(sscanf(buf,"%s %s", tmpstring, pars->runid)==2)             {IRUNID=1; irf++;}
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"opt")     !=NULL && IOPT==0)
            if(sscanf(buf,"%s %s", tmpstring, pars->opt)==2)               {IOPT=1; irf++;}
            else error_msg("Problem reading parameter file");

          if(strstr(buf,"hcons")   !=NULL && pars->hcons==-1)
            if(sscanf(buf,"%s %d",  tmpstring, &(pars->hcons))==2) irf++;
            else error_msg("Problem reading parameter file");
        }
      }
    }
  } else {
      printf("# WARNING ... no parameter file provided, using command line or defaults\n");
  }

/***** write info *****/

  if(IRUNID  ==0)        error_msg("No run-id given");
  if(IOPT    ==0)        error_msg("No opt string given");
  if(ITOPFILE==0)        error_msg("No topology filename given");
  if(IQMEFILE==0)        error_msg("No qme filename given");
  if(XYZFILE ==0)        error_msg("No xyz filename given");
  if(pars->n_conf == -1) error_msg("Number of conformers not defined");
  if(pars->n_qme == -1)  error_msg("Number of QM energies not defined");
  if(pars->n_qmf == -1)  error_msg("Number of QM forces not defined");

  if(pars->n_generations_max  == -1) {
    printf("# WARNING ... using default (100) for n_generations_max\n");
    pars->n_generations_max=100;
  }
  if(pars->n_individuals_max  == -1) {
    printf("# WARNING ... using default (200) for n_individuals_max\n");
    pars->n_individuals_max=200;
  }
  if(pars->initial_mutrate    == -1.0) {
    printf("# WARNING ... using default (0.1) for initial_mutrate\n");
    pars->initial_mutrate=0.01;
  }
  if(pars->codon_mutrate      == -1.0) {
    printf("# WARNING ... using default (0.006) for codon_mutrate\n");
    pars->codon_mutrate=0.006;
  }
  if(pars->individual_mutrate == -1.0) {
    printf("# WARNING ... using default (0.2) for individual_mutrate\n");
    pars->individual_mutrate=0.2;
  }
  if(pars->p_crossover        == -1.0) {
    printf("# WARNING ... using default (0.05) for p_crossover\n");
    pars->p_crossover=0.05;
  }
  if(pars->f_elite            == -1.0) {
    printf("# WARNING ... using default (0.2) for f_elite\n");
    pars->f_elite=0.2;
  }
  if(pars->f_lucky            == -1.0) {
    printf("# WARNING ... using default (0.05) for f_lucky\n");
    pars->f_lucky=0.05;
  }
  if(pars->rseed              == -1) {
    unsigned int localrseed = (unsigned) time(&t);
    printf("# WARNING ... using time (%u) for random seed\n", localrseed);
    pars->rseed=localrseed;
  }
  if(pars->vdw_cutoff         == -1.0) {
    printf("# WARNING ... using default (12.0 Angstrom) for vdw_cutoff\n");
    pars->vdw_cutoff         = 1.2;
  }

  printf("#\n# This is FOGAL version 21\n");
  printf("#\n");
  if(PFILE!=0) printf("# read %d parameters from command line and %d from %s\n", irc, irf, ARGV[PFILE]);
  else         printf("# read %d parameters from command line and zero from input parameter file\n", irc);
  printf("#\n");
  printf("# run id:                        %s\n",      pars->runid);
  printf("# opt string:                    %s\n",      pars->opt);
  printf("# H-bond contraints:             %d\n",      pars->hcons);

  printf("# Total number of conformers:    %12d\n",    pars->n_conf);
  printf("# Total number of generations:   %12d\n",    pars->n_generations_max);
  printf("# Total number of individuals:   %12d\n",    pars->n_individuals_max);

  printf("# Initial mutation rate:         %12.4lf\n", pars->initial_mutrate);
  printf("# Muation rate/individual:       %12.4lf\n", pars->individual_mutrate);
  printf("# Muation rate/codon:            %12.4lf\n", pars->codon_mutrate);
  printf("# Crossover probability:         %12.4lf\n", pars->p_crossover);
  printf("# Elite fraction:                %12.4lf\n", pars->f_elite);
  printf("# Lucky fraction:                %12.4lf\n", pars->f_lucky);

  printf("# Input topology file:           %s\n", pars->topology_filename);
  printf("# Input QME file:                %s\n", pars->qme_filename);
  printf("# Input coordinates file:        %s\n", pars->xyz_filename);
  printf("# Number of molecules:           %d\n", pars->n_mol);
  printf("# Number QM energies:            %d\n", pars->n_qme);
  printf("# Number QM forces:              %d\n", pars->n_qmf);

  printf("# Random seed:                   %12d\n",pars->rseed);
  printf("#\n");

  fflush(stdout);

  return(irc+irf);
}
