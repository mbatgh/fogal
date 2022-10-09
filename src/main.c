
/* ================================================================= */
/* === force field optimization with genetic algoriem, author MB === */
/* ================================================================= */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include "strucs.h"
#include "functions.h"

/* ===== global constants ============================================================= */

const double FRACDEPSMAX=2.0;
const double FRACDEPSMIN=0.5;
const double FRACDSIGMAX=1.25;
const double FRACDSIGMIN=0.75;

const double DANGMAX=5.0;
const double BONDFCMIN=1.5e+05;
const double BONDFCMAX=7.0e+05;
const double ANGLEFCMAX=900.0;
const double DIHEDRALFCMAX=40.0;
const double DIHEDRALFCMIN=0.0;
const double IMPROPERFCMAX=50.0;

const int XX=0;
const int YY=1;
const int ZZ=2;

/* ===== function declarations ============================================================= */

int get_parameters(int,char*[],t_parameters*);
int allocate_memory(t_parameters*,t_topology*,double***,double**,double**,double**,double**);
int read_topology(const t_parameters*,t_topology*);

int read_coords(t_parameters*,t_topology*,double**);
int read_qme(t_parameters*,t_topology*,double*,double*);

int initiate_genotypes(t_topology*,t_genotype**,t_genotype**,t_codonlimits*,t_parameters*);
int mutate(t_genotype*,t_codonlimits*,double);
int mate(t_genotype*,t_genotype*,t_genotype*,int,double);

int fitness (const t_parameters*,t_genotype*,t_topology*,t_codonlimits*,double**,double*,double*,double*,double*);
int write_top(const char*,const t_parameters*,t_topology*,t_genotype*, t_codonlimits* limits);
int copy_gene(t_genotype*,t_genotype*);
int setup_mm(const char*,const t_parameters*,t_genotype*,t_topology*,t_codonlimits*,double**,double*,double*,double*,double*);
int write_mme(const char*,const t_parameters*,t_genotype*,t_topology*,t_codonlimits*,double**,double*,double*,double*,double*);

double randn(double,double);
void error_msg(char*);


/* ===== main ============================================================================== */

int main(int ARGC, char* ARGV[]) {

/***** local variables *****/

  int i, i1, i2, p1, p2;
  int n_generation;
  int n_individual;
  int n_codons;
  int nread;
  int iret;
  int nmutpi;
  int tmpnmut;
  int nmut0=0;
  int nmut=0;
  int nrerun=0;
  int nmate1;
  int nmate2;
  int nrest;
  int reti;

  double frand;

  char tmpfilnam[256];

  int* partner_array;

  double* qme;
  double* w_qme;
  double* mme;
  double* mme0;

  double eltime;

  t_parameters    pars;
  t_topology      input_topology;
  t_genotype*     parent;
  t_genotype*     child;

  t_codonlimits   limits;

  double**        coords;

  struct timespec start, stop;

/***** record start time *****/

  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

/***** write version info *****/

  printf("\n# FOGAL version 21\n");

/***** get parameters from command line or file *****/

  reti = get_parameters(ARGC, ARGV, &pars);
  if(reti>0) printf("# returned from get_parameters, obtained %d parameter values\n", reti);
  else error_msg("Something went wrong reading input parameters");

/***** allocate memory *****/

  reti = allocate_memory(&pars,&input_topology,&coords,&qme,&w_qme,&mme,&mme0);
  if(reti<1) error_msg("Something went wrong allocating memory");

/***** read topology *****/

  reti = read_topology(&pars,&input_topology);
  if(reti!=0) error_msg("Something went wrong reading the topology");

/***** read coordinate file *****/

  reti = read_coords(&pars,&input_topology,coords);
  if(reti!=0) error_msg("Something went wrong reading the coordinates");

/***** read qm energies *****/

  reti = read_qme(&pars,&input_topology,qme,w_qme);
  if(reti!=0) error_msg("Problems reading QME file");

/***** initiate individuals/genotypes *****/

  n_codons = initiate_genotypes(&input_topology,&parent,&child,&limits,&pars);
  if(n_codons<1) error_msg("Problems initiating genotypes");

  printf("# initiated %d codons\n", n_codons);

/***** seed random number generator *****/

  srand(pars.rseed);
  printf("# Seeding RNG using %d\n", pars.rseed);

/***** determine fitness of input topology and write conformational energies qme vs gmx *****/

  sprintf(tmpfilnam,"%s-ef-init", pars.runid);
  reti = setup_mm(tmpfilnam,&pars,&(parent[0]),&(input_topology), &limits, coords, qme, w_qme, mme, mme0);

  printf("%12.2lf %16.1lf %10d %16.6lf\n", 0.0, 0.0, 0, parent[0].fitness);

  fflush(stdout);

/***** apply random mutations to the genome of each individual *****/

  for(n_individual=0; n_individual<pars.n_individuals_max; n_individual++) {
    nmut0+=mutate(&(parent[n_individual]),&limits,pars.initial_mutrate);
    nread = fitness(&pars, &(parent[n_individual]), &input_topology, &limits, coords, qme, w_qme, mme, mme0);
    parent[n_individual].index = n_individual;
    nrerun++;
  }

  qsort(parent, pars.n_individuals_max, sizeof(parent[0]), compare_genomes);

  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  eltime = (double)(stop.tv_sec - start.tv_sec)/60.0;

  printf("# Created %d parents (generation zero, initial mutation rate=%lf)\n", pars.n_individuals_max, pars.initial_mutrate);
  printf("# Optimal fitness is:\n");

  printf("%12.2lf %16.1lf %10d %16.6lf\n", 0.0, eltime, nrerun, parent[0].fitness);
  printf("#\n");

  fflush(stdout);

/***** allocate mem *****/

  partner_array = calloc(pars.n_individuals_max, sizeof(int));


/*###############################################################################/
/##### start main loop over generations #########################################/
/###############################################################################*/

  for(n_generation=1;n_generation<=pars.n_generations_max;n_generation++) {

    printf("# GENERATION %6d ===============================\n", n_generation);

/***** do sex, mutate, and procreate *****/

    nmate1=(int)(pars.f_elite*pars.n_individuals_max);
    nmate2=nmate1;
    nrest=pars.n_individuals_max-nmate1;

    for(i=0;i<nmate1;i++) partner_array[i]=i;

    for(i=0;i<nrest;i++) {
      frand = ((double)rand()/RAND_MAX);
      if(frand<pars.f_lucky) {
        i1=(int)(nrest*frand)+nmate1;
        partner_array[nmate2]=i1;
        nmate2++;
      }
    }

    iret=0;
    nmutpi=0;
    nmut=0;

    for(n_individual=0;n_individual<pars.n_individuals_max;n_individual++) {

      i1=(int)( (double)nmate2*(double)rand()/RAND_MAX );
      i2=(int)( (double)nmate2*(double)rand()/RAND_MAX );

      p1=partner_array[i1];
      p2=partner_array[i2];

      iret += mate(&(parent[p1]), &(parent[p2]), &(child[n_individual]), n_individual, pars.p_crossover);

      frand=((double)rand()/RAND_MAX);
      if(frand<pars.individual_mutrate) {
        tmpnmut=mutate(&(child[n_individual]),&limits,pars.codon_mutrate);
        if(tmpnmut>0) {
          nmutpi+=tmpnmut;
          nmut++;
        }
      }

      nread = fitness(&pars, &(child[n_individual]), &input_topology, &limits, coords, qme, w_qme, mme, mme0);

      nrerun++;
    }

    printf("# Mate %d pairings between individuals (<#co>=%.2lf) from a subpopulation of %d (%d lucky)\n",
            pars.n_individuals_max, (double)iret/pars.n_individuals_max, nmate2, nmate2-nmate1);
    printf("# Chose %d out of %d individuals for mutation\n", nmut, pars.n_individuals_max);
    if(nmut>0) printf("# Mutated on average %.2f out of %d codons in each of these individuals\n", (double)nmutpi/nmut, n_codons);
    else       printf("# Mutated on average %.2f out of %d codons in each of these individuals\n", 0.0, n_codons);

    qsort(child, pars.n_individuals_max, sizeof(child[0]), compare_genomes);

    for(i=0;i<pars.n_individuals_max;i++) nread=copy_gene(&(child[i]), &(parent[i]));

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
    eltime = (double)(stop.tv_sec - start.tv_sec)/60.0;

    printf("# generation        real-time     reruns Best-fitness:\n");
    printf("%12.2lf %16.1lf %10d %16.6lf\n", (double)n_generation, eltime, nrerun, parent[0].fitness);
    fflush(stdout);
/*
    sprintf(tmpfilnam,"%s-intermed-%06d.top", pars.runid, n_generation);
    nread=write_top(tmpfilnam,&pars,&(input_topology),&(parent[0]),&limits);
*/
  }

/*###############################################################################/
/##### end of main loop over generations ########################################/
/###############################################################################*/

/***** do housekeeping *****/

  if(pars.n_generations_max>0) {

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
    eltime = (double)(stop.tv_sec - start.tv_sec)/60.0;

    printf("# evolution of %d generations of %d individuals each with %d codons:\n",
           pars.n_generations_max, pars.n_individuals_max, parent[0].length);
    printf("# total real time elapsed for GA: %12.2lf minutes\n", eltime);

    sprintf(tmpfilnam,"%s-ga.top", pars.runid);
    nread=write_top(tmpfilnam,&pars,&(input_topology),&(parent[0]),&limits);

    sprintf(tmpfilnam,"%s-ef-ga", pars.runid);
    nread=write_mme(tmpfilnam,&(pars),&(parent[0]),&(input_topology),&limits,coords,qme,w_qme,mme,mme0);

    printf("# saved best topology in %s-ga.top and energies from QM vs gmx in e-%s-ga\n", pars.runid, pars.runid);
  }

  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  eltime = (double)(stop.tv_sec - start.tv_sec)/60.0;

  printf("# Evolution of %d generations of %d individuals each with %d codons\n",
         pars.n_generations_max, pars.n_individuals_max, parent[0].length);
  printf("# (%d fitness evaluations)\n",nrerun);
  printf("# total real time elapsed: %12.2lf mins\n", eltime);
  
  exit(1);
}
