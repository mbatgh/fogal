
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "strucs.h"
#include "functions.h"



/* =========================================================================================== */
void error_msg(char* msg) {

  fprintf(stderr,"\n%s\n\n", msg);
  fprintf(stderr,"usage: fogal [command-line parameters]\n\n");

  fprintf(stderr,"  -h          show this msg\n");

  fprintf(stderr,"  -z  string  parameter filename\n");
  fprintf(stderr,"  -o  string  type of parameters optimized [a]ngles,[d]ihedrals,[i]mpropers)\n");
  fprintf(stderr,"  -g  int     total number of generations\n");
  fprintf(stderr,"  -i  int     total number of individuals per generation\n");

  fprintf(stderr,"  -a  float   mutation rate applied to first generation population\n");
  fprintf(stderr,"  -m  float   mutation rate at the level of a single genome\n");
  fprintf(stderr,"  -u  float   mutation rate at the level of population\n");
  fprintf(stderr,"  -c  float   crossover probability (0-1)\n");
  fprintf(stderr,"  -f  float   fraction of individuals that always get a chance to mate\n");
  fprintf(stderr,"  -l  float   fraction of poor individuals that get a chance to mate\n");

  fprintf(stderr,"  -s  int     seed for random generator\n");

  fprintf(stderr,"  -w  float   cut-off radius for VdW interactions\n");

  fprintf(stderr,"  -p  string  name of input topology file\n");
  fprintf(stderr,"  -e  string  name of text file with QM energies\n");
  fprintf(stderr,"  -t  string  name of trajectory file with conformations\n");

  fprintf(stderr,"  -k  int     number of conformers\n");
  fprintf(stderr,"  -q  int     number of QM energies in efw file\n");
  fprintf(stderr,"  -n  int     number of molecules per conformer\n");
  fprintf(stderr,"  -v  int     number of forces in efw file\n");

  fprintf(stderr,"  -r  string  run id, for unique filenames\n");

  fprintf(stderr,"  -y          test: parameters are read and processed but optimization is not done\n");

  fprintf(stderr,"\n");

  exit(1);
}
        


/* =========================================================================================== */
int words(const char sentence[]) {

    int counted = 0;
    const char* it = sentence;
    int inword = 0;

    do switch(*it) {
        case '\0': 
        case ' ': case '\t': case '\n': case '\r':
            if (inword) { inword = 0; counted++; }
            break;
        default: inword = 1;
    } while(*it++);

    return counted;
}

/* =========================================================================================== */
double randn (double mu, double sigma) {

  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return(mu + sigma * (double) X1);
}


/* =========================================================================================== */
int compare_genomes(const void *a, const void *b){

    t_genotype* struct_a = (t_genotype*) a;
    t_genotype* struct_b = (t_genotype*) b;
    if (struct_a->fitness < struct_b->fitness) return -1;
    else if (struct_a->fitness == struct_b->fitness) return 0;
    else return 1;
}

/* =========================================================================================== */
int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}

/* =========================================================================================== */
void timeval_print(struct timeval *tv)
{
    char buffer[30];
    time_t curtime;

    printf("%ld.%06ld", tv->tv_sec, tv->tv_usec);
    curtime = tv->tv_sec;
    strftime(buffer, 30, "%m-%d-%Y  %T", localtime(&curtime));
    printf(" = %s.%06ld\n", buffer, tv->tv_usec);
}


