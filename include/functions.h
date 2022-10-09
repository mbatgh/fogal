

extern int get_parameters(int,char*[],t_parameters*);
extern int allocate_memory(t_parameters*,t_topology*,double***,double**,double**,double**,double**);
extern int read_topology(const t_parameters*,t_topology*);
extern int initiate_genotypes(t_topology*,t_genotype**,t_genotype**,t_codonlimits*,t_parameters*);
extern int write_top(const char*,const t_parameters*,t_topology*,t_genotype*, t_codonlimits* limits);
extern int mutate(t_genotype*,t_codonlimits*,double);
extern int copy_gene(t_genotype*,t_genotype*);
extern int mate(t_genotype*,t_genotype*,t_genotype*,int,double);
extern int fitness (const t_parameters*,t_genotype*,t_topology*,t_codonlimits*,double**,double*,double*,double*,double*);
extern int setup_mm(const char*,const t_parameters*,t_genotype*,t_topology*,t_codonlimits*,double**,double*,double*,double*,double*);
extern int read_qme(t_parameters*,t_topology*,double*,double*);
extern int read_coords(t_parameters*,t_topology*,double**);
extern int compare_genomes(const void*,const void*);
extern int write_mme(const char*,const t_parameters*,t_genotype*,t_topology*,t_codonlimits*,double**,double*,double*,double*,double*);


extern double randn(double,double);
extern int words(const char*);
extern void error_msg(char*);

extern void reset_input_mode(void);
extern void set_input_mode(void);
extern int ask_stop(void);
extern void presskey(void);
extern void mic_msg(char*);
extern char mic_getchar(char*);
extern void mabort(const char*);

extern double **dmat(int,int);
extern double *dvec(int);
extern float **fmat(int,int);
extern float *fvec(int);
extern int **imat(int,int);
extern char **cmat(int,int);
extern int *ivec(int);
extern char *cvec(int);
extern double **d0mat(int,int);
extern double *d0vec(int);
extern float **f0mat(int,int);
extern float *f0vec(int);
extern int **i0mat(int,int);
extern int *i0vec(int);
extern char **c0mat(int,int);
extern char *c0vec(int);
extern void waitsec(int);

extern const double BONDFCMIN;
extern const double BONDFCMAX;
extern const double ANGLEFCMAX;
extern const double DIHEDRALFCMAX;
extern const double DIHEDRALFCMIN;
extern const double IMPROPERFCMAX;

extern const double FRACDEPSMAX;
extern const double FRACDEPSMIN;
extern const double FRACDSIGMAX;
extern const double FRACDSIGMIN;
extern const double DANGMAX;

extern const int XX;
extern const int YY;
extern const int ZZ;

