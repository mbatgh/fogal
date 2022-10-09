

  typedef struct {

    int n_generations_max;
    int n_individuals_max;
    int n_eesimplex_max;

    double vdw_cutoff;
    double initial_mutrate;
    double codon_mutrate;
    double individual_mutrate;
    double p_crossover;
    double f_elite;
    double f_lucky;
    double stdeve;
    double stdevf;

    int n_mol;
    int n_qme;
    int n_qmf;
    int n_conf;
    int hcons;
    int rseed;

    char* topology_filename;
    char* qme_filename;
    char* xyz_filename;

    char* runid;
    char* opt;

  } t_parameters;


  typedef struct {

    int n_conf;
    int n_types;
    int n_atoms;
    int n_pairs;
    int n_bonds;
    int n_angles;
    int n_dihedrals;
    int n_impropers;
    int n_nb;
   
    char**  atomtype;
    double* sigma;
    double* epsilon;
    double* sigma0;
    double* epsilon0;
    int*    egrpidx;   
    int*    sgrpidx;

    char**  nbat1;
    char**  nbat2;
    int*    nbidx;
    char    resname[4];
    char**  type;
    char**  atomname;
    double* charge;
    double* mass;
    int*    ndex;
    int*    ish;

    int*    bondi;
    int*    bondj;
    double* bondr;
    double* bondk;
    int*    bgrpidx;
    int*    isconstrained;

    int*    pairi;
    int*    pairj;

    int*    anglei;
    int*    anglej;
    int*    anglek;
    double* angleeq;
    double* anglefc;
    int*    agrpidx;
    int*    cgrpidx;
    
    int*    dihi;
    int*    dihj;
    int*    dihk;
    int*    dihl;
    double* dihphase;
    double* dihfc;
    double* dihpn;
    int*    dgrpidx;   
    
    int*    impi;
    int*    impj;
    int*    impk;
    int*    impl;
    double* impphase;
    double* impfc;
    double* imppn;
    int*    igrpidx;   

    int*    neighbor1;
    int*    neighbor2;
    int     n_neighbors;

    double* nb_el_fac;
    double* nb_lj_fac;

    double* lj12;
    double* lj6;

  } t_topology;

  typedef struct {

    int length;
    double fitness;
    int index;
    double* codon;

  } t_genotype;

  typedef struct {

    double* lower;
    double* mean;
    double* upper;
    int* codonidx;
    int* topolidx;

  } t_codonlimits;

