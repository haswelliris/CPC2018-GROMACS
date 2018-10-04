#ifndef __83knmYT___TDK9dnsa9dw9NUE_TEW5g5SJ8__
#define __83knmYT___TDK9dnsa9dw9NUE_TEW5g5SJ8__

#define GMX_DOUBLE
#undef GMX_DOUBLE

/* With CPU kernels the i-cluster size is always 4 atoms.
 * With x86 SIMD the j-cluster size can be 2, 4 or 8, otherwise 4.
 */
#define NBNXN_CPU_CLUSTER_I_SIZE       4

#define NBNXN_CPU_CLUSTER_I_SIZE_2LOG  2

#define D_BOX_Z 1
#define D_BOX_Y 1
#define D_BOX_X 2
#define N_BOX_Z (2*D_BOX_Z+1)
#define N_BOX_Y (2*D_BOX_Y+1)
#define N_BOX_X (2*D_BOX_X+1)
#define N_IVEC  (N_BOX_Z*N_BOX_Y*N_BOX_X)
#define CENTRAL (N_IVEC/2)
#define SHIFTS  N_IVEC

#ifndef FALSE
#      define  FALSE   (0)
#endif
#ifndef TRUE
#      define  TRUE    (1)
#endif

/* To avoid NaN when excluded atoms are at zero distance, we add a small
 * number to r^2. NBNXN_AVOID_SING_R2_INC^-3 should fit in real.
 */
#ifndef GMX_DOUBLE
#define NBNXN_AVOID_SING_R2_INC  1.0e-12f
#else
/* The double prec. x86 SIMD kernels use a single prec. invsqrt, so > 1e-38 */
#define NBNXN_AVOID_SING_R2_INC  1.0e-36
#endif

typedef int gmx_bool;

typedef struct tMPI_Atomic
{
    volatile int value;
}
tMPI_Atomic_t;

typedef unsigned int gmx_bitmask_t;

typedef unsigned long t_excl;

typedef struct
{
    int             igeometry;    /* The type of list (atom, water, etc.)  */
    int             ielec;        /* Coulomb loop type index for kernels   */
    int             ielecmod;     /* Coulomb modifier (e.g. switch/shift)  */
    int             ivdw;         /* VdW loop type index for kernels       */
    int             ivdwmod;      /* VdW modifier (e.g. switch/shift)      */
    int             type;         /* Type of interaction, listed in
                                     gmx_nblist_interaction_type           */

    int             nri, maxnri;  /* Current/max number of i particles	   */
    int             nrj, maxnrj;  /* Current/max number of j particles	   */
    int *           iinr;         /* The i-elements                        */
    int *           iinr_end;     /* The end atom, only with enlistCG      */
    int *           gid;          /* Index in energy arrays                */
    int *           shift;        /* Shift vector index                    */
    int *           jindex;       /* Index in jjnr                         */
    int *           jjnr;         /* The j-atom list                       */
    int *           jjnr_end;     /* The end atom, only with enltypeCG     */
    char *          excl_fep;     /* Exclusions for FEP with Verlet scheme */
    t_excl *        excl;         /* Exclusions, only with enltypeCG       */

    /* We use separate pointers for kernels that compute both potential
     * and force (vf suffix), only potential (v) or only force (f)
     */
    void *          kernelptr_vf;
    void *          kernelptr_v;
    void *          kernelptr_f;

    /* Pad the list of neighbors for each i atom with "-1" entries up to the
     * simd_padding_width, if it is larger than 0. This is necessary for many
     * accelerated kernels using single-instruction multiple-data operations
     * internally.
     */
    int             simd_padding_width;

} t_nblist;

// =======================================================================================
//                                     GMX DEFs
// =======================================================================================


/*! \brief Double precision accuracy */
#define GMX_DOUBLE_EPS   2.2204460492503131e-16

/*! \brief Maximum double precision value - reduced 1 unit in last digit for MSVC */
#define GMX_DOUBLE_MAX   1.7976931348623157e+308

/*! \brief Minimum double precision value */
#define GMX_DOUBLE_MIN   2.2250738585072014e-308

/*! \brief Single precision accuracy */
#define GMX_FLOAT_EPS    1.19209290e-07F

/*! \brief Maximum single precision value - reduced 1 unit in last digit for MSVC */
#define GMX_FLOAT_MAX    3.40282346E+38F

/*! \brief Minimum single precision value */
#define GMX_FLOAT_MIN    1.175494351E-38F

#ifdef __PGI
/* The portland group x86 C/C++ compilers do not treat negative zero initializers
 * correctly, but "optimizes" them to positive zero, so we implement it explicitly.
 * These constructs are optimized to simple loads at compile time. If you want to
 * use them on other compilers those have to support gcc preprocessor extensions.
 * Note: These initializers might be sensitive to the endianness (which can
 * be different for byte and word order), so check that it works for your platform
 * and add a separate section if necessary before adding to the ifdef above.
 */
#    define GMX_DOUBLE_NEGZERO  ({ const union { int  di[2]; double d; } _gmx_dzero = {0, -2147483648}; _gmx_dzero.d; })
#    define GMX_FLOAT_NEGZERO   ({ const union { int  fi; float f; } _gmx_fzero = {-2147483648}; _gmx_fzero.f; })
#else
/*! \brief Negative zero in double */
#    define GMX_DOUBLE_NEGZERO  (-0.0)

/*! \brief Negative zero in float */
#    define GMX_FLOAT_NEGZERO   (-0.0f)
#endif

/*! \typedef real
 * \brief Precision-dependent \Gromacs floating-point type.
 */
/*! \def HAVE_REAL
 * \brief Used to check whether `real` is already defined.
 */
/*! \def GMX_MPI_REAL
 * \brief MPI data type for `real`.
 */
/*! \def GMX_REAL_EPS
 * \brief Accuracy for `real`.
 */
/*! \def GMX_REAL_MIN
 * \brief Smallest non-zero value for `real`.
 */
/*! \def GMX_REAL_MAX
 * \brief Largest finite value for `real`.
 */
/*! \def GMX_REAL_NEGZERO
 * \brief Negative zero for `real`.
 */
/*! \def gmx_real_fullprecision_pfmt
 * \brief Format string for full `real` precision.
 */
#ifdef GMX_DOUBLE

#ifndef HAVE_REAL
typedef double      real;
#define HAVE_REAL
#endif

#define GMX_MPI_REAL      MPI_DOUBLE
#define GMX_REAL_EPS      GMX_DOUBLE_EPS
#define GMX_REAL_MIN      GMX_DOUBLE_MIN
#define GMX_REAL_MAX      GMX_DOUBLE_MAX
#define GMX_REAL_NEGZERO  GMX_DOUBLE_NEGZERO
#define gmx_real_fullprecision_pfmt "%21.14e"

#define gmx_invsqrt(x) (1./sqrt(x))

#else /* GMX_DOUBLE */

#ifndef HAVE_REAL
typedef float           real;
#define HAVE_REAL
#endif

#define GMX_MPI_REAL      MPI_FLOAT
#define GMX_REAL_EPS      GMX_FLOAT_EPS
#define GMX_REAL_MIN      GMX_FLOAT_MIN
#define GMX_REAL_MAX      GMX_FLOAT_MAX
#define GMX_REAL_NEGZERO  GMX_FLOAT_NEGZERO
#define gmx_real_fullprecision_pfmt "%14.7e"

#define gmx_invsqrt(x) (1./sqrtf(x))

#endif /* GMX_DOUBLE */

#define XX      0 /* Defines for indexing in */
#define YY      1 /* vectors                 */
#define ZZ      2
#define DIM     3 /* Dimension of vectors    */

typedef real    rvec[DIM];

/* Used with force switching or a constant potential shift:
 * rsw       = max(r - r_switch, 0)
 * force/p   = r^-(p+1) + c2*rsw^2 + c3*rsw^3
 * potential = r^-p + c2/3*rsw^3 + c3/4*rsw^4 + cpot
 * With a constant potential shift c2 and c3 are both 0.
 */
typedef struct {
    real c2;
    real c3;
    real cpot;
} shift_consts_t;

/* Used with potential switching:
 * rsw        = max(r - r_switch, 0)
 * sw         = 1 + c3*rsw^3 + c4*rsw^4 + c5*rsw^5
 * dsw        = 3*c3*rsw^2 + 4*c4*rsw^3 + 5*c5*rsw^4
 * force      = force*dsw - potential*sw
 * potential *= sw
 */
typedef struct {
    real c3;
    real c4;
    real c5;
} switch_consts_t;

typedef struct {
    int             cutoff_scheme;

    /* VdW */
    int             vdwtype;
    int             vdw_modifier;
    real            rvdw;
    real            rvdw_switch;
    shift_consts_t  dispersion_shift;
    shift_consts_t  repulsion_shift;
    switch_consts_t vdw_switch;
    /* TODO: remove this variable, used for not modyfing the group kernels,
     * it is equal to -dispersion_shift->cpot
     */
    real sh_invrc6;

    /* type of electrostatics (defined in enums.h) */
    int  eeltype;
    int  coulomb_modifier;

    /* Coulomb */
    real rcoulomb;

    /* Cut-off */
    real rlist;
    real rlistlong;

    /* PME/Ewald */
    real ewaldcoeff_q;
    real ewaldcoeff_lj;
    int  ljpme_comb_rule; /* LJ combination rule for the LJ PME mesh part */
    real sh_ewald;        /* -sh_ewald is added to the direct space potential */
    real sh_lj_ewald;     /* sh_lj_ewald is added to the correction potential */

    /* Dielectric constant resp. multiplication factor for charges */
    real epsilon_r;
    real epsfac;

    /* Constants for reaction-field or plain cut-off */
    real epsilon_rf;
    real k_rf;
    real c_rf;

    /* Force/energy interpolation tables, linear in force, quadratic in V */
    real  tabq_scale;
    int   tabq_size;
    /* Coulomb force table, size of array is tabq_size (when used) */
    real *tabq_coul_F;
    /* Coulomb energy table, size of array is tabq_size (when used) */
    real *tabq_coul_V;
    /* Coulomb force+energy table, size of array is tabq_size*4,
       entry quadruplets are: F[i], F[i+1]-F[i], V[i], 0,
       this is used with single precision x86 SIMD for aligned loads */
    real *tabq_coul_FDV0;

    /* Vdw force table for LJ-PME, size of array is tabq_size (when used) */
    real *tabq_vdw_F;
    /* Vdw energy table for LJ-PME, size of array is tabq_size (when used) */
    real *tabq_vdw_V;
    /* Vdw force+energy table for LJ-PME, size of array is tabq_size*4, entry
       quadruplets are: F[i], F[i+1]-F[i], V[i], 0, this is used with
       single precision x86 SIMD for aligned loads */
    real *tabq_vdw_FDV0;

} interaction_const_t;

/* A buffer data structure of 64 bytes
 * to be placed at the beginning and end of structs
 * to avoid cache invalidation of the real contents
 * of the struct by writes to neighboring memory.
 */
typedef struct {
    int dummy[16];
} gmx_cache_protect_t;

/* Abstract type for pair searching data */
typedef struct nbnxn_search * nbnxn_search_t;

/* Function that should return a pointer *ptr to memory
 * of size nbytes.
 * Error handling should be done within this function.
 */
typedef void nbnxn_alloc_t (void **ptr, size_t nbytes);

/* Function that should free the memory pointed to by *ptr.
 * NULL should not be passed to this function.
 */
typedef void nbnxn_free_t (void *ptr);

/* This is the actual cluster-pair list j-entry.
 * cj is the j-cluster.
 * The interaction bits in excl are indexed i-major, j-minor.
 * The cj entries are sorted such that ones with exclusions come first.
 * This means that once a full mask (=NBNXN_INTERACTION_MASK_ALL)
 * is found, all subsequent j-entries in the i-entry also have full masks.
 */
typedef struct {
    int          cj;    /* The j-cluster                    */
    unsigned int excl;  /* The exclusion (interaction) bits */
} nbnxn_cj_t;

/* In nbnxn_ci_t the integer shift contains the shift in the lower 7 bits.
 * The upper bits contain information for non-bonded kernel optimization.
 * Simply calculating LJ and Coulomb for all pairs in a cluster pair is fine.
 * But three flags can be used to skip interactions, currently only for subc=0
 * !(shift & NBNXN_CI_DO_LJ(subc))   => we can skip LJ for all pairs
 * shift & NBNXN_CI_HALF_LJ(subc)    => we can skip LJ for the second half of i
 * !(shift & NBNXN_CI_DO_COUL(subc)) => we can skip Coulomb for all pairs
 */
#define NBNXN_CI_SHIFT          127
#define NBNXN_CI_DO_LJ(subc)    (1<<(7+3*(subc)))
#define NBNXN_CI_HALF_LJ(subc)  (1<<(8+3*(subc)))
#define NBNXN_CI_DO_COUL(subc)  (1<<(9+3*(subc)))

/* Simple pair-list i-unit */
typedef struct {
    int ci;             /* i-cluster             */
    int shift;          /* Shift vector index plus possible flags, see above */
    int cj_ind_start;   /* Start index into cj   */
    int cj_ind_end;     /* End index into cj     */
} nbnxn_ci_t;

/* Grouped pair-list i-unit */
typedef struct {
    int sci;            /* i-super-cluster       */
    int shift;          /* Shift vector index plus possible flags */
    int cj4_ind_start;  /* Start index into cj4  */
    int cj4_ind_end;    /* End index into cj4    */
} nbnxn_sci_t;

typedef struct {
    unsigned int imask;    /* The i-cluster interactions mask for 1 warp  */
    int          excl_ind; /* Index into the exclusion array for 1 warp   */
} nbnxn_im_ei_t;

typedef struct {
    int           cj[4];   /* The 4 j-clusters                            */
    nbnxn_im_ei_t imei[2]; /* The i-cluster mask data       for 2 warps   */
} nbnxn_cj4_t;

typedef struct {
    unsigned int pair[32]; /* Topology exclusion interaction bits for one warp,
                            * each unsigned has bitS for 4*8 i clusters
                            */
} nbnxn_excl_t;

typedef struct nbnxn_pairlist_t {
    gmx_cache_protect_t cp0;

    nbnxn_alloc_t      *alloc;
    nbnxn_free_t       *free;

    gmx_bool            bSimple;         /* Simple list has na_sc=na_s and uses cj   *
                                          * Complex list uses cj4                    */

    int                     na_ci;       /* The number of atoms per i-cluster        */
    int                     na_cj;       /* The number of atoms per j-cluster        */
    int                     na_sc;       /* The number of atoms per super cluster    */
    real                    rlist;       /* The radius for constructing the list     */
    int                     nci;         /* The number of i-clusters in the list     */
    nbnxn_ci_t             *ci;          /* The i-cluster list, size nci             */
    int                     ci_nalloc;   /* The allocation size of ci                */
    int                     nsci;        /* The number of i-super-clusters in the list */
    nbnxn_sci_t            *sci;         /* The i-super-cluster list                 */
    int                     sci_nalloc;  /* The allocation size of sci               */

    int                     ncj;         /* The number of j-clusters in the list     */
    nbnxn_cj_t             *cj;          /* The j-cluster list, size ncj             */
    int                     cj_nalloc;   /* The allocation size of cj                */

    int                     ncj4;        /* The total number of 4*j clusters         */
    nbnxn_cj4_t            *cj4;         /* The 4*j cluster list, size ncj4          */
    int                     cj4_nalloc;  /* The allocation size of cj4               */
    int                     nexcl;       /* The count for excl                       */
    nbnxn_excl_t           *excl;        /* Atom interaction bits (non-exclusions)   */
    int                     excl_nalloc; /* The allocation size for excl             */
    int                     nci_tot;     /* The total number of i clusters           */

    struct nbnxn_list_work *work;

    gmx_cache_protect_t     cp1;
} nbnxn_pairlist_t;

typedef struct {
    int                nnbl;        /* number of lists */
    nbnxn_pairlist_t **nbl;         /* lists */
    gmx_bool           bCombined;   /* TRUE if lists get combined into one (the 1st) */
    gmx_bool           bSimple;     /* TRUE if the list of of type "simple"
                                       (na_sc=na_s, no super-clusters used) */
    int                natpair_ljq; /* Total number of atom pairs for LJ+Q kernel */
    int                natpair_lj;  /* Total number of atom pairs for LJ kernel   */
    int                natpair_q;   /* Total number of atom pairs for Q kernel    */
    t_nblist         **nbl_fep;
} nbnxn_pairlist_set_t;

enum {
    nbatXYZ, nbatXYZQ, nbatX4, nbatX8
};

typedef struct {
    real *f;      /* f, size natoms*fstride                             */
    real *fshift; /* Shift force array, size SHIFTS*DIM                 */
    int   nV;     /* The size of *Vvdw and *Vc                          */
    real *Vvdw;   /* Temporary Van der Waals group energy storage       */
    real *Vc;     /* Temporary Coulomb group energy storage             */
    int   nVS;    /* The size of *VSvdw and *VSc                        */
    real *VSvdw;  /* Temporary SIMD Van der Waals group energy storage  */
    real *VSc;    /* Temporary SIMD Coulomb group energy storage        */
} nbnxn_atomdata_output_t;

/* Block size in atoms for the non-bonded thread force-buffer reduction,
 * should be a multiple of all cell and x86 SIMD sizes (i.e. 2, 4 and 8).
 * Should be small to reduce the reduction and zeroing cost,
 * but too small will result in overhead.
 * Currently the block size is NBNXN_BUFFERFLAG_SIZE*3*sizeof(real)=192 bytes.
 */
#ifdef GMX_DOUBLE
#define NBNXN_BUFFERFLAG_SIZE   8
#else
#define NBNXN_BUFFERFLAG_SIZE  16
#endif

/* We store the reduction flags as gmx_bitmask_t.
 * This limits the number of flags to BITMASK_SIZE.
 */
#define NBNXN_BUFFERFLAG_MAX_THREADS  (BITMASK_SIZE)

/* Flags for telling if threads write to force output buffers */
typedef struct {
    int               nflag;       /* The number of flag blocks                         */
    gmx_bitmask_t    *flag;        /* Bit i is set when thread i writes to a cell-block */
    int               flag_nalloc; /* Allocation size of cxy_flag                       */
} nbnxn_buffer_flags_t;

/* LJ combination rules: geometric, Lorentz-Berthelot, none */
enum {
    ljcrGEOM, ljcrLB, ljcrNONE, ljcrNR
};

typedef struct nbnxn_atomdata_t {
    nbnxn_alloc_t           *alloc;
    nbnxn_free_t            *free;
    int                      ntype;           /* The number of different atom types                 */
    real                    *nbfp;            /* Lennard-Jones 6*C6 and 12*C12 params, size ntype^2*2 */
    int                      comb_rule;       /* Combination rule, see enum above                   */
    real                    *nbfp_comb;       /* LJ parameter per atom type, size ntype*2           */
    real                    *nbfp_s4;         /* As nbfp, but with stride 4, size ntype^2*4. This
                                               * might suit 4-wide SIMD loads of two values (e.g.
                                               * two floats in single precision on x86).            */
    int                      natoms;          /* Number of atoms                                    */
    int                      natoms_local;    /* Number of local atoms                           */
    int                     *type;            /* Atom types                                         */
    real                    *lj_comb;         /* LJ parameters per atom for combining for pairs     */
    int                      XFormat;         /* The format of x (and q), enum                      */
    int                      FFormat;         /* The format of f, enum                              */
    real                    *q;               /* Charges, can be NULL if incorporated in x          */
    int                      na_c;            /* The number of atoms per cluster                    */
    int                      nenergrp;        /* The number of energy groups                        */
    int                      neg_2log;        /* Log2 of nenergrp                                   */
    int                     *energrp;         /* The energy groups per cluster, can be NULL         */
    gmx_bool                 bDynamicBox;     /* Do we need to update shift_vec every step?    */
    rvec                    *shift_vec;       /* Shift vectors, copied from t_forcerec              */
    int                      xstride;         /* stride for a coordinate in x (usually 3 or 4)      */
    int                      fstride;         /* stride for a coordinate in f (usually 3 or 4)      */
    real                    *x;               /* x and possibly q, size natoms*xstride              */

    /* j-atom minus i-atom index for generating self and Newton exclusions
     * cluster-cluster pairs of the diagonal, for 4xn and 2xnn kernels.
     */
    real                    *simd_4xn_diagonal_j_minus_i;
    real                    *simd_2xnn_diagonal_j_minus_i;
    /* Filters for topology exclusion masks for the SIMD kernels.
     * filter2 is the same as filter1, but with each element duplicated.
     */
    unsigned int            *simd_exclusion_filter1;
    unsigned int            *simd_exclusion_filter2;
    real                    *simd_interaction_array; /* Array of masks needed for exclusions on QPX */
    int                      nout;                   /* The number of force arrays                         */
    nbnxn_atomdata_output_t *out;                    /* Output data structures               */
    int                      nalloc;                 /* Allocation size of all arrays (for x/f *x/fstride) */
    gmx_bool                 bUseBufferFlags;        /* Use the flags or operate on all atoms     */
    nbnxn_buffer_flags_t     buffer_flags;           /* Flags for buffer zeroing+reduc.  */
    gmx_bool                 bUseTreeReduce;         /* Use tree for force reduction */
    tMPI_Atomic_t           *syncStep;               /* Synchronization step for tree reduce */
} nbnxn_atomdata_t;

#endif
