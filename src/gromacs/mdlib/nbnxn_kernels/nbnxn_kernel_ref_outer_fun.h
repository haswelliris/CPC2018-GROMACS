// custom define
#define GMX_DOUBLE

#ifdef GMX_DOUBLE
	typedef double real;
#else /* GMX_DOUBLE */
	typedef float real;
#endif

// end of custom define

//-----------------------------------------------------------------------------------

// outer define

/* We could use nbat->xstride and nbat->fstride, but macros might be faster */
#define NBNXN_CPU_CLUSTER_I_SIZE	   4
#define UNROLLI	NBNXN_CPU_CLUSTER_I_SIZE
#define UNROLLJ	NBNXN_CPU_CLUSTER_I_SIZE

#define X_STRIDE   3
#define F_STRIDE   3
/* Local i-atom buffer strides */
#define XI_STRIDE  3
#define FI_STRIDE  3

#define CALC_SHIFTFORCES

//end of outer define

// ------------------------------------------------------------------------------------

// other define

#define D_BOX_Z 1
#define D_BOX_Y 1
#define D_BOX_X 2
#define N_BOX_Z (2*D_BOX_Z+1)
#define N_BOX_Y (2*D_BOX_Y+1)
#define N_BOX_X (2*D_BOX_X+1)
#define N_IVEC  (N_BOX_Z*N_BOX_Y*N_BOX_X)
#define CENTRAL (N_IVEC/2)


#define DIM     3 /* Dimension of vectors    */

#define FI_STRIDE  3

#ifndef GMX_DOUBLE
#define NBNXN_AVOID_SING_R2_INC  1.0e-12f
#else
/* The double prec. x86 SIMD kernels use a single prec. invsqrt, so > 1e-38 */
#define NBNXN_AVOID_SING_R2_INC  1.0e-36
#endif

#define NBNXN_CI_SHIFT          127

#define XX      0 /* Defines for indexing in */
#define YY      1 /* vectors                 */
#define ZZ      2

//-------------------------------------------------------------------------------------

enum {para_CALC_COUL_RF, para_CALC_COUL_TAB, para_CALC_ENERGIES, para_ENERGY_GROUPS, 
	para_LJ_CUT, para_LJ_EWALD, para_LJ_EWALD_COMB_GEOM, para_LJ_EWALD_COMB_LB, 
	para_LJ_FORCE_SWITCH, para_LJ_POT_SWITCH, para_VDW_CUTOFF_CHECK, 
	para_EXCL_FORCES, para_count
};

// ------------------------------------------------------------------------------------

// sub define for nbnxn_pairlist_t

typedef struct {
	int dummy[16];
} gmx_cache_protect_t;

typedef void nbnxn_alloc_t (void **ptr, size_t nbytes);

typedef void nbnxn_free_t (void *ptr);

typedef int gmx_bool;

/* Simple pair-list i-unit */
typedef struct {
	int ci;			 /* i-cluster			 */
	int shift;		  /* Shift vector index plus possible flags, see above */
	int cj_ind_start;   /* Start index into cj   */
	int cj_ind_end;	 /* End index into cj	 */
} nbnxn_ci_t;

/* Grouped pair-list i-unit */
typedef struct {
	int sci;			/* i-super-cluster	   */
	int shift;		  /* Shift vector index plus possible flags */
	int cj4_ind_start;  /* Start index into cj4  */
	int cj4_ind_end;	/* End index into cj4	*/
} nbnxn_sci_t;

typedef struct {
	int		  cj;	/* The j-cluster					*/
	unsigned int excl;  /* The exclusion (interaction) bits */
} nbnxn_cj_t;

typedef struct {
	int		   cj[4];   /* The 4 j-clusters							*/
	nbnxn_im_ei_t imei[2]; /* The i-cluster mask data	   for 2 warps   */
} nbnxn_cj4_t;

typedef struct {
	unsigned int pair[32]; /* Topology exclusion interaction bits for one warp,
							* each unsigned has bitS for 4*8 i clusters
							*/
} nbnxn_excl_t;

// sub sub define of nbnxn_list_work

/* Bounding box for a nbnxn atom cluster */
#define NNBSBB_C		 4
typedef struct {
	float lower[NNBSBB_C];
	float upper[NNBSBB_C];
} nbnxn_bb_t;

// end sub sub define of nbnxn_list_work

/* Working data for the actual i-supercell during pair search */
typedef struct nbnxn_list_work {
	gmx_cache_protect_t	 cp0;	/* Protect cache between threads			   */

	nbnxn_bb_t			 *bb_ci;  /* The bounding boxes, pbc shifted, for each cluster */
	float				  *pbb_ci; /* As bb_ci, but in xxxx packed format			   */
	real				   *x_ci;   /* The coordinates, pbc shifted, for each atom	   */
// #ifdef GMX_NBNXN_SIMD
//	 nbnxn_x_ci_simd_4xn_t  *x_ci_simd_4xn;
//	 nbnxn_x_ci_simd_2xnn_t *x_ci_simd_2xnn;
// #endif
	int					 cj_ind;		  /* The current cj_ind index for the current list	 */
	int					 cj4_init;		/* The first unitialized cj4 block				   */

	float				  *d2;			  /* Bounding box distance work array				  */

	nbnxn_cj_t			 *cj;			  /* The j-cell list								   */
	int					 cj_nalloc;	   /* Allocation size of cj							 */

	int					 ncj_noq;		 /* Nr. of cluster pairs without Coul for flop count  */
	int					 ncj_hlj;		 /* Nr. of cluster pairs with 1/2 LJ for flop count   */

	int					*sort;			/* Sort index					*/
	int					 sort_nalloc;	 /* Allocation size of sort	   */

	nbnxn_sci_t			*sci_sort;		/* Second sci array, for sorting */
	int					 sci_sort_nalloc; /* Allocation size of sci_sort   */

	gmx_cache_protect_t	 cp1;			 /* Protect cache between threads			   */
} nbnxn_list_work_t;

// end sub define for nbnxn_pairlist_t

typedef struct nbnxn_pairlist_t {
	gmx_cache_protect_t cp0;

	nbnxn_alloc_t	  *alloc;
	nbnxn_free_t	   *free;

	gmx_bool			bSimple;		 /* Simple list has na_sc=na_s and uses cj   *
										  * Complex list uses cj4					*/

	int					 na_ci;	   /* The number of atoms per i-cluster		*/
	int					 na_cj;	   /* The number of atoms per j-cluster		*/
	int					 na_sc;	   /* The number of atoms per super cluster	*/
	real					rlist;	   /* The radius for constructing the list	 */
	int					 nci;		 /* The number of i-clusters in the list	 */
	nbnxn_ci_t			 *ci;		  /* The i-cluster list, size nci			 */
	int					 ci_nalloc;   /* The allocation size of ci				*/
	int					 nsci;		/* The number of i-super-clusters in the list */
	nbnxn_sci_t			*sci;		 /* The i-super-cluster list				 */
	int					 sci_nalloc;  /* The allocation size of sci			   */

	int					 ncj;		 /* The number of j-clusters in the list	 */
	nbnxn_cj_t			 *cj;		  /* The j-cluster list, size ncj			 */
	int					 cj_nalloc;   /* The allocation size of cj				*/

	int					 ncj4;		/* The total number of 4*j clusters		 */
	nbnxn_cj4_t			*cj4;		 /* The 4*j cluster list, size ncj4		  */
	int					 cj4_nalloc;  /* The allocation size of cj4			   */
	int					 nexcl;	   /* The count for excl					   */
	nbnxn_excl_t		   *excl;		/* Atom interaction bits (non-exclusions)   */
	int					 excl_nalloc; /* The allocation size for excl			 */
	int					 nci_tot;	 /* The total number of i clusters		   */

	struct nbnxn_list_work *work;

	gmx_cache_protect_t	 cp1;
} nbnxn_pairlist_t;

// -------------------------------------------------------------------------------------

// sub define of 

typedef struct {
	real *f;	  /* f, size natoms*fstride							 */
	real *fshift; /* Shift force array, size SHIFTS*DIM				 */
	int   nV;	 /* The size of *Vvdw and *Vc						  */
	real *Vvdw;   /* Temporary Van der Waals group energy storage	   */
	real *Vc;	 /* Temporary Coulomb group energy storage			 */
	int   nVS;	/* The size of *VSvdw and *VSc						*/
	real *VSvdw;  /* Temporary SIMD Van der Waals group energy storage  */
	real *VSc;	/* Temporary SIMD Coulomb group energy storage		*/
} nbnxn_atomdata_output_t;

// sub sub define of nbnxn_buffer_flags_t

// end sub sub define of nbnxn_buffer_flags_t

/* Flags for telling if threads write to force output buffers */
typedef struct {
	int			   nflag;	   /* The number of flag blocks						 */
	unsigned __int64	*flag;		/* Bit i is set when thread i writes to a cell-block */
	int			   flag_nalloc; /* Allocation size of cxy_flag					   */
} nbnxn_buffer_flags_t;

typedef struct tMPI_Atomic
{
	int value; /**< The atomic value.*/
}
tMPI_Atomic_t;

// end sub define of 

typedef struct nbnxn_atomdata_t {
	nbnxn_alloc_t		   *alloc;
	nbnxn_free_t			*free;
	int					  ntype;		   /* The number of different atom types				 */
	real					*nbfp;			/* Lennard-Jones 6*C6 and 12*C12 params, size ntype^2*2 */
	int					  comb_rule;	   /* Combination rule, see enum above				   */
	real					*nbfp_comb;	   /* LJ parameter per atom type, size ntype*2		   */
	real					*nbfp_s4;		 /* As nbfp, but with stride 4, size ntype^2*4. This
											   * might suit 4-wide SIMD loads of two values (e.g.
											   * two floats in single precision on x86).			*/
	int					  natoms;		  /* Number of atoms									*/
	int					  natoms_local;	/* Number of local atoms						   */
	int					 *type;			/* Atom types										 */
	real					*lj_comb;		 /* LJ parameters per atom for combining for pairs	 */
	int					  XFormat;		 /* The format of x (and q), enum					  */
	int					  FFormat;		 /* The format of f, enum							  */
	real					*q;			   /* Charges, can be NULL if incorporated in x		  */
	int					  na_c;			/* The number of atoms per cluster					*/
	int					  nenergrp;		/* The number of energy groups						*/
	int					  neg_2log;		/* Log2 of nenergrp								   */
	int					 *energrp;		 /* The energy groups per cluster, can be NULL		 */
	gmx_bool				 bDynamicBox;	 /* Do we need to update shift_vec every step?	*/
	rvec					*shift_vec;	   /* Shift vectors, copied from t_forcerec			  */
	int					  xstride;		 /* stride for a coordinate in x (usually 3 or 4)	  */
	int					  fstride;		 /* stride for a coordinate in f (usually 3 or 4)	  */
	real					*x;			   /* x and possibly q, size natoms*xstride			  */

	/* j-atom minus i-atom index for generating self and Newton exclusions
	 * cluster-cluster pairs of the diagonal, for 4xn and 2xnn kernels.
	 */
	real					*simd_4xn_diagonal_j_minus_i;
	real					*simd_2xnn_diagonal_j_minus_i;
	/* Filters for topology exclusion masks for the SIMD kernels.
	 * filter2 is the same as filter1, but with each element duplicated.
	 */
	unsigned int			*simd_exclusion_filter1;
	unsigned int			*simd_exclusion_filter2;
	real					*simd_interaction_array; /* Array of masks needed for exclusions on QPX */
	int					  nout;				   /* The number of force arrays						 */
	nbnxn_atomdata_output_t *out;					/* Output data structures			   */
	int					  nalloc;				 /* Allocation size of all arrays (for x/f *x/fstride) */
	gmx_bool				 bUseBufferFlags;		/* Use the flags or operate on all atoms	 */
	nbnxn_buffer_flags_t	 buffer_flags;		   /* Flags for buffer zeroing+reduc.  */
	gmx_bool				 bUseTreeReduce;		 /* Use tree for force reduction */
	tMPI_Atomic_t		   *syncStep;			   /* Synchronization step for tree reduce */
} nbnxn_atomdata_t;

// ---------------------------------------------------------------------------------------------------

// sub define of interaction_const_t

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

// end sub define of interaction_const_t

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

// -----------------------------------------------------------------------------

typedef real    rvec[DIM];