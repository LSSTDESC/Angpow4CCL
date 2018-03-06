#ifndef ANGPOW_CCL_PRIVATE_SEEN
#define ANGPOW_CCL_PRIVATE_SEEN



#ifdef __cplusplus
extern "C" {
#endif
//delcaration of the C external interface functions of Angpow and related typedefs

// these declarations are here to substitute the corresponding CCL declaration to compile agpow_ccl.cc file without reference to CCL  
#ifdef INFILE_ANGPOWCCL_CC   
//declaration of structures
typedef struct s_splpar SplPar; // correspond to SplPar
typedef struct s_cosmology ccl_cosmology; // correspond to ccl_cosmology
typedef struct s_cltracer CCL_ClTracer; // correspond to CCL_ClTracer
typedef struct s_clworkspace CCL_ClWorkspace; // correspond to CCL_ClWorksapce
#endif

struct angpow_cltracer { int has_rsd; int has_magnification; double chimin; double chimax; double zmin; double zmax; SplPar * spl_bz; SplPar * spl_nz;};
struct angpow_clworkspace { int lmax; int l_limber; int l_linstep; double l_logstep; int * l_arr;};

//declarations of functions
typedef double (* Spline_eval)(double , SplPar *); // correspond to ccl_spline_eval
typedef double (* Comoving_radial_distance)(ccl_cosmology *, double , int *); // correspond to ccl_comoving_radial_distance
typedef double (* Scale_factor_of_chi)(ccl_cosmology *, double , int *); // correspond to ccl_scale_factor_of_chi
typedef double (* Growth_rate)(ccl_cosmology *, double , int *); // correspond to ccl_growth_rate
typedef double (* Nonlin_matter_power)(ccl_cosmology *, double , double , int *); // correspond to ccl_non_linear_matter_power

void ccl_angular_cls_angpow_default(Spline_eval spl, Comoving_radial_distance chi,
				    Scale_factor_of_chi a, Growth_rate f, Nonlin_matter_power P,
				    ccl_cosmology *ccl_cosmo, CCL_ClWorkspace *w,
				    CCL_ClTracer *clt1,CCL_ClTracer *clt2,
				    struct angpow_cltracer *a_clt1, struct angpow_cltracer *a_clt2, struct angpow_clworkspace *a_w,	    
				    char * message, double *cl_out,int * status);

#ifdef __cplusplus
}
#endif

  
#endif// ANGPOW_CCL_PRIVATE_SEEN
