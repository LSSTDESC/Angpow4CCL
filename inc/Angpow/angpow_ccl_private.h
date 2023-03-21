#ifndef ANGPOW_CCL_PRIVATE_SEEN
#define ANGPOW_CCL_PRIVATE_SEEN

#ifdef __cplusplus
extern "C" {
#endif
//delcaration of the C external interface functions of Angpow and related typedefs

// these declarations are here to substitute the corresponding CCL declaration to compile agpow_ccl.cc file without reference to CCL  
#ifdef INFILE_ANGPOWCCL_CC   
//declaration of structures
typedef struct s_f1d ccl_f1d_t; // correspond to ccl_f1d_t
typedef struct s_f2d ccl_f2d_t; // correspond to ccl_f2d_t (type of psp)
typedef struct s_cosmology ccl_cosmology; // correspond to ccl_cosmology
typedef struct s_cltracer ccl_cl_tracer_t; // correspond to ccl_cl_tracer_t
#endif

struct tracer_range { double chimin; double chimax; double zmin; double zmax;};

//declarations of functions
typedef double (* F1d_eval)(ccl_f1d_t *, double);  // correspond to ccl_f1d_t_eval
typedef double (* F2d_eval)(ccl_f2d_t *, double, double, void *, int *);  // correspond to ccl_f2d_t_eval
typedef double (* Tracer_get_kernel)(ccl_cl_tracer_t *, double, int *);  //correspond to ccl_cl_tracer_t_get_kernel
typedef double (* Tracer_get_transfer)(ccl_cl_tracer_t *, double, double, int *);  //correspond to ccl_cl_tracer_t_get_transfer
typedef double (* Comoving_radial_distance)(ccl_cosmology *, double , int *); // correspond to ccl_comoving_radial_distance
typedef double (* Scale_factor_of_chi)(ccl_cosmology *, double , int *); // correspond to ccl_scale_factor_of_chi

void ccl_angular_cls_angpow_default(F1d_eval f1d_eval, F2d_eval f2d_eval,
					Comoving_radial_distance chi_of_a, Scale_factor_of_chi a_of_chi,
					Tracer_get_kernel get_kernel, Tracer_get_transfer transfer,
					ccl_f2d_t *psp, ccl_cosmology *ccl_cosmo, ccl_cl_tracer_t *clt1, ccl_cl_tracer_t *clt2,
				    struct tracer_range *clt1_range, struct tracer_range *clt2_range, int *ell_values, int n_ells,
				    char *message, double *cl_out, int *status);

#ifdef __cplusplus
}
#endif

#endif// ANGPOW_CCL_PRIVATE_SEEN
