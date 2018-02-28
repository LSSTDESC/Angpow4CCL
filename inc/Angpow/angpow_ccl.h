/** @file */
#ifdef __cplusplus
extern "C" {
#endif
#include "ccl_cls.h"
#include "ccl_power.h"
#include "ccl_background.h"
#include "ccl_error.h"
#include "ccl_utils.h"
#include "ccl_params.h"
#ifdef __cplusplus
}
#endif //__cplusplus
 
#pragma once


void ccl_angular_cls_angpow(ccl_cosmology *ccl_cosmo,CCL_ClWorkspace *w,
				   CCL_ClTracer *clt1,CCL_ClTracer *clt2,
				   double *cl_out,int * status);

{
  ccl_angular_cls_angpow_default(ccl_spline_eval, ccl_comoving_radial_distance,
			    ccl_scale_factor_of_chi, ccl_growth_rate,
			    ccl_nonlin_matter_power,
			    ccl_cosmo, w, clt1, clt2, cl_out, status)

}
#ifdef __cplusplus
}
#endif
