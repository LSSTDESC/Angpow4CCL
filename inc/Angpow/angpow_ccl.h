#ifndef ANGPOW_CCL_SEEN
#define ANGPOW_CCL_SEEN

#include "ccl_cls.h"
#include "ccl_power.h"
#include "ccl_background.h"
#include "ccl_error.h"
#include "ccl_utils.h"
#include "ccl_params.h"
#include "Angpow/angpow_ccl_private.h"

//#pragma once

#ifdef HAVE_ANGPOW

//this is the public Angpow interface function to be used within CCL
void ccl_angular_cls_angpow(ccl_cosmology *ccl_cosmo,CCL_ClWorkspace *w,
				   CCL_ClTracer *clt1,CCL_ClTracer *clt2,
				   double *cl_out,int * status)
{

  struct angpow_cltracer a_clt1;
  struct angpow_cltracer a_clt2;
  struct angpow_clworkspace a_w;
  a_clt1.has_rsd = clt1->has_rsd;  a_clt1.has_magnification = clt1->has_magnification;
  a_clt1.chimin = clt1->chimin;    a_clt1.chimax = clt1->chimax;
  a_clt1.zmin = clt1->zmin;        a_clt1.zmax = clt1->zmax;
  a_clt1.spl_bz = clt1->spl_bz;    a_clt1.spl_nz = clt1->spl_nz;

  a_clt2.has_rsd = clt2->has_rsd;  a_clt2.has_magnification = clt2->has_magnification;
  a_clt2.chimin = clt2->chimin;    a_clt2.chimax = clt2->chimax;
  a_clt2.zmin = clt2->zmin;        a_clt2.zmax = clt2->zmax;
  a_clt2.spl_bz = clt2->spl_bz;    a_clt2.spl_nz = clt2->spl_nz;
  
  a_w.lmax=w->lmax; a_w.l_limber=w->l_limber; a_w.l_linstep = w->l_linstep;
  a_w.l_logstep = w->l_logstep; a_w.l_arr = w->l_arr;
  
  ccl_angular_cls_angpow_default(ccl_spline_eval, ccl_comoving_radial_distance,
				 ccl_scale_factor_of_chi, ccl_growth_rate,
				 ccl_nonlin_matter_power,
				 ccl_cosmo, w, clt1, clt2,
				 &a_clt1, &a_clt2, &a_w,    
				 ccl_cosmo->status_message, cl_out, status);

}

#else

//this is the public Angpow interface function to be used within CCL
void ccl_angular_cls_angpow(ccl_cosmology *ccl_cosmo,CCL_ClWorkspace *w,
				   CCL_ClTracer *clt1,CCL_ClTracer *clt2,
				   double *cl_out,int * status)
{
  printf("ccl_angular_cls_angpow: HAVE_ANGPOW has not been defined then Angpow is not included in CCL\n");
}

#endif /* HAVE_ANGPOW */
  
#ifdef __cplusplus
}
#endif

#endif /*ANGPOW_CCL_SEEN*/
