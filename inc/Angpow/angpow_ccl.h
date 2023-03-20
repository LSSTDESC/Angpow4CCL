#ifndef ANGPOW_CCL_SEEN
#define ANGPOW_CCL_SEEN

#include "ccl.h"

#include <math.h>

#include "Angpow/angpow_ccl_private.h"

//#pragma once

#ifdef HAVE_ANGPOW

ccl_f1d_t *new_redshift_kernel(ccl_f1d_t *kernel, ccl_cosmology *cosmo, int *status)
{
  if (kernel == NULL)
    return NULL;

  ccl_f1d_t *kn_z = NULL;
  size_t n = kernel->spline->size;
  double *z, *wz;

  z = malloc(sizeof(double) * n);
  if (z == NULL)
    *status = CCL_ERROR_MEMORY;
  wz = malloc(sizeof(double) * n);
  if (wz == NULL)
    *status = CCL_ERROR_MEMORY;

  if (*status == 0)
  {
    for (int iz = 0; iz < n; iz++)
    {
      double a = ccl_scale_factor_of_chi(cosmo, kernel->spline->x[iz], status);
      double hz = cosmo->params.h * ccl_h_over_h0(cosmo, a, status) / ccl_constants.CLIGHT_HMPC;

      z[iz] = 1. / a - 1.;
      wz[iz] = kernel->spline->y[iz] / hz; // W(z) = W(chi) * dchi/dz = W(chi)/H(z)
    }
    kn_z = ccl_f1d_t_new(n, z, wz, kernel->y0, kernel->yf, kernel->extrap_lo_type, kernel->extrap_hi_type, status);
  }

  if (z != NULL)
    free(z);
  if (wz != NULL)
    free(wz);

  return kn_z;
}

ccl_cl_tracer_t *new_redshift_tracer(ccl_cl_tracer_t *tracer, ccl_cosmology *cosmo, int *status)
{
  ccl_cl_tracer_t *tr_z = malloc(sizeof(ccl_cl_tracer_t));
  if (tr_z == NULL)
    *status = CCL_ERROR_MEMORY;

  if (*status == 0)
  {
    tr_z->der_angles = tracer->der_angles;
    tr_z->der_bessel = tracer->der_bessel;
    tr_z->transfer = tracer->transfer;
    tr_z->chi_min = 1. / ccl_scale_factor_of_chi(cosmo, tracer->chi_min, status) - 1.;
    tr_z->chi_max = 1. / ccl_scale_factor_of_chi(cosmo, tracer->chi_max, status) - 1.;
    tr_z->kernel = NULL;
    if (tracer->kernel != NULL)
    {
      tr_z->kernel = new_redshift_kernel(tracer->kernel, cosmo, status);
    }
  }

  if (*status)
  {
    ccl_cl_tracer_t_free(tr_z);
    tr_z = NULL;
  }
  return tr_z;
}

double redshift_tracer_get_kernel(ccl_cl_tracer_t *tr, double z, int *status) {
  if (tr != NULL && tr->kernel != NULL)
    return ccl_f1d_t_eval(tr->kernel, z);
  else
    return 1.0;
}

void set_tracer_range(struct tracer_range *range, ccl_cl_tracer_t *tr_chi, ccl_cl_tracer_t *tr_z)
{
  range->chimin = tr_chi->chi_min;
  range->chimax = tr_chi->chi_max;
  range->zmin = tr_z->chi_min;
  range->zmax = tr_z->chi_max;
}

// this is the public Angpow interface function to be used within CCL
void ccl_angular_cls_angpow(ccl_cosmology *ccl_cosmo,
				   ccl_cl_tracer_t *clt1, ccl_cl_tracer_t *clt2, ccl_f2d_t *psp,
				   int nl_out, int *l_out, double *cl_out, int *status)
{
  struct tracer_range clt1_range, clt2_range;
  ccl_cl_tracer_t *clt1_z, *clt2_z;

  clt1_z = new_redshift_tracer(clt1, ccl_cosmo, status);
  clt2_z = new_redshift_tracer(clt2, ccl_cosmo, status);
  if (*status == 0)
  {
    set_tracer_range(&clt1_range, clt1, clt1_z);
    set_tracer_range(&clt2_range, clt2, clt2_z);

    ccl_angular_cls_angpow_default(
        ccl_f1d_t_eval, ccl_f2d_t_eval,
        ccl_comoving_radial_distance, ccl_scale_factor_of_chi,
        redshift_tracer_get_kernel, ccl_cl_tracer_t_get_transfer,
        psp, ccl_cosmo, clt1_z, clt2_z, &clt1_range, &clt2_range, l_out, nl_out,
        ccl_cosmo->status_message, cl_out, status);
  }

  if (clt1_z != NULL)
  {
    clt1_z->transfer = NULL; // prevent the transfer object of origin tracer being freed
    ccl_cl_tracer_t_free(clt1_z);
  }
  if (clt2_z != NULL)
  {
    clt2_z->transfer = NULL;
    ccl_cl_tracer_t_free(clt2_z);
  }
}

#else

//this is the public Angpow interface function to be used within CCL
void ccl_angular_cls_angpow(ccl_cosmology *ccl_cosmo,
				   ccl_cl_tracer_t *clt1, ccl_cl_tracer_t *clt2, ccl_f2d_t *psp,
				   int nl_out, int *l_out, double *cl_out, int *status)
{
  printf("ccl_angular_cls_angpow: HAVE_ANGPOW has not been defined then Angpow is not included in CCL\n");
}

#endif /* HAVE_ANGPOW */
  
#ifdef __cplusplus
}
#endif

#endif /*ANGPOW_CCL_SEEN*/
