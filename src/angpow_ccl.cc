#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <numeric>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_integration.h"
#include "Angpow/angpow_tools.h"
#include "Angpow/angpow_parameters.h"
#include "Angpow/angpow_pk2cl.h"
#include "Angpow/angpow_powspec_base.h"
#include "Angpow/angpow_cosmo_base.h"
#include "Angpow/angpow_radial.h"
#include "Angpow/angpow_radial_base.h"
#include "Angpow/angpow_clbase.h"
#include "Angpow/angpow_ctheta.h"
#include "Angpow/angpow_exceptions.h"  //exceptions
#include "Angpow/angpow_integrand_base.h"

#define INFILE_ANGPOWCCL_CC 
#include "Angpow/angpow_ccl_private.h"

#define CCL_ERROR_ANGPOW 1042


//! Here CCL is passed to Angpow base classes
namespace Angpow {
  //! Selection window W(z) with spline input from CCL
  class RadSelectCCL : public RadSelectBase
  {
  public:
    RadSelectCCL(Tracer_get_kernel get_kernel, ccl_cl_tracer_t *clt, double x0, double xf) :
      RadSelectBase(x0, xf), get_kernel_(get_kernel), clt_(clt) {}

    virtual r_8 operator()(r_8 z) const {
      int status = 0;
      return get_kernel_(clt_, z, &status);
    }

  private:
    Tracer_get_kernel get_kernel_;
    ccl_cl_tracer_t *clt_;
  }; // RadSplineSelect

  //! This class import the cosmology from CCL to make the conversion z <-> r(z)
  //! and then to compute the cut-off in the integrals
  class CosmoCoordCCL : public CosmoCoordBase {
  public:
    //! Ctor
  CosmoCoordCCL(Comoving_radial_distance chi, Scale_factor_of_chi a, ccl_cosmology * cosmo)
    : chi_(chi), a_(a), ccl_cosmo_(cosmo) 
    {
    }//Ctor
    
    //! Dtor
    virtual ~CosmoCoordCCL() {}
    
    //! r(z): radial comoving distance Mpc
    inline virtual double r(double z) const {
      int status=0;
      return chi_(ccl_cosmo_, 1./(1.+z), &status);
    }
    //! z(r): the inverse of radial comoving distance (Mpc)
    inline virtual double z(double r) const {
      int status=0;
      return 1./a_(ccl_cosmo_,r,&status)-1.;
    }
    
  protected:
    Comoving_radial_distance chi_;
    Scale_factor_of_chi a_;
    ccl_cosmology* ccl_cosmo_;  //!< access to CCL cosmology
  };// CosmoCoordCCL


  //! Base class of half the integrand function
  //! Basically Angpow does three integrals of the form
  //!    C_ell = int dz1 int dz2 int dk f1(ell,k,z1)*f2(ell,k,z2)
  //! This class set the f1 and f2 functions from ccl_cl_tracer_c class
  class IntegrandCCL : public IntegrandBase {
  public:
    IntegrandCCL(Comoving_radial_distance chi, F2d_eval f2d_eval, Tracer_get_transfer get_transfer,
                 ccl_f2d_t *psp, ccl_cl_tracer_t *clt, ccl_cosmology *cosmo, int ell = 0, r_8 z = 0.) :
                 chi_(chi), f2d_eval_(f2d_eval), get_transfer_(get_transfer), psp_(psp), clt_(clt), cosmo_(cosmo), ell_(ell), z_(z)
    {
      Init(ell, z);
    }
    //! Initialize the function f(ell,k,z) at a given ell and z
    //! to be integrated over k
    void Init(int ell, r_8 z) {
      int status = 0;
      ell_ = ell;
      z_ = z;
      R_ = chi_(cosmo_, 1.0 / (1 + z), &status);
      jlR_ = new JBess1(ell_, R_);
    }
    virtual ~IntegrandCCL() {}
    //! Return f(ell,k,z) for a given k
    //! (ell and z must be initialized before)
    virtual r_8 operator()(r_8 k) const {
      int status = 0;
      r_8 Pk = f2d_eval_(psp_, log(k), 1./(1+z_), cosmo_, &status);
      r_8 jlRk = (*jlR_)(k);
      r_8 transfer = get_transfer_(clt_, log(k), 1./(1.+z_), &status); // transfer function defined in tracer object
      return (k * sqrt(fabs(Pk)) * transfer * jlRk);
    }
    //! Clone function for OpenMP integration
    virtual IntegrandCCL* clone() const {
      return new IntegrandCCL(static_cast<const IntegrandCCL&>(*this));
    }
    virtual void ExplicitDestroy() {
      if(jlR_) delete jlR_;
    }
  private:
    Comoving_radial_distance chi_; // function to get comoving radial distance from redshift
    F2d_eval f2d_eval_;  // function to evaluate ccl power spectrum
    Tracer_get_transfer get_transfer_;  // function to evaluate transfer function
    ccl_f2d_t *psp_; // ccl power spectrum object
    ccl_cl_tracer_t *clt_;  // no ownership
    ccl_cosmology* cosmo_;  //no ownership
    int ell_;  // multipole ell
    r_8 z_;   // redshift z
    r_8 R_;   // radial comoving distance r(z)
    JBess1* jlR_;  // j_ell(k*R)
    
    //Minimal copy to allow Main operator(int, r_8, r_8) to work
    //JEC 22/4/17 use cloning of PowerSpectrum
    IntegrandCCL(const IntegrandCCL& copy) : chi_(copy.chi_), f2d_eval_(copy.f2d_eval_), get_transfer_(copy.get_transfer_),
					     psp_(copy.psp_), clt_(copy.clt_), cosmo_(copy.cosmo_), ell_(0), z_(0), jlR_(0){} 

  };//IntegrandCCL

  //! This class represents the output Cls from Angpow
  //! Allows passing the ell values to be computed explicitly
  class CloutCCL : public Clbase
  {
  public:
    CloutCCL(std::vector<int> ells) : Clbase(0), cl_map_updated_(false) {
      ells_ = ells;
      sort(ells_.begin(), ells_.end());

      cls_.resize(ells_.size());
      for (int i = 0; i < ells_.size(); i++) {
          cls_[i].first = ells_[i];
      }

      Lmax_ = ells_.back();
      ellsAll_.resize(Lmax_);
      std::iota(std::begin(ellsAll_), std::end(ellsAll_), 0); // 0, 1,..., Lmax-1
      maximal_ = (ells_ == ellsAll_);
    }

    Acl &operator[](int index_l) {
      cl_map_updated_ = false;
      return Clbase::operator[](index_l);
    }

    void Interpolate() {
      Clbase::Interpolate();
      updateClMap();
    }

    void updateClMap() {
      for (auto &cl : cls_) {
        cl_map_[cl.first] = cl.second;
      }
      cl_map_updated_ = true;
    }

    double getClOfEll(int ell) {
      if (!cl_map_updated_) {
        updateClMap();
      }
      return cl_map_.at(ell);
    }

  private:
    std::map<int, double> cl_map_;
    bool cl_map_updated_;
  }; //CloutCCL
}//end namespace

void ccl_angular_cls_angpow_default(F1d_eval f1d_eval, F2d_eval f2d_eval,
                                    Comoving_radial_distance chi_of_a, Scale_factor_of_chi a_of_chi,
                                    Tracer_get_kernel get_kernel, Tracer_get_transfer get_transfer,
                                    ccl_f2d_t *psp, ccl_cosmology *ccl_cosmo, ccl_cl_tracer_t *clt1, ccl_cl_tracer_t *clt2,
                                    struct tracer_range *clt1_range, struct tracer_range *clt2_range, int *ell_values, int n_ells,
                                    char *message, double *cl_out, int *status)
{
  // Initialize the Angpow parameters
  int chebyshev_order_1 = 9; // polynoms of order 2^{N}
  int chebyshev_order_2 = 9;
  int nroots = 200;

  // Initialize the Cl with parameters to select the ell set which is interpolated after the processing
  std::vector<int> ells(ell_values, ell_values + n_ells);
  Angpow::CloutCCL clout(ells);

  int l_max_use = clout.InitialElls().back();
  double cl_kmax = 3.14 * l_max_use / std::min(clt1_range->chimin, clt2_range->chimin);
  int nsamp_z_1 = (int)(std::max((double)15, 0.124 * cl_kmax * (clt1_range->chimin + clt1_range->chimax) - 0.76 * l_max_use));
  int nsamp_z_2 = (int)(std::max((double)15, 0.124 * cl_kmax * (clt2_range->chimin + clt2_range->chimax) - 0.76 * l_max_use));

  // Initialize the radial selection windows W(z)
  Angpow::RadSelectCCL Z1win(get_kernel, clt1, clt1_range->zmin, clt1_range->zmax);
  Angpow::RadSelectCCL Z2win(get_kernel, clt2, clt2_range->zmin, clt2_range->zmax);
  // The cosmological distance tool to make the conversion z <-> r(z)
  Angpow::CosmoCoordCCL cosmo(chi_of_a, a_of_chi, ccl_cosmo);
  // Initilaie the two integrand functions f(ell,k,z)
  Angpow::IntegrandCCL int1(chi_of_a, f2d_eval, get_transfer, psp, clt1, ccl_cosmo);
  Angpow::IntegrandCCL int2(chi_of_a, f2d_eval, get_transfer, psp, clt2, ccl_cosmo);

  // Main class to compute Cl with Angpow
  Angpow::Pk2Cl pk2cl; // Default: the user parameters are used in the Constructor
  pk2cl.SetOrdFunc(chebyshev_order_1, chebyshev_order_2);
  pk2cl.SetRadOrder(nsamp_z_1, nsamp_z_2);
  pk2cl.SetKmax(cl_kmax);
  pk2cl.SetNRootPerInt(nroots);
  pk2cl.Compute(int1, int2, cosmo, &Z1win, &Z2win, clout[clout.Size() - 1].first + 1, clout);

  // Pass the Clbase class values (ell and C_ell) to the output spline
  for (int index_l = 0; index_l < n_ells; index_l++)
    cl_out[index_l] = clout.getClOfEll(ell_values[index_l]);
}
