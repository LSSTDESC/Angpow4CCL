#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_integration.h"
//#include "Angpow/angpow_ccl.h"
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
  class RadSplineSelect : public RadSelectBase {
  public:
  RadSplineSelect(Spline_eval f, SplPar * spl,double x0,double xf): 
    RadSelectBase(x0, xf), f_(f), spl_(spl) {}
    
    
    virtual r_8 operator()(r_8 z) const {
      return f_(z,spl_);
    }
    
  private:
    Spline_eval f_;
    SplPar * spl_;
  };//RadSplineSelect
  

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
      return chi_(ccl_cosmo_,1./(1+z),&status);
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
  //! This class set the f1 and f2 functions from CCL_ClTracer class
  //! for galaxy counts :
  //!    f(ell,k,z) = k*sqrt(P(k,z))*[ b(z)*j_ell(r(z)*k) + f(z)*j"_ell(r(z)*k) ]
  class IntegrandCCL : public IntegrandBase {
  public:
  IntegrandCCL(Comoving_radial_distance chi, Growth_rate f, Nonlin_matter_power P, Spline_eval spl,
	        angpow_cltracer* clt, ccl_cosmology* cosmo,int ell=0, r_8 z=0.):
    chi_(chi), f_(f), P_(P), spl_(spl), clt_(clt), cosmo_(cosmo), ell_(ell), z_(z) {
      Init(ell,z);
    }
    //! Initialize the function f(ell,k,z) at a given ell and z
    //! to be integrated over k
    void Init(int ell, r_8 z){
      int status=0;
      ell_=ell; z_=z;
      R_=chi_(cosmo_,1.0/(1+z),&status);
      jlR_ = new JBess1(ell_,R_);
      if(clt_->has_rsd) {
	jlp1R_ = new JBess1(ell_+1,R_);
	// WARNING: here we want to store dlnD/dln(+1z) = - dlnD/dlna
	fz_= - f_(cosmo_,1.0/(1+z), &status);
      }
      bz_ = spl_(z,clt_->spl_bz);
    }
    virtual ~IntegrandCCL() {}
    //! Return f(ell,k,z) for a given k
    //! (ell and z must be initialized before) 
    virtual r_8 operator()(r_8 k) const {
      int status=0;
      r_8 Pk=P_(cosmo_,k,1./(1+z_),&status);
      r_8 x = k*R_;
      r_8 jlRk = (*jlR_)(k);
      r_8 delta = bz_*jlRk; // density term with bias
      if(clt_->has_rsd){ // RSD term
	r_8 jlRksecond = 0.;
	if(x<1e-40) { // compute second derivative j"_ell(r(z)*k)
	  if(ell_==0) {
	    jlRksecond = -1./3. + x*x/10.;
	  } else if(ell_==2) {
	    jlRksecond = 2./15. - 2*x*x/35.;
	  } else {
	    jlRksecond = 0.;
	  }
	} else {
	  jlRksecond = 2.*(*jlp1R_)(k)/x + (ell_*(ell_-1.)/(x*x) - 1.)*jlRk;
	}
	delta += fz_*jlRksecond;
      }
      return(k*sqrt(fabs(Pk))*delta);
    }
    //! Clone function for OpenMP integration
    virtual IntegrandCCL* clone() const {
      return new IntegrandCCL(static_cast<const IntegrandCCL&>(*this));
    }
    virtual void ExplicitDestroy() {
      if(jlR_) delete jlR_;
      if(jlp1R_) delete jlp1R_;
    }
  private:
    Comoving_radial_distance chi_; // function to get comoving radial distance from redshift
    Growth_rate f_; // function to get growth rate from redshift
    Nonlin_matter_power P_; // funciton to get matter power spectrum
    Spline_eval spl_;
    angpow_cltracer* clt_;  //no ownership
    ccl_cosmology* cosmo_;  //no ownership
    int ell_;  // multipole ell
    r_8 z_;   // redshift z
    r_8 R_;   // radial comoving distance r(z)
    r_8 fz_;  // growth rate f(z)
    r_8 bz_;  // bias b(z)
    JBess1* jlR_;  // j_ell(k*R)
    JBess1* jlp1R_;   // j_(ell+1)(k*R)
    
    //Minimal copy to allow Main operator(int, r_8, r_8) to work
    //JEC 22/4/17 use cloning of PowerSpectrum
    IntegrandCCL(const IntegrandCCL& copy) : clt_(copy.clt_), chi_(copy.chi_),f_(copy.f_),
					     P_(copy.P_),spl_(copy.spl_),
					     cosmo_(copy.cosmo_), ell_(0), z_(0), jlR_(0), jlp1R_(0){} 

  };//IntegrandCCL


}//end namespace



void ccl_angular_cls_angpow_default(Spline_eval spl, Comoving_radial_distance chi,
			    Scale_factor_of_chi a, Growth_rate f, Nonlin_matter_power P,
			    ccl_cosmology *ccl_cosmo, CCL_ClWorkspace *w,
			    CCL_ClTracer *clt1,CCL_ClTracer *clt2,
			    struct angpow_cltracer *a_clt1, struct angpow_cltracer *a_clt2,
			    struct angpow_clworkspace *a_w,	    
			    char * message, double *cl_out,int * status)
{
  // Initialize the Angpow parameters
  int chebyshev_order_1=9; // polynoms of order 2^{N}
  int chebyshev_order_2=9;
  int nroots=200; 
  int l_max_use=std::min(a_w->l_limber,a_w->lmax);
  double cl_kmax= 3.14*l_max_use/std::min(a_clt1->chimin,a_clt2->chimin); 
  int nsamp_z_1=(int)(std::max((double)15,0.124*cl_kmax*(a_clt1->chimin+a_clt1->chimax)-0.76*l_max_use));
  int nsamp_z_2=(int)(std::max((double)15,0.124*cl_kmax*(a_clt2->chimin+a_clt2->chimax)-0.76*l_max_use));
  //benchmark prints
  /*
  printf("lmax= %d l_limber= %d\n", l_max_use,a_w->l_limber);
  printf("chebyshevorder1: %d chebyshevorder2: %d\n",chebyshev_order_1,chebyshev_order_2);
  printf("radialorder1: %d radialorder2: %d\n", nsamp_z_1,nsamp_z_2);
  printf("kmax= %.lg\n", cl_kmax);
  printf("nroots= %d\n",nroots);
  printf(" %.lg %.lg\n", (double)a_clt1->chimax,(double)a_clt1->chimin);//,a_w->dchi);
  printf(" %.lg %.lg\n", (double)a_clt2->chimax,(double)a_clt2->chimin);//,a_w->dchi);
  */
  // Initialize the radial selection windows W(z)
  Angpow::RadSplineSelect Z1win(spl,a_clt1->spl_nz,a_clt1->zmin,a_clt1->zmax);
  Angpow::RadSplineSelect Z2win(spl,a_clt2->spl_nz,a_clt2->zmin,a_clt2->zmax);
 
  // The cosmological distance tool to make the conversion z <-> r(z)
  Angpow::CosmoCoordCCL cosmo(chi,a,ccl_cosmo);
  // Initilaie the two integrand functions f(ell,k,z)
  Angpow::IntegrandCCL int1(chi,f,P,spl,a_clt1,ccl_cosmo);
  Angpow::IntegrandCCL int2(chi,f,P,spl,a_clt2,ccl_cosmo);

  // Initialize the Cl with parameters to select the ell set which is interpolated after the processing
  Angpow::Clbase clout(l_max_use+1,a_w->l_linstep,a_w->l_logstep);
  // Check Angpow's ells match those of CCL (maybe we could to pass those ells explicitly to Angpow)
  for(int index_l=0; index_l<clout.Size(); index_l++) {
    if(clout[index_l].first!=a_w->l_arr[index_l]) {
      *status=CCL_ERROR_ANGPOW;
      strcpy(message,"ccl_cls.c: ccl_angular_cls_angpow_default(); "
	     "ell-bins defined in angpow don't match those of CCL\n");
      return;
    }
  }

  // Main class to compute Cl with Angpow
  Angpow::Pk2Cl pk2cl; //Default: the user parameters are used in the Constructor 
  pk2cl.SetOrdFunc(chebyshev_order_1,chebyshev_order_2);
  pk2cl.SetRadOrder(nsamp_z_1,nsamp_z_2);
  pk2cl.SetKmax(cl_kmax);
  pk2cl.SetNRootPerInt(nroots); 
  pk2cl.Compute(int1,int2,cosmo,&Z1win,&Z2win,clout[clout.Size()-1].first+1,clout);

  // Pass the Clbase class values (ell and C_ell) to the output spline
  for(int index_l=0; index_l<clout.Size(); index_l++)
    cl_out[index_l]=clout[index_l].second;
}


  
