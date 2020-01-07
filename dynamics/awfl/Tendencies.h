
#ifndef _TENDENCIES_H_
#define _TENDENCIES_H_

#include "const.h"
#include "Parallel.h"
#include "SArray.h"
#include "Domain.h"
#include "Exchange.h"
#include "WenoLimiter.h"
#include "AderDT.h"

class Tendencies {

  // Stores reconstruction samplts of the state vector at the left and right limits of the
  // cell interface. The cell of index i is bounded by interfaces of index i and i+1
  // For the second index, the "left" limit is index 0, and the "right" limit is index 1
  // dims(numState+1,2,nz+1,ny+1,nx+1)
  realArr stateLimits;

  // tord GLL weights normalized such that they sum to 1
  SArray<real,tord> gllWts;

  // Transforms either:
  // ord stencil cell averages to tord GLL points (without WENO limiting); OR
  // ord coefficients to tord GLL points (with WENO limiting)
  SArray<real,ord,tord> to_gll;

  // Transforms ord stencil averages into tord GLL points
  SArray<real,ord,tord> s2g;

  // Transforms (in the x-direction) either:
  // ord stencil cell averages to tord GLL points of the spatial derivative (without WENO limiting); OR
  // ord coefficients to tord GLL points of the spatial derivative (with WENO limiting)
  SArray<real,ord,tord> to_derivX_gll;

  // Transforms (in the y-direction) either:
  // ord stencil cell averages to tord GLL points of the spatial derivative (without WENO limiting); OR
  // ord coefficients to tord GLL points of the spatial derivative (with WENO limiting)
  SArray<real,ord,tord> to_derivY_gll;

  // Transforms (in the z-direction) either:
  // ord stencil cell averages to tord GLL points of the spatial derivative (without WENO limiting); OR
  // ord coefficients to tord GLL points of the spatial derivative (with WENO limiting)
  SArray<real,ord,tord> to_derivZ_gll;

  // Transform ord stencil averages into tord GLL points of the spatial derivative (x-direction)
  SArray<real,ord,tord> s2d2gX;

  // Transform ord stencil averages into tord GLL points of the spatial derivative (y-direction)
  SArray<real,ord,tord> s2d2gY;

  // Transform ord stencil averages into tord GLL points of the spatial derivative (z-direction)
  SArray<real,ord,tord> s2d2gZ;

  // In the x-direction, transforms tord GLL points into tord GLL points of the spatial derivative
  SArray<real,tord,tord> aderDerivX;

  // In the y-direction, transforms tord GLL points into tord GLL points of the spatial derivative
  SArray<real,tord,tord> aderDerivY;

  // In the z-direction, transforms tord GLL points into tord GLL points of the spatial derivative
  SArray<real,tord,tord> aderDerivZ;

  // Transforms stencil cell averages into coefficients for the candidate polynomials of WENO limiting
  SArray<real,ord,ord,ord> wenoRecon;

  // The ideal weights of the WENO candidate polynomials
  SArray<real,hs+2> wenoIdl;

  // The handicapping of the Total Variation of the high-order polynomial in WENO limiting
  real wenoSigma;

public :



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Spatially reconstruct ord stencil cell averages into tord GLL points with WENO limiting
  // 
  // INPUTS
  //   stencil: ord real values (cell averages centered about the cell being reconstructed)
  //   doWeno: 1 if we're doing WENO limiting, not 1 otherwise
  //   wenoRecon: Transformation matrices for biased stencil WENO Reconstruction
  //   to_gll: Transforms ord stencil values to tord GLL points (no WENO) or ord coefficients to tord
  //           GLL points (WENO)
  //   wenoIdl: Ideal weights for WENO candidate polynomials
  //   wenoSigma: Handicap factor for the Total Variation of the high-order polynomial
  // 
  // OUTPUTS: 
  //   gll: The reconstructed GLL points values
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  YAKL_INLINE void reconStencil(SArray<real,ord> const &stencil, SArray<real,tord> &gll, int const doWeno,
                                SArray<real,ord,ord,ord> const &wenoRecon, SArray<real,ord,tord> const &to_gll,
                                SArray<real,hs+2> const &wenoIdl, real wenoSigma) {
    SArray<real,ord> coefs;
    if (doWeno) {
      compute_weno_coefs(wenoRecon,stencil,coefs,wenoIdl,wenoSigma);
    } else {
      for (int ii=0; ii<ord; ii++) {
        coefs(ii) = stencil(ii);
      }
    }

    for (int ii=0; ii<tord; ii++) {
      gll(ii) = 0.;
      for (int s=0; s<ord; s++) {
        gll(ii) += to_gll(s,ii) * coefs(s);
      }
    }
  }



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Allocate arrays, pre-compute transformation arrays, WENO stuff, and GLL weights
  // 
  // INPUTS
  //   dom: The Domain object
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  void initialize(Domain const &dom);



  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute time tendencies for rho, u, v, w, and theta using ADER Differential Transform time
  // stepping in the x-direction
  // 
  // INPUTS
  //   state: The state vector: rho,u,v,w,theta. dims  (numState,nz+2*hs,ny+2*hs,nx+2*hs)
  //   dom: The Domain class object   
  //   par: The Parallel class object
  // 
  // OUTPUTS
  //   state: The state vector unchanged except for halos 
  //   exch: The Exchange class object for halo and edge exchanges
  //   tend: The time tendencies of the state vector   (numstate,nz,ny,nx)
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  void compEulerTend_X(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend);



  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute time tendencies for rho, u, v, w, and theta using ADER Differential Transform time
  // stepping in the x-direction
  // 
  // INPUTS
  //   state: The state vector: rho,u,v,w,theta. dims  (numState,nz+2*hs,ny+2*hs,nx+2*hs)
  //   dom: The Domain class object   
  //   par: The Parallel class object
  // 
  // OUTPUTS
  //   state: The state vector unchanged except for halos 
  //   exch: The Exchange class object for halo and edge exchanges
  //   tend: The time tendencies of the state vector   (numstate,nz,ny,nx)
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  void compEulerTend_Y(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend);



  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute time tendencies for rho, u, v, w, and theta using ADER Differential Transform time
  // stepping in the x-direction
  // 
  // INPUTS
  //   state: The state vector: rho,u,v,w,theta. dims  (numState,nz+2*hs,ny+2*hs,nx+2*hs)
  //   dom: The Domain class object   
  // 
  // OUTPUTS
  //   state: The state vector unchanged except for halos 
  //   tend: The time tendencies of the state vector   (numstate,nz,ny,nx)
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  void compEulerTend_Z(realArr &state, Domain const &dom, realArr &tend);



  ////////////////////////////////////////////////////////////////////////////////////////////
  // Compute the source term tendencies
  // 
  // INPUTS
  //   state: The fluid state vector, rho,u,v,w,theta; dims(numState,nz+2*hs,ny+2*hs,nx+2*hs)
  //   dom: The Domain class object
  // 
  // OUTPUTS
  //   tend: State tendency; dims(numState,nz,ny,nx)
  ////////////////////////////////////////////////////////////////////////////////////////////
  void compEulerTend_S(realArr const &state, Domain const &dom, realArr &tend);



  ////////////////////////////////////////////////////////////////////////////////////////////
  // Fill in vertical halo cells with boundary forcing data. vertical velocity is zero, and
  // all other variables replicate the last domain value to mimick zero gradient free-slip
  // 
  // INPUTS
  //   state: The fluid state vector, rho,u,v,w,theta; dims(numState,nz+2*hs,ny+2*hs,nx+2*hs)
  //   dom: The Domain class object
  // 
  // OUTPUTS
  //   state: fluid state with halos filled with boundary forcing data
  ////////////////////////////////////////////////////////////////////////////////////////////
  void stateBoundariesZ(realArr &state, Domain const &dom);



  ////////////////////////////////////////////////////////////////////////////////////////////
  // Fill in stateLimits vertical domain edge data with boundary forcing
  // Vertical velocity is zero, and all other variables replicate the last domain value to 
  // mimick zero gradient free-slip
  // 
  // INPUTS
  //   state: stateLimits for rho,u,v,w,theta,pressure ;  dims(numState+1,2,nz+1,ny+1,nx+1)
  //   dom: The Domain class object
  // 
  // OUTPUTS
  //   state: stateLimits with vertical domain edge data filled in
  ////////////////////////////////////////////////////////////////////////////////////////////
  void edgeBoundariesZ(realArr &stateLimits, Domain const &dom);



  // void compStrakaTend(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend);

};

#endif

