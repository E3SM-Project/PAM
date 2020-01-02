
#ifndef _ADERDT_H_
#define _ADERDT_H_

#include "const.h"
#include "SArray.h"


///////////////////////////////////////////////////////////////////////////////////////////
// Computes tord-1 time derivatives of rho, u, v, w, p, utend, vtend, and wtend using
// Differential Transforms in the time dimension using the x-direction flux Jacobian.
// This uses vector form for density to save calculations and Jacobian form for everything
// else. 
// 
// INPUTS
//   state: state values at tord GLL points stored in state(:,0,:). This routine expects
//          full density and pressure, not perturbations. dims are (var,time,space)
//   deriv: spatial derivative values at tord GLL points stored in deriv(:,0,:). drho and
//          dp can be perturbation or full differentials. dims are (var,time,space)
//   deriv_mat: Matrix that transforms tord GLL points into the spatial derivative stored
//              at the same tord GLL points.
// 
// OUTPUTS
//   state: 0th- to (tord-1)th-order time derivatives of the state values
//   deriv: 0th- to (tord-1)th-order time derivatives of the state spatial derivatives
//   utend: 0th- to (tord-1)th-order time derivatives of the u-tendency (RHS)
//          dims are (time,space)
//   vtend: 0th- to (tord-1)th-order time derivatives of the v-tendency (RHS)
//          dims are (time,space)
//   wtend: 0th- to (tord-1)th-order time derivatives of the w-tendency (RHS)
//          dims are (time,space)
///////////////////////////////////////////////////////////////////////////////////////////
YAKL_INLINE void diffTransformEulerX( SArray<real,numState,tord,tord> &state, 
                                      SArray<real,numState,tord,tord> &deriv,
                                      SArray<real,numState,tord,tord> &tend ,
                                      SArray<real,tord,tord> const &deriv_mat ) {
  SArray<real,tord> r, cs2or, cs2ot;
  for (int ii=0; ii<tord; ii++) {
    r(ii) = state(idR,0,ii);
    real t = state(idT,0,ii);
    real p = C0*pow(r(ii)*t,GAMMA);
    cs2or(ii) = GAMMA*p/(r(ii)*r(ii));
    cs2ot(ii) = GAMMA*p/(r(ii)*t    );
  }
  // tend will be used to store u*dr, u*du, u*dv, u*dw, u*dt

  // Compute the zeroth-order DTs of the intermediate functions and fluxes
  for (int ii=0; ii<tord; ii++) {
    real u  = state(idU,0,ii);
    
    // Initialize the 0th-order time DTs (i.e., the values)
    tend(idR,0,ii) = u*deriv(idR,0,ii);
    tend(idU,0,ii) = u*deriv(idU,0,ii);
    tend(idV,0,ii) = u*deriv(idV,0,ii);
    tend(idW,0,ii) = u*deriv(idW,0,ii);
    tend(idT,0,ii) = u*deriv(idT,0,ii);
  } // ii-loop

  // Loop over the time derivatives, computing the (kt+1)th time DTs in each iteration
  for (int kt=0; kt<tord-1; kt++) {

    // Compute (kt+1)th DT of u, v, w, and t (Jacobian form)
    for (int ii=0; ii<tord; ii++) {
      state(idR,kt+1,ii) = -( tend(idR,kt,ii) + r(ii)*deriv(idU,kt,ii)                                  )/(kt+1);  // r
      state(idU,kt+1,ii) = -( tend(idU,kt,ii) + cs2or(ii)*deriv(idR,kt,ii) + cs2ot(ii)*deriv(idT,kt,ii) )/(kt+1);  // u
      state(idV,kt+1,ii) = -( tend(idV,kt,ii)                                                           )/(kt+1);  // v
      state(idW,kt+1,ii) = -( tend(idW,kt,ii)                                                           )/(kt+1);  // w
      state(idT,kt+1,ii) = -( tend(idT,kt,ii)                                                           )/(kt+1);  // p
    }

    // Compute spatial derivative of the (kt+1)th DTs of the state
    for (int l=0; l<numState; l++) {
      for (int ii=0; ii<tord; ii++) {
        real d_dx = 0;
        // Matrix-vector multiply against the spatial differentiation matrix
        for (int s=0; s<tord; s++) {
          d_dx += deriv_mat(s,ii) * state(l,kt+1,s);
        }
        deriv(l,kt+1,ii) = d_dx;
      }
    }

    // Compute the (kt+1)th DT of all temporary variables
    for (int ii=0; ii<tord; ii++) {
      real tot_tmp_u_dr  = 0;
      real tot_tmp_u_du  = 0;
      real tot_tmp_u_dv  = 0;
      real tot_tmp_u_dw  = 0;
      real tot_tmp_u_dt  = 0;
      for (int rt=0; rt<=kt+1; rt++) {
        tot_tmp_u_dr  += state(idU,rt,ii) * deriv(idR,kt+1-rt,ii);
        tot_tmp_u_du  += state(idU,rt,ii) * deriv(idU,kt+1-rt,ii);
        tot_tmp_u_dv  += state(idU,rt,ii) * deriv(idV,kt+1-rt,ii);
        tot_tmp_u_dw  += state(idU,rt,ii) * deriv(idW,kt+1-rt,ii);
        tot_tmp_u_dt  += state(idU,rt,ii) * deriv(idT,kt+1-rt,ii);
      }
      tend(idR,kt+1,ii) = tot_tmp_u_dr;
      tend(idU,kt+1,ii) = tot_tmp_u_du;
      tend(idV,kt+1,ii) = tot_tmp_u_dv;
      tend(idW,kt+1,ii) = tot_tmp_u_dw;
      tend(idT,kt+1,ii) = tot_tmp_u_dt;
    } // ii-loop
  } // kt-loop

  // Transform utend, vtend, and wtend into the RHS of the wind equations
  for (int kt=0; kt<tord; kt++) {
    for (int ii=0; ii<tord; ii++) {
      tend(idR,kt,ii) = -( tend(idR,kt,ii) + r(ii)*deriv(idU,kt,ii)                                  );
      tend(idU,kt,ii) = -( tend(idU,kt,ii) + cs2or(ii)*deriv(idR,kt,ii) + cs2ot(ii)*deriv(idT,kt,ii) );
      tend(idV,kt,ii) = -( tend(idV,kt,ii)                                                           );
      tend(idW,kt,ii) = -( tend(idW,kt,ii)                                                           );
      tend(idT,kt,ii) = -( tend(idT,kt,ii)                                                           );
    }
  }
}



///////////////////////////////////////////////////////////////////////////////////////////
// Computes tord-1 time derivatives of rho, u, v, w, p, utend, vtend, and wtend using
// Differential Transforms in the time dimension using the y-direction flux Jacobian
// This uses vector form for density to save calculations and Jacobian form for everything
// else. 
// 
// INPUTS
//   state: state values at tord GLL points stored in state(:,0,:). This routine expects
//          full density and pressure, not perturbations. dims are (var,time,space)
//   deriv: spatial derivative values at tord GLL points stored in deriv(:,0,:). drho and
//          dp can be perturbation or full differentials. dims are (var,time,space)
//   deriv_mat: Matrix that transforms tord GLL points into the spatial derivative stored
//              at the same tord GLL points.
// 
// OUTPUTS
//   state: 0th- to (tord-1)th-order time derivatives of the state values
//   deriv: 0th- to (tord-1)th-order time derivatives of the state spatial derivatives
//   utend: 0th- to (tord-1)th-order time derivatives of the u-tendency (RHS)
//          dims are (time,space)
//   vtend: 0th- to (tord-1)th-order time derivatives of the v-tendency (RHS)
//          dims are (time,space)
//   wtend: 0th- to (tord-1)th-order time derivatives of the w-tendency (RHS)
//          dims are (time,space)
///////////////////////////////////////////////////////////////////////////////////////////
YAKL_INLINE void diffTransformEulerY( SArray<real,numState,tord,tord> &state, 
                                      SArray<real,numState,tord,tord> &deriv,
                                      SArray<real,numState,tord,tord> &tend ,
                                      SArray<real,tord,tord> const &deriv_mat ) {
  SArray<real,tord> r, cs2or, cs2ot;
  for (int ii=0; ii<tord; ii++) {
    r(ii) = state(idR,0,ii);
    real t = state(idT,0,ii);
    real p = C0*pow(r(ii)*t,GAMMA);
    cs2or(ii) = GAMMA*p/(r(ii)*r(ii));
    cs2ot(ii) = GAMMA*p/(r(ii)*t    );
  }
  // tend will be used to store v*dr, v*du, v*dv, v*dw, v*dt

  // Compute the zeroth-order DTs of the intermediate functions and fluxes
  for (int ii=0; ii<tord; ii++) {
    real v  = state(idV,0,ii);
    
    // Initialize the 0th-order time DTs (i.e., the values)
    tend(idR,0,ii) = v*deriv(idR,0,ii);
    tend(idU,0,ii) = v*deriv(idU,0,ii);
    tend(idV,0,ii) = v*deriv(idV,0,ii);
    tend(idW,0,ii) = v*deriv(idW,0,ii);
    tend(idT,0,ii) = v*deriv(idT,0,ii);
  } // ii-loop

  // Loop over the time derivatives, computing the (kt+1)th time DTs in each iteration
  for (int kt=0; kt<tord-1; kt++) {

    // Compute (kt+1)th DT of u, v, w, and t (Jacobian form)
    for (int ii=0; ii<tord; ii++) {
      state(idR,kt+1,ii) = -( tend(idR,kt,ii) + r(ii)*deriv(idV,kt,ii)                                  )/(kt+1);  // r
      state(idU,kt+1,ii) = -( tend(idU,kt,ii)                                                           )/(kt+1);  // u
      state(idV,kt+1,ii) = -( tend(idV,kt,ii) + cs2or(ii)*deriv(idR,kt,ii) + cs2ot(ii)*deriv(idT,kt,ii) )/(kt+1);  // v
      state(idW,kt+1,ii) = -( tend(idW,kt,ii)                                                           )/(kt+1);  // w
      state(idT,kt+1,ii) = -( tend(idT,kt,ii)                                                           )/(kt+1);  // p
    }

    // Compute spatial derivative of the (kt+1)th DTs of the state
    for (int l=0; l<numState; l++) {
      for (int ii=0; ii<tord; ii++) {
        real d_dy = 0;
        // Matrix-vector multiply against the spatial differentiation matrix
        for (int s=0; s<tord; s++) {
          d_dy += deriv_mat(s,ii) * state(l,kt+1,s);
        }
        deriv(l,kt+1,ii) = d_dy;
      }
    }

    // Compute the (kt+1)th DT of all temporary variables
    for (int ii=0; ii<tord; ii++) {
      real tot_tmp_v_dr  = 0;
      real tot_tmp_v_du  = 0;
      real tot_tmp_v_dv  = 0;
      real tot_tmp_v_dw  = 0;
      real tot_tmp_v_dt  = 0;
      for (int rt=0; rt<=kt+1; rt++) {
        tot_tmp_v_dr += state(idV,rt,ii) * deriv(idR,kt+1-rt,ii);
        tot_tmp_v_du += state(idV,rt,ii) * deriv(idU,kt+1-rt,ii);
        tot_tmp_v_dv += state(idV,rt,ii) * deriv(idV,kt+1-rt,ii);
        tot_tmp_v_dw += state(idV,rt,ii) * deriv(idW,kt+1-rt,ii);
        tot_tmp_v_dt += state(idV,rt,ii) * deriv(idT,kt+1-rt,ii);
      }
      tend(idR,kt+1,ii) = tot_tmp_v_dr;
      tend(idU,kt+1,ii) = tot_tmp_v_du;
      tend(idV,kt+1,ii) = tot_tmp_v_dv;
      tend(idW,kt+1,ii) = tot_tmp_v_dw;
      tend(idT,kt+1,ii) = tot_tmp_v_dt;
    } // ii-loop
  } // kt-loop

  // Transform utend, vtend, and wtend into the RHS of the wind equations
  for (int kt=0; kt<tord; kt++) {
    for (int ii=0; ii<tord; ii++) {
      tend(idR,kt,ii) = -( tend(idR,kt,ii) + r(ii)*deriv(idV,kt,ii)                                  );
      tend(idU,kt,ii) = -( tend(idU,kt,ii)                                                           );
      tend(idV,kt,ii) = -( tend(idV,kt,ii) + cs2or(ii)*deriv(idR,kt,ii) + cs2ot(ii)*deriv(idT,kt,ii) );
      tend(idW,kt,ii) = -( tend(idW,kt,ii)                                                           );
      tend(idT,kt,ii) = -( tend(idT,kt,ii)                                                           );
    }
  }
}



///////////////////////////////////////////////////////////////////////////////////////////
// Computes tord-1 time derivatives of rho, u, v, w, p, utend, vtend, and wtend using
// Differential Transforms in the time dimension using the z-direction flux Jacobian
// This uses vector form for density to save calculations and Jacobian form for everything
// else. 
// 
// INPUTS
//   state: state values at tord GLL points stored in state(:,0,:). This routine expects
//          full density and pressure, not perturbations. dims are (var,time,space)
//   deriv: spatial derivative values at tord GLL points stored in deriv(:,0,:). dp must
//          be the spatial derivative of perturbation pressure, not full pressure.
//   deriv_mat: Matrix that transforms tord GLL points into the spatial derivative stored
//              at the same tord GLL points.
//   dph: hydrostatic pressure spatial derivative in the vertical direction at tord GLL
//        points. Used to add hydrostatic balance in pressure advection term
// 
// OUTPUTS
//   state: 0th- to (tord-1)th-order time derivatives of the state values
//   deriv: 0th- to (tord-1)th-order time derivatives of the state spatial derivatives
//   utend: 0th- to (tord-1)th-order time derivatives of the u-tendency (RHS)
//          dims are (time,space)
//   vtend: 0th- to (tord-1)th-order time derivatives of the v-tendency (RHS)
//          dims are (time,space)
//   wtend: 0th- to (tord-1)th-order time derivatives of the w-tendency (RHS)
//          dims are (time,space)
///////////////////////////////////////////////////////////////////////////////////////////
YAKL_INLINE void diffTransformEulerZ( SArray<real,numState,tord,tord> &state, 
                                      SArray<real,numState,tord,tord> &deriv,
                                      SArray<real         ,tord,tord> &utend,
                                      SArray<real         ,tord,tord> &vtend,
                                      SArray<real         ,tord,tord> &wtend,
                                      SArray<real              ,tord> &dph,
                                      SArray<real,tord,tord> const &deriv_mat ) {
  SArray<real,tord,tord> tmp_r_w;     // r*w
  SArray<real,tord,tord> tmp_w_dp;    // w*dp/dz
  SArray<real,tord,tord> tmp_rr_dp;   // (1/rho)*dp'/dz
  SArray<real,tord,tord> tmp_p_dw;    // p*dw/dz
  // utend will be used to store w*du
  // vtend will be used to store w*dv
  // wtend will be used to store w*dw

  // The term (dp/dz)/rho involves division. Because of this, the time DTs
  // must be zeroed out to kill the term that would otherwise look like recursion
  for (int kt=1; kt<tord; kt++) {
    for (int ii=0; ii<tord; ii++) {
      tmp_rr_dp(kt,ii) = 0;
    }
  }

  // Compute the zeroth-order DTs of the intermediate functions and fluxes
  for (int ii=0; ii<tord; ii++) {
    // state values
    real r  = state(idR,0,ii);
    real u  = state(idU,0,ii);
    real v  = state(idV,0,ii);
    real w  = state(idW,0,ii);
    real p  = state(idT,0,ii);
    // state derivatives
    real du = deriv(idU,0,ii);
    real dv = deriv(idV,0,ii);
    real dw = deriv(idW,0,ii);
    real dp = deriv(idT,0,ii);
    
    // Initialize the 0th-order time DTs (i.e., the values)
    tmp_r_w  (0,ii) = r*w;
    utend    (0,ii) = w*du;
    vtend    (0,ii) = w*dv;
    wtend    (0,ii) = w*dw;
    tmp_w_dp (0,ii) = w*dp;
    tmp_rr_dp(0,ii) = dp/r;
    tmp_p_dw (0,ii) = p*dw;
  } // ii-loop

  // Loop over the time derivatives, computing the (kt+1)th time DTs in each iteration
  for (int kt=0; kt<tord-1; kt++) {

    // Compute (kt+1)th DT of u, v, w, and p
    for (int ii=0; ii<tord; ii++) {
      state(idU,kt+1,ii) = -(utend   (kt,ii)                        )/(kt+1);  // u
      state(idV,kt+1,ii) = -(vtend   (kt,ii)                        )/(kt+1);  // v
      state(idW,kt+1,ii) = -(wtend   (kt,ii) +      tmp_rr_dp(kt,ii))/(kt+1);  // w
      state(idT,kt+1,ii) = -(tmp_w_dp(kt,ii) + GAMMA*tmp_p_dw(kt,ii))/(kt+1);  // p
    }

    // Compute (kt+1)th DT of rho
    for (int ii=0; ii<tord; ii++) {
      real drflux_dz  = 0;
      // Matrix-vector multiply against the spatial differentiation matrix
      for (int s=0; s<tord; s++) {
        drflux_dz += deriv_mat(s,ii) * tmp_r_w(kt,s);
      }
      state(idR,kt+1,ii) = -drflux_dz/(kt+1);
    }

    // Compute spatial derivative of the (kt+1)th DTs of the state
    for (int l=0; l<numState; l++) {
      for (int ii=0; ii<tord; ii++) {
        real d_dz = 0;
        // Matrix-vector multiply against the spatial differentiation matrix
        for (int s=0; s<tord; s++) {
          d_dz += deriv_mat(s,ii) * state(l,kt+1,s);
        }
        deriv(l,kt+1,ii) = d_dz;
      }
    }

    // Compute the (kt+1)th DT of all temporary variables
    // Nearly all of these are of the form f*g
    // Except for (dp/dz)/rho, which is of the form f/g
    for (int ii=0; ii<tord; ii++) {
      real tot_tmp_r_w   = 0;
      real tot_tmp_w_du  = 0;
      real tot_tmp_w_dv  = 0;
      real tot_tmp_w_dw  = 0;
      real tot_tmp_w_dp  = 0;
      real tot_tmp_p_dw  = 0;
      real tot_tmp_rr_dp = 0;
      for (int rt=0; rt<=kt+1; rt++) {
        tot_tmp_r_w   += state(idR,rt,ii) * state(idW,kt+1-rt,ii);
        tot_tmp_w_du  += state(idW,rt,ii) * deriv(idU,kt+1-rt,ii);
        tot_tmp_w_dv  += state(idW,rt,ii) * deriv(idV,kt+1-rt,ii);
        tot_tmp_w_dw  += state(idW,rt,ii) * deriv(idW,kt+1-rt,ii);
        tot_tmp_w_dp  += state(idW,rt,ii) * deriv(idT,kt+1-rt,ii);
        tot_tmp_p_dw  += state(idT,rt,ii) * deriv(idW,kt+1-rt,ii);
        tot_tmp_rr_dp += deriv(idT,rt,ii) - state(idR,rt,ii) * tmp_rr_dp(kt+1-rt,ii);
      }
      tmp_r_w  (kt+1,ii) = tot_tmp_r_w ;
      utend    (kt+1,ii) = tot_tmp_w_du;
      vtend    (kt+1,ii) = tot_tmp_w_dv;
      wtend    (kt+1,ii) = tot_tmp_w_dw;
      tmp_w_dp (kt+1,ii) = tot_tmp_w_dp;
      tmp_p_dw (kt+1,ii) = tot_tmp_p_dw;
      tmp_rr_dp(kt+1,ii) = tot_tmp_rr_dp / state(idR,0,ii);
    } // ii-loop
  } // kt-loop

  // Transform utend, vtend, and wtend into the RHS of the wind equations
  for (int kt=0; kt<tord; kt++) {
    for (int ii=0; ii<tord; ii++) {
      utend(kt,ii) = -utend(kt,ii);
      vtend(kt,ii) = -vtend(kt,ii);
      wtend(kt,ii) = -wtend(kt,ii) - tmp_rr_dp(kt,ii);
    }
  }
}



////////////////////////////////////////////////////////////////////////////////////
// Compute the time average from tord-1 time derivatives at tord GLL points in
// space, and store into the 0th index in time for an array of numState variables
// 
// INPUTS
//   dts: tord-1 time derivatives of numState varaibles at tord GLL points in space
//   dom: Domain class object (needed for time step size)
// OUTPUTS
//   dts: time-average of numState variables at tord GLL points in space stored in
//        the 0th time index
////////////////////////////////////////////////////////////////////////////////////
YAKL_INLINE void timeAvg( SArray<real,numState,tord,tord> &dts , Domain const &dom ) {
  real dtmult = dom.dt;
  for (int kt=1; kt<tord; kt++) {
    for (int l=0; l<numState; l++) {
      for (int ii=0; ii<tord; ii++) {
        dts(l,0,ii) += dts(l,kt,ii) * dtmult / (kt+1);
      }
    }
    dtmult *= dom.dt;
  }
}



////////////////////////////////////////////////////////////////////////////////////
// Compute the time average from tord-1 time derivatives at tord GLL points in
// space, and store into the 0th index in time for a single varaible
// 
// INPUTS
//   dts: tord-1 time derivatives of one varaible at tord GLL points in space
//   dom: Domain class object (needed for time step size)
// OUTPUTS
//   dts: time-average of one varaible at tord GLL points in space stored in
//        the 0th time index
////////////////////////////////////////////////////////////////////////////////////
YAKL_INLINE void timeAvg( SArray<real,tord,tord> &dts , Domain const &dom ) {
  real dtmult = dom.dt;
  for (int kt=1; kt<tord; kt++) {
    for (int ii=0; ii<tord; ii++) {
      dts(0,ii) += dts(kt,ii) * dtmult / (kt+1);
    }
    dtmult *= dom.dt;
  }
}


#endif
