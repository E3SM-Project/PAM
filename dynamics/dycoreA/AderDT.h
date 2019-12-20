
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
                                      SArray<real         ,tord,tord> &utend,
                                      SArray<real         ,tord,tord> &vtend,
                                      SArray<real         ,tord,tord> &wtend,
                                      SArray<real,tord,tord> const &deriv_mat ) {
  SArray<real,tord,tord> tmp_r_u;     // r*u
  SArray<real,tord,tord> tmp_u_dp;    // u*dp/dx
  SArray<real,tord,tord> tmp_rr_dp;   // (1/rho)*dp/dx
  SArray<real,tord,tord> tmp_p_du;    // p*du/dx
  // utend will be used to store u*du
  // vtend will be used to store u*dv
  // wtend will be used to store u*dw

  // The term (dp/dx)/rho involves division. Because of this, the time DTs
  // must be zeroed out to kill the term that would cause a term to depends on itself
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
    tmp_r_u  (0,ii) = r*u;
    utend    (0,ii) = u*du;
    vtend    (0,ii) = u*dv;
    wtend    (0,ii) = u*dw;
    tmp_u_dp (0,ii) = u*dp;
    tmp_rr_dp(0,ii) = dp/r;
    tmp_p_du (0,ii) = p*du;
  } // ii-loop

  // Loop over the time derivatives, computing the (kt+1)th time DTs in each iteration
  for (int kt=0; kt<tord-1; kt++) {

    // Compute (kt+1)th DT of u, v, w, and p (Jacobian form)
    for (int ii=0; ii<tord; ii++) {
      state(idU,kt+1,ii) = -(utend   (kt,ii) +      tmp_rr_dp(kt,ii))/(kt+1);  // u
      state(idV,kt+1,ii) = -(vtend   (kt,ii)                        )/(kt+1);  // v
      state(idW,kt+1,ii) = -(wtend   (kt,ii)                        )/(kt+1);  // w
      state(idT,kt+1,ii) = -(tmp_u_dp(kt,ii) + GAMMA*tmp_p_du(kt,ii))/(kt+1);  // p
    }

    // Compute (kt+1)th DT of rho (vector form)
    for (int ii=0; ii<tord; ii++) {
      real drflux_dx = 0;
      // Matrix-vector multiply against the spatial differentiation matrix
      for (int s=0; s<tord; s++) {
        drflux_dx += deriv_mat(s,ii) * tmp_r_u(kt,s);
      }
      state(idR,kt+1,ii) = -drflux_dx/(kt+1);
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
    // Nearly all of these use the non-linear transform form: f*g
    // Except for (dp/dx)/rho, which is of the form: f/g
    for (int ii=0; ii<tord; ii++) {
      real tot_tmp_r_u   = 0;
      real tot_tmp_u_du  = 0;
      real tot_tmp_u_dv  = 0;
      real tot_tmp_u_dw  = 0;
      real tot_tmp_u_dp  = 0;
      real tot_tmp_p_du  = 0;
      real tot_tmp_rr_dp = 0;
      for (int rt=0; rt<=kt+1; rt++) {
        tot_tmp_r_u   += state(idR,rt,ii) * state(idU,kt+1-rt,ii);
        tot_tmp_u_du  += state(idU,rt,ii) * deriv(idU,kt+1-rt,ii);
        tot_tmp_u_dv  += state(idU,rt,ii) * deriv(idV,kt+1-rt,ii);
        tot_tmp_u_dw  += state(idU,rt,ii) * deriv(idW,kt+1-rt,ii);
        tot_tmp_u_dp  += state(idU,rt,ii) * deriv(idT,kt+1-rt,ii);
        tot_tmp_p_du  += state(idT,rt,ii) * deriv(idU,kt+1-rt,ii);
        tot_tmp_rr_dp += deriv(idT,rt,ii) - state(idR,rt,ii) * tmp_rr_dp(kt+1-rt,ii);
      }
      tmp_r_u  (kt+1,ii) = tot_tmp_r_u ;
      utend    (kt+1,ii) = tot_tmp_u_du;
      vtend    (kt+1,ii) = tot_tmp_u_dv;
      wtend    (kt+1,ii) = tot_tmp_u_dw;
      tmp_u_dp (kt+1,ii) = tot_tmp_u_dp;
      tmp_p_du (kt+1,ii) = tot_tmp_p_du;
      tmp_rr_dp(kt+1,ii) = tot_tmp_rr_dp / state(idR,0,ii);
    } // ii-loop
  } // kt-loop

  // Transform utend, vtend, and wtend into the RHS of the wind equations
  for (int kt=0; kt<tord; kt++) {
    for (int ii=0; ii<tord; ii++) {
      utend(kt,ii) = -utend(kt,ii) - tmp_rr_dp(kt,ii);
      vtend(kt,ii) = -vtend(kt,ii);
      wtend(kt,ii) = -wtend(kt,ii);
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
                                      SArray<real         ,tord,tord> &utend,
                                      SArray<real         ,tord,tord> &vtend,
                                      SArray<real         ,tord,tord> &wtend,
                                      SArray<real,tord,tord> const &deriv_mat ) {
  SArray<real,tord,tord> tmp_r_v;     // r*v
  SArray<real,tord,tord> tmp_v_dp;    // v*dp/dy
  SArray<real,tord,tord> tmp_rr_dp;   // (1/rho)*dp/dy
  SArray<real,tord,tord> tmp_p_dv;    // p*dv/dy
  // utend will be used to store v*du
  // vtend will be used to store v*dv
  // wtend will be used to store v*dw

  // The term (dp/dy)/rho involves division. Because of this, the time DTs
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
    tmp_r_v  (0,ii) = r*v;
    utend    (0,ii) = v*du;
    vtend    (0,ii) = v*dv;
    wtend    (0,ii) = v*dw;
    tmp_v_dp (0,ii) = v*dp;
    tmp_rr_dp(0,ii) = dp/r;
    tmp_p_dv (0,ii) = p*dv;
  } // ii-loop

  // Loop over the time derivatives, computing the (kt+1)th time DTs in each iteration
  for (int kt=0; kt<tord-1; kt++) {

    // Compute (kt+1)th DT of u, v, w, and p
    for (int ii=0; ii<tord; ii++) {
      state(idU,kt+1,ii) = -(utend   (kt,ii)                        )/(kt+1);  // u
      state(idV,kt+1,ii) = -(vtend   (kt,ii) +      tmp_rr_dp(kt,ii))/(kt+1);  // v
      state(idW,kt+1,ii) = -(wtend   (kt,ii)                        )/(kt+1);  // w
      state(idT,kt+1,ii) = -(tmp_v_dp(kt,ii) + GAMMA*tmp_p_dv(kt,ii))/(kt+1);  // p
    }

    // Compute (kt+1)th DT of rho
    for (int ii=0; ii<tord; ii++) {
      real drflux_dy  = 0;
      // Matrix-vector multiply against the spatial differentiation matrix
      for (int s=0; s<tord; s++) {
        drflux_dy += deriv_mat(s,ii) * tmp_r_v(kt,s);
      }
      state(idR,kt+1,ii) = -drflux_dy/(kt+1);
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
    // Nearly all of these are of the form f*g
    // Except for (dp/dy)/rho, which is of the form f/g
    for (int ii=0; ii<tord; ii++) {
      real tot_tmp_r_v   = 0;
      real tot_tmp_v_du  = 0;
      real tot_tmp_v_dv  = 0;
      real tot_tmp_v_dw  = 0;
      real tot_tmp_v_dp  = 0;
      real tot_tmp_p_dv  = 0;
      real tot_tmp_rr_dp = 0;
      for (int rt=0; rt<=kt+1; rt++) {
        tot_tmp_r_v   += state(idR,rt,ii) * state(idV,kt+1-rt,ii);
        tot_tmp_v_du  += state(idV,rt,ii) * deriv(idU,kt+1-rt,ii);
        tot_tmp_v_dv  += state(idV,rt,ii) * deriv(idV,kt+1-rt,ii);
        tot_tmp_v_dw  += state(idV,rt,ii) * deriv(idW,kt+1-rt,ii);
        tot_tmp_v_dp  += state(idV,rt,ii) * deriv(idT,kt+1-rt,ii);
        tot_tmp_p_dv  += state(idT,rt,ii) * deriv(idV,kt+1-rt,ii);
        tot_tmp_rr_dp += deriv(idT,rt,ii) - state(idR,rt,ii) * tmp_rr_dp(kt+1-rt,ii);
      }
      tmp_r_v  (kt+1,ii) = tot_tmp_r_v ;
      utend    (kt+1,ii) = tot_tmp_v_du;
      vtend    (kt+1,ii) = tot_tmp_v_dv;
      wtend    (kt+1,ii) = tot_tmp_v_dw;
      tmp_v_dp (kt+1,ii) = tot_tmp_v_dp;
      tmp_p_dv (kt+1,ii) = tot_tmp_p_dv;
      tmp_rr_dp(kt+1,ii) = tot_tmp_rr_dp / state(idR,0,ii);
    } // ii-loop
  } // kt-loop

  // Transform utend, vtend, and wtend into the RHS of the wind equations
  for (int kt=0; kt<tord; kt++) {
    for (int ii=0; ii<tord; ii++) {
      utend(kt,ii) = -utend(kt,ii);
      vtend(kt,ii) = -vtend(kt,ii) - tmp_rr_dp(kt,ii);
      wtend(kt,ii) = -wtend(kt,ii);
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
    tmp_w_dp (0,ii) = w*(dp+dph(ii));
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
