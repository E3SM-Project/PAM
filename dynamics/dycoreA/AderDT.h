
#ifndef _ADERDT_H_
#define _ADERDT_H_

#include "const.h"
#include "SArray.h"


///////////////////////////////////////////////////////////////////////////////////////////
// Computes tord-1 time derivatives at tord spatial GLL points located thoughout the cell
// of the state values, state spatial derivatives, and time tendencies of the state using
// Differential Transforms in the time dimension using the x-direction flux Jacobian.
// 
// Each vector contains rho, u, v, w, and theta
// 
// Acoustic dynamics are linearized in the time dimension (but not in space) because they
// can be handled less accurately without much impact overall. This reduces local storage
// requirements and computations.
// 
// INPUTS
//   state: state values at tord GLL points stored in state(:,0,:). This routine expects
//          full density and theta, not perturbations. dims are (var,time,space)
//   deriv: spatial derivative values at tord GLL points stored in deriv(:,0,:). drho and
//          dp can be perturbation or full derivatives. dims are (var,time,space)
//   deriv_mat: Matrix that transforms tord GLL points into the spatial derivative stored
//              at the same tord GLL points.
// 
// OUTPUTS
//   state: 0th- to (tord-1)th-order time derivatives of the state values
//   deriv: 0th- to (tord-1)th-order time derivatives of the state spatial derivatives
//   tend:  0th- to (tord-1)th-order time derivatives of the state time tendencies
//          dims are (var,time,space)
///////////////////////////////////////////////////////////////////////////////////////////
YAKL_INLINE void diffTransformEulerX( SArray<real,numState,tord,tord> &state, 
                                      SArray<real,numState,tord,tord> &deriv,
                                      SArray<real,numState,tord,tord> &tend ,
                                      SArray<real,tord,tord> const &deriv_mat ) {
  // tend will be used to store u*dr, u*du, u*dv, u*dw, u*dt
  SArray<real,tord> r, cs2or, cs2ot; // density, cs^2/rho, cs^2/theta
  // Precompute density, cs^2/rho, and cs^2/theta
  for (int ii=0; ii<tord; ii++) {
    r(ii) = state(idR,0,ii);
    real t = state(idT,0,ii);
    real p = C0*pow(r(ii)*t,GAMMA);
    cs2or(ii) = GAMMA*p/(r(ii)*r(ii));
    cs2ot(ii) = GAMMA*p/(r(ii)*t    );
  }

  // Compute the zeroth-order DTs of the intermediate functions and fluxes
  for (int ii=0; ii<tord; ii++) {
    real u = state(idU,0,ii);
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
      state(idT,kt+1,ii) = -( tend(idT,kt,ii)                                                           )/(kt+1);  // t
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
        tot_tmp_u_dr += state(idU,rt,ii) * deriv(idR,kt+1-rt,ii);
        tot_tmp_u_du += state(idU,rt,ii) * deriv(idU,kt+1-rt,ii);
        tot_tmp_u_dv += state(idU,rt,ii) * deriv(idV,kt+1-rt,ii);
        tot_tmp_u_dw += state(idU,rt,ii) * deriv(idW,kt+1-rt,ii);
        tot_tmp_u_dt += state(idU,rt,ii) * deriv(idT,kt+1-rt,ii);
      }
      tend(idR,kt+1,ii) = tot_tmp_u_dr;
      tend(idU,kt+1,ii) = tot_tmp_u_du;
      tend(idV,kt+1,ii) = tot_tmp_u_dv;
      tend(idW,kt+1,ii) = tot_tmp_u_dw;
      tend(idT,kt+1,ii) = tot_tmp_u_dt;
    } // ii-loop
  } // kt-loop

  // Transform utend, vtend, and wtend into the RHS of the wind equations
  // These are needed for the local tencencies in high-order flux-difference splitting
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
// Computes tord-1 time derivatives at tord spatial GLL points located thoughout the cell
// of the state values, state spatial derivatives, and time tendencies of the state using
// Differential Transforms in the time dimension using the y-direction flux Jacobian.
// 
// Each vector contains rho, u, v, w, and theta
// 
// Acoustic dynamics are linearized in the time dimension (but not in space) because they
// can be handled less accurately without much impact overall. This reduces local storage
// requirements and computations.
// 
// INPUTS
//   state: state values at tord GLL points stored in state(:,0,:). This routine expects
//          full density and theta, not perturbations. dims are (var,time,space)
//   deriv: spatial derivative values at tord GLL points stored in deriv(:,0,:). drho and
//          dp can be perturbation or full derivatives. dims are (var,time,space)
//   deriv_mat: Matrix that transforms tord GLL points into the spatial derivative stored
//              at the same tord GLL points.
// 
// OUTPUTS
//   state: 0th- to (tord-1)th-order time derivatives of the state values
//   deriv: 0th- to (tord-1)th-order time derivatives of the state spatial derivatives
//   tend:  0th- to (tord-1)th-order time derivatives of the state time tendencies
//          dims are (var,time,space)
///////////////////////////////////////////////////////////////////////////////////////////
YAKL_INLINE void diffTransformEulerY( SArray<real,numState,tord,tord> &state, 
                                      SArray<real,numState,tord,tord> &deriv,
                                      SArray<real,numState,tord,tord> &tend ,
                                      SArray<real,tord,tord> const &deriv_mat ) {
  // tend will be used to store v*dr, v*du, v*dv, v*dw, v*dt
  SArray<real,tord> r, cs2or, cs2ot; // density, cs^2/rho, cs^2/theta
  // Precompute density, cs^2/rho, and cs^2/theta
  for (int ii=0; ii<tord; ii++) {
    r(ii) = state(idR,0,ii);
    real t = state(idT,0,ii);
    real p = C0*pow(r(ii)*t,GAMMA);
    cs2or(ii) = GAMMA*p/(r(ii)*r(ii));
    cs2ot(ii) = GAMMA*p/(r(ii)*t    );
  }

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
      state(idT,kt+1,ii) = -( tend(idT,kt,ii)                                                           )/(kt+1);  // t
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
  // These are needed for the local tencencies in high-order flux-difference splitting
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
// Computes tord-1 time derivatives at tord spatial GLL points located thoughout the cell
// of the state values, state spatial derivatives, and time tendencies of the state using
// Differential Transforms in the time dimension using the z-direction flux Jacobian.
// 
// state vector contains rho, u, v, w, theta, and pressure
// deriv vector contains rho, u, v, w, theta, and perturbation pressure
// tend  vector is for   rho, u, v, w, and theta
// 
// Acoustic dynamics are linearized in the time dimension (but not in space) because they
// can be handled less accurately without much impact overall. This reduces local storage
// requirements and computations.
// 
// INPUTS
//   state: state values at tord GLL points stored in state(:,0,:). This routine expects
//          full density and theta, not perturbations. dims are (var,time,space)
//   deriv: spatial derivative values at tord GLL points stored in deriv(:,0,:).
//          dims are (var,time,space)
//   deriv_mat: Matrix that transforms tord GLL points into the spatial derivative stored
//              at the same tord GLL points.
// 
// OUTPUTS
//   state: 0th- to (tord-1)th-order time derivatives of the state values
//   deriv: 0th- to (tord-1)th-order time derivatives of the state spatial derivatives
//   tend:  0th- to (tord-1)th-order time derivatives of the state time tendencies
//          dims are (var,time,space)
///////////////////////////////////////////////////////////////////////////////////////////
YAKL_INLINE void diffTransformEulerZ( SArray<real,numState+1,tord,tord> &state, 
                                      SArray<real,numState+1,tord,tord> &deriv,
                                      SArray<real,numState  ,tord,tord> &tend ,
                                      SArray<real                ,tord> &rh   ,
                                      SArray<real,tord,tord> const &deriv_mat ) {
  // tend will be used to store w*dr, w*du, w*dv, w*dw, w*dt
  SArray<real,tord> r, cs2, rcs2ot; // density, cs^2, and rho*cs^2/theta
  // Precompute density, cs^2, and rho*cs^2/theta
  for (int ii=0; ii<tord; ii++) {
    r(ii) = state(idR,0,ii);
    real t = state(idT,0,ii);
    real p = state(idP,0,ii);
    cs2(ii) = GAMMA*p/r(ii);
    rcs2ot(ii) = GAMMA*p/t;
  }

  // Compute the zeroth-order DTs of the intermediate functions and fluxes
  for (int ii=0; ii<tord; ii++) {
    real w = state(idW,0,ii);
    tend(idR,0,ii) = w*deriv(idR,0,ii);
    tend(idU,0,ii) = w*deriv(idU,0,ii);
    tend(idV,0,ii) = w*deriv(idV,0,ii);
    tend(idW,0,ii) = w*deriv(idW,0,ii);
    tend(idT,0,ii) = w*deriv(idT,0,ii);
  } // ii-loop

  // Loop over the time derivatives, computing the (kt+1)th time DTs in each iteration
  for (int kt=0; kt<tord-1; kt++) {

    // Compute (kt+1)th DT of u, v, w, and t (Jacobian form)
    for (int ii=0; ii<tord; ii++) {
      state(idR,kt+1,ii) = -( tend(idR,kt,ii) + r(ii)*deriv(idW,kt,ii) )/(kt+1);  // r
      state(idU,kt+1,ii) = -( tend(idU,kt,ii)                          )/(kt+1);  // u
      state(idV,kt+1,ii) = -( tend(idV,kt,ii)                          )/(kt+1);  // v
      state(idW,kt+1,ii) = -( tend(idW,kt,ii) + deriv(idP,kt,ii)/r(ii) )/(kt+1);  // w
      state(idT,kt+1,ii) = -( tend(idT,kt,ii)                          )/(kt+1);  // t
      state(idP,kt+1,ii) = rcs2ot(ii)*state(idT,kt+1,ii) + cs2(ii)*state(idR,kt+1,ii);    // p
      // if (kt == 0) {
      //   state(idW,kt+1,ii) -= GRAV * ( state(idR,0,ii) - rh(ii) ) / state(idR,0,ii);
      // }
    }

    // Compute spatial derivative of the (kt+1)th DTs of the state
    for (int l=0; l<numState+1; l++) {
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
    for (int ii=0; ii<tord; ii++) {
      real tot_tmp_w_dr = 0;
      real tot_tmp_w_du = 0;
      real tot_tmp_w_dv = 0;
      real tot_tmp_w_dw = 0;
      real tot_tmp_w_dt = 0;
      for (int rt=0; rt<=kt+1; rt++) {
        tot_tmp_w_dr += state(idW,rt,ii) * deriv(idR,kt+1-rt,ii);
        tot_tmp_w_du += state(idW,rt,ii) * deriv(idU,kt+1-rt,ii);
        tot_tmp_w_dv += state(idW,rt,ii) * deriv(idV,kt+1-rt,ii);
        tot_tmp_w_dw += state(idW,rt,ii) * deriv(idW,kt+1-rt,ii);
        tot_tmp_w_dt += state(idW,rt,ii) * deriv(idT,kt+1-rt,ii);
      }
      tend(idR,kt+1,ii) = tot_tmp_w_dr;
      tend(idU,kt+1,ii) = tot_tmp_w_du;
      tend(idV,kt+1,ii) = tot_tmp_w_dv;
      tend(idW,kt+1,ii) = tot_tmp_w_dw;
      tend(idT,kt+1,ii) = tot_tmp_w_dt;
    } // ii-loop
  } // kt-loop

  // Transform utend, vtend, and wtend into the RHS of the wind equations
  // These are needed for the local tencencies in high-order flux-difference splitting
  for (int kt=0; kt<tord; kt++) {
    for (int ii=0; ii<tord; ii++) {
      tend(idR,kt,ii) = -( tend(idR,kt,ii) + r(ii)*deriv(idW,kt,ii) );
      tend(idU,kt,ii) = -( tend(idU,kt,ii)                          );
      tend(idV,kt,ii) = -( tend(idV,kt,ii)                          );
      tend(idW,kt,ii) = -( tend(idW,kt,ii) + deriv(idP,kt,ii)/r(ii) );
      tend(idT,kt,ii) = -( tend(idT,kt,ii)                          );
      // if (kt == 0) {
      //   tend(idW,0,ii) -= GRAV * ( state(idR,0,ii) - rh(ii) ) / state(idR,0,ii);
      // }
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
// space, and store into the 0th index in time for an array of numState variables
// 
// INPUTS
//   dts: tord-1 time derivatives of numState varaibles at tord GLL points in space
//   dom: Domain class object (needed for time step size)
// OUTPUTS
//   dts: time-average of numState variables at tord GLL points in space stored in
//        the 0th time index
////////////////////////////////////////////////////////////////////////////////////
YAKL_INLINE void timeAvg( SArray<real,numState+1,tord,tord> &dts , Domain const &dom ) {
  real dtmult = dom.dt;
  for (int kt=1; kt<tord; kt++) {
    for (int l=0; l<numState+1; l++) {
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
