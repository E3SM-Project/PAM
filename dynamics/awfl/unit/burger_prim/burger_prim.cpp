
#include "TransformMatrices.h"
#include "WenoLimiter.h"


real func(real xloc) {
  return (cos(2*M_PI*xloc) + 1.5) / 2;
}



real L1( realHost1d const &state_lo , realHost1d const &state_hi ) {
  int nx_hi = state_hi.dimension[0];
  int nx_lo = state_lo.dimension[0];
  int factor = nx_hi / nx_lo;
  realHost1d state_interp("state_interp",nx_lo);
  for (int i=0; i < nx_lo; i++) {
    real avg = 0;
    for (int ii=0; ii < factor; ii++) {
      avg += state_hi(i*factor + ii);
    }
    state_interp(i) = avg / factor;
  }

  real l1_numer = 0;
  real l1_denom = 0;
  for (int i=0; i < nx_lo; i++) {
    l1_numer += abs( state_interp(i) - state_lo(i) );
    l1_denom += abs( state_interp(i) );
  }

  return l1_numer / l1_denom;
}



template <unsigned nAder, unsigned ngll>
void compute_timeAvg( SArray<real,2,nAder,ngll> &dts , real dt ) {
  real dtmult = dt;
  for (int kt=1; kt<nAder; kt++) {
    for (int ii=0; ii<ngll; ii++) {
      dts(0,ii) += dts(kt,ii) * dtmult / (kt+1._fp);
    }
    dtmult *= dt;
  }
}



void print_state( realHost1d const &state ) {
  for (int i=0; i < state.dimension[0]; i++) {
    std::cout << std::scientific << state(i) << "\n";
  }
}



// Conservative approach
template <int ord>
realHost1d simulate_cons( int nx , real sim_time ) {
  int constexpr hs = (ord-1)/2;

  // Setup transformation matrices and WENO
  SArray<real,1,ord> gllPts_ord;
  SArray<real,1,ord> gllWts_ord;
  SArray<real,2,ord,ord> s2c;
  SArray<real,2,ord,ord> c2g;
  SArray<real,2,ord,ord> g2c;
  SArray<real,2,ord,ord> c2d;

  TransformMatrices::coefs_to_deriv(c2d);
  TransformMatrices::sten_to_coefs (s2c);
  TransformMatrices::coefs_to_gll  (c2g);
  TransformMatrices::gll_to_coefs  (g2c);

  TransformMatrices::get_gll_points (gllPts_ord);
  TransformMatrices::get_gll_weights(gllWts_ord);

  real dx = 1./nx;
  real dt = 0.4 * dx / 1.25;
  realHost1d state("state",nx+2*hs);
  realHost1d flux ("flux" ,nx+1);
  real etime = 0;
  for (int i=0; i < nx; i++) {
    state(hs+i) = 0;
    for (int ii=0; ii < ord; ii++) {
      real xloc = (i+0.5)*dx + gllPts_ord(ii)*dx;
      state(hs+i) += func(xloc) * gllWts_ord(ii);
    }
  }
  while (etime < sim_time) {
    if (sim_time - etime < dt) dt = sim_time - etime;

    // periodic boundaries
    for (int i=0; i < hs; i++) {
      state(i) = state(nx+i);
      state(hs+nx+i) = state(hs+i);
    }

    // Compute state at cell boundaries
    for (int i=0; i < nx; i++) {
      // Form stencil
      SArray<real,1,ord> stencil;
      for (int ii=0; ii < ord; ii++) { stencil(ii) = state(i+ii); }
      // Reconstruct
      SArray<real,1,ord> gll = c2g * s2c * stencil;
      SArray<real,2,ord,ord> state_DTs;
      SArray<real,2,ord,ord> flux_DTs;
      for (int ii=0; ii < ord; ii++) {
        state_DTs(0,ii) = gll(ii);
        flux_DTs(0,ii) = gll(ii) * gll(ii) / 2;
      }
      SArray<real,1,ord> next_deriv;
      SArray<real,1,ord> flux_DTs_curr;
      // Ader
      for (int kt=0; kt < ord-1; kt++) {
        for (int ii=0; ii < ord; ii++) { flux_DTs_curr(ii) = flux_DTs(kt,ii); }
        next_deriv = c2g * c2d * g2c * flux_DTs_curr;
        for (int ii=0; ii < ord; ii++) { state_DTs(kt+1,ii) = - next_deriv(ii) / dx / (kt+1.); }
        for (int ii=0; ii < ord; ii++) {
          real tot = 0;
          for (int rt=0; rt <= kt+1; rt++) {
            tot += state_DTs(rt,ii) * state_DTs(kt+1-rt,ii);
          }
          flux_DTs(kt+1,ii) = tot / 2;
        }
      }
      compute_timeAvg( flux_DTs , dt );
      flux(i+1) = flux_DTs(0,ord-1);
    }

    // periodic boundaries
    flux(0) = flux(nx);
    
    // Apply tendencies
    for (int i=0; i < nx; i++) {
      state(hs+i) = state(hs+i) - dt / dx * (flux(i+1) - flux(i));
    }
    etime += dt;
  }

  realHost1d state_ret("state_ret",nx);
  for (int i=0; i < nx; i++) { state_ret(i) = state(hs+i); }
  return state_ret;
}




// The most accurate primitive equation approach
// Conservative DTs, then exact central integral, then exact fwaves
template <int ord>
realHost1d simulate_prim1( int nx , real sim_time ) {
  int constexpr hs = (ord-1)/2;

  // Setup transformation matrices and WENO
  SArray<real,1,ord> gllPts_ord;
  SArray<real,1,ord> gllWts_ord;
  SArray<real,2,ord,ord> s2c;
  SArray<real,2,ord,ord> c2g;
  SArray<real,2,ord,ord> g2c;
  SArray<real,2,ord,ord> c2d;

  TransformMatrices::coefs_to_deriv(c2d);
  TransformMatrices::sten_to_coefs (s2c);
  TransformMatrices::coefs_to_gll  (c2g);
  TransformMatrices::gll_to_coefs  (g2c);

  TransformMatrices::get_gll_points (gllPts_ord);
  TransformMatrices::get_gll_weights(gllWts_ord);

  real dx = 1./nx;
  real dt = 0.4 * dx / 1.25;
  realHost1d state ("state" ,nx+2*hs);
  realHost2d fwaves("fwaves",2,nx+1);
  realHost1d tend  ("tend"  ,nx);
  real etime = 0;
  for (int i=0; i < nx; i++) {
    state(hs+i) = 0;
    for (int ii=0; ii < ord; ii++) {
      real xloc = (i+0.5)*dx + gllPts_ord(ii)*dx;
      state(hs+i) += func(xloc) * gllWts_ord(ii);
    }
  }
  while (etime < sim_time) {
    if (sim_time - etime < dt) dt = sim_time - etime;

    // periodic boundaries
    for (int i=0; i < hs; i++) {
      state(i) = state(nx+i);
      state(hs+nx+i) = state(hs+i);
    }

    // Compute state at cell boundaries
    for (int i=0; i < nx; i++) {
      // Form stencil
      SArray<real,1,ord> stencil;
      for (int ii=0; ii < ord; ii++) { stencil(ii) = state(i+ii); }
      // Reconstruct
      SArray<real,1,ord> u  = c2g       * s2c * stencil;
      SArray<real,1,ord> du = c2g * c2d * s2c * stencil;
      for (int ii=0; ii < ord; ii++) { du(ii) /= dx; }
      SArray<real,2,ord,ord> u_DTs;
      SArray<real,2,ord,ord> du_DTs;
      SArray<real,2,ord,ord> u_du_DTs;

      for (int ii=0; ii < ord; ii++) {
        u_DTs   (0,ii) = u (ii);
        du_DTs  (0,ii) = du(ii);
        u_du_DTs(0,ii) = u(ii)*du(ii);
      }

      // Ader
      for (int kt=0; kt < ord-1; kt++) {
        for (int ii=0; ii < ord; ii++) {
          u_DTs(kt+1,ii) = - u_du_DTs(kt,ii) / (kt+1.);
        }
        for (int ii=0; ii < ord; ii++) { u(ii) = u_DTs(kt+1,ii); }
        du = c2g * c2d * g2c * u;
        for (int ii=0; ii < ord; ii++) { du_DTs(kt+1,ii) = du(ii)/dx; }
        for (int ii=0; ii < ord; ii++) {
          real tot = 0;
          for (int rt=0; rt <= kt+1; rt++) {
            tot += u_DTs(rt,ii) * du_DTs(kt+1-rt,ii);
          }
          u_du_DTs(kt+1,ii) = tot;
        }
      }
      compute_timeAvg( u_DTs    , dt );
      compute_timeAvg( u_du_DTs , dt );

      tend(i) = 0;
      for (int ii=0; ii < ord; ii++) {
        tend(i) += -gllWts_ord(ii) * u_du_DTs(0,ii);
      }

      fwaves(1,i  ) = u_DTs(0,0    );
      fwaves(0,i+1) = u_DTs(0,ord-1);
    }

    // periodic boundaries
    fwaves(0,0 ) = fwaves(0,nx);
    fwaves(1,nx) = fwaves(1,0 );

    // Compute f-waves
    for (int i=0; i < nx+1; i++) {
      real um = fwaves(0,i);
      real up = fwaves(1,i);
      real u = 0.5 * (um + up);
      real du = up - um;
      fwaves(0,i) = 0;
      fwaves(1,i) = u*du;
    }
    
    // Apply tendencies
    for (int i=0; i < nx; i++) {
      tend(i) += - (fwaves(1,i) + fwaves(0,i+1)) / dx;
      state(hs+i) = state(hs+i) + dt * tend(i);
    }
    etime += dt;
  }

  realHost1d state_ret("state_ret",nx);
  for (int i=0; i < nx; i++) { state_ret(i) = state(hs+i); }
  return state_ret;
}




// Primitive equations with boundary terms in central integral using upwind values
template <int ord>
realHost1d simulate_prim2( int nx , real sim_time ) {
  int constexpr hs = (ord-1)/2;

  // Setup transformation matrices and WENO
  SArray<real,1,ord> gllPts_ord;
  SArray<real,1,ord> gllWts_ord;
  SArray<real,2,ord,ord> s2c;
  SArray<real,2,ord,ord> c2g;
  SArray<real,2,ord,ord> g2c;
  SArray<real,2,ord,ord> c2d;

  TransformMatrices::coefs_to_deriv(c2d);
  TransformMatrices::sten_to_coefs (s2c);
  TransformMatrices::coefs_to_gll  (c2g);
  TransformMatrices::gll_to_coefs  (g2c);

  TransformMatrices::get_gll_points (gllPts_ord);
  TransformMatrices::get_gll_weights(gllWts_ord);

  real dx = 1./nx;
  real dt = 0.4 * dx / 1.25;
  realHost1d state ("state" ,nx+2*hs);
  realHost2d fwaves("fwaves",2,nx+1);
  realHost2d ctend ("ctend" ,2,nx+1);
  realHost1d tend  ("tend"  ,nx);
  real etime = 0;
  for (int i=0; i < nx; i++) {
    state(hs+i) = 0;
    for (int ii=0; ii < ord; ii++) {
      real xloc = (i+0.5)*dx + gllPts_ord(ii)*dx;
      state(hs+i) += func(xloc) * gllWts_ord(ii);
    }
  }
  while (etime < sim_time) {
    if (sim_time - etime < dt) dt = sim_time - etime;

    // periodic boundaries
    for (int i=0; i < hs; i++) {
      state(i) = state(nx+i);
      state(hs+nx+i) = state(hs+i);
    }

    // Compute state at cell boundaries
    for (int i=0; i < nx; i++) {
      // Form stencil
      SArray<real,1,ord> stencil;
      for (int ii=0; ii < ord; ii++) { stencil(ii) = state(i+ii); }
      // Reconstruct
      SArray<real,1,ord> u  = c2g       * s2c * stencil;
      SArray<real,1,ord> du = c2g * c2d * s2c * stencil;
      for (int ii=0; ii < ord; ii++) { du(ii) /= dx; }
      SArray<real,2,ord,ord> u_DTs;
      SArray<real,2,ord,ord> du_DTs;
      SArray<real,2,ord,ord> u_du_DTs;

      for (int ii=0; ii < ord; ii++) {
        u_DTs   (0,ii) = u (ii);
        du_DTs  (0,ii) = du(ii);
        u_du_DTs(0,ii) = u(ii)*du(ii);
      }

      // Ader
      for (int kt=0; kt < ord-1; kt++) {
        for (int ii=0; ii < ord; ii++) {
          u_DTs(kt+1,ii) = - u_du_DTs(kt,ii) / (kt+1.);
        }
        for (int ii=0; ii < ord; ii++) { u(ii) = u_DTs(kt+1,ii); }
        du = c2g * c2d * g2c * u;
        for (int ii=0; ii < ord; ii++) { du_DTs(kt+1,ii) = du(ii)/dx; }
        for (int ii=0; ii < ord; ii++) {
          real tot = 0;
          for (int rt=0; rt <= kt+1; rt++) {
            tot += u_DTs(rt,ii) * du_DTs(kt+1-rt,ii);
          }
          u_du_DTs(kt+1,ii) = tot;
        }
      }
      compute_timeAvg( u_DTs    , dt );
      compute_timeAvg( u_du_DTs , dt );

      tend(i) = 0;
      for (int ii=1; ii < ord-1; ii++) {
        tend(i) += -gllWts_ord(ii) * u_du_DTs(0,ii);
      }

      fwaves(1,i  ) = u_DTs(0,0    );
      fwaves(0,i+1) = u_DTs(0,ord-1);

      ctend (1,i  ) = -u_du_DTs(0,0    );
      ctend (0,i+1) = -u_du_DTs(0,ord-1);
    }

    // periodic boundaries
    fwaves(0,0 ) = fwaves(0,nx);
    fwaves(1,nx) = fwaves(1,0 );
    ctend (0,0 ) = ctend (0,nx);
    ctend (1,nx) = ctend (1,0 );

    // Compute f-waves
    for (int i=0; i < nx+1; i++) {
      // Compute fwaves
      real um = fwaves(0,i);
      real up = fwaves(1,i);
      real u = 0.5 * (um + up);
      real du = up - um;
      fwaves(0,i) = 0;
      fwaves(1,i) = u*du;

      // Update left and right quadrature using upwind values
      if (i > 0 ) tend(i-1) += ctend(0,i) * gllWts_ord(0    );
      if (i < nx) tend(i  ) += ctend(0,i) * gllWts_ord(ord-1);
    }
    
    // Apply tendencies
    for (int i=0; i < nx; i++) {
      tend(i) += - (fwaves(1,i) + fwaves(0,i+1)) / dx;
      state(hs+i) = state(hs+i) + dt * tend(i);
    }
    etime += dt;
  }

  realHost1d state_ret("state_ret",nx);
  for (int i=0; i < nx; i++) { state_ret(i) = state(hs+i); }
  return state_ret;
}




// Primitive equations, except I average out the Jacobian in the central term
// and pull it out of the integral expression.
template <int ord>
realHost1d simulate_prim3( int nx , real sim_time ) {
  int constexpr hs = (ord-1)/2;

  // Setup transformation matrices and WENO
  SArray<real,1,ord> gllPts_ord;
  SArray<real,1,ord> gllWts_ord;
  SArray<real,2,ord,ord> s2c;
  SArray<real,2,ord,ord> c2g;
  SArray<real,2,ord,ord> g2c;
  SArray<real,2,ord,ord> c2d;

  TransformMatrices::coefs_to_deriv(c2d);
  TransformMatrices::sten_to_coefs (s2c);
  TransformMatrices::coefs_to_gll  (c2g);
  TransformMatrices::gll_to_coefs  (g2c);

  TransformMatrices::get_gll_points (gllPts_ord);
  TransformMatrices::get_gll_weights(gllWts_ord);

  real dx = 1./nx;
  real dt = 0.4 * dx / 1.25;
  realHost1d state ("state" ,nx+2*hs);
  realHost2d fwaves("fwaves",2,nx+1);
  realHost1d tend  ("tend"  ,nx);
  real etime = 0;
  for (int i=0; i < nx; i++) {
    state(hs+i) = 0;
    for (int ii=0; ii < ord; ii++) {
      real xloc = (i+0.5)*dx + gllPts_ord(ii)*dx;
      state(hs+i) += func(xloc) * gllWts_ord(ii);
    }
  }
  while (etime < sim_time) {
    if (sim_time - etime < dt) dt = sim_time - etime;

    // periodic boundaries
    for (int i=0; i < hs; i++) {
      state(i) = state(nx+i);
      state(hs+nx+i) = state(hs+i);
    }

    // Compute state at cell boundaries
    for (int i=0; i < nx; i++) {
      // Form stencil
      SArray<real,1,ord> stencil;
      for (int ii=0; ii < ord; ii++) { stencil(ii) = state(i+ii); }
      // Reconstruct
      SArray<real,1,ord> u  = c2g       * s2c * stencil;
      SArray<real,1,ord> du = c2g * c2d * s2c * stencil;
      for (int ii=0; ii < ord; ii++) { du(ii) /= dx; }
      SArray<real,2,ord,ord> u_DTs;
      SArray<real,2,ord,ord> du_DTs;
      SArray<real,2,ord,ord> u_du_DTs;

      for (int ii=0; ii < ord; ii++) {
        u_DTs   (0,ii) = u (ii);
        du_DTs  (0,ii) = du(ii);
        u_du_DTs(0,ii) = u(ii)*du(ii);
      }

      // Ader
      for (int kt=0; kt < ord-1; kt++) {
        for (int ii=0; ii < ord; ii++) {
          u_DTs(kt+1,ii) = - u_du_DTs(kt,ii) / (kt+1.);
        }
        for (int ii=0; ii < ord; ii++) { u(ii) = u_DTs(kt+1,ii); }
        du = c2g * c2d * g2c * u;
        for (int ii=0; ii < ord; ii++) { du_DTs(kt+1,ii) = du(ii)/dx; }
        for (int ii=0; ii < ord; ii++) {
          real tot = 0;
          for (int rt=0; rt <= kt+1; rt++) {
            tot += u_DTs(rt,ii) * du_DTs(kt+1-rt,ii);
          }
          u_du_DTs(kt+1,ii) = tot;
        }
      }
      compute_timeAvg( u_DTs    , dt );
      compute_timeAvg( u_du_DTs , dt );

      real uavg = 0;
      for (int ii=0; ii < ord; ii++) {
        uavg += gllWts_ord(ii) * u_DTs(0,ii);
      }

      tend(i) = - uavg * (u_DTs(0,ord-1) - u_DTs(0,0)) / dx;

      fwaves(1,i  ) = u_DTs(0,0    );
      fwaves(0,i+1) = u_DTs(0,ord-1);
    }

    // periodic boundaries
    fwaves(0,0 ) = fwaves(0,nx);
    fwaves(1,nx) = fwaves(1,0 );

    // Compute f-waves
    for (int i=0; i < nx+1; i++) {
      real um = fwaves(0,i);
      real up = fwaves(1,i);
      real u = 0.5 * (um + up);
      real du = up - um;
      fwaves(0,i) = 0;
      fwaves(1,i) = u*du;
    }
    
    // Apply tendencies
    for (int i=0; i < nx; i++) {
      tend(i) += - (fwaves(1,i) + fwaves(0,i+1)) / dx;
      state(hs+i) = state(hs+i) + dt * tend(i);
    }
    etime += dt;
  }

  realHost1d state_ret("state_ret",nx);
  for (int i=0; i < nx; i++) { state_ret(i) = state(hs+i); }
  return state_ret;
}




int main() {
  real constexpr sim_time = 0.1;

  realHost1d state_hi = simulate_cons<5>( 1000 , sim_time );

  real l1_100 = L1( simulate_prim3<5>( 100 , sim_time ) , state_hi );
  real l1_50  = L1( simulate_prim3<5>( 50  , sim_time ) , state_hi );

  std::cout << log(l1_100 / l1_50) / log(50./100.) << "\n";
  // print_state( simulate_prim1<5>( 100 , sim_time ) );

  std::cout << L1( simulate_cons <5>( 50 , sim_time ) , state_hi ) << "\n";
  std::cout << L1( simulate_prim3<5>( 50 , sim_time ) , state_hi ) << "\n";
}

