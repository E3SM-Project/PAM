
#include "riemann_rho_theta_full.h"

namespace awfl {

  void riemann_rho_theta_full_x( pam::PamCoupler const &coupler ,
                                 real6d const & state_limits    ,
                                 real5d const &state_flux       ,
                                 real6d const & tracer_limits   ,
                                 real5d const &tracer_flux      ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    int constexpr idR = 0;  // density perturbation
    int constexpr idU = 1;  // u
    int constexpr idV = 2;  // v
    int constexpr idW = 3;  // w
    int constexpr idT = 4;  // potential temperature perturbation
    auto num_state   = state_limits .extent(0);
    auto num_tracers = tracer_limits.extent(0);
    auto nz          = coupler.get_nz     ();
    auto ny          = coupler.get_ny     ();
    auto nx          = coupler.get_nx     ();
    auto nens        = coupler.get_nens   ();
    auto Rd          = coupler.get_R_d    ();
    auto cp          = coupler.get_cp_d   ();
    auto p0          = coupler.get_p0     ();
    auto gamma       = coupler.get_gamma_d();
    auto kappa       = Rd/cp;
    auto C0          = pow( Rd * pow( p0 , -kappa ) , gamma );
    parallel_for( "Spatial.h X Riemann" , SimpleBounds<4>(nz,ny,nx+1,nens) ,
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      if (i == 0 ) {
        for (int l=0; l < num_state  ; l++) { state_limits (l,0,k,j,0 ,iens) = state_limits (l,0,k,j,nx,iens); }
        for (int l=0; l < num_tracers; l++) { tracer_limits(l,0,k,j,0 ,iens) = tracer_limits(l,0,k,j,nx,iens); }
      }
      if (i == nx) {
        for (int l=0; l < num_state  ; l++) { state_limits (l,1,k,j,nx,iens) = state_limits (l,1,k,j,0 ,iens); }
        for (int l=0; l < num_tracers; l++) { tracer_limits(l,1,k,j,nx,iens) = tracer_limits(l,1,k,j,0 ,iens); }
      }
      // Get left and right state
      real r_L = state_limits(idR,0,k,j,i,iens)    ;   real r_R = state_limits(idR,1,k,j,i,iens)    ;
      real u_L = state_limits(idU,0,k,j,i,iens)/r_L;   real u_R = state_limits(idU,1,k,j,i,iens)/r_R;
      real v_L = state_limits(idV,0,k,j,i,iens)/r_L;   real v_R = state_limits(idV,1,k,j,i,iens)/r_R;
      real w_L = state_limits(idW,0,k,j,i,iens)/r_L;   real w_R = state_limits(idW,1,k,j,i,iens)/r_R;
      real t_L = state_limits(idT,0,k,j,i,iens)/r_L;   real t_R = state_limits(idT,1,k,j,i,iens)/r_R;
      // Compute average state
      real r = 0.5_fp * (r_L + r_R);
      real u = 0.5_fp * (u_L + u_R);
      real v = 0.5_fp * (v_L + v_R);
      real w = 0.5_fp * (w_L + w_R);
      real t = 0.5_fp * (t_L + t_R);
      real p = C0 * pow(r*t,gamma);
      real cs2 = gamma*p/r;
      real cs  = sqrt(cs2);

      // COMPUTE UPWIND STATE FLUXES
      // Get left and right fluxes
      real q1_L = state_limits(idR,0,k,j,i,iens);   real q1_R = state_limits(idR,1,k,j,i,iens);
      real q2_L = state_limits(idU,0,k,j,i,iens);   real q2_R = state_limits(idU,1,k,j,i,iens);
      real q3_L = state_limits(idV,0,k,j,i,iens);   real q3_R = state_limits(idV,1,k,j,i,iens);
      real q4_L = state_limits(idW,0,k,j,i,iens);   real q4_R = state_limits(idW,1,k,j,i,iens);
      real q5_L = state_limits(idT,0,k,j,i,iens);   real q5_R = state_limits(idT,1,k,j,i,iens);
      // Compute upwind characteristics
      // Waves 1-3, velocity: u
      real w1, w2, w3;
      if (u > 0) {
        w1 = q1_L - q5_L/t;
        w2 = q3_L - v*q5_L/t;
        w3 = q4_L - w*q5_L/t;
      } else {
        w1 = q1_R - q5_R/t;
        w2 = q3_R - v*q5_R/t;
        w3 = q4_R - w*q5_R/t;
      }
      // Wave 5, velocity: u-cs
      real w5 =  u*q1_R/(2*cs) - q2_R/(2*cs) + q5_R/(2*t);
      // Wave 6, velocity: u+cs
      real w6 = -u*q1_L/(2*cs) + q2_L/(2*cs) + q5_L/(2*t);
      // Use right eigenmatrix to compute upwind flux
      real q1 = w1 + w5 + w6;
      real q2 = u*w1 + (u-cs)*w5 + (u+cs)*w6;
      real q3 = w2 + v*w5 + v*w6;
      real q4 = w3 + w*w5 + w*w6;
      real q5 =      t*w5 + t*w6;

      state_flux(idR,k,j,i,iens) = q2;
      state_flux(idU,k,j,i,iens) = q2*q2/q1 + C0*pow(q5,gamma);
      state_flux(idV,k,j,i,iens) = q2*q3/q1;
      state_flux(idW,k,j,i,iens) = q2*q4/q1;
      state_flux(idT,k,j,i,iens) = q2*q5/q1;

      // COMPUTE UPWIND TRACER FLUXES
      // Handle it one tracer at a time
      for (int tr=0; tr < num_tracers; tr++) {
        if (u > 0) {
          tracer_flux(tr,k,j,i,iens) = q2 * tracer_limits(tr,0,k,j,i,iens) / r_L;
        } else {
          tracer_flux(tr,k,j,i,iens) = q2 * tracer_limits(tr,1,k,j,i,iens) / r_R;
        }
      }
    });
  }

}


