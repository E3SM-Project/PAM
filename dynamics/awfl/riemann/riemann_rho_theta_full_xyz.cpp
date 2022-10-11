
#include "riemann_rho_theta_full.h"

namespace awfl {

  void riemann_rho_theta_full_xyz( pam::PamCoupler const &coupler  ,
                                   real6d const & state_limits_x   ,
                                   real5d const &state_flux_x      ,
                                   real6d const & tracer_limits_x  ,
                                   real5d const &tracer_flux_x     ,
                                   real6d const & state_limits_y   ,
                                   real5d const &state_flux_y      ,
                                   real6d const & tracer_limits_y  ,
                                   real5d const &tracer_flux_y     ,
                                   real6d const & state_limits_z   ,
                                   real5d const &state_flux_z      ,
                                   real6d const & tracer_limits_z  ,
                                   real5d const &tracer_flux_z     ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    int constexpr idR = 0;  // density perturbation
    int constexpr idU = 1;  // u
    int constexpr idV = 2;  // v
    int constexpr idW = 3;  // w
    int constexpr idT = 4;  // potential temperature perturbation
    auto num_state   = state_limits_x .extent(0);
    auto num_tracers = tracer_limits_x.extent(0);
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
    bool sim2d       = ny == 1;
    parallel_for( SimpleBounds<4>(nz+1,ny+1,nx+1,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      ////////////////////////////////
      // x-direction
      ////////////////////////////////
      if (j < ny && k < nz) {
        if (i == 0 ) {
          for (int l=0; l < num_state  ; l++) { state_limits_x (l,0,k,j,0 ,iens) = state_limits_x (l,0,k,j,nx,iens); }
          for (int l=0; l < num_tracers; l++) { tracer_limits_x(l,0,k,j,0 ,iens) = tracer_limits_x(l,0,k,j,nx,iens); }
        }
        if (i == nx) {
          for (int l=0; l < num_state  ; l++) { state_limits_x (l,1,k,j,nx,iens) = state_limits_x (l,1,k,j,0 ,iens); }
          for (int l=0; l < num_tracers; l++) { tracer_limits_x(l,1,k,j,nx,iens) = tracer_limits_x(l,1,k,j,0 ,iens); }
        }
        // Get left and right state
        real r_L = state_limits_x(idR,0,k,j,i,iens)    ;   real r_R = state_limits_x(idR,1,k,j,i,iens)    ;
        real u_L = state_limits_x(idU,0,k,j,i,iens)/r_L;   real u_R = state_limits_x(idU,1,k,j,i,iens)/r_R;
        real v_L = state_limits_x(idV,0,k,j,i,iens)/r_L;   real v_R = state_limits_x(idV,1,k,j,i,iens)/r_R;
        real w_L = state_limits_x(idW,0,k,j,i,iens)/r_L;   real w_R = state_limits_x(idW,1,k,j,i,iens)/r_R;
        real t_L = state_limits_x(idT,0,k,j,i,iens)/r_L;   real t_R = state_limits_x(idT,1,k,j,i,iens)/r_R;
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
        real q1_L = state_limits_x(idR,0,k,j,i,iens);   real q1_R = state_limits_x(idR,1,k,j,i,iens);
        real q2_L = state_limits_x(idU,0,k,j,i,iens);   real q2_R = state_limits_x(idU,1,k,j,i,iens);
        real q3_L = state_limits_x(idV,0,k,j,i,iens);   real q3_R = state_limits_x(idV,1,k,j,i,iens);
        real q4_L = state_limits_x(idW,0,k,j,i,iens);   real q4_R = state_limits_x(idW,1,k,j,i,iens);
        real q5_L = state_limits_x(idT,0,k,j,i,iens);   real q5_R = state_limits_x(idT,1,k,j,i,iens);
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

        state_flux_x(idR,k,j,i,iens) = q2;
        state_flux_x(idU,k,j,i,iens) = q2*q2/q1 + C0*pow(q5,gamma);
        state_flux_x(idV,k,j,i,iens) = q2*q3/q1;
        state_flux_x(idW,k,j,i,iens) = q2*q4/q1;
        state_flux_x(idT,k,j,i,iens) = q2*q5/q1;

        // COMPUTE UPWIND TRACER FLUXES
        // Handle it one tracer at a time
        for (int tr=0; tr < num_tracers; tr++) {
          if (u > 0) {
            tracer_flux_x(tr,k,j,i,iens) = q2 * tracer_limits_x(tr,0,k,j,i,iens) / r_L;
          } else {
            tracer_flux_x(tr,k,j,i,iens) = q2 * tracer_limits_x(tr,1,k,j,i,iens) / r_R;
          }
        }
      }

      ////////////////////////////////
      // y-direction
      ////////////////////////////////
      if (i < nx && k < nz && (! sim2d)) {
        if (j == 0 ) {
          for (int l=0; l < num_state  ; l++) { state_limits_y (l,0,k,0 ,i,iens) = state_limits_y (l,0,k,ny,i,iens); }
          for (int l=0; l < num_tracers; l++) { tracer_limits_y(l,0,k,0 ,i,iens) = tracer_limits_y(l,0,k,ny,i,iens); }
        }
        if (j == ny) {
          for (int l=0; l < num_state  ; l++) { state_limits_y (l,1,k,ny,i,iens) = state_limits_y (l,1,k,0 ,i,iens); }
          for (int l=0; l < num_tracers; l++) { tracer_limits_y(l,1,k,ny,i,iens) = tracer_limits_y(l,1,k,0 ,i,iens); }
        }
        // Get left and right state
        real r_L = state_limits_y(idR,0,k,j,i,iens)    ;   real r_R = state_limits_y(idR,1,k,j,i,iens)    ;
        real u_L = state_limits_y(idU,0,k,j,i,iens)/r_L;   real u_R = state_limits_y(idU,1,k,j,i,iens)/r_R;
        real v_L = state_limits_y(idV,0,k,j,i,iens)/r_L;   real v_R = state_limits_y(idV,1,k,j,i,iens)/r_R;
        real w_L = state_limits_y(idW,0,k,j,i,iens)/r_L;   real w_R = state_limits_y(idW,1,k,j,i,iens)/r_R;
        real t_L = state_limits_y(idT,0,k,j,i,iens)/r_L;   real t_R = state_limits_y(idT,1,k,j,i,iens)/r_R;
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
        real q1_L = state_limits_y(idR,0,k,j,i,iens);   real q1_R = state_limits_y(idR,1,k,j,i,iens);
        real q2_L = state_limits_y(idU,0,k,j,i,iens);   real q2_R = state_limits_y(idU,1,k,j,i,iens);
        real q3_L = state_limits_y(idV,0,k,j,i,iens);   real q3_R = state_limits_y(idV,1,k,j,i,iens);
        real q4_L = state_limits_y(idW,0,k,j,i,iens);   real q4_R = state_limits_y(idW,1,k,j,i,iens);
        real q5_L = state_limits_y(idT,0,k,j,i,iens);   real q5_R = state_limits_y(idT,1,k,j,i,iens);
        // Compute upwind characteristics
        // Waves 1-3, velocity: v
        real w1, w2, w3;
        if (v > 0) {
          w1 = q1_L - q5_L/t;
          w2 = q2_L - u*q5_L/t;
          w3 = q4_L - w*q5_L/t;
        } else {
          w1 = q1_R - q5_R/t;
          w2 = q2_R - u*q5_R/t;
          w3 = q4_R - w*q5_R/t;
        }
        // Wave 5, velocity: v-cs
        real w5 =  v*q1_R/(2*cs) - q3_R/(2*cs) + q5_R/(2*t);
        // Wave 6, velocity: v+cs
        real w6 = -v*q1_L/(2*cs) + q3_L/(2*cs) + q5_L/(2*t);
        // Use right eigenmatrix to compute upwind flux
        real q1 = w1 + w5 + w6;
        real q2 = w2 + u*w5 + u*w6;
        real q3 = v*w1 + (v-cs)*w5 + (v+cs)*w6;
        real q4 = w3 + w*w5 + w*w6;
        real q5 =      t*w5 + t*w6;

        state_flux_y(idR,k,j,i,iens) = q3;
        state_flux_y(idU,k,j,i,iens) = q3*q2/q1;
        state_flux_y(idV,k,j,i,iens) = q3*q3/q1 + C0*pow(q5,gamma);
        state_flux_y(idW,k,j,i,iens) = q3*q4/q1;
        state_flux_y(idT,k,j,i,iens) = q3*q5/q1;

        // COMPUTE UPWIND TRACER FLUXES
        // Handle it one tracer at a time
        for (int tr=0; tr < num_tracers; tr++) {
          if (v > 0) {
            tracer_flux_y(tr,k,j,i,iens) = q3 * tracer_limits_y(tr,0,k,j,i,iens) / r_L;
          } else {
            tracer_flux_y(tr,k,j,i,iens) = q3 * tracer_limits_y(tr,1,k,j,i,iens) / r_R;
          }
        }
      }

      ////////////////////////////////
      // z-direction
      ////////////////////////////////
      if (i < nx && j < ny) {
        if (k == 0) {
          for (int l = 0; l < num_state  ; l++) { state_limits_z (l,0,0 ,j,i,iens) = state_limits_z (l,1,0 ,j,i,iens); }
          for (int l = 0; l < num_tracers; l++) { tracer_limits_z(l,0,0 ,j,i,iens) = tracer_limits_z(l,1,0 ,j,i,iens); }
        }
        if (k == nz) {
          for (int l = 0; l < num_state  ; l++) { state_limits_z (l,1,nz,j,i,iens) = state_limits_z (l,0,nz,j,i,iens); }
          for (int l = 0; l < num_tracers; l++) { tracer_limits_z(l,1,nz,j,i,iens) = tracer_limits_z(l,0,nz,j,i,iens); }
        }
        // Get left and right state
        real r_L = state_limits_z(idR,0,k,j,i,iens)    ;   real r_R = state_limits_z(idR,1,k,j,i,iens)    ;
        real u_L = state_limits_z(idU,0,k,j,i,iens)/r_L;   real u_R = state_limits_z(idU,1,k,j,i,iens)/r_R;
        real v_L = state_limits_z(idV,0,k,j,i,iens)/r_L;   real v_R = state_limits_z(idV,1,k,j,i,iens)/r_R;
        real w_L = state_limits_z(idW,0,k,j,i,iens)/r_L;   real w_R = state_limits_z(idW,1,k,j,i,iens)/r_R;
        real t_L = state_limits_z(idT,0,k,j,i,iens)/r_L;   real t_R = state_limits_z(idT,1,k,j,i,iens)/r_R;
        // Compute average state
        real r = 0.5_fp * (r_L + r_R);
        real u = 0.5_fp * (u_L + u_R);
        real v = 0.5_fp * (v_L + v_R);
        real w = 0.5_fp * (w_L + w_R);
        real t = 0.5_fp * (t_L + t_R);
        real p = C0 * pow(r*t,gamma);
        real cs2 = gamma*p/r;
        real cs  = sqrt(cs2);
        // Get left and right fluxes
        real q1_L = state_limits_z(idR,0,k,j,i,iens);   real q1_R = state_limits_z(idR,1,k,j,i,iens);
        real q2_L = state_limits_z(idU,0,k,j,i,iens);   real q2_R = state_limits_z(idU,1,k,j,i,iens);
        real q3_L = state_limits_z(idV,0,k,j,i,iens);   real q3_R = state_limits_z(idV,1,k,j,i,iens);
        real q4_L = state_limits_z(idW,0,k,j,i,iens);   real q4_R = state_limits_z(idW,1,k,j,i,iens);
        real q5_L = state_limits_z(idT,0,k,j,i,iens);   real q5_R = state_limits_z(idT,1,k,j,i,iens);
        // Compute upwind characteristics
        // Waves 1-3, velocity: w
        real w1, w2, w3;
        if (w > 0) {
          w1 = q1_L - q5_L/t;
          w2 = q2_L - u*q5_L/t;
          w3 = q3_L - v*q5_L/t;
        } else {
          w1 = q1_R - q5_R/t;
          w2 = q2_R - u*q5_R/t;
          w3 = q3_R - v*q5_R/t;
        }
        // Wave 5, velocity: w-cs
        real w5 =  w*q1_R/(2*cs) - q4_R/(2*cs) + q5_R/(2*t);
        // Wave 6, velocity: w+cs
        real w6 = -w*q1_L/(2*cs) + q4_L/(2*cs) + q5_L/(2*t);
        // Use right eigenmatrix to compute upwind flux
        real q1 = w1 + w5 + w6;
        real q2 = w2 + u*w5 + u*w6;
        real q3 = w3 + v*w5 + v*w6;
        real q4 = w*w1 + (w-cs)*w5 + (w+cs)*w6;
        real q5 =      t*w5 + t*w6;

        state_flux_z(idR,k,j,i,iens) = q4;
        state_flux_z(idU,k,j,i,iens) = q4*q2/q1;
        state_flux_z(idV,k,j,i,iens) = q4*q3/q1;
        state_flux_z(idW,k,j,i,iens) = q4*q4/q1 + C0*pow(q5,gamma);
        state_flux_z(idT,k,j,i,iens) = q4*q5/q1;

        // COMPUTE UPWIND TRACER FLUXES
        // Handle it one tracer at a time
        for (int tr=0; tr < num_tracers; tr++) {
          if (w > 0) {
            tracer_flux_z(tr,k,j,i,iens) = q4 * tracer_limits_z(tr,0,k,j,i,iens) / r_L;
          } else {
            tracer_flux_z(tr,k,j,i,iens) = q4 * tracer_limits_z(tr,1,k,j,i,iens) / r_R;
          }
        }
      }

    });
  }

}


