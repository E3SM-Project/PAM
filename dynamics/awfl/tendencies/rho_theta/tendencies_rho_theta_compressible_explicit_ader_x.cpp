
#include "tendencies_rho_theta_compressible_explicit_ader.h"


namespace awfl {
namespace tendencies_rho_theta {
namespace compressible_explicit_ader {


  void compute_tendencies_x( pam::PamCoupler const &coupler ,
                             real5d const &state   , real5d const &state_tend  ,
                             real5d const &tracers , real5d const &tracer_tend ,
                             awfl::Recon const &recon , real &dt ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    using awfl::tendencies_rho_theta::num_state  ;
    using awfl::tendencies_rho_theta::max_tracers;
    using awfl::tendencies_rho_theta::idR;  // density perturbation
    using awfl::tendencies_rho_theta::idU;  // rho*u
    using awfl::tendencies_rho_theta::idV;  // rho*v
    using awfl::tendencies_rho_theta::idW;  // rho*w
    using awfl::tendencies_rho_theta::idT;  // rho*potential temperature perturbation
    int static constexpr nAder = ngll;

    auto wrapx = YAKL_LAMBDA (int i, int ii, int nx) {
      int ret = i+ii;
      if (ret < hs+0   ) ret += nx;
      if (ret > hs+nx-1) ret -= nx;
      return ret;
    };

    auto num_tracers = coupler.get_num_tracers();
    auto nz          = coupler.get_nz     ();
    auto ny          = coupler.get_ny     ();
    auto nx          = coupler.get_nx     ();
    auto dx          = coupler.get_dx     ();
    auto nens        = coupler.get_nens   ();
    auto Rd          = coupler.get_R_d    ();
    auto cp          = coupler.get_cp_d   ();
    auto p0          = coupler.get_p0     ();
    auto gamma       = coupler.get_gamma_d();
    auto kappa       = Rd/cp;
    auto C0          = pow( Rd * pow( p0 , -kappa ) , gamma );
    auto tracer_pos  = coupler.get_tracer_positivity_array();
    auto sim2d = ny == 1;

    YAKL_SCOPE( weno_scalars            , recon.weno_scalars           );
    YAKL_SCOPE( weno_winds              , recon.weno_winds             );
    YAKL_SCOPE( c2g                     , recon.coefs_to_gll           );
    YAKL_SCOPE( s2g                     , recon.sten_to_gll            );
    YAKL_SCOPE( s2c                     , recon.sten_to_coefs          );
    YAKL_SCOPE( weno_recon_lower        , recon.weno_recon_lower       );
    YAKL_SCOPE( idl                     , recon.idl                    );
    YAKL_SCOPE( sigma                   , recon.sigma                  );
    YAKL_SCOPE( derivMatrix             , recon.derivMatrix            );

    real6d state_limits ("state_limits" ,num_state  ,2,nz,ny,nx+1,nens);
    real6d tracer_limits("tracer_limits",num_tracers,2,nz,ny,nx+1,nens);

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( "Spatial.h X recon" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      SArray<real,2,nAder,ngll> r_DTs , ru_DTs;

      { // State
        SArray<real,2,nAder,ngll> rv_DTs , rw_DTs , rt_DTs, ruu_DTs , ruv_DTs , ruw_DTs , rut_DTs , rt_gamma_DTs;

        { // Recon
          SArray<real,1,ord>  stencil;
          SArray<real,1,ngll> gll;

          // Density
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idR,hs+k,hs+j,wrapx(i,ii,nx),iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
          // Add hydrostasis back on
          for (int ii=0; ii < ngll; ii++) { r_DTs(0,ii) = gll(ii); }

          // u
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idU,hs+k,hs+j,wrapx(i,ii,nx),iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
          for (int ii=0; ii < ngll; ii++) { ru_DTs(0,ii) = gll(ii); }

          // v
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idV,hs+k,hs+j,wrapx(i,ii,nx),iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
          for (int ii=0; ii < ngll; ii++) { rv_DTs(0,ii) = gll(ii); }

          // w
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idW,hs+k,hs+j,wrapx(i,ii,nx),iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
          for (int ii=0; ii < ngll; ii++) { rw_DTs(0,ii) = gll(ii); }

          // theta
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idT,hs+k,hs+j,wrapx(i,ii,nx),iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
          // Add hydrostasis back on
          for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) = gll(ii); }
        }

        for (int ii=0; ii < ngll; ii++) {
          real r = r_DTs (0,ii);
          real u = ru_DTs(0,ii) / r;
          real v = rv_DTs(0,ii) / r;
          real w = rw_DTs(0,ii) / r;
          real t = rt_DTs(0,ii) / r;
          ruu_DTs     (0,ii) = r*u*u;
          ruv_DTs     (0,ii) = r*u*v;
          ruw_DTs     (0,ii) = r*u*w;
          rut_DTs     (0,ii) = r*u*t;
          rt_gamma_DTs(0,ii) = pow(r*t,gamma);
        }

        if (nAder > 1) {
          awfl::diffTransformEulerConsX( r_DTs , ru_DTs , rv_DTs , rw_DTs , rt_DTs , ruu_DTs , ruv_DTs , ruw_DTs ,
                                         rut_DTs , rt_gamma_DTs , derivMatrix , C0 , gamma , dx );
        }

        SArray<real,1,ngll> r_tavg, ru_tavg;
        awfl::compute_timeAvg( r_DTs  , r_tavg  , dt );
        awfl::compute_timeAvg( ru_DTs , ru_tavg , dt );
        awfl::compute_timeAvg( rv_DTs           , dt );
        awfl::compute_timeAvg( rw_DTs           , dt );
        awfl::compute_timeAvg( rt_DTs           , dt );

        // Left interface
        state_limits(idR,1,k,j,i  ,iens) = r_tavg  (0     );
        state_limits(idU,1,k,j,i  ,iens) = ru_tavg (0     );
        state_limits(idV,1,k,j,i  ,iens) = rv_DTs(0,0     );
        state_limits(idW,1,k,j,i  ,iens) = rw_DTs(0,0     );
        state_limits(idT,1,k,j,i  ,iens) = rt_DTs(0,0     );
        // Right interface
        state_limits(idR,0,k,j,i+1,iens) = r_tavg  (ngll-1);
        state_limits(idU,0,k,j,i+1,iens) = ru_tavg (ngll-1);
        state_limits(idV,0,k,j,i+1,iens) = rv_DTs(0,ngll-1);
        state_limits(idW,0,k,j,i+1,iens) = rw_DTs(0,ngll-1);
        state_limits(idT,0,k,j,i+1,iens) = rt_DTs(0,ngll-1);
      }

      { // Tracers
        for (int tr=0; tr < num_tracers; tr++) {
          SArray<real,2,nAder,ngll> rt_DTs, rut_DTs;

          { // Recon
            SArray<real,1,ord>  stencil;
            SArray<real,1,ngll> gll;

            for (int ii=0; ii < ord; ii++) { stencil(ii) = tracers(tr,hs+k,hs+j,wrapx(i,ii,nx),iens); }
            awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
            if (tracer_pos(tr)) {
              for (int ii=0; ii < ngll; ii++) { gll(ii) = max( 0._fp , gll(ii) ); }
            }
            for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) = gll(ii); }
          }

          for (int ii=0; ii < ngll; ii++) {
            rut_DTs(0,ii) = rt_DTs(0,ii) * ru_DTs(0,ii) / r_DTs(0,ii);
          }

          if (nAder > 1) {
            awfl::diffTransformTracer( r_DTs , ru_DTs , rt_DTs , rut_DTs , derivMatrix , dx );
          }

          awfl::compute_timeAvg( rt_DTs  , dt );

          if (tracer_pos(tr)) {
            for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) = max( 0._fp , rt_DTs(0,ii) ); }
          }

          tracer_limits(tr,1,k,j,i  ,iens) = rt_DTs(0,0     ); // Left interface
          tracer_limits(tr,0,k,j,i+1,iens) = rt_DTs(0,ngll-1); // Right interface
        }
      }
    });

    real5d state_flux ("state_flux" ,num_state  ,nz,ny,nx+1,nens);
    real5d tracer_flux("tracer_flux",num_tracers,nz,ny,nx+1,nens);

    awfl::riemann_rho_theta_full_x( coupler , state_limits , state_flux , tracer_limits , tracer_flux );

    state_limits  = real6d();
    tracer_limits = real6d();

    awfl::fct_positivity_x( coupler , tracers , tracer_flux , tracer_pos , dt );

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( "Spatial.h X tendencies" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l = 0; l < num_state; l++) {
        if (sim2d && l == idV) {
          state_tend(l,k,j,i,iens) = 0;
        } else {
          state_tend(l,k,j,i,iens) = - ( state_flux(l,k,j,i+1,iens) - state_flux(l,k,j,i,iens) ) / dx;
        }
      }
      for (int l = 0; l < num_tracers; l++) {
        // Compute tracer tendency
        tracer_tend(l,k,j,i,iens) = - ( tracer_flux(l,k,j,i+1,iens) - tracer_flux(l,k,j,i,iens) ) / dx;
      }
    });
  }


}
}
}


