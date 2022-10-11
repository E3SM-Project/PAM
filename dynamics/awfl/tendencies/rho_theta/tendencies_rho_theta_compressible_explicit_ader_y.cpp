
#include "tendencies_rho_theta_compressible_explicit_ader.h"


namespace awfl {
namespace tendencies_rho_theta {
namespace compressible_explicit_ader {


  void compute_tendencies_y( pam::PamCoupler const &coupler ,
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

    auto wrapy = YAKL_LAMBDA (int j, int jj, int ny) {
      int ret = j+jj;
      if (ret < hs+0   ) ret += ny;
      if (ret > hs+ny-1) ret -= ny;
      return ret;
    };

    auto num_tracers = coupler.get_num_tracers();
    auto nz          = coupler.get_nz     ();
    auto ny          = coupler.get_ny     ();
    auto nx          = coupler.get_nx     ();
    auto dy          = coupler.get_dy     ();
    auto nens        = coupler.get_nens   ();
    auto Rd          = coupler.get_R_d    ();
    auto cp          = coupler.get_cp_d   ();
    auto p0          = coupler.get_p0     ();
    auto gamma       = coupler.get_gamma_d();
    auto kappa       = Rd/cp;
    auto C0          = pow( Rd * pow( p0 , -kappa ) , gamma );
    auto tracer_pos  = coupler.get_tracer_positivity_array();
    auto sim2d = ny == 1;

    if (sim2d) {
      state_tend  = 0;
      tracer_tend = 0;
      return;
    }

    YAKL_SCOPE( weno_scalars            , recon.weno_scalars           );
    YAKL_SCOPE( weno_winds              , recon.weno_winds             );
    YAKL_SCOPE( c2g                     , recon.coefs_to_gll           );
    YAKL_SCOPE( s2g                     , recon.sten_to_gll            );
    YAKL_SCOPE( s2c                     , recon.sten_to_coefs          );
    YAKL_SCOPE( weno_recon_lower        , recon.weno_recon_lower       );
    YAKL_SCOPE( idl                     , recon.idl                    );
    YAKL_SCOPE( sigma                   , recon.sigma                  );
    YAKL_SCOPE( derivMatrix             , recon.derivMatrix            );

    real6d state_limits ("state_limits" ,num_state  ,2,nz,ny+1,nx,nens);
    real6d tracer_limits("tracer_limits",num_tracers,2,nz,ny+1,nx,nens);

    // Loop through all cells, reconstruct in y-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( "Spatial.h Y recon" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      SArray<real,2,nAder,ngll> r_DTs , rv_DTs;

      { // State
        SArray<real,2,nAder,ngll> ru_DTs , rw_DTs , rt_DTs, rvu_DTs , rvv_DTs , rvw_DTs , rvt_DTs , rt_gamma_DTs;

        { // Recon
          SArray<real,1,ord>  stencil;
          SArray<real,1,ngll> gll;

          // Density
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idR,hs+k,wrapy(j,jj,ny),hs+i,iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
          // Add hydrostasis back on
          for (int jj=0; jj < ngll; jj++) { r_DTs(0,jj) = gll(jj); }

          // u
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idU,hs+k,wrapy(j,jj,ny),hs+i,iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
          for (int jj=0; jj < ngll; jj++) { ru_DTs(0,jj) = gll(jj); }

          // v
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idV,hs+k,wrapy(j,jj,ny),hs+i,iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
          for (int jj=0; jj < ngll; jj++) { rv_DTs(0,jj) = gll(jj); }

          // w
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idW,hs+k,wrapy(j,jj,ny),hs+i,iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
          for (int jj=0; jj < ngll; jj++) { rw_DTs(0,jj) = gll(jj); }

          // theta
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idT,hs+k,wrapy(j,jj,ny),hs+i,iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
          // Add hydrostasis back on
          for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) = gll(jj); }
        }

        for (int jj=0; jj < ngll; jj++) {
          real r = r_DTs (0,jj);
          real u = ru_DTs(0,jj) / r;
          real v = rv_DTs(0,jj) / r;
          real w = rw_DTs(0,jj) / r;
          real t = rt_DTs(0,jj) / r;
          rvu_DTs     (0,jj) = r*v*u;
          rvv_DTs     (0,jj) = r*v*v;
          rvw_DTs     (0,jj) = r*v*w;
          rvt_DTs     (0,jj) = r*v*t;
          rt_gamma_DTs(0,jj) = pow(r*t,gamma);
        }

        if (nAder > 1) {
          awfl::diffTransformEulerConsY( r_DTs , ru_DTs , rv_DTs , rw_DTs , rt_DTs , rvu_DTs , rvv_DTs , rvw_DTs ,
                                         rvt_DTs , rt_gamma_DTs , derivMatrix , C0 , gamma , dy );
        }

        SArray<real,1,ngll> r_tavg, rv_tavg;
        awfl::compute_timeAvg( r_DTs  , r_tavg  , dt );
        awfl::compute_timeAvg( ru_DTs           , dt );
        awfl::compute_timeAvg( rv_DTs , rv_tavg , dt );
        awfl::compute_timeAvg( rw_DTs           , dt );
        awfl::compute_timeAvg( rt_DTs           , dt );

        // Left interface
        state_limits(idR,1,k,j  ,i,iens) = r_tavg  (0     );
        state_limits(idU,1,k,j  ,i,iens) = ru_DTs(0,0     );
        state_limits(idV,1,k,j  ,i,iens) = rv_tavg (0     );
        state_limits(idW,1,k,j  ,i,iens) = rw_DTs(0,0     );
        state_limits(idT,1,k,j  ,i,iens) = rt_DTs(0,0     );
        // Right interface
        state_limits(idR,0,k,j+1,i,iens) = r_tavg  (ngll-1);
        state_limits(idU,0,k,j+1,i,iens) = ru_DTs(0,ngll-1);
        state_limits(idV,0,k,j+1,i,iens) = rv_tavg (ngll-1);
        state_limits(idW,0,k,j+1,i,iens) = rw_DTs(0,ngll-1);
        state_limits(idT,0,k,j+1,i,iens) = rt_DTs(0,ngll-1);
      }

      { // Tracers
        for (int tr=0; tr < num_tracers; tr++) {
          SArray<real,2,nAder,ngll> rt_DTs, rvt_DTs;
          
          { // Recon
            SArray<real,1,ord>  stencil;
            SArray<real,1,ngll> gll;

            for (int jj=0; jj < ord; jj++) { stencil(jj) = tracers(tr,hs+k,wrapy(j,jj,ny),hs+i,iens); }
            awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
            if (tracer_pos(tr)) {
              for (int jj=0; jj < ngll; jj++) { gll(jj) = max( 0._fp , gll(jj) ); }
            }
            for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) = gll(jj); }
          }

          for (int jj=0; jj < ngll; jj++) {
            rvt_DTs(0,jj) = rt_DTs(0,jj) * rv_DTs(0,jj) / r_DTs(0,jj);
          }

          if (nAder > 1) {
            awfl::diffTransformTracer( r_DTs , rv_DTs , rt_DTs , rvt_DTs , derivMatrix , dy );
          }

          awfl::compute_timeAvg( rt_DTs  , dt );

          if (tracer_pos(tr)) {
            for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) = max( 0._fp , rt_DTs(0,jj) ); }
          }

          tracer_limits(tr,1,k,j  ,i,iens) = rt_DTs (0,0     ); // Left interface
          tracer_limits(tr,0,k,j+1,i,iens) = rt_DTs (0,ngll-1); // Right interface
        }
      }
    });

    real5d state_flux ("state_flux" ,num_state  ,nz,ny+1,nx,nens);
    real5d tracer_flux("tracer_flux",num_tracers,nz,ny+1,nx,nens);

    awfl::riemann_rho_theta_full_y( coupler , state_limits , state_flux , tracer_limits , tracer_flux );

    state_limits  = real6d();
    tracer_limits = real6d();

    awfl::fct_positivity_y( coupler, tracers, tracer_flux , tracer_pos, dt);

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( "Spatial.h Y tendendies" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l=0; l < num_state; l++) {
        state_tend(l,k,j,i,iens) = - ( state_flux(l,k,j+1,i,iens) - state_flux(l,k,j,i,iens) ) / dy;
      }
      for (int l=0; l < num_tracers; l++) {
        // Compute the tracer tendency
        tracer_tend(l,k,j,i,iens) = - ( tracer_flux(l,k,j+1,i,iens) - tracer_flux(l,k,j,i,iens) ) / dy;
      }
    });
  }


}
}
}


