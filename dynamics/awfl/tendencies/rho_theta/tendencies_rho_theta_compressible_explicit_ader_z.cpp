
#include "tendencies_rho_theta_compressible_explicit_ader.h"


namespace awfl {
namespace tendencies_rho_theta {
namespace compressible_explicit_ader {


  void compute_tendencies_z( pam::PamCoupler const &coupler ,
                             real5d const &state   , real5d const &state_tend       ,
                             real5d const &tracers , real5d const &tracer_tend      ,
                             awfl::Recon const &recon ,
                             awfl::tendencies_rho_theta::Hydrostasis const &hydrostasis ,
                             real &dt ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    using awfl::tendencies_rho_theta::num_state  ;
    using awfl::tendencies_rho_theta::max_tracers;
    using awfl::tendencies_rho_theta::idR;  // density perturbation
    using awfl::tendencies_rho_theta::idU;  // rho*u
    using awfl::tendencies_rho_theta::idV;  // rho*v
    using awfl::tendencies_rho_theta::idW;  // rho*w
    using awfl::tendencies_rho_theta::idT;  // rho*potential temperature perturbation
    using awfl::tendencies_rho_theta::Hydrostasis;
    using awfl::Recon;
    using awfl::tendencies_rho_theta::Hydrostasis;
    int static constexpr nAder = ngll;

    auto wrapz = YAKL_LAMBDA (int k, int kk, int nz) {
      int ret = k+kk;
      if (ret < hs+0   ) ret = hs+0;
      if (ret > hs+nz-1) ret = hs+nz-1;
      return ret;
    };

    auto num_tracers = coupler.get_num_tracers();
    auto nz          = coupler.get_nz  ();
    auto ny          = coupler.get_ny  ();
    auto nx          = coupler.get_nx  ();
    auto dz          = coupler.get_data_manager_readonly().get<real const,2>("vertical_cell_dz");
    auto nens        = coupler.get_nens();
    auto grav        = coupler.get_grav   ();
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
    YAKL_SCOPE( idl                     , recon.idl                    );
    YAKL_SCOPE( sigma                   , recon.sigma                  );
    YAKL_SCOPE( derivMatrix             , recon.derivMatrix            );
    YAKL_SCOPE( vert_sten_to_gll        , recon.vert_sten_to_gll       );
    YAKL_SCOPE( vert_sten_to_coefs      , recon.vert_sten_to_coefs     );
    YAKL_SCOPE( vert_weno_recon_lower   , recon.vert_weno_recon_lower  );
    YAKL_SCOPE( hyDensSten              , hydrostasis.hyDensSten       );
    YAKL_SCOPE( hyDensThetaSten         , hydrostasis.hyDensThetaSten  );
    YAKL_SCOPE( hyDensGLL               , hydrostasis.hyDensGLL        );
    YAKL_SCOPE( hyDensThetaGLL          , hydrostasis.hyDensThetaGLL   );

    SArray<real,1,ngll> gllWts_ngll, gllPts_ngll;
    TransformMatrices::get_gll_points (gllPts_ngll);
    TransformMatrices::get_gll_weights(gllWts_ngll);

    real6d state_limits ("state_limits" ,num_state  ,2,nz+1,ny,nx,nens);
    real6d tracer_limits("tracer_limits",num_tracers,2,nz+1,ny,nx,nens);

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( "Spatial.h Z recon" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      SArray<real,2,ord,ngll>       s2g_loc;
      SArray<real,2,ord,ord>        s2c_loc;
      SArray<real,3,hs+1,hs+1,hs+1> weno_recon_lower_loc;
      for (int jj=0; jj < ord; jj++) {
        for (int ii=0; ii < ngll; ii++) {
          s2g_loc(jj,ii) = vert_sten_to_gll(k,jj,ii,iens);
        }
      }
      for (int jj=0; jj < ord; jj++) {
        for (int ii=0; ii < ord; ii++) {
          s2c_loc(jj,ii) = vert_sten_to_coefs(k,jj,ii,iens);
        }
      }
      for (int kk=0; kk < hs+1; kk++) {
        for (int jj=0; jj < hs+1; jj++) {
          for (int ii=0; ii < hs+1; ii++) {
            weno_recon_lower_loc(kk,jj,ii) = vert_weno_recon_lower(k,kk,jj,ii,iens);
          }
        }
      }

      SArray<real,2,nAder,ngll> r_DTs , rw_DTs;

      { // State
        SArray<real,2,nAder,ngll> ru_DTs , rv_DTs , rt_DTs, rwu_DTs , rwv_DTs , rww_DTs , rwt_DTs , rt_gamma_DTs;
        { // Recon
          SArray<real,1,ord>  stencil;
          SArray<real,1,ngll> gll;

          // Density
          for (int kk=0; kk < ord; kk++) {
            if (k+kk < hs || k+kk > hs+nz-1) {
              stencil(kk) = state(idR,hs+k,hs+j,hs+i,iens) - hyDensSten(k,hs,iens);
            } else {
              stencil(kk) = state(idR,k+kk,hs+j,hs+i,iens) - hyDensSten(k,kk,iens);
            }
          }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc ,
                                        idl , sigma , weno_scalars );
          // Add hydrostasis back on
          for (int kk=0; kk < ngll; kk++) { r_DTs(0,kk) = gll(kk) + hyDensGLL(k,kk,iens); }

          // u values and derivatives
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idU,wrapz(k,kk,nz),hs+j,hs+i,iens) /
                                                         state(idR,wrapz(k,kk,nz),hs+j,hs+i,iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc ,
                                        idl , sigma , weno_winds );
          for (int kk=0; kk < ngll; kk++) { ru_DTs(0,kk) = gll(kk) * r_DTs(0,kk); }

          // v
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idV,wrapz(k,kk,nz),hs+j,hs+i,iens) /
                                                         state(idR,wrapz(k,kk,nz),hs+j,hs+i,iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc ,
                                        idl , sigma , weno_winds );
          for (int kk=0; kk < ngll; kk++) { rv_DTs(0,kk) = gll(kk) * r_DTs(0,kk); }

          // w
          for (int kk=0; kk < ord; kk++) {
            stencil(kk) = state(idW,wrapz(k,kk,nz),hs+j,hs+i,iens) /
                          state(idR,wrapz(k,kk,nz),hs+j,hs+i,iens);
            if (k+kk > hs+nz-1 || k+kk < hs) stencil(kk) = 0;
          }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc ,
                                        idl , sigma , weno_winds );
          if (k == nz-1) gll(ngll-1) = 0;
          if (k == 0   ) gll(0     ) = 0;
          for (int kk=0; kk < ngll; kk++) { rw_DTs(0,kk) = gll(kk) * r_DTs(0,kk); }

          // rho*theta
          for (int kk=0; kk < ord; kk++) {
            if (k+kk < hs || k+kk > hs+nz-1) {
              stencil(kk) = state(idT,hs+k,hs+j,hs+i,iens) - hyDensThetaSten(k,hs,iens);
            } else {
              stencil(kk) = state(idT,k+kk,hs+j,hs+i,iens) - hyDensThetaSten(k,kk,iens);
            }
          }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc ,
                                        idl , sigma , weno_scalars );
          // Add hydrostasis back on
          for (int kk=0; kk < ngll; kk++) { rt_DTs(0,kk) = gll(kk) + hyDensThetaGLL(k,kk,iens); }
        }

        for (int kk=0; kk < ngll; kk++) {
          real r = r_DTs (0,kk);
          real u = ru_DTs(0,kk) / r;
          real v = rv_DTs(0,kk) / r;
          real w = rw_DTs(0,kk) / r;
          real t = rt_DTs(0,kk) / r;
          rwu_DTs     (0,kk) = r*w*u;
          rwv_DTs     (0,kk) = r*w*v;
          rww_DTs     (0,kk) = r*w*w;
          rwt_DTs     (0,kk) = r*w*t;
          rt_gamma_DTs(0,kk) = pow(r*t,gamma);
        }

        if (nAder > 1) {
          awfl::diffTransformEulerConsZ( r_DTs , ru_DTs , rv_DTs , rw_DTs , rt_DTs , rwu_DTs , rwv_DTs , rww_DTs ,
                                         rwt_DTs , rt_gamma_DTs , derivMatrix , C0 , gamma ,
                                         grav , k , dz(k,iens) , nz );
        }

        SArray<real,1,ngll> r_tavg, rw_tavg;
        awfl::compute_timeAvg( r_DTs  , r_tavg  , dt );
        awfl::compute_timeAvg( ru_DTs           , dt );
        awfl::compute_timeAvg( rv_DTs           , dt );
        awfl::compute_timeAvg( rw_DTs , rw_tavg , dt );
        awfl::compute_timeAvg( rt_DTs           , dt );

        // Left interface
        state_limits(idR,1,k  ,j,i,iens) = r_tavg  (0     );
        state_limits(idU,1,k  ,j,i,iens) = ru_DTs(0,0     );
        state_limits(idV,1,k  ,j,i,iens) = rv_DTs(0,0     );
        state_limits(idW,1,k  ,j,i,iens) = rw_tavg (0     );
        state_limits(idT,1,k  ,j,i,iens) = rt_DTs(0,0     );
        // Right interface
        state_limits(idR,0,k+1,j,i,iens) = r_tavg  (ngll-1);
        state_limits(idU,0,k+1,j,i,iens) = ru_DTs(0,ngll-1);
        state_limits(idV,0,k+1,j,i,iens) = rv_DTs(0,ngll-1);
        state_limits(idW,0,k+1,j,i,iens) = rw_tavg (ngll-1);
        state_limits(idT,0,k+1,j,i,iens) = rt_DTs(0,ngll-1);

        real ravg = 0;
        for (int kk=0; kk < ngll; kk++) {
          ravg += r_tavg(kk) * gllWts_ngll(kk);
        }
        state_tend(idR,k,j,i,iens) = 0;
        state_tend(idU,k,j,i,iens) = 0;
        state_tend(idV,k,j,i,iens) = 0;
        state_tend(idW,k,j,i,iens) = -grav * ravg;
        state_tend(idT,k,j,i,iens) = 0;
      }

      { // Tracers
        for (int tr=0; tr < num_tracers; tr++) {
          SArray<real,2,nAder,ngll> rt_DTs;  // Density * tracer
          SArray<real,2,nAder,ngll> rwt_DTs; // Density * wwind * tracer
          { // Recon
            SArray<real,1,ord>  stencil;
            SArray<real,1,ngll> gll;

            for (int kk=0; kk < ord; kk++) { stencil(kk) = tracers(tr ,wrapz(k,kk,nz),hs+j,hs+i,iens) /
                                                           state  (idR,wrapz(k,kk,nz),hs+j,hs+i,iens); }
            awfl::reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc ,
                                          idl , sigma , weno_scalars );
            for (int kk=0; kk < ngll; kk++) { gll(kk) *= r_DTs(0,kk); }
            if (tracer_pos(tr)) {
              for (int kk=0; kk < ngll; kk++) { gll(kk) = max( 0._fp , gll(kk) ); }
            }
            for (int kk=0; kk < ngll; kk++) { rt_DTs(0,kk) = gll(kk); }
          }

          for (int kk=0; kk < ngll; kk++) {
            rwt_DTs(0,kk) = rt_DTs(0,kk) * rw_DTs(0,kk) / r_DTs(0,kk);
          }

          if (nAder > 1) {
            awfl::diffTransformTracer( r_DTs , rw_DTs , rt_DTs , rwt_DTs , derivMatrix , dz(k,iens) );
          }

          awfl::compute_timeAvg( rt_DTs  , dt );

          if (tracer_pos(tr)) {
            for (int kk=0; kk < ngll; kk++) { rt_DTs(0,kk) = max( 0._fp , rt_DTs(0,kk) ); }
          }

          if (k == nz-1) rwt_DTs(0,ngll-1) = 0;
          if (k == 0   ) rwt_DTs(0,0     ) = 0;

          tracer_limits(tr,1,k  ,j,i,iens) = rt_DTs (0,0     ); // Left interface
          tracer_limits(tr,0,k+1,j,i,iens) = rt_DTs (0,ngll-1); // Right interface
        }
      }
    });

    real5d state_flux ("state_flux" ,num_state  ,nz+1,ny,nx,nens);
    real5d tracer_flux("tracer_flux",num_tracers,nz+1,ny,nx,nens);

    awfl::riemann_rho_theta_full_z(coupler, state_limits , state_flux , tracer_limits , tracer_flux );

    state_limits  = real6d();
    tracer_limits = real6d();

    awfl::fct_positivity_z( coupler , tracers , tracer_flux , tracer_pos , dt );

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( "Spatial.h Z tendencies" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l=0; l < num_state; l++) {
        if (sim2d && l == idV) {
          state_tend(l,k,j,i,iens) = 0;
        } else {
          state_tend(l,k,j,i,iens) += - ( state_flux(l,k+1,j,i,iens) - state_flux(l,k,j,i,iens) ) / dz(k,iens);
        }
      }
      for (int l=0; l < num_tracers; l++) {
        // Compute tracer tendency
        tracer_tend(l,k,j,i,iens) = - ( tracer_flux(l,k+1,j,i,iens) - tracer_flux(l,k,j,i,iens) ) / dz(k,iens);
      }
    });
  }


}
}
}


