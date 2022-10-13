
#include "tendencies_rho_theta_compressible_explicit_semidiscrete.h"


namespace awfl {
namespace tendencies_rho_theta {
namespace compressible_explicit_semidiscrete {


  void compute_tendencies_xyz( pam::PamCoupler const &coupler ,
                               realConst5d state   , real5d const &state_tend       ,
                               realConst5d tracers , real5d const &tracer_tend      ,
                               realConst5d tracers_start ,
                               awfl::Recon const &recon ,
                               awfl::tendencies_rho_theta::Hydrostasis const &hydrostasis ,
                               real dt ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    using awfl::tendencies_rho_theta::num_state  ;
    using awfl::tendencies_rho_theta::max_tracers;
    using awfl::tendencies_rho_theta::idR;  // density
    using awfl::tendencies_rho_theta::idU;  // rho*u
    using awfl::tendencies_rho_theta::idV;  // rho*v
    using awfl::tendencies_rho_theta::idW;  // rho*w
    using awfl::tendencies_rho_theta::idT;  // rho*potential temperature
    using awfl::reconstruct_gll_values;

    auto wrapx = YAKL_LAMBDA (int i, int ii, int nx) {
      int ret = i+ii;
      if (ret < hs+0   ) ret += nx;
      if (ret > hs+nx-1) ret -= nx;
      return ret;
    };
    auto wrapy = YAKL_LAMBDA (int j, int jj, int ny) {
      int ret = j+jj;
      if (ret < hs+0   ) ret += ny;
      if (ret > hs+ny-1) ret -= ny;
      return ret;
    };
    auto wrapz = YAKL_LAMBDA (int k, int kk, int nz) {
      int ret = k+kk;
      if (ret < hs+0   ) ret = hs+0;
      if (ret > hs+nz-1) ret = hs+nz-1;
      return ret;
    };

    auto num_tracers = coupler.get_num_tracers();
    auto nz          = coupler.get_nz     ();
    auto ny          = coupler.get_ny     ();
    auto nx          = coupler.get_nx     ();
    auto dx          = coupler.get_dx     ();
    auto dy          = coupler.get_dy     ();
    auto dz          = coupler.get_data_manager_readonly().get<real const,2>("vertical_cell_dz");
    auto nens        = coupler.get_nens   ();
    auto Rd          = coupler.get_R_d    ();
    auto cp          = coupler.get_cp_d   ();
    auto p0          = coupler.get_p0     ();
    auto gamma       = coupler.get_gamma_d();
    auto grav        = coupler.get_grav   ();
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
    YAKL_SCOPE( vert_sten_to_gll        , recon.vert_sten_to_gll       );
    YAKL_SCOPE( vert_sten_to_coefs      , recon.vert_sten_to_coefs     );
    YAKL_SCOPE( vert_weno_recon_lower   , recon.vert_weno_recon_lower  );
    YAKL_SCOPE( hyDensSten              , hydrostasis.hyDensSten       );
    YAKL_SCOPE( hyDensThetaSten         , hydrostasis.hyDensThetaSten  );
    YAKL_SCOPE( hyDensGLL               , hydrostasis.hyDensGLL        );
    YAKL_SCOPE( hyDensThetaGLL          , hydrostasis.hyDensThetaGLL   );

    real6d state_limits_x ("state_limits_x" ,num_state  ,2,nz,ny,nx+1,nens);
    real6d tracer_limits_x("tracer_limits_x",num_tracers,2,nz,ny,nx+1,nens);
    real6d state_limits_y ("state_limits_y" ,num_state  ,2,nz,ny+1,nx,nens);
    real6d tracer_limits_y("tracer_limits_y",num_tracers,2,nz,ny+1,nx,nens);
    real6d state_limits_z ("state_limits_z" ,num_state  ,2,nz+1,ny,nx,nens);
    real6d tracer_limits_z("tracer_limits_z",num_tracers,2,nz+1,ny,nx,nens);

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( "Spatial.h X recon" , SimpleBounds<4>(nz,ny,nx,nens) ,
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

      SArray<real,1,ord>  stencil;
      SArray<real,1,ngll> gll;

      //////////////
      // Density
      //////////////
      // x direction
      for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idR,hs+k,hs+j,wrapx(i,ii,nx),iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
      state_limits_x(idR,1,k,j,i  ,iens) = gll(0     );
      state_limits_x(idR,0,k,j,i+1,iens) = gll(ngll-1);
      // y direction
      if (! sim2d) {
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idR,hs+k,wrapy(j,jj,ny),hs+i,iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
        state_limits_y(idR,1,k,j  ,i,iens) = gll(0     );
        state_limits_y(idR,0,k,j+1,i,iens) = gll(ngll-1);
      }
      // z direction
      for (int kk=0; kk < ord; kk++) {
        if (k+kk < hs || k+kk > hs+nz-1) { stencil(kk) = state(idR,hs+k,hs+j,hs+i,iens) - hyDensSten(k,hs,iens); }
        else                             { stencil(kk) = state(idR,k+kk,hs+j,hs+i,iens) - hyDensSten(k,kk,iens); }
      }
      reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc , idl , sigma , weno_scalars );
      state_limits_z(idR,1,k  ,j,i,iens) = gll(0     ) + hyDensGLL(k,0     ,iens);
      state_limits_z(idR,0,k+1,j,i,iens) = gll(ngll-1) + hyDensGLL(k,ngll-1,iens);

      //////////////
      // u
      //////////////
      // x direction
      for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idU,hs+k,hs+j,wrapx(i,ii,nx),iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
      state_limits_x(idU,1,k,j,i  ,iens) = gll(0     );
      state_limits_x(idU,0,k,j,i+1,iens) = gll(ngll-1);
      // y direction
      if (! sim2d) {
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idU,hs+k,wrapy(j,jj,ny),hs+i,iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
        state_limits_y(idU,1,k,j  ,i,iens) = gll(0     );
        state_limits_y(idU,0,k,j+1,i,iens) = gll(ngll-1);
      }
      // z direction
      for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idU,wrapz(k,kk,nz),hs+j,hs+i,iens) /
                                                     state(idR,wrapz(k,kk,nz),hs+j,hs+i,iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc , idl , sigma , weno_winds );
      state_limits_z(idU,1,k  ,j,i,iens) = gll(0     ) * state_limits_z(idR,1,k  ,j,i,iens);
      state_limits_z(idU,0,k+1,j,i,iens) = gll(ngll-1) * state_limits_z(idR,0,k+1,j,i,iens);

      //////////////
      // v
      //////////////
      if (! sim2d) {
        // x direction
        for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idV,hs+k,hs+j,wrapx(i,ii,nx),iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
        state_limits_x(idV,1,k,j,i  ,iens) = gll(0     );
        state_limits_x(idV,0,k,j,i+1,iens) = gll(ngll-1);
        // y direction
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idV,hs+k,wrapy(j,jj,ny),hs+i,iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
        state_limits_y(idV,1,k,j  ,i,iens) = gll(0     );
        state_limits_y(idV,0,k,j+1,i,iens) = gll(ngll-1);
        // z direction
        for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idV,wrapz(k,kk,nz),hs+j,hs+i,iens) /
                                                       state(idR,wrapz(k,kk,nz),hs+j,hs+i,iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc , idl , sigma , weno_winds );
        state_limits_z(idV,1,k  ,j,i,iens) = gll(0     ) * state_limits_z(idR,1,k  ,j,i,iens);
        state_limits_z(idV,0,k+1,j,i,iens) = gll(ngll-1) * state_limits_z(idR,0,k+1,j,i,iens);
      }

      //////////////
      // w
      //////////////
      // x direction
      for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idW,hs+k,hs+j,wrapx(i,ii,nx),iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
      state_limits_x(idW,1,k,j,i  ,iens) = gll(0     );
      state_limits_x(idW,0,k,j,i+1,iens) = gll(ngll-1);
      // y direction
      if (! sim2d) {
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idW,hs+k,wrapy(j,jj,ny),hs+i,iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
        state_limits_y(idW,1,k,j  ,i,iens) = gll(0     );
        state_limits_y(idW,0,k,j+1,i,iens) = gll(ngll-1);
      }
      // z direction
      for (int kk=0; kk < ord; kk++) {
        stencil(kk) = state(idW,wrapz(k,kk,nz),hs+j,hs+i,iens) /
                      state(idR,wrapz(k,kk,nz),hs+j,hs+i,iens);
        if (k+kk > hs+nz-1 || k+kk < hs) stencil(kk) = 0;
      }
      reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc , idl , sigma , weno_winds );
      state_limits_z(idW,1,k  ,j,i,iens) = gll(0     ) * state_limits_z(idR,1,k  ,j,i,iens);
      state_limits_z(idW,0,k+1,j,i,iens) = gll(ngll-1) * state_limits_z(idR,0,k+1,j,i,iens);
      if (k == 0   ) state_limits_z(idW,1,k  ,j,i,iens) = 0;
      if (k == nz-1) state_limits_z(idW,0,k+1,j,i,iens) = 0;

      //////////////
      // theta
      //////////////
      // x direction
      for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idT,hs+k,hs+j,wrapx(i,ii,nx),iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
      state_limits_x(idT,1,k,j,i  ,iens) = gll(0     );
      state_limits_x(idT,0,k,j,i+1,iens) = gll(ngll-1);
      // y direction
      if (! sim2d) {
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idT,hs+k,wrapy(j,jj,ny),hs+i,iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
        state_limits_y(idT,1,k,j  ,i,iens) = gll(0     );
        state_limits_y(idT,0,k,j+1,i,iens) = gll(ngll-1);
      }
      // z direction
      for (int kk=0; kk < ord; kk++) {
        if (k+kk < hs || k+kk > hs+nz-1) { stencil(kk) = state(idT,hs+k,hs+j,hs+i,iens) - hyDensThetaSten(k,hs,iens); }
        else                             { stencil(kk) = state(idT,k+kk,hs+j,hs+i,iens) - hyDensThetaSten(k,kk,iens); }
      }
      reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc , idl , sigma , weno_scalars );
      state_limits_z(idT,1,k  ,j,i,iens) = gll(0     ) + hyDensThetaGLL(k,0     ,iens);
      state_limits_z(idT,0,k+1,j,i,iens) = gll(ngll-1) + hyDensThetaGLL(k,ngll-1,iens);

      for (int tr=0; tr < num_tracers; tr++) {
        // x direction
        for (int ii=0; ii < ord; ii++) { stencil(ii) = tracers(tr,hs+k,hs+j,wrapx(i,ii,nx),iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
        if (tracer_pos(tr)) { for (int ii=0; ii < ngll; ii++) { gll(ii) = max( 0._fp , gll(ii) ); } }
        tracer_limits_x(tr,1,k,j,i  ,iens) = gll(0     ); // Left interface
        tracer_limits_x(tr,0,k,j,i+1,iens) = gll(ngll-1); // Right interface
        // y direction
        if (! sim2d) {
          for (int jj=0; jj < ord; jj++) { stencil(jj) = tracers(tr,hs+k,wrapy(j,jj,ny),hs+i,iens); }
          reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
          if (tracer_pos(tr)) { for (int jj=0; jj < ngll; jj++) { gll(jj) = max( 0._fp , gll(jj) ); } }
          tracer_limits_y(tr,1,k,j  ,i,iens) = gll(0     ); // Left interface
          tracer_limits_y(tr,0,k,j+1,i,iens) = gll(ngll-1); // Right interface
        }
        // z direction
        for (int kk=0; kk < ord; kk++) { stencil(kk) = tracers(tr ,wrapz(k,kk,nz),hs+j,hs+i,iens) /
                                                       state  (idR,wrapz(k,kk,nz),hs+j,hs+i,iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc , idl , sigma , weno_scalars );
        if (tracer_pos(tr)) { for (int kk=0; kk < ngll; kk++) { gll(kk) = max( 0._fp , gll(kk) ); } }
        tracer_limits_z(tr,1,k  ,j,i,iens) = gll(0     ) * state_limits_z(idR,1,k  ,j,i,iens);
        tracer_limits_z(tr,0,k+1,j,i,iens) = gll(ngll-1) * state_limits_z(idR,0,k+1,j,i,iens);
      }
    });

    real5d state_flux_x ("state_flux_x" ,num_state  ,nz,ny,nx+1,nens);
    real5d tracer_flux_x("tracer_flux_x",num_tracers,nz,ny,nx+1,nens);
    real5d state_flux_y ("state_flux_y" ,num_state  ,nz,ny+1,nx,nens);
    real5d tracer_flux_y("tracer_flux_y",num_tracers,nz,ny+1,nx,nens);
    real5d state_flux_z ("state_flux_z" ,num_state  ,nz+1,ny,nx,nens);
    real5d tracer_flux_z("tracer_flux_z",num_tracers,nz+1,ny,nx,nens);

    awfl::riemann_rho_theta_full_xyz( coupler , state_limits_x , state_flux_x , tracer_limits_x , tracer_flux_x ,
                                                state_limits_y , state_flux_y , tracer_limits_y , tracer_flux_y ,
                                                state_limits_z , state_flux_z , tracer_limits_z , tracer_flux_z );

    state_limits_x  = real6d();
    tracer_limits_x = real6d();
    state_limits_y  = real6d();
    tracer_limits_y = real6d();
    state_limits_z  = real6d();
    tracer_limits_z = real6d();

    awfl::fct_positivity_xyz( coupler , tracers_start , tracer_flux_x , tracer_flux_y , tracer_flux_z , tracer_pos , dt );

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( "Spatial.h X tendencies" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l = 0; l < num_state; l++) {
        if (sim2d && l == idV) {
          state_tend(l,k,j,i,iens) = 0;
        } else {
          state_tend(l,k,j,i,iens)  = - ( state_flux_x(l,k,j,i+1,iens) - state_flux_x(l,k,j,i,iens) ) / dx;
          if (! sim2d) state_tend(l,k,j,i,iens) += - ( state_flux_y(l,k,j+1,i,iens) - state_flux_y(l,k,j,i,iens) ) / dy;
          state_tend(l,k,j,i,iens) += - ( state_flux_z(l,k+1,j,i,iens) - state_flux_z(l,k,j,i,iens) ) / dz(k,iens);
        }
      }
      state_tend(idW,k,j,i,iens) -= state(idR,hs+k,hs+j,hs+i,iens) * grav;
      for (int l = 0; l < num_tracers; l++) {
        // Compute tracer tendency
        tracer_tend(l,k,j,i,iens)  = - ( tracer_flux_x(l,k,j,i+1,iens) - tracer_flux_x(l,k,j,i,iens) ) / dx;
        if (! sim2d) tracer_tend(l,k,j,i,iens) += - ( tracer_flux_y(l,k,j+1,i,iens) - tracer_flux_y(l,k,j,i,iens) ) / dy;
        tracer_tend(l,k,j,i,iens) += - ( tracer_flux_z(l,k+1,j,i,iens) - tracer_flux_z(l,k,j,i,iens) ) / dz(k,iens);
      }
    });
  }


}
}
}


