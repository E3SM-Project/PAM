
#include "fct_positivity.h"

namespace awfl {


  void fct_positivity_x( pam::PamCoupler const &coupler, realConst5d tracers, real5d const &tracer_flux ,
                         boolConst1d tracer_pos, real dt) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    auto num_tracers = tracers.extent(0);
    auto nz          = coupler.get_nz  ();
    auto ny          = coupler.get_ny  ();
    auto nx          = coupler.get_nx  ();
    auto nens        = coupler.get_nens();
    auto dx          = coupler.get_dx  ();
    parallel_for( "Spatial.h X FCT" , SimpleBounds<5>(num_tracers,nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int tr, int k, int j, int i, int iens) {
      if (tracer_pos(tr)) {
        real my_mass = tracers(tr,hs+k,hs+j,hs+i,iens);
        real mass_out = dt*(std::max(0._fp,tracer_flux(tr,k,j,i+1,iens)) -
                            std::min(0._fp,tracer_flux(tr,k,j,i  ,iens)))/dx;
        if (mass_out > my_mass) {
          // Limit the fluxes involved in removing mass from this cell
          // Make sure periodic boundary fluxes match when changed
          if (tracer_flux(tr,k,j,i+1,iens) > 0) {
            tracer_flux(tr,k,j,i+1,iens) *= my_mass / mass_out;
            if (i+1 == nx) tracer_flux(tr,k,j,0 ,iens) = tracer_flux(tr,k,j,i+1,iens);
          }
          if (tracer_flux(tr,k,j,i  ,iens) < 0) {
            tracer_flux(tr,k,j,i  ,iens) *= my_mass / mass_out;
            if (i   == 0 ) tracer_flux(tr,k,j,nx,iens) = tracer_flux(tr,k,j,i  ,iens);
          }
        }
      }
    });
  }



  void fct_positivity_y( pam::PamCoupler const &coupler, realConst5d tracers, real5d const &tracer_flux ,
                         boolConst1d tracer_pos, real dt) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    auto num_tracers = tracers.extent(0);
    auto nz          = coupler.get_nz  ();
    auto ny          = coupler.get_ny  ();
    auto nx          = coupler.get_nx  ();
    auto nens        = coupler.get_nens();
    auto dy          = coupler.get_dy  ();
    parallel_for( "Spatial.h Y FCT" , SimpleBounds<5>(num_tracers,nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int tr, int k, int j, int i, int iens) {
      if (tracer_pos(tr)) {
        real my_mass = tracers(tr,hs+k,hs+j,hs+i,iens);
        real mass_out = dt*(std::max(0._fp,tracer_flux(tr,k,j+1,i,iens)) -
                            std::min(0._fp,tracer_flux(tr,k,j  ,i,iens)))/dy;
        if (mass_out > my_mass) {
          // Limit the fluxes involved in removing mass from this cell
          // Make sure periodic boundary fluxes match when changed
          if (tracer_flux(tr,k,j+1,i,iens) > 0) {
            tracer_flux(tr,k,j+1,i,iens) *= my_mass / mass_out;
            if (j+1 == ny) tracer_flux(tr,k,0 ,i,iens) = tracer_flux(tr,k,j+1,i,iens);
          }
          if (tracer_flux(tr,k,j  ,i,iens) < 0) {
            tracer_flux(tr,k,j  ,i,iens) *= my_mass / mass_out;
            if (j   == 0 ) tracer_flux(tr,k,ny,i,iens) = tracer_flux(tr,k,j  ,i,iens);
          }
        }
      }
    });
  }



  void fct_positivity_z( pam::PamCoupler const &coupler, realConst5d tracers, real5d const &tracer_flux ,
                         boolConst1d tracer_pos, real dt) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    auto num_tracers = tracers.extent(0);
    auto nz          = coupler.get_nz  ();
    auto ny          = coupler.get_ny  ();
    auto nx          = coupler.get_nx  ();
    auto nens        = coupler.get_nens();
    auto dz          = coupler.get_data_manager_readonly().get<real const,2>("vertical_cell_dz");
    parallel_for( "Spatial.h Z FCT" , SimpleBounds<5>(num_tracers,nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int tr, int k, int j, int i, int iens) {
      if (tracer_pos(tr)) {
        real my_mass = tracers(tr,hs+k,hs+j,hs+i,iens);
        real mass_out = dt*(std::max(0._fp,tracer_flux(tr,k+1,j,i,iens)) -
                            std::min(0._fp,tracer_flux(tr,k  ,j,i,iens)))/dz(k,iens);
        if (mass_out > my_mass) {
          // Limit the fluxes involved in removing mass from this cell
          if (tracer_flux(tr,k+1,j,i,iens) > 0) tracer_flux(tr,k+1,j,i,iens) *= my_mass / mass_out;
          if (tracer_flux(tr,k  ,j,i,iens) < 0) tracer_flux(tr,k  ,j,i,iens) *= my_mass / mass_out;
        }
      }
    });
  }



  void fct_positivity_xyz( pam::PamCoupler const &coupler ,
                           realConst5d tracers            ,
                           real5d const &tracer_flux_x    ,
                           real5d const &tracer_flux_y    ,
                           real5d const &tracer_flux_z    ,
                           boolConst1d tracer_pos, real dt) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    auto num_tracers = tracers.extent(0);
    auto nz          = coupler.get_nz  ();
    auto ny          = coupler.get_ny  ();
    auto nx          = coupler.get_nx  ();
    auto nens        = coupler.get_nens();
    auto dx          = coupler.get_dx  ();
    auto dy          = coupler.get_dy  ();
    auto dz          = coupler.get_data_manager_readonly().get<real const,2>("vertical_cell_dz");
    parallel_for( "Spatial.h Z FCT" , SimpleBounds<5>(num_tracers,nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int tr, int k, int j, int i, int iens) {
      if (tracer_pos(tr)) {
        real my_mass = tracers(tr,hs+k,hs+j,hs+i,iens);
        // x-direction
        real mass_out  = dt*(std::max(0._fp,tracer_flux_x(tr,k,j,i+1,iens)) -
                             std::min(0._fp,tracer_flux_x(tr,k,j,i  ,iens)))/dx;
        // y-direction
        mass_out += dt*(std::max(0._fp,tracer_flux_y(tr,k,j+1,i,iens)) -
                        std::min(0._fp,tracer_flux_y(tr,k,j  ,i,iens)))/dy;
        // z-direction
        mass_out += dt*(std::max(0._fp,tracer_flux_z(tr,k+1,j,i,iens)) -
                        std::min(0._fp,tracer_flux_z(tr,k  ,j,i,iens)))/dz(k,iens);
        if (mass_out > my_mass) {
          // x-direction
          if (tracer_flux_x(tr,k,j,i+1,iens) > 0) {
            tracer_flux_x(tr,k,j,i+1,iens) *= my_mass / mass_out;
            if (i+1 == nx) tracer_flux_x(tr,k,j,0 ,iens) = tracer_flux_x(tr,k,j,i+1,iens);
          }
          if (tracer_flux_x(tr,k,j,i  ,iens) < 0) {
            tracer_flux_x(tr,k,j,i  ,iens) *= my_mass / mass_out;
            if (i   == 0 ) tracer_flux_x(tr,k,j,nx,iens) = tracer_flux_x(tr,k,j,i  ,iens);
          }
          // y-direction
          if (tracer_flux_y(tr,k,j+1,i,iens) > 0) {
            tracer_flux_y(tr,k,j+1,i,iens) *= my_mass / mass_out;
            if (j+1 == ny) tracer_flux_y(tr,k,0 ,i,iens) = tracer_flux_y(tr,k,j+1,i,iens);
          }
          if (tracer_flux_y(tr,k,j  ,i,iens) < 0) {
            tracer_flux_y(tr,k,j  ,i,iens) *= my_mass / mass_out;
            if (j   == 0 ) tracer_flux_y(tr,k,ny,i,iens) = tracer_flux_y(tr,k,j  ,i,iens);
          }
          // z-direction
          if (tracer_flux_z(tr,k+1,j,i,iens) > 0) tracer_flux_z(tr,k+1,j,i,iens) *= my_mass / mass_out;
          if (tracer_flux_z(tr,k  ,j,i,iens) < 0) tracer_flux_z(tr,k  ,j,i,iens) *= my_mass / mass_out;
        }
      }
    });
  }


}


