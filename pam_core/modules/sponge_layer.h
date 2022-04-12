
#pragma once

#include "pam_const.h"
#include "pam_coupler.h"

namespace modules {

  inline void sponge_layer( PamCoupler &coupler , real dt ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    int nz   = coupler.get_nz  ();
    int ny   = coupler.get_ny  ();
    int nx   = coupler.get_nx  ();
    int nens = coupler.get_nens();

    int num_layers = 10;  // Number of model top vertical layers that participate in the sponge relaxation

    int WFLD = 3; // fourth entry into "fields" is the "w velocity" field. Set the havg to zero for WFLD

    // Get a list of tracer names for retrieval
    std::vector<std::string> tracer_names = coupler.get_tracer_names();
    int num_tracers = coupler.get_num_tracers();

    // Allocate MultiField of all state and tracer havg variables, since we're doing the same operation on each
    pam::MultiField<real,2> havg_fields;
    havg_fields.add_field( real2d("havg_rho_d",nz,nens) );
    havg_fields.add_field( real2d("havg_uvel" ,nz,nens) );
    havg_fields.add_field( real2d("havg_vvel" ,nz,nens) );
    havg_fields.add_field( real2d("havg_wvel" ,nz,nens) );
    havg_fields.add_field( real2d("havg_temp" ,nz,nens) );
    for (int tr=0; tr < num_tracers; tr++) {
      char const * name = (std::string("havg_")+tracer_names[tr]).c_str();
      havg_fields.add_field( real2d(name,nz,nens) );
    }

    auto &dm = coupler.get_data_manager_readwrite();

    // Create MultiField of all state and tracer full variables, since we're doing the same operation on each
    pam::MultiField<real,4> full_fields;
    full_fields.add_field( dm.get<real,4>("density_dry") );
    full_fields.add_field( dm.get<real,4>("uvel"       ) );
    full_fields.add_field( dm.get<real,4>("vvel"       ) );
    full_fields.add_field( dm.get<real,4>("wvel"       ) );
    full_fields.add_field( dm.get<real,4>("temp"       ) );
    for (int tr=0; tr < num_tracers; tr++) {
      full_fields.add_field( dm.get<real,4>(tracer_names[tr]) );
    }

    int num_fields = havg_fields.get_num_fields();

    real r_nx_ny = 1._fp / (nx*ny);

    // Zero out the havg_fields
    parallel_for( Bounds<3>(num_fields,nz,nens) , YAKL_LAMBDA (int ifld, int k, int iens) {
      havg_fields(ifld,k,iens) = 0;
    });

    // Compute the horizontal average for each vertical level (that we use for the sponge layer) and ensemble
    parallel_for( Bounds<5>(num_fields,num_layers,ny,nx,nens) , YAKL_LAMBDA (int ifld, int kloc, int j, int i, int iens) {
      int k = nz - 1 - kloc;
      if (ifld != WFLD) yakl::atomicAdd( havg_fields(ifld,k,iens) , full_fields(ifld,k,j,i,iens) * r_nx_ny );
    });

    auto zint = dm.get<real const,2>("vertical_interface_height");
    auto zmid = dm.get<real const,2>("vertical_midpoint_height" );

    real constexpr time_scale = 60;  // strength of each application is dt / time_scale  (same as SAM's tau_min)
    real time_factor = dt / time_scale;

    // use a cosine relaxation in space:  ((cos(pi*rel_dist)+1)/2)^2
    parallel_for( Bounds<5>(num_fields,num_layers,ny,nx,nens) , YAKL_LAMBDA (int ifld, int kloc, int j, int i, int iens) {
      int k = nz - 1 - kloc;
      real rel_dist = ( zint(nz,iens) - zmid(k,iens) ) / ( zint(nz,iens) - zmid(nz-1-(num_layers-1),iens) );
      real space_factor = ( cos(M_PI*rel_dist) + 1 ) / 2;
      real factor = space_factor * time_factor;
      full_fields(ifld,k,j,i,iens) += ( havg_fields(ifld,k,iens) - full_fields(ifld,k,j,i,iens) ) * factor;
    });
  }

}


