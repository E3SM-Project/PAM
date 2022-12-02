
#pragma once

#include "pam_coupler.h"

namespace modules {

  void time_average_init( pam::PamCoupler &coupler , std::vector<std::string> var_names ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    
    int nz   = coupler.get_nz  ();
    int ny   = coupler.get_ny  ();
    int nx   = coupler.get_nx  ();
    int nens = coupler.get_nens();

    auto &dm = coupler.get_data_manager_readwrite();

    MultiField<real,1> fields;

    int num_vars = var_names.size();
    int max_dim = 0;

    for (int i=0; i < num_vars; i++) {
      auto var_name = var_names[i];
      auto tavg_name = var_name + std::string("_time_average");
      dm.register_and_allocate<real>( tavg_name , "" , dm.get_shape(var_name) );
      auto tavg_var = dm.get_collapsed<real>( tavg_name );
      fields.add_field( tavg_var );
      max_dim = max( tavg_var.size() , max_dim );
    }

    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(num_vars,max_dim) , YAKL_LAMBDA (int l, int i) {
      if (i < fields(l).size()) fields(l,i) = 0;
    });
  }


  void time_average_accumulate( pam::PamCoupler &coupler , std::vector<std::string> var_names ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    
    int nz   = coupler.get_nz  ();
    int ny   = coupler.get_ny  ();
    int nx   = coupler.get_nx  ();
    int nens = coupler.get_nens();

    auto &dm = coupler.get_data_manager_readwrite();

    auto crm_dt = coupler.get_option<real>("crm_dt");
    auto gcm_dt = coupler.get_option<real>("gcm_physics_dt");

    MultiField<real,1> fields_var;
    MultiField<real,1> fields_tavg;

    int num_vars = var_names.size();
    int max_dim = 0;

    for (int i=0; i < num_vars; i++) {
      auto var_name = var_names[i];
      auto tavg_name = var_name + std::string("_time_average");
      fields_var .add_field( dm.get_collapsed<real>( var_name  ) );
      fields_tavg.add_field( dm.get_collapsed<real>( tavg_name ) );
      max_dim = max( fields_var(i).size() , max_dim );
    }

    real factor = crm_dt / gcm_dt
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(num_vars,max_dim) , YAKL_LAMBDA (int l, int i) {
      if (i < fields(l).size()) fields_tavg(l,i) += fields_var(l,i) * factor;
    });
  }

}


