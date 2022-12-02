
#pragma once

#include "pam_coupler.h"

namespace modules {


  void horizontal_average( pam::PamCoupler &coupler , std::vector<std::string> var_names ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    
    int nz   = coupler.get_nz  ();
    int ny   = coupler.get_ny  ();
    int nx   = coupler.get_nx  ();
    int nens = coupler.get_nens();

    auto &dm = coupler.get_data_manager_readwrite();

    MultiField<real,3> fields_var;
    MultiField<real,2> fields_havg;

    int num_vars = var_names.size();
    int max_dim = 0;

    for (int i=0; i < num_vars; i++) {
      auto var_name = var_names[i];
      auto havg_name = var_name + std::string("_horizontal_average");
      auto shape = dm.get_shape(var_name);
      int k, ncol, nens;
      if (shape.size() == 1) endrun("ERROR: Cannot horizontally average a 1-D variable");
      if (shape.size() == 2) { k = 1;   ncol = shape[0];   nens = shape[1]; }
      if (shape.size() == 3) { k = shape[0];   ncol = shape[1];   nens = shape[2]; }
      if (shape.size() == 4) { k = shape[0];   ncol = shape[1]*shape[2];   nens = shape[3]; }
      dm.register_and_allocate<real>( havg_name , "" , {k,nens} );
      fields_var .add_field( dm.get_collapsed<real>(var_name).reshape(k,ncol,nens) );
      fields_havg.add_field( dm.get<real,2>( havg_name ) );
    }

    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(num_vars,nz,nens) , YAKL_LAMBDA (int l, int k, int iens) {
      int ncol = fields_var(l).extent(1);
      fields_havg(l,k,iens) = 0;
      r_ncol = 1._fp / ncol;
      for (int i=0; i < ncol; i++) { fields_havg(l,k,iens) += fields_var(l,k,i,iens) * r_ncol; }
    });
  }


}


