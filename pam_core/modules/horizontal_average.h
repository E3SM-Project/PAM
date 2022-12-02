
#pragma once

#include "pam_coupler.h"

namespace modules {


  // Computes the horizontal averaging over columns of a list of variable names and stores each in the coupler as
  // varname_horizontal_average. All variables are assumed to have an "nens" dimension at the end.
  // You'll input a list of {string,bool} pairs, where the bool tells us if the variable has a vertical dimension.
  // 
  // If a variable has *no* vertical dimension, it is interpreted as follows:
  // 1-D: Not allowed. This implies only a scalar for each CRM
  // 2-D: ncol,nens
  // 3-D: ny,nx,nens
  // 4-D: Not allowed. Only two horizontal dimensions makes sense
  // 
  // If a variable has a vertical dimension, it is interpreted as follows:
  // 1-D: Not allowed. This implies only a scalar for each CRM.
  // 2-D: Not allowed. This implies nz,nens, which does not have horizontal dimensions
  // 3-D: nz,ncol,nens
  // 4-D: nz,ny,nx,nens
  // 5-D: Not allowed. Only two horizontal dimensions makes sense
  void horizontal_average( pam::PamCoupler &coupler ,
                           std::vector<std::tuple<std::string,bool>> var_list ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    
    int nz   = coupler.get_nz  ();
    int nens = coupler.get_nens();

    auto &dm = coupler.get_data_manager_readwrite();

    MultiField<real,3> fields_var;
    MultiField<real,2> fields_havg;

    int num_vars = var_list.size();
    int max_nz = 1;

    for (int i=0; i < num_vars; i++) {
      auto var_name         = var_list[i].get<0>();
      auto has_vertical_dim = var_list[i].get<1>();
      auto havg_name = var_name + std::string("_horizontal_average");
      auto shape = dm.get_shape(var_name);
      int nz, ncol, nens;
      if (shape[shape.size()-1] != nens) endrun("ERROR: Last dimension must be nens");
      if (has_vertical_dim) {
        if (shape.size() == 1) endrun("ERROR: Cannot horizontally average a 1-D variable");
        if (shape.size() == 2) endrun("ERROR: Cannot horizontally average a nz,nens variable");
        if (shape.size() == 3) { nz = shape[0];   ncol = shape[1];   nens = shape[2]; }
        if (shape.size() == 4) { nz = shape[0];   ncol = shape[1]*shape[2];   nens = shape[3]; }
      } else {
        if (shape.size() == 1) endrun("ERROR: Cannot horizontally average a 1-D variable");
        if (shape.size() == 2) { nz = 1;   ncol = shape[0];   nens = shape[1]; }
        if (shape.size() == 3) { nz = 1;   ncol = shape[0]*shape[1];   nens = shape[2]; }
        if (shape.size() == 4) endrun("ERROR: Only two horizontal dimensions allowed");
      }
      dm.register_and_allocate<real>( havg_name , "" , {nz,nens} );
      fields_var .add_field( dm.get_collapsed<real>(var_name).reshape(nz,ncol,nens) );
      fields_havg.add_field( dm.get<real,2>( havg_name ) );
      max_nz = max(nz,max_nz);
    }

    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(num_vars,max_nz,nens) , YAKL_LAMBDA (int l, int k, int iens) {
      int nz   = fields_var(l).extent(0);
      int ncol = fields_var(l).extent(1);
      if (k < nz) {
        fields_havg(l,k,iens) = 0;
        r_ncol = 1._fp / ncol;
        for (int i=0; i < ncol; i++) { fields_havg(l,k,iens) += fields_var(l,k,i,iens) * r_ncol; }
      }
    });
  }


}


