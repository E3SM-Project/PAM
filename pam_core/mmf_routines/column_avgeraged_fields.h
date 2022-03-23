
#pragma once

#include "pam_const.h"
#include "pam_coupler.h"


// Templated on the number of dimensions, N, in each field and the type, T, of the fields being column-averaged
// It is expected that nz is the first dimension, and nens is the last dimension for each field
// It is expected that the fields all have the same dimension sizes.
// The returned array always has 3 dimensions: { names.size() , vertical_size , nens }
// N must be > 2 for this routine to be meaningful: nz, nens, and the "column" dimensions
template <class T, int N>
inline real3d column_averaged_fields( PamCoupler &coupler , std::vector<std::string> names ) {
  auto tmp = coupler.dm.get<T,N>( names[0] );

  int num_fields = names.size();
  int nz = tmp.dimensions[0];
  int ncol = 1;
  for (int i=1; i < N-1; i++) { ncol *= tmp.dimension[i]; } // dims {1,...,N-2} are column dimensions
  int nens = coupler.dm.get_dimension_size("nens");

  // Aggregate the fields requested by the caller
  // Resize each to level,column,ensemble dimensions
  MultiField<T,3> fields;
  for (int i = 0; i < num_fields; i++) {
    auto fld = coupler.dm.get<T,N>( names[i] ).resize<3>({nz,ncol,nens});
    fields.add_field( fld );
  }

  real3d column_averaged("column_averaged",num_fields,nz,nens);
  memset( column_averaged , 0._fp );

  using yakl::atomicAdd;

  real r_ncol = 1._fp / ncol;
  parallel_for( SimpleBounds<4>(num_fields,nz,ncol,nens) , YAKL_LAMBDA (int l, int k, int i, int iens) {
    atomicAdd( column_averaged(l,k,iens) , fields(l,k,i,iens) * r_ncol );
  });

  return column_averaged;
}


