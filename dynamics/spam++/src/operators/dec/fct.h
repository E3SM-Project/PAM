#ifndef _FCT_H_
#define _FCT_H_

#include "common.h"

template<uint ndofs> void YAKL_INLINE calculate_edgeflux(SArray<real,ndofs,ndims> &edgeflux, SArray<real,ndofs,ndims> const &recon, SArray<real,ndims> const &flux)
{
  for (int l=0; l<ndofs; l++) {
  for (int k=0; k<ndims; k++) {
    edgeflux(l,k) = recon(l,k) * flux(k);
  }}
}

yakl::parallel_for("ComputeEdgeFlux", topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
SArray<real,nqdofs,ndims> edgeflux;
SArray<real,nqdofs,ndims> recon;
SArray<real,ndims> flux;
for (int d=0; d<ndims; d++) {
flux(d) = auxiliary_vars.fields_arr[UVAR].data(d, k+ks, j+js, i+is);
for (int l=0; l<nqdofs; l++) {
  recon(l,d) = auxiliary_vars.fields_arr[QRECONVAR].data(l+d*nqdofs, k+ks, j+js, i+is);
}}
calculate_edgeflux<nqdofs>(edgeflux, recon, flux);

for (int d=0; d<ndims; d++) {
for (int l=0; l<nqdofs; l++) {
auxiliary_vars.fields_arr[EDGEFLUXVAR].data(l+d*nqdofs, k+ks, j+js, i+is) = edgeflux(l,d);
}}
});





template<uint ndofs> void YAKL_INLINE calculate_Mf(SArray<real,ndofs> &Mf, SArray<real,ndofs,ndims,2> const &edgeflux, real dt)
{
    real eps = 1.0e-8;
        for (int l=0; l<ndofs; l++) {
          Mf(l) = 0.;
          for (int k=0; k<ndims; k++) {
            Mf(l) += mymax(edgeflux(l,k,1),0.0) - mymin(edgeflux(l,k,0), 0.0);
          }
            Mf(l) *= dt;
            Mf(l) += eps;
        }
}
yakl::parallel_for("ComputeMf", topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
SArray<real,nqdofs> Mf;
SArray<real,nqdofs,ndims,2> edgeflux;
for (int d=0; d<ndims; d++) {
    for (int l=0; l<nqdofs; l++) {
      for (int m=0; m<2; m++) {
        if (d==0) {edgeflux(l,d,m) = auxiliary_vars.fields_arr[EDGEFLUXVAR].data(l+d*nqdofs, k+ks, j+js, i+is+m); }
        if (d==1) {edgeflux(l,d,m) = auxiliary_vars.fields_arr[EDGEFLUXVAR].data(l+d*nqdofs, k+ks, j+js+m, i+is); }
        if (d==2) {edgeflux(l,d,m) = auxiliary_vars.fields_arr[EDGEFLUXVAR].data(l+d*nqdofs, k+ks+m, j+js, i+is); }
  }}}
    calculate_Mf<nqdofs>(Mf, edgeflux, dt);
    for (int l=0; l<nqdofs; l++) { auxiliary_vars.fields_arr[MFVAR].data(l, k+ks, j+js, i+is) = Mf(l); }
    });



template<uint ndofs> void YAKL_INLINE calculate_phi(SArray<real,ndofs,ndims> &Phi, SArray<real,ndofs,ndims,2> const &q, SArray<real,ndofs,ndims,2> const &Mf, SArray<real,ndofs,ndims> const &edgeflux)
{

    real upwind_param;

        for (int l=0; l<ndofs; l++) {
          for (int k=0; k<ndims; k++) {
            upwind_param = copysign(1.0, edgeflux(l,k));
            upwind_param = 0.5*(upwind_param + fabs(upwind_param));
            Phi(l,k) = mymin(1., q(l,k,1)/Mf(l,k,1)) * (1. - upwind_param) + mymin(1., q(l,k,0)/Mf(l,k,0)) * upwind_param;
          }
    }
}


yakl::parallel_for("ComputePhi", topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
SArray<real,nqdofs,ndims> Phi;
SArray<real,nqdofs,ndims,2> const q;
SArray<real,nqdofs,ndims,2> const Mf;
SArray<real,nqdofs,ndims> edgeflux;
for (int d=0; d<ndims; d++) {
    for (int l=0; l<nqdofs; l++) {
      edgeflux(l,d) = auxiliary_vars.fields_arr[EDGEFLUXVAR].data(l+d*nqdofs, k+ks, j+js, i+is);
      for (int m=0; m<2; m++) {
        if (d==0) {
          Mf(l,d,m) = auxiliary_vars.fields_arr[MFVAR].data(l, k+ks, j+js, i+is+m-1);
          q(l,d,m) = x.fields_arr[QVAR].data(l, k+ks, j+js, i+is+m-1);
        }
        if (d==1) {
          Mf(l,d,m) = auxiliary_vars.fields_arr[MFVAR].data(l, k+ks, j+js+m-1, i+is);
          q(l,d,m) = x.fields_arr[QVAR].data(l, k+ks, j+js+m-1, i+is);
        }
        if (d==2) {
          Mf(l,d,m) = auxiliary_vars.fields_arr[MFVAR].data(l, k+ks+m-1, j+js, i+is);
          q(l,d,m) = x.fields_arr[QVAR].data(l, k+ks+m-1, j+js, i+is);
        }
  }}}
calculate_phi<nqdofs>(Phi, q, Mf , edgeflux);
for (int d=0; d<ndims; d++) {
    for (int l=0; l<nqdofs; l++) {
  auxiliary_vars.fields_arr[PHIVAR].data(l+d*nqdofs, k+ks, j+js, i+is) = Phi(l,d);
    }}
});

#endif
