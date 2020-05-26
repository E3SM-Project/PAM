#ifndef _FCT_H_
#define _FCT_H_

#include "common.h"



template<uint ndofs> void YAKL_INLINE calculate_edgeflux(SArray<real,ndofs,ndims> &edgeflux, SArray<real,ndofs,ndims> const &recon, SArray<real,ndims> const &flux)
{
  for (int l=0; l<ndofs; l++) {
  for (int d=0; d<ndims; d++) {
    edgeflux(l,d) = recon(l,d) * flux(d);
  }}
}

template<uint ndofs> YAKL_INLINE void compute_edgefluxes(realArr edgefluxvar, const realArr reconvar, const realArr U, int is, int js, int ks, int i, int j, int k)
{
SArray<real,ndofs,ndims> edgeflux;
SArray<real,ndofs,ndims> recon;
SArray<real,ndims> flux;
for (int d=0; d<ndims; d++) {
flux(d) = U(d, k+ks, j+js, i+is);
for (int l=0; l<ndofs; l++) {
  recon(l,d) = reconvar(l+d*ndofs, k+ks, j+js, i+is);
}}
calculate_edgeflux<ndofs>(edgeflux, recon, flux);

for (int d=0; d<ndims; d++) {
for (int l=0; l<ndofs; l++) {
edgefluxvar(l+d*ndofs, k+ks, j+js, i+is) = edgeflux(l,d);
}}
}





template<uint ndofs> void YAKL_INLINE calculate_Mf(SArray<real,ndofs> &Mf, SArray<real,ndofs,ndims,2> const &edgeflux, real dt)
{
    real eps = 1.0e-8;
        for (int l=0; l<ndofs; l++) {
          Mf(l) = 0.;
          for (int d=0; d<ndims; d++) {
            Mf(l) += mymax(edgeflux(l,d,1),0.0) - mymin(edgeflux(l,d,0), 0.0);
          }
            Mf(l) *= dt;
            Mf(l) += eps;
        }
}

template<uint ndofs> YAKL_INLINE void compute_Mf(realArr Mfvar, const realArr edgefluxvar, real dt, int is, int js, int ks, int i, int j, int k)
{
SArray<real,ndofs> Mf;
SArray<real,ndofs,ndims,2> edgeflux;
for (int d=0; d<ndims; d++) {
    for (int l=0; l<ndofs; l++) {
      for (int m=0; m<2; m++) {
        if (d==0) {edgeflux(l,d,m) = edgefluxvar(l+d*ndofs, k+ks, j+js, i+is+m); }
        if (d==1) {edgeflux(l,d,m) = edgefluxvar(l+d*ndofs, k+ks, j+js+m, i+is); }
        if (d==2) {edgeflux(l,d,m) = edgefluxvar(l+d*ndofs, k+ks+m, j+js, i+is); }
  }}}
    calculate_Mf<ndofs>(Mf, edgeflux, dt);
    for (int l=0; l<ndofs; l++) { Mfvar(l, k+ks, j+js, i+is) = Mf(l); }
}


template<uint ndofs> void YAKL_INLINE calculate_phi(SArray<real,ndofs,ndims> &Phi, SArray<real,ndofs,ndims,2> const &q, SArray<real,ndofs,ndims,2> const &Mf, SArray<real,ndofs,ndims> const &edgeflux)
{

    real upwind_param;

        for (int l=0; l<ndofs; l++) {
          for (int d=0; d<ndims; d++) {
            upwind_param = copysign(1.0, edgeflux(l,d));
            upwind_param = 0.5*(upwind_param + fabs(upwind_param));
            Phi(l,d) = mymin(1., q(l,d,1)/Mf(l,d,1)) * (1. - upwind_param) + mymin(1., q(l,d,0)/Mf(l,d,0)) * upwind_param;
          }
    }
}

template<uint ndofs> YAKL_INLINE void compute_Phi(realArr Phivar, const realArr edgefluxvar, const realArr Mfvar, const realArr qvar, int is, int js, int ks, int i, int j, int k)
{
SArray<real,ndofs,ndims> Phi;
SArray<real,ndofs,ndims,2> const q;
SArray<real,ndofs,ndims,2> const Mf;
SArray<real,ndofs,ndims> edgeflux;
for (int d=0; d<ndims; d++) {
    for (int l=0; l<ndofs; l++) {
      edgeflux(l,d) = edgefluxvar(l+d*ndofs, k+ks, j+js, i+is);
      for (int m=0; m<2; m++) {
        if (d==0) {
          Mf(l,d,m) = Mfvar(l, k+ks, j+js, i+is+m-1);
          q(l,d,m) = qvar(l, k+ks, j+js, i+is+m-1);
        }
        if (d==1) {
          Mf(l,d,m) = Mfvar(l, k+ks, j+js+m-1, i+is);
          q(l,d,m) = qvar(l, k+ks, j+js+m-1, i+is);
        }
        if (d==2) {
          Mf(l,d,m) = Mfvar(l, k+ks+m-1, j+js, i+is);
          q(l,d,m) = qvar(l, k+ks+m-1, j+js, i+is);
        }
  }}}
calculate_phi<ndofs>(Phi, q, Mf , edgeflux);
for (int d=0; d<ndims; d++) {
    for (int l=0; l<ndofs; l++) {
  Phivar(l+d*ndofs, k+ks, j+js, i+is) = Phi(l,d);
    }}
}
#endif
