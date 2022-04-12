#pragma once

#include "common.h"



template<uint ndofs> void YAKL_INLINE calculate_edgeflux(SArray<real,2,ndofs,ndims> &edgeflux, SArray<real,2,ndofs,ndims> const &recon, SArray<real,1,ndims> const &flux)
{
  for (int l=0; l<ndofs; l++) {
  for (int d=0; d<ndims; d++) {
    edgeflux(l,d) = recon(l,d) * flux(d);
  }}
}

template<uint ndofs> YAKL_INLINE void compute_edgefluxes(real5d edgefluxvar, const real5d reconvar, const real5d U, int is, int js, int ks, int i, int j, int k, int n)
{
SArray<real,2,ndofs,ndims> edgeflux;
SArray<real,2,ndofs,ndims> recon;
SArray<real,1,ndims> flux;
for (int d=0; d<ndims; d++) {
flux(d) = U(d, k+ks, j+js, i+is, n);
for (int l=0; l<ndofs; l++) {
  recon(l,d) = reconvar(l+d*ndofs, k+ks, j+js, i+is, n);
}}
calculate_edgeflux<ndofs>(edgeflux, recon, flux);

for (int d=0; d<ndims; d++) {
for (int l=0; l<ndofs; l++) {
edgefluxvar(l+d*ndofs, k+ks, j+js, i+is, n) = edgeflux(l,d);
}}
}

template<uint ndofs> void YAKL_INLINE calculate_vertedgeflux(SArray<real,1,ndofs> &edgeflux, SArray<real,1,ndofs> const &recon, real const &flux)
{
  for (int l=0; l<ndofs; l++) {
    edgeflux(l) = recon(l) * flux;
  }
}

template<uint ndofs> YAKL_INLINE void compute_vertedgefluxes(real5d vertedgefluxvar, const real5d vertreconvar, const real5d UW, int is, int js, int ks, int i, int j, int k, int n)
{
  SArray<real,1,ndofs> edgeflux;
  SArray<real,1,ndofs> recon;
  real flux;
  flux = UW(0, k+ks, j+js, i+is, n);
  for (int l=0; l<ndofs; l++) {
    recon(l) = vertreconvar(l, k+ks, j+js, i+is, n);
  }
  calculate_vertedgeflux<ndofs>(edgeflux, recon, flux);

  for (int l=0; l<ndofs; l++) {
  vertedgefluxvar(l, k+ks, j+js, i+is, n) = edgeflux(l);
  }
}


template<uint ndofs> void YAKL_INLINE calculate_Mf(SArray<real,1,ndofs> &Mf, SArray<real,3,ndofs,ndims,2> const &edgeflux, real dt)
{
    real eps = 1.0e-8_fp;
        for (int l=0; l<ndofs; l++) {
          Mf(l) = 0.;
          for (int d=0; d<ndims; d++) {
            Mf(l) += mymax(edgeflux(l,d,1),0.0) - mymin(edgeflux(l,d,0), 0.0);
          }
            Mf(l) *= dt;
            Mf(l) += eps;
        }
}

template<uint ndofs> void YAKL_INLINE calculate_Mfext(SArray<real,1,ndofs> &Mf, SArray<real,3,ndofs,ndims,2> const &edgeflux, SArray<real,2,ndofs,2> const &vertedgeflux, real dt)
{
    real eps = 1.0e-8_fp;
        for (int l=0; l<ndofs; l++) {
          Mf(l) = 0.;
          for (int d=0; d<ndims; d++) {
            Mf(l) += mymax(edgeflux(l,d,1),0.0) - mymin(edgeflux(l,d,0), 0.0);
          }
          Mf(l) += mymax(vertedgeflux(l,1),0.0) - mymin(vertedgeflux(l,0), 0.0);
            Mf(l) *= dt;
            Mf(l) += eps;
        }
}

template<uint ndofs> YAKL_INLINE void compute_Mf(real5d Mfvar, const real5d edgefluxvar, real dt, int is, int js, int ks, int i, int j, int k, int n)
{
SArray<real,1,ndofs> Mf;
SArray<real,3,ndofs,ndims,2> edgeflux;
for (int l=0; l<ndofs; l++) {
  for (int m=0; m<2; m++) {
    for (int d=0; d<ndims; d++) {
        if (d==0) {edgeflux(l,d,m) = edgefluxvar(l+d*ndofs, k+ks, j+js, i+is+m, n); }
        if (d==1) {edgeflux(l,d,m) = edgefluxvar(l+d*ndofs, k+ks, j+js+m, i+is, n); }
  }}}
    calculate_Mf<ndofs>(Mf, edgeflux, dt);
    for (int l=0; l<ndofs; l++) { Mfvar(l, k+ks, j+js, i+is, n) = Mf(l); }
}

template<uint ndofs> YAKL_INLINE void compute_Mfext(real5d Mfvar, const real5d edgefluxvar, const real5d vertedgefluxvar, real dt, int is, int js, int ks, int i, int j, int k, int n)
{
SArray<real,1,ndofs> Mf;
SArray<real,3,ndofs,ndims,2> edgeflux;
SArray<real,2,ndofs,2> vertedgeflux;
for (int l=0; l<ndofs; l++) {
  for (int m=0; m<2; m++) {
    for (int d=0; d<ndims; d++) {
        if (d==0) {edgeflux(l,d,m) = edgefluxvar(l+d*ndofs, k+ks, j+js, i+is+m, n); }
        if (d==1) {edgeflux(l,d,m) = edgefluxvar(l+d*ndofs, k+ks, j+js+m, i+is, n); }
  }
vertedgeflux(l,m) = vertedgefluxvar(l, k+ks+m, j+js, i+is, n);
}}
    calculate_Mfext<ndofs>(Mf, edgeflux, vertedgeflux, dt);
    for (int l=0; l<ndofs; l++) { Mfvar(l, k+ks, j+js, i+is, n) = Mf(l); }
}



template<uint ndofs> void YAKL_INLINE calculate_phi(SArray<real,2,ndofs,ndims> &Phi, SArray<real,3,ndofs,ndims,2> const &q, SArray<real,3,ndofs,ndims,2> const &Mf, SArray<real,2,ndofs,ndims> const &edgeflux)
{

    real upwind_param;

        for (int l=0; l<ndofs; l++) {
          for (int d=0; d<ndims; d++) {
            upwind_param = copysign(1.0, edgeflux(l,d));
            upwind_param = 0.5_fp*(upwind_param + fabs(upwind_param));
            Phi(l,d) = mymin(1._fp, q(l,d,1)/Mf(l,d,1)) * (1._fp - upwind_param) + mymin(1._fp, q(l,d,0)/Mf(l,d,0)) * upwind_param;
          }
    }
}

template<uint ndofs> void YAKL_INLINE calculate_phivert(SArray<real,1,ndofs> &Phivert, SArray<real,2,ndofs,2> const &q, SArray<real,2,ndofs,2> const &Mf, SArray<real,1,ndofs> const &vertedgeflux)
{
    real upwind_param;
        for (int l=0; l<ndofs; l++) {
            upwind_param = copysign(1.0, vertedgeflux(l));
            upwind_param = 0.5_fp*(upwind_param + fabs(upwind_param));
            Phivert(l) = mymin(1._fp, q(l,1)/Mf(l,1)) * (1. - upwind_param) + mymin(1._fp, q(l,0)/Mf(l,0)) * upwind_param;
    }
}

template<uint ndofs> YAKL_INLINE void compute_Phi(real5d Phivar, const real5d edgefluxvar, const real5d Mfvar, const real5d qvar, int is, int js, int ks, int i, int j, int k, int n)
{
SArray<real,2,ndofs,ndims> Phi;
SArray<real,3,ndofs,ndims,2> const q;
SArray<real,3,ndofs,ndims,2> const Mf;
SArray<real,2,ndofs,ndims> edgeflux;

for (int d=0; d<ndims; d++) {
    for (int l=0; l<ndofs; l++) {
      edgeflux(l,d) = edgefluxvar(l+d*ndofs, k+ks, j+js, i+is, n);
      for (int m=0; m<2; m++) {
        if (d==0) {
          Mf(l,d,m) = Mfvar(l, k+ks, j+js, i+is+m-1, n);
          q(l,d,m) = qvar(l, k+ks, j+js, i+is+m-1, n);
        }
        if (d==1) {
          Mf(l,d,m) = Mfvar(l, k+ks, j+js+m-1, i+is, n);
          q(l,d,m) = qvar(l, k+ks, j+js+m-1, i+is, n);
        }

  }}}
calculate_phi<ndofs>(Phi, q, Mf , edgeflux);
for (int d=0; d<ndims; d++) {
    for (int l=0; l<ndofs; l++) {
  Phivar(l+d*ndofs, k+ks, j+js, i+is, n) = Phi(l,d);
    }}
}


template<uint ndofs> YAKL_INLINE void compute_Phivert(real5d Phivertvar, const real5d vertedgefluxvar, const real5d Mfvar, const real5d qvar, int is, int js, int ks, int i, int j, int k, int n)
{
SArray<real,1,ndofs> Phivert;
SArray<real,2,ndofs,2> const qvert;
SArray<real,2,ndofs,2> const Mfvert;
SArray<real,1,ndofs> vertedgeflux;

for (int l=0; l<ndofs; l++) {
vertedgeflux(l) = vertedgefluxvar(l, k+ks, j+js, i+is, n);
for (int m=0; m<2; m++) {
 Mfvert(l,m) = Mfvar(l, k+ks+m-1, j+js, i+is, n);
 qvert(l,m) = qvar(l, k+ks+m-1, j+js, i+is, n);
}
}
calculate_phivert<ndofs>(Phivert, qvert, Mfvert , vertedgeflux);
for (int l=0; l<ndofs; l++) {
Phivertvar(l, k+ks, j+js, i+is, n) = Phivert(l);
}
}
