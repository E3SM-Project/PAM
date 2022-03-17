#ifndef _FCT_H_
#define _FCT_H_

#include "common.h"



template<uint ndofs> void YAKL_INLINE calculate_edgeflux(SArray<real,2,ndofs,ndims> &edgeflux, SArray<real,2,ndofs,ndims> const &recon, SArray<real,1,ndims> const &flux)
{
  for (int l=0; l<ndofs; l++) {
  for (int d=0; d<ndims; d++) {
    edgeflux(l,d) = recon(l,d) * flux(d);
  }}
}

template<uint ndofs> YAKL_INLINE void compute_edgefluxes(real4d edgefluxvar, const real4d reconvar, const real4d U, int is, int js, int ks, int i, int j, int k)
{
SArray<real,2,ndofs,ndims> edgeflux;
SArray<real,2,ndofs,ndims> recon;
SArray<real,1,ndims> flux;
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

template<uint ndofs> void YAKL_INLINE calculate_vertedgeflux(SArray<real,1,ndofs> &edgeflux, SArray<real,1,ndofs> const &recon, real const &flux)
{
  for (int l=0; l<ndofs; l++) {
    edgeflux(l) = recon(l) * flux;
  }
}

template<uint ndofs> YAKL_INLINE void compute_vertedgefluxes(real4d vertedgefluxvar, const real4d vertreconvar, const real4d UW, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,ndofs> edgeflux;
  SArray<real,1,ndofs> recon;
  real flux;
  flux = UW(0, k+ks, j+js, i+is);
  for (int l=0; l<ndofs; l++) {
    recon(l) = vertreconvar(l, k+ks, j+js, i+is);
  }
  calculate_vertedgeflux<ndofs>(edgeflux, recon, flux);

  for (int l=0; l<ndofs; l++) {
  vertedgefluxvar(l, k+ks, j+js, i+is) = edgeflux(l);
  }
}


template<uint ndofs> void YAKL_INLINE calculate_Mf(SArray<real,1,ndofs> &Mf, SArray<real,3,ndofs,ndims,2> const &edgeflux, real dt)
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

template<uint ndofs> void YAKL_INLINE calculate_Mfext(SArray<real,1,ndofs> &Mf, SArray<real,3,ndofs,ndims,2> const &edgeflux, SArray<real,2,ndofs,2> const &vertedgeflux, real dt)
{
    real eps = 1.0e-8;
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

template<uint ndofs> YAKL_INLINE void compute_Mf(real4d Mfvar, const real4d edgefluxvar, real dt, int is, int js, int ks, int i, int j, int k)
{
SArray<real,1,ndofs> Mf;
SArray<real,3,ndofs,ndims,2> edgeflux;
for (int l=0; l<ndofs; l++) {
  for (int m=0; m<2; m++) {
    for (int d=0; d<ndims; d++) {
        if (d==0) {edgeflux(l,d,m) = edgefluxvar(l+d*ndofs, k+ks, j+js, i+is+m); }
        if (d==1) {edgeflux(l,d,m) = edgefluxvar(l+d*ndofs, k+ks, j+js+m, i+is); }
  }}}
    calculate_Mf<ndofs>(Mf, edgeflux, dt);
    for (int l=0; l<ndofs; l++) { Mfvar(l, k+ks, j+js, i+is) = Mf(l); }
}

template<uint ndofs> YAKL_INLINE void compute_Mfext(real4d Mfvar, const real4d edgefluxvar, const real4d vertedgefluxvar, real dt, int is, int js, int ks, int i, int j, int k)
{
SArray<real,1,ndofs> Mf;
SArray<real,3,ndofs,ndims,2> edgeflux;
SArray<real,2,ndofs,2> vertedgeflux;
for (int l=0; l<ndofs; l++) {
  for (int m=0; m<2; m++) {
    for (int d=0; d<ndims; d++) {
        if (d==0) {edgeflux(l,d,m) = edgefluxvar(l+d*ndofs, k+ks, j+js, i+is+m); }
        if (d==1) {edgeflux(l,d,m) = edgefluxvar(l+d*ndofs, k+ks, j+js+m, i+is); }
  }
vertedgeflux(l,m) = vertedgefluxvar(l, k+ks+m, j+js, i+is);
}}
    calculate_Mfext<ndofs>(Mf, edgeflux, vertedgeflux, dt);
    for (int l=0; l<ndofs; l++) { Mfvar(l, k+ks, j+js, i+is) = Mf(l); }
}



template<uint ndofs> void YAKL_INLINE calculate_phi(SArray<real,2,ndofs,ndims> &Phi, SArray<real,3,ndofs,ndims,2> const &q, SArray<real,3,ndofs,ndims,2> const &Mf, SArray<real,2,ndofs,ndims> const &edgeflux)
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

template<uint ndofs> void YAKL_INLINE calculate_phivert(SArray<real,1,ndofs> &Phivert, SArray<real,2,ndofs,2> const &q, SArray<real,2,ndofs,2> const &Mf, SArray<real,1,ndofs> const &vertedgeflux)
{
    real upwind_param;
        for (int l=0; l<ndofs; l++) {
            upwind_param = copysign(1.0, vertedgeflux(l));
            upwind_param = 0.5*(upwind_param + fabs(upwind_param));
            Phivert(l) = mymin(1., q(l,1)/Mf(l,1)) * (1. - upwind_param) + mymin(1., q(l,0)/Mf(l,0)) * upwind_param;
    }
}

template<uint ndofs> YAKL_INLINE void compute_Phi(real4d Phivar, const real4d edgefluxvar, const real4d Mfvar, const real4d qvar, int is, int js, int ks, int i, int j, int k)
{
SArray<real,2,ndofs,ndims> Phi;
SArray<real,3,ndofs,ndims,2> const q;
SArray<real,3,ndofs,ndims,2> const Mf;
SArray<real,2,ndofs,ndims> edgeflux;

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

  }}}
calculate_phi<ndofs>(Phi, q, Mf , edgeflux);
for (int d=0; d<ndims; d++) {
    for (int l=0; l<ndofs; l++) {
  Phivar(l+d*ndofs, k+ks, j+js, i+is) = Phi(l,d);
    }}
}


template<uint ndofs> YAKL_INLINE void compute_Phivert(real4d Phivertvar, const real4d vertedgefluxvar, const real4d Mfvar, const real4d qvar, int is, int js, int ks, int i, int j, int k)
{
SArray<real,1,ndofs> Phivert;
SArray<real,2,ndofs,2> const qvert;
SArray<real,2,ndofs,2> const Mfvert;
SArray<real,1,ndofs> vertedgeflux;

for (int l=0; l<ndofs; l++) {
vertedgeflux(l) = vertedgefluxvar(l, k+ks, j+js, i+is);
for (int m=0; m<2; m++) {
 Mfvert(l,m) = Mfvar(l, k+ks+m-1, j+js, i+is);
 qvert(l,m) = qvar(l, k+ks+m-1, j+js, i+is);
}
}
calculate_phivert<ndofs>(Phivert, qvert, Mfvert , vertedgeflux);
for (int l=0; l<ndofs; l++) {
Phivertvar(l, k+ks, j+js, i+is) = Phivert(l);
}
}

#endif
