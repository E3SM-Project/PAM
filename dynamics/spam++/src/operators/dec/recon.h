
#ifndef _RECON_H_
#define _RECON_H_

#include "common.h"
#include "weno_recon.h"
#include "cfv_recon.h"
#include "weno_func_recon.h"

template<uint ndofs> void YAKL_INLINE centered_recon(SArray<real,ndofs,ndims> &recon, SArray<real,ndofs,ndims,2> const &edgerecon) {

    for (int l=0; l<ndofs; l++) {
      for (int k=0; k<ndims; k++) {
      recon(l,k) = 0.5* (edgerecon(l,k,1) + edgerecon(l,k,0));
    }
  }
}

template<uint ndofs> void YAKL_INLINE upwind_recon(SArray<real,ndofs,ndims> &recon, SArray<real,ndofs,ndims,2> const &edgerecon, SArray<real,ndims> const &flux) {

    real upwind_param;
    for (int l=0; l<ndofs; l++) {
      for (int k=0; k<ndims; k++) {
      upwind_param = copysign(1.0, flux(k));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      recon(l,k) = edgerecon(l,k,1) * (1. - upwind_param) + edgerecon(l,k,0) * upwind_param;
    }
  }
}







template <uint ndofs, RECONSTRUCTION_TYPE recontype, uint ord, uint tord=2, uint hs=(ord-1)/2> void YAKL_INLINE compute_primal_edge_recon(
  realArr edgereconvar, const realArr var, int is, int js, int ks, int i, int j, int k,
  SArray<real,ord,ord,ord> const &wenoRecon, SArray<real,ord,tord> const &to_gll, SArray<real,hs+2> const &wenoIdl, real wenoSigma)
  {

SArray<real,ndofs,ndims,ord> stencil;
SArray<real,ndofs,ndims,2> edgerecon;

for (int p=0; p<ord; p++) {
for (int l=0; l<ndofs; l++) {
  for (int d=0; d<ndims; d++) {
  if (d==0) {stencil(l,d,p) = var(l, k+ks, j+js, i+is+p-hs);}
  if (d==1) {stencil(l,d,p) = var(l, k+ks, j+js+p-hs, i+is);}
  if (d==2) {stencil(l,d,p) = var(l, k+ks+p-hs, j+js, i+is);}
}}}

if (recontype == RECONSTRUCTION_TYPE::CFV) { cfv<ndofs> (edgerecon, stencil); }
if (recontype == RECONSTRUCTION_TYPE::WENO) { weno<ndofs> (edgerecon, stencil); }
if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) { weno_func<ndofs, ord> (edgerecon, stencil, wenoRecon, to_gll, wenoIdl, wenoSigma); }

for (int d=0; d<ndims; d++) {
for (int l=0; l<ndofs; l++) {
edgereconvar(l*ndofs+2*d, k+ks, j+js, i+is) = edgerecon(l,d,0);
edgereconvar(l*ndofs+2*d+1, k+ks, j+js, i+is) = edgerecon(l,d,1);
}}

}





template <uint ndofs, RECONSTRUCTION_TYPE recontype, uint ord, uint tord=2, uint hs=(ord-1)/2> void YAKL_INLINE compute_dual_edge_recon(
  realArr edgereconvar, const realArr var, int is, int js, int ks, int i, int j, int k,
  SArray<real,ord,ord,ord> const &wenoRecon, SArray<real,ord,tord> const &to_gll, SArray<real,hs+2> const &wenoIdl, real wenoSigma)
  {

SArray<real,ndofs,ndims,ord> stencil;
SArray<real,ndofs,ndims,2> edgerecon;

for (int p=0; p<ord; p++) {
for (int l=0; l<ndofs; l++) {
  for (int d=ndims-1; d>=0; d--) {
  if (d==ndims-1) {stencil(l,d,p) = var(l, k+ks, j+js, i+is+p-hs);}
  if (d==ndims-2) {stencil(l,d,p) = var(l, k+ks, j+js+p-hs, i+is);}
  if (d==ndims-3) {stencil(l,d,p) = var(l, k+ks+p-hs, j+js, i+is);}
}}}

if (dual_reconstruction_type == RECONSTRUCTION_TYPE::CFV) { cfv<ndofs> (edgerecon, stencil); }
if (dual_reconstruction_type == RECONSTRUCTION_TYPE::WENO) { weno<ndofs> (edgerecon, stencil); }
if (dual_reconstruction_type == RECONSTRUCTION_TYPE::WENOFUNC) { weno_func<ndofs,ord> (edgerecon, stencil, wenoRecon, to_gll, wenoIdl, wenoSigma); }

for (int d=ndims-1; d>=0; d--) {
for (int l=0; l<ndofs; l++) {
edgereconvar(l*ndofs+2*d, k+ks, j+js, i+is) = edgerecon(l,d,0);
edgereconvar(l*ndofs+2*d+1, k+ks, j+js, i+is) = edgerecon(l,d,1);
}}

}



template <uint ndofs, RECONSTRUCTION_TYPE recontype> void YAKL_INLINE compute_primal_recon(
  realArr reconvar, const realArr edgereconvar, const realArr U, int is, int js, int ks, int i, int j, int k)
  {

SArray<real,ndofs,ndims> recon;
SArray<real,ndims> uvar;
SArray<real,ndofs,ndims,2> edgerecon;

for (int d=0; d<ndims; d++) {
  uvar(d) = U(d, k+ks, j+js, i+is);
for (int l=0; l<ndofs; l++) {
  if (d==0)
  {
  edgerecon(l,d,0) = edgereconvar(l*ndofs+2*d+1, k+ks, j+js, i+is-1);
  edgerecon(l,d,1) = edgereconvar(l*ndofs+2*d, k+ks, j+js, i+is);
}
if (d==1)
{
edgerecon(l,d,0) = edgereconvar(l*ndofs+2*d+1, k+ks, j+js-1, i+is);
edgerecon(l,d,1) = edgereconvar(l*ndofs+2*d, k+ks, j+js, i+is);
}
if (d==2)
{
edgerecon(l,d,0) = edgereconvar(l*ndofs+2*d+1, k+ks-1, j+js, i+is);
edgerecon(l,d,1) = edgereconvar(l*ndofs+2*d, k+ks, j+js, i+is);
}
}}

if (recontype == RECONSTRUCTION_TYPE::CFV) { centered_recon<ndofs>(recon, edgerecon); }
if (recontype == RECONSTRUCTION_TYPE::WENO) { upwind_recon<ndofs>(recon, edgerecon, uvar); }
if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) { upwind_recon<ndofs>(recon, edgerecon, uvar); }

for (int d=0; d<ndims; d++) {
for (int l=0; l<ndofs; l++) {
reconvar(l+d*ndofs, k+ks, j+js, i+is) = recon(l,d);
}}

}

template <uint ndofs, RECONSTRUCTION_TYPE recontype> void YAKL_INLINE compute_dual_recon(
  realArr reconvar, const realArr edgereconvar, const realArr UT, int is, int js, int ks, int i, int j, int k)
  {
    SArray<real,ndofs,ndims> recon;
    SArray<real,ndims> uvar;
    SArray<real,ndofs,ndims,2> edgerecon;

    for (int d=ndims-1; d>=0; d--) {
      uvar(d) = UT(d, k+ks, j+js, i+is);
      if (ndims ==2 && d==ndims-2) {uvar(d) = -uvar(d);} //corrects for "twist" in 2D

    for (int l=0; l<ndofs; l++) {
      if (d==ndims-1)
      {
      edgerecon(l,d,0) = edgereconvar(l*ndofs+2*d+1, k+ks, j+js, i+is);
      edgerecon(l,d,1) = edgereconvar(l*ndofs+2*d, k+ks, j+js, i+is+1);
    }
    if (d==ndims-2)
    {
    edgerecon(l,d,0) = edgereconvar(l*ndofs+2*d+1, k+ks, j+js, i+is);
    edgerecon(l,d,1) = edgereconvar(l*ndofs+2*d, k+ks, j+js+1, i+is);
  }
  if (d==ndims-3)
  {
  edgerecon(l,d,0) = edgereconvar(l*ndofs+2*d+1, k+ks, j+js, i+is);
  edgerecon(l,d,1) = edgereconvar(l*ndofs+2*d, k+ks+1, j+js, i+is);
  }
    }}

    if (recontype == RECONSTRUCTION_TYPE::CFV) { centered_recon<ndofs>(recon, edgerecon); }
    if (recontype == RECONSTRUCTION_TYPE::WENO) { upwind_recon<ndofs>(recon, edgerecon, uvar); }
    if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) { upwind_recon<ndofs>(recon, edgerecon, uvar); }

    for (int d=ndims-1; d>=0; d--) {
    for (int l=0; l<ndofs; l++) {
    reconvar(l+d*ndofs, k+ks, j+js, i+is) = recon(l,d);
  }}
  }
#endif
