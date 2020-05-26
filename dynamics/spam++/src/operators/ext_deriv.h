#ifndef _EXT_DERIV_H_
#define _EXT_DERIV_H_

#include "common.h"


template<uint ndofs> void YAKL_INLINE wDbar2( SArray<real,ndofs> &var, SArray<real,ndofs,ndims,2> const &recon, SArray<real,ndims,2> const &flux ) {

  for (int l=0; l<ndofs; l++) {
    var(l) = 0.;
    for (int d=0; d<ndims; d++) {
      var(l) += flux(d,1) * recon(l,d,1) - flux(d,0) * recon(l,d,0);
    }
  }
}

template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> YAKL_INLINE void compute_wDbar2(realArr tendvar, const realArr reconvar, const realArr U, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,ndofs> tend;
  SArray<real,ndofs,ndims,2> recon;
  SArray<real,ndims,2> flux;

  for (int d=0; d<ndims; d++) {
    for (int m=0; m<2; m++) {
    if (d==0) { flux(d,m) = U(d, k+ks, j+js, i+is+m);}
    if (d==1) { flux(d,m) = U(d, k+ks, j+js+m, i+is);}
    if (d==2) { flux(d,m) = U(d, k+ks+m, j+js, i+is);}
  for (int l=0; l<ndofs; l++) {
    if (d==0) { recon(l,d,m) = reconvar(l+d*ndofs, k+ks, j+js, i+is+m);}
    if (d==1) { recon(l,d,m) = reconvar(l+d*ndofs, k+ks, j+js+m, i+is);}
    if (d==2) { recon(l,d,m) = reconvar(l+d*ndofs, k+ks+m, j+js, i+is);}
  }}}

  wDbar2<ndofs>( tend, recon, flux );

  for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) = tend(l);}

  if (addmode == ADD_MODE::REPLACE) {for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) = tend(l);}}
  if (addmode == ADD_MODE::ADD) {for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) += tend(l);}}
}

template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> YAKL_INLINE void compute_wDbar2_fct(realArr tendvar, const realArr reconvar, const realArr phivar, const realArr U, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,ndofs> tend;
  SArray<real,ndofs,ndims,2> recon;
  SArray<real,ndims,2> flux;

  for (int d=0; d<ndims; d++) {
    for (int m=0; m<2; m++) {
    if (d==0) { flux(d,m) = U(d, k+ks, j+js, i+is+m);}
    if (d==1) { flux(d,m) = U(d, k+ks, j+js+m, i+is);}
    if (d==2) { flux(d,m) = U(d, k+ks+m, j+js, i+is);}
  for (int l=0; l<ndofs; l++) {
    if (d==0) { recon(l,d,m) = reconvar(l+d*ndofs, k+ks, j+js, i+is+m) * phivar(l+d*ndofs, k+ks, j+js, i+is+m);}
    if (d==1) { recon(l,d,m) = reconvar(l+d*ndofs, k+ks, j+js+m, i+is) * phivar(l+d*ndofs, k+ks, j+js+m, i+is);}
    if (d==2) { recon(l,d,m) = reconvar(l+d*ndofs, k+ks+m, j+js, i+is) * phivar(l+d*ndofs, k+ks+m, j+js, i+is);}
  }}}

  wDbar2<ndofs>( tend, recon, flux );

if (addmode == ADD_MODE::REPLACE) {for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) = tend(l);}}
if (addmode == ADD_MODE::ADD) {for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) += tend(l);}}

}

template<uint ndofs> void YAKL_INLINE wD1( SArray<real,ndofs,ndims> &var, SArray<real,ndofs,ndims> const &recon, SArray<real,ndofs,ndims,2> const &dens) {

  for (int l=0; l<ndofs; l++) {
    for (int d=0; d<ndims; d++) {
      var(l,d) = recon(l,d) * (dens(l,d,1) - dens(l,d,0));
    }
  }
}

template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_wD1(realArr tendvar, const realArr reconvar, const realArr densvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,ndofs,ndims> tend;
  SArray<real,ndofs,ndims> recon;
  SArray<real,ndofs,ndims,2> dens;
  for (int l=0; l<ndofs; l++) {
    for (int d=0; d<ndims; d++) {
      recon(l,d) = reconvar(l+d*ndofs,k+ks,j+js,i+is);
      dens(l,d,1) = densvar(l,k+ks,j+js,i+is);
      if (d==0) {dens(l,d,0) = densvar(l,k+ks,j+js,i+is-1);}
      if (d==1) {dens(l,d,0) = densvar(l,k+ks,j+js-1,i+is);}
      if (d==2) {dens(l,d,0) = densvar(l,k+ks-1,j+js,i+is);}
    }}
  wD1<ndofs>( tend, recon, dens );
  for (int l=0; l<ndofs; l++) {
    for (int d=0; d<ndims; d++) {
      if (addmode == ADD_MODE::REPLACE) {tendvar(l+d*ndofs,k+ks,j+js,i+is) = tend(l,d);}
      if (addmode == ADD_MODE::ADD) {tendvar(l+d*ndofs,k+ks,j+js,i+is) += tend(l,d);}
    }}
}

template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_wD1_fct(realArr tendvar, const realArr reconvar, const realArr Phivar, const realArr densvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,ndofs,ndims> tend;
  SArray<real,ndofs,ndims> recon;
  SArray<real,ndofs,ndims,2> dens;
  for (int l=0; l<ndofs; l++) {
    for (int d=0; d<ndims; d++) {
      recon(l,d) = reconvar(l+d*ndofs,k+ks,j+js,i+is) * Phivar(l+d*ndofs,k+ks,j+js,i+is);
      dens(l,d,1) = densvar(l,k+ks,j+js,i+is);
      if (d==0) {dens(l,d,0) = densvar(l,k+ks,j+js,i+is-1);}
      if (d==1) {dens(l,d,0) = densvar(l,k+ks,j+js-1,i+is);}
      if (d==2) {dens(l,d,0) = densvar(l,k+ks-1,j+js,i+is);}
    }}
  wD1<ndofs>( tend, recon, dens );
  for (int l=0; l<ndofs; l++) {
    for (int d=0; d<ndims; d++) {
      if (addmode == ADD_MODE::REPLACE) {tendvar(l+d*ndofs,k+ks,j+js,i+is) = tend(l,d);}
      if (addmode == ADD_MODE::ADD) {tendvar(l+d*ndofs,k+ks,j+js,i+is) += tend(l,d);}
    }}
}

// CAN WE GENERALIZE THIS TO 2D/3D?...maybe...

template<uint ndofs> void YAKL_INLINE D2( SArray<real,ndofs> &var, SArray<real,ndofs,4> const &flux) {

      for (int l=0; l<ndofs; l++) {

      var(l) = (flux(l,0) - flux(l,1) - flux(l,2) + flux(l,3));
  }

}

template<uint ndofs> void YAKL_INLINE compute_D2(SArray<real,ndofs> &tend, const realArr fluxvar, int is, int js, int ks, int i, int j, int k)
{
SArray<real,ndofs,4> flux;

for (int l=0; l<ndofs; l++) {
  flux(l,0) = fluxvar(l+1*ndofs, k+ks, j+js, i+is);
  flux(l,1) = fluxvar(l+0*ndofs, k+ks, j+js, i+is);
  flux(l,2) = fluxvar(l+1*ndofs, k+ks, j+js, i+is-1);
  flux(l,3) = fluxvar(l+0*ndofs, k+ks, j+js-1, i+is);
}

D2<ndofs>( tend, flux );

}

template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_D2(realArr tendvar, const realArr fluxvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,ndofs> tend;
  compute_D2<ndofs> (tend, fluxvar, is, js, ks, i, j, k);
  if (addmode == ADD_MODE::REPLACE) {for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) = tend(l);}}
  if (addmode == ADD_MODE::ADD) {for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) += tend(l);}}

}


#endif
