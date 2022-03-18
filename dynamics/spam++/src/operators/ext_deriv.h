#ifndef _EXT_DERIV_H_
#define _EXT_DERIV_H_

#include "common.h"


template<uint ndofs> void YAKL_INLINE wDbar2( SArray<real,1,ndofs> &var, SArray<real,3,ndofs,ndims,2> const &recon, SArray<real,2,ndims,2> const &flux ) {

  for (int l=0; l<ndofs; l++) {
    var(l) = 0.;
    for (int d=0; d<ndims; d++) {
      var(l) += flux(d,1) * recon(l,d,1) - flux(d,0) * recon(l,d,0);
    }
  }
}

template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> YAKL_INLINE void compute_wDbar2(real4d tendvar, const real4d reconvar, const real4d U, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,ndofs> tend;
  SArray<real,3,ndofs,ndims,2> recon;
  SArray<real,1,ndims,2> flux;

  for (int d=0; d<ndims; d++) {
    for (int m=0; m<2; m++) {
    if (d==0) { flux(d,m) = U(d, k+ks, j+js, i+is+m);}
    if (d==1) { flux(d,m) = U(d, k+ks, j+js+m, i+is);}
    //if (d==2) { flux(d,m) = U(d, k+ks+m, j+js, i+is);}
  for (int l=0; l<ndofs; l++) {
    if (d==0) { recon(l,d,m) = reconvar(l+d*ndofs, k+ks, j+js, i+is+m);}
    if (d==1) { recon(l,d,m) = reconvar(l+d*ndofs, k+ks, j+js+m, i+is);}
    //if (d==2) { recon(l,d,m) = reconvar(l+d*ndofs, k+ks+m, j+js, i+is);}
  }}}

  wDbar2<ndofs>( tend, recon, flux );

  for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) = tend(l);}

  if (addmode == ADD_MODE::REPLACE) {for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) = tend(l);}}
  if (addmode == ADD_MODE::ADD) {for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) += tend(l);}}
}

template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> YAKL_INLINE void compute_wDbar2_fct(real4d tendvar, const real4d reconvar, const real4d phivar, const real4d U, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,ndofs> tend;
  SArray<real,3,ndofs,ndims,2> recon;
  SArray<real,2,ndims,2> flux;

  for (int d=0; d<ndims; d++) {
    for (int m=0; m<2; m++) {
    if (d==0) { flux(d,m) = U(d, k+ks, j+js, i+is+m);}
    if (d==1) { flux(d,m) = U(d, k+ks, j+js+m, i+is);}
    //if (d==2) { flux(d,m) = U(d, k+ks+m, j+js, i+is);}
  for (int l=0; l<ndofs; l++) {
    if (d==0) { recon(l,d,m) = reconvar(l+d*ndofs, k+ks, j+js, i+is+m) * phivar(l+d*ndofs, k+ks, j+js, i+is+m);}
    if (d==1) { recon(l,d,m) = reconvar(l+d*ndofs, k+ks, j+js+m, i+is) * phivar(l+d*ndofs, k+ks, j+js+m, i+is);}
    //if (d==2) { recon(l,d,m) = reconvar(l+d*ndofs, k+ks+m, j+js, i+is) * phivar(l+d*ndofs, k+ks+m, j+js, i+is);}
  }}}

  wDbar2<ndofs>( tend, recon, flux );

if (addmode == ADD_MODE::REPLACE) {for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) = tend(l);}}
if (addmode == ADD_MODE::ADD) {for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) += tend(l);}}

}

template<uint ndofs> void YAKL_INLINE wDvbar( SArray<real,1,ndofs> &var, SArray<real,2,ndofs,2> const &recon, SArray<real,1,2> const &flux ) {
  for (int l=0; l<ndofs; l++) {
    var(l) = flux(1) * recon(l,1) - flux(0) * recon(l,0);
  }
}

template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> YAKL_INLINE void compute_wDvbar_fct(real4d tendvar, const real4d vertreconvar, const real4d vertphivar, const real4d UW, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,ndofs> tend;
  SArray<real,2,ndofs,2> recon;
  SArray<real,1,2> flux;
  for (int m=0; m<2; m++) {
    flux(m) = UW(0, k+ks+m, j+js, i+is);
    for (int l=0; l<ndofs; l++) {
      recon(l,m) = vertreconvar(l, k+ks+m, j+js, i+is) * vertphivar(l, k+ks+m, j+js, i+is);
    }}

wDvbar<ndofs>( tend, recon, flux );

if (addmode == ADD_MODE::REPLACE) {for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) = tend(l);}}
if (addmode == ADD_MODE::ADD) {for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) += tend(l);}}
}

template<uint ndofs> void YAKL_INLINE wD1( SArray<real,1,ndims> &var, SArray<real,2,ndofs,ndims> const &recon, SArray<real,3,ndofs,ndims,2> const &dens) {

  for (int d=0; d<ndims; d++) {
    var(d) = 0.0;
  for (int l=0; l<ndofs; l++) {
      var(d) += recon(l,d) * (dens(l,d,1) - dens(l,d,0));
    }
  }
}

template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_wD1(real4d tendvar, const real4d reconvar, const real4d densvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,ndims> tend;
  SArray<real,2,ndofs,ndims> recon;
  SArray<real,3,ndofs,ndims,2> dens;
  for (int l=0; l<ndofs; l++) {
    for (int d=0; d<ndims; d++) {
      recon(l,d) = reconvar(l+d*ndofs,k+ks,j+js,i+is);
      dens(l,d,1) = densvar(l,k+ks,j+js,i+is);
      if (d==0) {dens(l,d,0) = densvar(l,k+ks,j+js,i+is-1);}
      if (d==1) {dens(l,d,0) = densvar(l,k+ks,j+js-1,i+is);}
      //if (d==2) {dens(l,d,0) = densvar(l,k+ks-1,j+js,i+is);}
    }}
  wD1<ndofs>( tend, recon, dens );
    for (int d=0; d<ndims; d++) {
      if (addmode == ADD_MODE::REPLACE) {tendvar(d,k+ks,j+js,i+is) = tend(d);}
      if (addmode == ADD_MODE::ADD) {tendvar(d,k+ks,j+js,i+is) += tend(d);}
    }
}

template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_wD1_fct(real4d tendvar, const real4d reconvar, const real4d Phivar, const real4d densvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,ndims> tend;
  SArray<real,2,ndofs,ndims> recon;
  SArray<real,3,ndofs,ndims,2> dens;
  for (int l=0; l<ndofs; l++) {
    for (int d=0; d<ndims; d++) {
      recon(l,d) = reconvar(l+d*ndofs,k+ks,j+js,i+is) * Phivar(l+d*ndofs,k+ks,j+js,i+is);
      dens(l,d,1) = densvar(l,k+ks,j+js,i+is);
      if (d==0) {dens(l,d,0) = densvar(l,k+ks,j+js,i+is-1);}
      if (d==1) {dens(l,d,0) = densvar(l,k+ks,j+js-1,i+is);}
      //if (d==2) {dens(l,d,0) = densvar(l,k+ks-1,j+js,i+is);}
    }}
  wD1<ndofs>( tend, recon, dens );
    for (int d=0; d<ndims; d++) {
      if (addmode == ADD_MODE::REPLACE) {tendvar(d,k+ks,j+js,i+is) = tend(d);}
      if (addmode == ADD_MODE::ADD) {tendvar(d,k+ks,j+js,i+is) += tend(d);}
    }
}



template<uint ndofs> void YAKL_INLINE wDv( real &var, SArray<real,1,ndofs> const &recon, SArray<real,2,ndofs,2> const &dens) {

    var = 0.0;
  for (int l=0; l<ndofs; l++) {
      var += recon(l) * (dens(l,1) - dens(l,0));
    }
  }

  template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_wDv_fct(real4d tendvar, const real4d vertreconvar, const real4d Phivertvar, const real4d densvar, int is, int js, int ks, int i, int j, int k)
  {
    real tend;
    SArray<real,1,ndofs> recon;
    SArray<real,2,ndofs,2> dens;
    for (int l=0; l<ndofs; l++) {
      //Need to add 1 to k here because UW has an extra dof at the bottom
        recon(l) = vertreconvar(l,k+ks+1,j+js,i+is) * Phivertvar(l,k+ks+1,j+js,i+is);
        dens(l,1) = densvar(l,k+ks+1,j+js,i+is);
        dens(l,0) = densvar(l,k+ks,j+js,i+is);
      }
    wDv<ndofs>( tend, recon, dens );
        if (addmode == ADD_MODE::REPLACE) {tendvar(0,k+ks,j+js,i+is) = tend;}
        if (addmode == ADD_MODE::ADD) {tendvar(0,k+ks,j+js,i+is) += tend;}
  }


template<uint ndofs> void YAKL_INLINE D2( SArray<real,1,ndofs> &var, SArray<real,2,ndofs,4> const &flux) {

      for (int l=0; l<ndofs; l++) {

      var(l) = (flux(l,0) - flux(l,1) - flux(l,2) + flux(l,3));
  }

}

template<uint ndofs> void YAKL_INLINE compute_D2(SArray<real,1,ndofs> &tend, const real4d fluxvar, int is, int js, int ks, int i, int j, int k)
{
SArray<real,2,ndofs,4> flux;

for (int l=0; l<ndofs; l++) {
  flux(l,0) = fluxvar(l+1*ndofs, k+ks, j+js, i+is);
  flux(l,1) = fluxvar(l+0*ndofs, k+ks, j+js, i+is);
  flux(l,2) = fluxvar(l+1*ndofs, k+ks, j+js, i+is-1);
  flux(l,3) = fluxvar(l+0*ndofs, k+ks, j+js-1, i+is);
}

D2<ndofs>( tend, flux );

}

template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_D2(real4d tendvar, const real4d fluxvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,ndofs> tend;
  compute_D2<ndofs> (tend, fluxvar, is, js, ks, i, j, k);
  if (addmode == ADD_MODE::REPLACE) {for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) = tend(l);}}
  if (addmode == ADD_MODE::ADD) {for (int l=0; l<ndofs; l++) { tendvar(l, k+ks, j+js, i+is) += tend(l);}}

}






template<uint ndofs> void YAKL_INLINE compute_Dxz(SArray<real,1,ndofs> &tend, const real4d v, const real4d w, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,4> flux;
  for (int l=0; l<ndofs; l++) {
    flux(0) = v(l, k+ks, j+js, i+is); //v1 +
    flux(1) = v(l, k+ks+1, j+js, i+is); //v1 -
    flux(2) = w(l, k+ks, j+js, i+is); // w + 
    flux(3) = w(l, k+ks, j+js, i+is-1); //w -
    tend(l) = (flux(0) - flux(1) + flux(2) - flux(3));
  }
}
template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Dxz(real4d tendvar, const real4d v, const real4d w, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,ndofs> tend;
  compute_Dxz<ndofs>(tend, v, w, is, js, ks, i, j, k);
  for (int l=0; l<ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {tendvar(l, k+ks, j+js, i+is) = tend(l);}
    if (addmode == ADD_MODE::ADD) {tendvar(l, k+ks, j+js, i+is) += tend(l);}
  }
}



#endif
