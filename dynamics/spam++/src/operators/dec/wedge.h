

#ifndef _WEDGE_H_
#define _WEDGE_H_

#include "common.h"
#include <cmath>

// Q
  template<uint ndofs> void YAKL_INLINE Q2D( SArray<real,ndofs,2> &vel, SArray<real,ndofs,2,5> const &recon, SArray<real,2,4> const &flux) {

      for (int l=0; l<ndofs; l++) {

      //x-dir
      vel(l,0) = -1./8. * (flux(0,0)*recon(l,0,0) + flux(0,1)*recon(l,0,1) + flux(0,2)*recon(l,0,2) + flux(0,3)*recon(l,0,3) +
                           flux(0,0)*recon(l,0,4) + flux(0,1)*recon(l,0,4) + flux(0,2)*recon(l,0,4) + flux(0,3)*recon(l,0,4));

      //y-dir
      vel(l,1) =  1./8. * (flux(1,0)*recon(l,1,0) + flux(1,1)*recon(l,1,1) + flux(1,2)*recon(l,1,2) + flux(1,3)*recon(l,1,3) +
                           flux(1,0)*recon(l,1,4) + flux(1,1)*recon(l,1,4) + flux(1,2)*recon(l,1,4) + flux(1,3)*recon(l,1,4));
  }
}


template<uint ndofs> void YAKL_INLINE Q2D_nonEC( SArray<real,ndofs,2> &vel, SArray<real,ndofs,2> const &recon, SArray<real,2,4> const &flux) {
  for (int l=0; l<ndofs; l++) {
  //x-dir
  vel(l,0) = -1./4. * recon(l,0) * (flux(0,0) + flux(0,1) + flux(0,2) + flux(0,3));
  //y-dir
  vel(l,1) = 1./4. * recon(l,1) * (flux(1,0) + flux(1,1) + flux(1,2) + flux(1,3));
}

}


template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Q_EC(realArr qflux, const realArr reconvar, const realArr Uvar, int is, int js, int ks, int i, int j, int k)
{

SArray<real,ndofs,2> vel;
SArray<real,2,4> flux;
flux(0,0) = Uvar(1, k+ks, j+js, i+is);
flux(0,1) = Uvar(1, k+ks, j+js, i+is-1);
flux(0,2) = Uvar(1, k+ks, j+js+1, i+is);
flux(0,3) = Uvar(1, k+ks, j+js+1, i+is-1);

flux(1,0) = Uvar(0, k+ks, j+js, i+is);
flux(1,1) = Uvar(0, k+ks, j+js, i+is+1);
flux(1,2) = Uvar(0, k+ks, j+js-1, i+is);
flux(1,3) = Uvar(0, k+ks, j+js-1, i+is+1);

   SArray<real,ndofs,2,5> recon;
   for (int l=0; l<ndofs; l++) {
     recon(l,0,0) = reconvar(l+1*ndofs, k+ks, j+js, i+is);
     recon(l,0,1) = reconvar(l+1*ndofs, k+ks, j+js, i+is-1);
     recon(l,0,2) = reconvar(l+1*ndofs, k+ks, j+js+1, i+is);
     recon(l,0,3) = reconvar(l+1*ndofs, k+ks, j+js+1, i+is-1);
     recon(l,0,4) = reconvar(l+0*ndofs, k+ks, j+js, i+is);

     recon(l,1,0) = reconvar(l+0*ndofs, k+ks, j+js, i+is);
     recon(l,1,1) = reconvar(l+0*ndofs, k+ks, j+js, i+is+1);
     recon(l,1,2) = reconvar(l+0*ndofs, k+ks, j+js-1, i+is);
     recon(l,1,3) = reconvar(l+0*ndofs, k+ks, j+js-1, i+is+1);
     recon(l,1,4) = reconvar(l+1*ndofs, k+ks, j+js, i+is);
}
   Q2D<ndofs>( vel, recon, flux);

   if (addmode == ADD_MODE::REPLACE){
       for (int l=0; l<ndofs; l++) {
         qflux(l+0*ndofs, k+ks, j+js, i+is) = vel(l,0) ;
         qflux(l+1*ndofs, k+ks, j+js, i+is) = vel(l,1) ;
       }
     }
   if (addmode == ADD_MODE::ADD){
       for (int l=0; l<ndofs; l++) {
         qflux(l+0*ndofs, k+ks, j+js, i+is) += vel(l,0) ;
         qflux(l+1*ndofs, k+ks, j+js, i+is) += vel(l,1) ;
       }
     }

}


template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Q_nonEC(realArr qflux, const realArr reconvar, const realArr Uvar, int is, int js, int ks, int i, int j, int k)
{

  SArray<real,ndofs,2> vel;
  SArray<real,2,4> flux;
  flux(0,0) = Uvar(1, k+ks, j+js, i+is);
  flux(0,1) = Uvar(1, k+ks, j+js, i+is-1);
  flux(0,2) = Uvar(1, k+ks, j+js+1, i+is);
  flux(0,3) = Uvar(1, k+ks, j+js+1, i+is-1);

  flux(1,0) = Uvar(0, k+ks, j+js, i+is);
  flux(1,1) = Uvar(0, k+ks, j+js, i+is+1);
  flux(1,2) = Uvar(0, k+ks, j+js-1, i+is);
  flux(1,3) = Uvar(0, k+ks, j+js-1, i+is+1);

    SArray<real,ndofs,2> recon;
    for (int l=0; l<ndofs; l++) {
      recon(l,0) = reconvar(l+0*ndofs, k+ks, j+js, i+is);
      recon(l,1) = reconvar(l+1*ndofs, k+ks, j+js, i+is);
  }
    Q2D_nonEC<ndofs>( vel, recon, flux);

if (addmode == ADD_MODE::REPLACE){
    for (int l=0; l<ndofs; l++) {
      qflux(l+0*ndofs, k+ks, j+js, i+is) = vel(l,0) ;
      qflux(l+1*ndofs, k+ks, j+js, i+is) = vel(l,1) ;
    }
  }
if (addmode == ADD_MODE::ADD){
    for (int l=0; l<ndofs; l++) {
      qflux(l+0*ndofs, k+ks, j+js, i+is) += vel(l,0) ;
      qflux(l+1*ndofs, k+ks, j+js, i+is) += vel(l,1) ;
    }
  }

}

// W
void YAKL_INLINE W2D( SArray<real,2> &vel, SArray<real,2,4> const &flux) {
        //x-dir
        vel(0) = -1./4. * (flux(0,0) + flux(0,1) + flux(0,2) + flux(0,3));
        //y-dir
        vel(1) = 1./4. * (flux(1,0) + flux(1,1) + flux(1,2) + flux(1,3));
  }



void YAKL_INLINE compute_W(realArr UTvar, realArr Uvar, int is, int js, int ks, int i, int j, int k)
{

SArray<real,2> ut;
SArray<real,2,4> flux;
flux(0,0) = Uvar(1, k+ks, j+js, i+is);
flux(0,1) = Uvar(1, k+ks, j+js, i+is-1);
flux(0,2) = Uvar(1, k+ks, j+js+1, i+is);
flux(0,3) = Uvar(1, k+ks, j+js+1, i+is-1);

flux(1,0) = Uvar(0, k+ks, j+js, i+is);
flux(1,1) = Uvar(0, k+ks, j+js, i+is+1);
flux(1,2) = Uvar(0, k+ks, j+js-1, i+is);
flux(1,3) = Uvar(0, k+ks, j+js-1, i+is+1);

W2D(ut, flux);
  UTvar(0, k+ks, j+js, i+is) = ut(0) ;
  UTvar(1, k+ks, j+js, i+is) = ut(1) ;
}



// R
void YAKL_INLINE R( SArray<real,1> &vard, SArray<real,2> const &varp) {
  vard(0) = 1./2. * (varp(0) + varp(1));
}
void YAKL_INLINE R( SArray<real,1> &vard, SArray<real,4> const &varp) {
  vard(0) = 1./4. * (varp(0) + varp(1) + varp(2) + varp(3));
}
void YAKL_INLINE R( SArray<real,1> &vard, SArray<real,8> const &varp) {
  vard(0) = 1./8. * (varp(0) + varp(1) + varp(2) + varp(3) + varp(4) + varp(5) + varp(6) + varp(7));
}

template<uint dof> void YAKL_INLINE compute_R(real &var, realArr densvar, int is, int js, int ks, int i, int j, int k)
{
SArray<real,1> dualdens;
  if (ndims==1)
  {
    SArray<real,2> dens;
    dens(0) = densvar(dof, k+ks, j+js, i+is);
    dens(1) = densvar(dof, k+ks, j+js, i+is-1);
    R(dualdens, dens);
  }
  if (ndims==2)
  {
    SArray<real,4> dens;
    dens(0) = densvar(dof, k+ks, j+js, i+is);
    dens(1) = densvar(dof, k+ks, j+js, i+is-1);
    dens(2) = densvar(dof, k+ks, j+js-1, i+is);
    dens(3) = densvar(dof, k+ks, j+js-1, i+is-1);
    R(dualdens, dens);
  }
  if (ndims==3)
  {
    SArray<real,8> dens;
    dens(0) = densvar(dof, k+ks, j+js, i+is);
    dens(1) = densvar(dof, k+ks, j+js, i+is-1);
    dens(2) = densvar(dof, k+ks, j+js-1, i+is);
    dens(3) = densvar(dof, k+ks, j+js-1, i+is-1);
    dens(4) = densvar(dof, k+ks-1, j+js, i+is);
    dens(5) = densvar(dof, k+ks-1, j+js, i+is-1);
    dens(6) = densvar(dof, k+ks-1, j+js-1, i+is);
    dens(7) = densvar(dof, k+ks-1, j+js-1, i+is-1);
    R(dualdens, dens);
  }
var = dualdens(0);
}

template<uint dof> void YAKL_INLINE compute_R(realArr var, realArr densvar, int is, int js, int ks, int i, int j, int k)
{
real var2;
compute_R<dof> (var2,densvar,is,js,ks,i,j,k);
var(dof, k+ks, j+js, i+is) = var2;
}


//phi
void YAKL_INLINE phi( SArray<real,ndims> &vare, SArray<real,ndims,2> const &var0) {
  for(int d=0; d<ndims; d++)
  {
    vare(d) = 1./2. * (var0(d,0) + var0(d,1));
  }
  }

  template<uint dof> void YAKL_INLINE compute_phi(realArr var, realArr densvar, int is, int js, int ks, int i, int j, int k)
  {
  SArray<real,ndims> xe;
  SArray<real,ndims,2> x0;
  for (int d=0; d<ndims; d++)
  {
    if (d == 0)
    {
      x0(d,0) = densvar(dof, k+ks, j+js, i+is);
      x0(d,1) = densvar(dof, k+ks, j+js, i+is-1);
    }
    if (d == 1)
    {
      x0(d,0) = densvar(dof, k+ks, j+js, i+is);
      x0(d,1) = densvar(dof, k+ks, j+js-1, i+is);
    }
    if (d == 2)
    {
      x0(d,0) = densvar(dof, k+ks, j+js, i+is);
      x0(d,1) = densvar(dof, k+ks-1, j+js, i+is);
    }
  }
  phi(xe, x0);
  for (int d=0; d<ndims; d++) { var(d, k+ks, j+js, i+is) = xe(d);}

  }

//phiT


void YAKL_INLINE phiT( SArray<real,1> &ke, SArray<real,ndims,2> const &u, SArray<real,ndims,2> const &v) {
  ke(0) = 0.0;
  for (int d=0; d<ndims; d++)
  {
    ke(0) += v(d,0)*u(d,0) + v(d,1)*u(d,1);
  }
  ke(0) *= 1./2.;
  }

  void YAKL_INLINE compute_phiT(realArr var, realArr uvar, realArr vvar, int is, int js, int ks, int i, int j, int k)
  {
    SArray<real,1> ke;
    SArray<real,ndims,2> u;
    SArray<real,ndims,2> v;
    for(int d=0; d<ndims; d++)
    {
      if (d==0)
      {
      u(d,0) = uvar(d,k+ks,j+js,i+is);
      u(d,1) = uvar(d,k+ks,j+js,i+is+1);
      v(d,0) = vvar(d,k+ks,j+js,i+is);
      v(d,1) = vvar(d,k+ks,j+js,i+is+1);
      }
      if (d==1)
      {
      u(d,0) = uvar(d,k+ks,j+js,i+is);
      u(d,1) = uvar(d,k+ks,j+js+1,i+is);
      v(d,0) = vvar(d,k+ks,j+js,i+is);
      v(d,1) = vvar(d,k+ks,j+js+1,i+is);
      }
      if (d==2)
      {
      u(d,0) = uvar(d,k+ks,j+js,i+is);
      u(d,1) = uvar(d,k+ks+1,j+js,i+is);
      v(d,0) = vvar(d,k+ks,j+js,i+is);
      v(d,1) = vvar(d,k+ks+1,j+js,i+is);
      }
    }
    phiT(ke, u, v);
    var(0, k+ks, j+js, i+is) = ke(0);
  }

  #endif
