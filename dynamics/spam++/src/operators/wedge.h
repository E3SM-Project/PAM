

#ifndef _WEDGE_H_
#define _WEDGE_H_

#include "common.h"
#include <cmath>
#include <iostream>

// Q
  template<uint ndofs> void YAKL_INLINE Q2D( SArray<real,2,ndofs,2> &vel, SArray<real,3,ndofs,2,5> const &recon, SArray<real,2,2,4> const &flux) {

      for (int l=0; l<ndofs; l++) {
      //x-dir
      vel(l,0) = -1./8. * (flux(0,0)*recon(l,0,0) + flux(0,1)*recon(l,0,1) + flux(0,2)*recon(l,0,2) + flux(0,3)*recon(l,0,3) +
                           flux(0,0)*recon(l,0,4) + flux(0,1)*recon(l,0,4) + flux(0,2)*recon(l,0,4) + flux(0,3)*recon(l,0,4));

      //y-dir
      vel(l,1) =  1./8. * (flux(1,0)*recon(l,1,0) + flux(1,1)*recon(l,1,1) + flux(1,2)*recon(l,1,2) + flux(1,3)*recon(l,1,3) +
                           flux(1,0)*recon(l,1,4) + flux(1,1)*recon(l,1,4) + flux(1,2)*recon(l,1,4) + flux(1,3)*recon(l,1,4));
  }
}


template<uint ndofs> void YAKL_INLINE Q2D_nonEC( SArray<real,2,ndofs,2> &vel, SArray<real,2,ndofs,2> const &recon, SArray<real,2,2,4> const &flux) {
  for (int l=0; l<ndofs; l++) {
  //x-dir
  vel(l,0) = -1./4. * recon(l,0) * (flux(0,0) + flux(0,1) + flux(0,2) + flux(0,3));
  //y-dir
  vel(l,1) = 1./4. * recon(l,1) * (flux(1,0) + flux(1,1) + flux(1,2) + flux(1,3));
}

}


template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Q_EC(real4d qflux, const real4d reconvar, const real4d Uvar, int is, int js, int ks, int i, int j, int k)
{

SArray<real,2,ndofs,2> vel;
SArray<real,2,2,4> flux;
flux(0,0) = Uvar(1, k+ks, j+js, i+is);
flux(0,1) = Uvar(1, k+ks, j+js, i+is-1);
flux(0,2) = Uvar(1, k+ks, j+js+1, i+is);
flux(0,3) = Uvar(1, k+ks, j+js+1, i+is-1);

flux(1,0) = Uvar(0, k+ks, j+js, i+is);
flux(1,1) = Uvar(0, k+ks, j+js, i+is+1);
flux(1,2) = Uvar(0, k+ks, j+js-1, i+is);
flux(1,3) = Uvar(0, k+ks, j+js-1, i+is+1);

   SArray<real,3,ndofs,2,5> recon;
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


template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Q_nonEC(real4d qflux, const real4d reconvar, const real4d Uvar, int is, int js, int ks, int i, int j, int k)
{

  SArray<real,2,ndofs,2> vel;
  SArray<real,2,2,4> flux;
  flux(0,0) = Uvar(1, k+ks, j+js, i+is);
  flux(0,1) = Uvar(1, k+ks, j+js, i+is-1);
  flux(0,2) = Uvar(1, k+ks, j+js+1, i+is);
  flux(0,3) = Uvar(1, k+ks, j+js+1, i+is-1);

  flux(1,0) = Uvar(0, k+ks, j+js, i+is);
  flux(1,1) = Uvar(0, k+ks, j+js, i+is+1);
  flux(1,2) = Uvar(0, k+ks, j+js-1, i+is);
  flux(1,3) = Uvar(0, k+ks, j+js-1, i+is+1);

    SArray<real,2,ndofs,2> recon;
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

template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Qxz_w_EC(real4d qflux, const real4d reconvar, const real4d vertreconvar, const real4d Uvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,4> flux;
  SArray<real,1,4> recon;
  flux(0) = Uvar(0, k+ks, j+js, i+is);
  flux(1) = Uvar(0, k+ks, j+js, i+is+1);
  flux(2) = Uvar(0, k+ks+1, j+js, i+is);
  flux(3) = Uvar(0, k+ks+1, j+js, i+is+1);
  for (int l=0; l<ndofs; l++) {
  recon(0) = (vertreconvar(l, k+ks, j+js, i+is) + reconvar(l, k+ks, j+js, i+is)) /2.;
  recon(1) = (vertreconvar(l, k+ks, j+js, i+is+1)+ reconvar(l, k+ks, j+js, i+is)) /2.;
  recon(2) = (vertreconvar(l, k+ks+1, j+js, i+is)+ reconvar(l, k+ks, j+js, i+is)) /2.;
  recon(3) = (vertreconvar(l, k+ks+1, j+js, i+is+1)+ reconvar(l, k+ks, j+js, i+is)) /2.;
  if (addmode == ADD_MODE::REPLACE){ qflux(l, k+ks, j+js, i+is) = 1./4. * (flux(0)*recon(0) + flux(1)*recon(1) + flux(2)*recon(2) + flux(3)*recon(3)); }
  if (addmode == ADD_MODE::ADD){ qflux(l, k+ks, j+js, i+is) += 1./4. * (flux(0)*recon(0) + flux(1)*recon(1) + flux(2)*recon(2) + flux(3)*recon(3)); }
  }
}
template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Qxz_w_nonEC(real4d qflux, const real4d reconvar, const real4d Uvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,4> flux;
  flux(0) = Uvar(0, k+ks, j+js, i+is);
  flux(1) = Uvar(0, k+ks, j+js, i+is+1);
  flux(2) = Uvar(0, k+ks+1, j+js, i+is);
  flux(3) = Uvar(0, k+ks+1, j+js, i+is+1);
  if (addmode == ADD_MODE::REPLACE) {
  for (int l=0; l<ndofs; l++) { qflux(l, k+ks, j+js, i+is) = 1./4. * (flux(0) + flux(1) + flux(2) + flux(3)) * reconvar(l, k+ks, j+js, i+is); }}
  if (addmode == ADD_MODE::ADD){
  for (int l=0; l<ndofs; l++) { qflux(l, k+ks, j+js, i+is) += 1./4. * (flux(0) + flux(1) + flux(2) + flux(3)) * reconvar(l, k+ks, j+js, i+is); }}
}
template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Qxz_w_EC_top(real4d qflux, const real4d reconvar, const real4d vertreconvar, const real4d Uvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,2> flux;
  SArray<real,1,2> recon;
  flux(0) = Uvar(0, k+ks, j+js, i+is);
  flux(1) = Uvar(0, k+ks, j+js, i+is+1);
  for (int l=0; l<ndofs; l++) {
  recon(0) = (vertreconvar(l, k+ks, j+js, i+is) + reconvar(l, k+ks, j+js, i+is)) /2.;
  recon(1) = (vertreconvar(l, k+ks, j+js, i+is+1)+ reconvar(l, k+ks, j+js, i+is)) /2.;
  if (addmode == ADD_MODE::REPLACE){ qflux(l, k+ks, j+js, i+is) = 1./4. * (flux(0)*recon(0) + flux(1)*recon(1)); }
  if (addmode == ADD_MODE::ADD){ qflux(l, k+ks, j+js, i+is) += 1./4. * (flux(0)*recon(0) + flux(1)*recon(1)); }
  }
}
template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Qxz_w_nonEC_top(real4d qflux, const real4d reconvar, const real4d Uvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,2> flux;
  flux(0) = Uvar(0, k+ks, j+js, i+is);
  flux(1) = Uvar(0, k+ks, j+js, i+is+1);
  if (addmode == ADD_MODE::REPLACE) {
  for (int l=0; l<ndofs; l++) { qflux(l, k+ks, j+js, i+is) = 1./4. * (flux(0) + flux(1)) * reconvar(l, k+ks, j+js, i+is); }}
  if (addmode == ADD_MODE::ADD){
  for (int l=0; l<ndofs; l++) { qflux(l, k+ks, j+js, i+is) += 1./4. * (flux(0) + flux(1)) * reconvar(l, k+ks, j+js, i+is); }}
}
template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Qxz_w_EC_bottom(real4d qflux, const real4d reconvar, const real4d vertreconvar, const real4d Uvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,2> flux;
  SArray<real,1,2> recon;
  flux(0) = Uvar(0, k+ks+1, j+js, i+is);
  flux(1) = Uvar(0, k+ks+1, j+js, i+is+1);
  for (int l=0; l<ndofs; l++) {
  recon(0) = (vertreconvar(l, k+ks+1, j+js, i+is)+ reconvar(l, k+ks, j+js, i+is)) /2.;
  recon(1) = (vertreconvar(l, k+ks+1, j+js, i+is+1)+ reconvar(l, k+ks, j+js, i+is)) /2.;
  if (addmode == ADD_MODE::REPLACE){ qflux(l, k+ks, j+js, i+is) = 1./4. * (flux(0)*recon(0) + flux(1)*recon(1)); }
  if (addmode == ADD_MODE::ADD){ qflux(l, k+ks, j+js, i+is) += 1./4. * (flux(0)*recon(0) + flux(1)*recon(1)); }
  }
}
template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Qxz_w_nonEC_bottom(real4d qflux, const real4d reconvar, const real4d Uvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,2> flux;
  flux(0) = Uvar(0, k+ks+1, j+js, i+is);
  flux(1) = Uvar(0, k+ks+1, j+js, i+is+1);
  if (addmode == ADD_MODE::REPLACE) {
  for (int l=0; l<ndofs; l++) { qflux(l, k+ks, j+js, i+is) = 1./4. * (flux(0) + flux(1)) * reconvar(l, k+ks, j+js, i+is); }}
  if (addmode == ADD_MODE::ADD){
  for (int l=0; l<ndofs; l++) { qflux(l, k+ks, j+js, i+is) += 1./4. * (flux(0) + flux(1)) * reconvar(l, k+ks, j+js, i+is); }}
}





template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Qxz_u_EC(real4d qvertflux, const real4d reconvar, const real4d vertreconvar, const real4d UWvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,4> flux;
  SArray<real,1,4> recon;
  flux(0) = UWvar(0, k+ks, j+js, i+is);
  flux(1) = UWvar(0, k+ks, j+js, i+is-1);
  flux(2) = UWvar(0, k+ks+1, j+js, i+is);
  flux(3) = UWvar(0, k+ks+1, j+js, i+is-1);
  for (int l=0; l<ndofs; l++) {
  //Have to subtract 1 in k here because UW has an extra dof at boundary compared to v!
  recon(0) = (reconvar(l, k+ks-1, j+js, i+is) + vertreconvar(l, k+ks, j+js, i+is)) /2.;
  recon(1) = (reconvar(l, k+ks-1, j+js, i+is-1) + vertreconvar(l, k+ks, j+js, i+is)) /2.;
  recon(2) = (reconvar(l, k+ks, j+js, i+is) + vertreconvar(l, k+ks, j+js, i+is)) /2.;
  recon(3) = (reconvar(l, k+ks, j+js, i+is-1)+ vertreconvar(l, k+ks, j+js, i+is)) /2.;
  //Added the minus sign here
  if (addmode == ADD_MODE::REPLACE){qvertflux(l, k+ks, j+js, i+is) = -1./4. * (flux(0)*recon(0) + flux(1)*recon(1) + flux(2)*recon(2) + flux(3)*recon(3)); }
  if (addmode == ADD_MODE::ADD){qvertflux(l, k+ks, j+js, i+is) += -1./4. * (flux(0)*recon(0) + flux(1)*recon(1) + flux(2)*recon(2) + flux(3)*recon(3)); }
}
}
template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Qxz_u_nonEC(real4d qvertflux, const real4d vertreconvar, const real4d UWvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,4> flux;
  flux(0) = UWvar(0, k+ks, j+js, i+is);
  flux(1) = UWvar(0, k+ks, j+js, i+is-1);
  flux(2) = UWvar(0, k+ks+1, j+js, i+is);
  flux(3) = UWvar(0, k+ks+1, j+js, i+is-1);
  //Added the minus sign here
  if (addmode == ADD_MODE::REPLACE){
  for (int l=0; l<ndofs; l++) { qvertflux(l, k+ks, j+js, i+is) = -1./4. * (flux(0) + flux(1) + flux(2) + flux(3)) * vertreconvar(l, k+ks, j+js, i+is); }}
  if (addmode == ADD_MODE::ADD){
  for (int l=0; l<ndofs; l++) { qvertflux(l, k+ks, j+js, i+is) += -1./4. * (flux(0) + flux(1) + flux(2) + flux(3)) * vertreconvar(l, k+ks, j+js, i+is); }}
}


template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Qxz_u_top(real4d qvertflux, const real4d vertreconvar, const real4d UWvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,2> flux;
  flux(0) = UWvar(0, k+ks+1, j+js, i+is);
  flux(1) = UWvar(0, k+ks+1, j+js, i+is-1);
  //Added the minus sign here
  if (addmode == ADD_MODE::REPLACE){
  for (int l=0; l<ndofs; l++) { qvertflux(l, k+ks, j+js, i+is) = -1./2. * (flux(0) + flux(1)) * vertreconvar(l, k+ks, j+js, i+is); }}
  if (addmode == ADD_MODE::ADD){
  for (int l=0; l<ndofs; l++) { qvertflux(l, k+ks, j+js, i+is) += -1./2. * (flux(0) + flux(1)) * vertreconvar(l, k+ks, j+js, i+is); }}
}
template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Qxz_u_bottom(real4d qvertflux, const real4d vertreconvar, const real4d UWvar, int is, int js, int ks, int i, int j, int k)
{
SArray<real,1,2> flux;
flux(0) = UWvar(0, k+ks, j+js, i+is);
flux(1) = UWvar(0, k+ks, j+js, i+is-1);
//Added the minus sign here
if (addmode == ADD_MODE::REPLACE){
for (int l=0; l<ndofs; l++) { qvertflux(l, k+ks, j+js, i+is) = -1./2. * (flux(0) + flux(1)) * vertreconvar(l, k+ks, j+js, i+is); }}
if (addmode == ADD_MODE::ADD){
for (int l=0; l<ndofs; l++) { qvertflux(l, k+ks, j+js, i+is) += -1./2. * (flux(0) + flux(1)) * vertreconvar(l, k+ks, j+js, i+is); }}
}
template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Qxz_u_nonEC_top(real4d qvertflux, const real4d vertreconvar, const real4d UWvar, int is, int js, int ks, int i, int j, int k)
{
  compute_Qxz_u_top<ndofs, addmode>(qvertflux, vertreconvar, UWvar, is, js, ks, i, j, k);
}
template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Qxz_u_EC_top(real4d qvertflux, const real4d reconvar, const real4d vertreconvar, const real4d UWvar, int is, int js, int ks, int i, int j, int k)
{
 compute_Qxz_u_top<ndofs, addmode>(qvertflux, vertreconvar, UWvar, is, js, ks, i, j, k);
}
template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Qxz_u_nonEC_bottom(real4d qvertflux, const real4d vertreconvar, const real4d UWvar, int is, int js, int ks, int i, int j, int k)
{
  compute_Qxz_u_bottom<ndofs, addmode>(qvertflux, vertreconvar, UWvar, is, js, ks, i, j, k);
}
template<uint ndofs, ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_Qxz_u_EC_bottom(real4d qvertflux, const real4d reconvar, const real4d vertreconvar, const real4d UWvar, int is, int js, int ks, int i, int j, int k)
{
 compute_Qxz_u_bottom<ndofs, addmode>(qvertflux, vertreconvar, UWvar, is, js, ks, i, j, k);
}

// W
void YAKL_INLINE W2D( SArray<real,1,2> &vel, SArray<real,2,2,4> const &flux) {
        //x-dir
        vel(0) = -1./4. * (flux(0,0) + flux(0,1) + flux(0,2) + flux(0,3));
        //y-dir
        vel(1) = 1./4. * (flux(1,0) + flux(1,1) + flux(1,2) + flux(1,3));
  }



void YAKL_INLINE compute_W(real4d UTvar, const real4d Uvar, int is, int js, int ks, int i, int j, int k)
{

SArray<real,1,2> ut;
SArray<real,2,2,4> flux;
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

// Wxz
void YAKL_INLINE compute_Wxz_u(real4d VTvar, const real4d UWvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,4> flux;
  flux(0) = UWvar(0, k+ks, j+js, i+is);
  flux(1) = UWvar(0, k+ks, j+js, i+is-1);
  flux(2) = UWvar(0, k+ks+1, j+js, i+is);
  flux(3) = UWvar(0, k+ks+1, j+js, i+is-1);
  //Added the minus sign here
  VTvar(0, k+ks, j+js, i+is) = -1./4. * (flux(0) + flux(1) + flux(2) + flux(3));
}
void YAKL_INLINE compute_Wxz_u_top(real4d VTvar, const real4d UWvar, int is, int js, int ks, int i, int j, int k)
{
SArray<real,1,2> flux;
flux(0) = UWvar(0, k+ks+1, j+js, i+is) ;
flux(1) = UWvar(0, k+ks+1, j+js, i+is-1);
//Added the minus sign here
VTvar(0, k+ks, j+js, i+is) = -1./2. * (flux(0) + flux(1));
}
void YAKL_INLINE compute_Wxz_u_bottom(real4d VTvar, const real4d UWvar, int is, int js, int ks, int i, int j, int k)
{
SArray<real,1,4> flux;
flux(0) = UWvar(0, k+ks, j+js, i+is);
flux(1) = UWvar(0, k+ks, j+js, i+is-1);
//Added the minus sign here
VTvar(0, k+ks, j+js, i+is) = -1./2. * (flux(0) + flux(1));
}

void YAKL_INLINE compute_Wxz_w(real4d WTvar, const real4d Uvar, int is, int js, int ks, int i, int j, int k)
{
SArray<real,1,4> flux;
flux(0) = Uvar(0, k+ks, j+js, i+is);
flux(1) = Uvar(0, k+ks, j+js, i+is+1);
flux(2) = Uvar(0, k+ks+1, j+js, i+is);
flux(3) = Uvar(0, k+ks+1, j+js, i+is+1);
WTvar(0, k+ks, j+js, i+is) = 1./4. * (flux(0) + flux(1) + flux(2) + flux(3));
}
void YAKL_INLINE compute_Wxz_w_top(real4d WTvar, const real4d Uvar, int is, int js, int ks, int i, int j, int k)
{
SArray<real,1,2> flux;
flux(0) = Uvar(0, k+ks, j+js, i+is);
flux(1) = Uvar(0, k+ks, j+js, i+is+1);
WTvar(0, k+ks, j+js, i+is) = 1./4. * (flux(0) + flux(1));
}
void YAKL_INLINE compute_Wxz_w_bottom(real4d WTvar, const real4d Uvar, int is, int js, int ks, int i, int j, int k)
{
SArray<real,1,2> flux;
flux(0) = Uvar(0, k+ks+1, j+js, i+is);
flux(1) = Uvar(0, k+ks+1, j+js, i+is+1);
WTvar(0, k+ks, j+js, i+is) = 1./4. * (flux(0) + flux(1));
}






// R
void YAKL_INLINE R( SArray<real,1,1> &vard, SArray<real,1,2> const &varp) {
  vard(0) = 1./2. * (varp(0) + varp(1));
}
void YAKL_INLINE R( SArray<real,1,1> &vard, SArray<real,1,4> const &varp) {
  vard(0) = 1./4. * (varp(0) + varp(1) + varp(2) + varp(3));
}
void YAKL_INLINE R( SArray<real,1,1> &vard, SArray<real,1,8> const &varp) {
  vard(0) = 1./8. * (varp(0) + varp(1) + varp(2) + varp(3) + varp(4) + varp(5) + varp(6) + varp(7));
}

template<uint dof> void YAKL_INLINE compute_R(real &var, const real4d densvar, int is, int js, int ks, int i, int j, int k)
{
SArray<real,1,1> dualdens;
  if (ndims==1)
  {
    SArray<real,1,2> dens;
    dens(0) = densvar(dof, k+ks, j+js, i+is);
    dens(1) = densvar(dof, k+ks, j+js, i+is-1);
    R(dualdens, dens);
  }
  if (ndims==2)
  {
    SArray<real,1,4> dens;
    dens(0) = densvar(dof, k+ks, j+js, i+is);
    dens(1) = densvar(dof, k+ks, j+js, i+is-1);
    dens(2) = densvar(dof, k+ks, j+js-1, i+is);
    dens(3) = densvar(dof, k+ks, j+js-1, i+is-1);
    R(dualdens, dens);
  }
  if (ndims==3)
  {
    SArray<real,1,8> dens;
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

template<uint dof> void YAKL_INLINE compute_R(real4d var, const real4d densvar, int is, int js, int ks, int i, int j, int k)
{
real var2;
compute_R<dof> (var2,densvar,is,js,ks,i,j,k);
var(dof, k+ks, j+js, i+is) = var2;
}


void YAKL_INLINE Rbnd( SArray<real,1,1> &vard, SArray<real,1,4> const &varp) {
  vard(0) = 1./4. * (varp(0) + varp(1)) + 1./2.*(varp(2) + varp(3));
}












//phi
void YAKL_INLINE phi( SArray<real,1,ndims> &vare, SArray<real,2,ndims,2> const &var0) {
  for(int d=0; d<ndims; d++)
  {
    vare(d) = 1./2. * (var0(d,0) + var0(d,1));
  }
  }

 void YAKL_INLINE compute_phi(real4d var, const real4d densvar, int is, int js, int ks, int i, int j, int k)
  {
  SArray<real,1,ndims> xe;
  SArray<real,2,ndims,2> x0;
  for (int d=0; d<ndims; d++)
  {
    if (d == 0)
    {
      x0(d,0) = densvar(0, k+ks, j+js, i+is);
      x0(d,1) = densvar(0, k+ks, j+js, i+is-1);
    }
    if (d == 1)
    {
      x0(d,0) = densvar(0, k+ks, j+js, i+is);
      x0(d,1) = densvar(0, k+ks, j+js-1, i+is);
    }
    if (d == 2)
    {
      x0(d,0) = densvar(0, k+ks, j+js, i+is);
      x0(d,1) = densvar(0, k+ks-1, j+js, i+is);
    }
  }
  phi(xe, x0);
  for (int d=0; d<ndims; d++) { var(d, k+ks, j+js, i+is) = xe(d);}

  }

  real YAKL_INLINE phiW( SArray<real,1,2> const &var0) 
{
  return 1./2. * (var0(0) + var0(1));
}

void YAKL_INLINE compute_phiW(real4d var, const real4d densvar, int is, int js, int ks, int i, int j, int k)
{
SArray<real,1,2> x0;
x0(0) = densvar(0, k+ks, j+js, i+is);
x0(1) = densvar(0, k+ks-1, j+js, i+is);
var(0, k+ks, j+js, i+is) = phiW(x0);
}





//phiT
void YAKL_INLINE phiT( SArray<real,1,1> &ke, SArray<real,2,ndims,2> const &u, SArray<real,2,ndims,2> const &v) {
  ke(0) = 0.0;
  for (int d=0; d<ndims; d++)
  {
    ke(0) += v(d,0)*u(d,0) + v(d,1)*u(d,1);
  }
  ke(0) *= 1./2.;
  }

  template<ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_phiT(real4d var, const real4d uvar, const real4d vvar, int is, int js, int ks, int i, int j, int k)
  {
    SArray<real,1,1> ke;
    SArray<real,2,ndims,2> u;
    SArray<real,2,ndims,2> v;
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
    if (addmode == ADD_MODE::REPLACE) {var(0, k+ks, j+js, i+is) = ke(0);}
    if (addmode == ADD_MODE::ADD) {var(0, k+ks, j+js, i+is) = ke(0);}
  }
  void YAKL_INLINE compute_phiT(SArray<real,1,1> &var, SArray<real,2,ndims,2> const &u, const real4d vvar, int is, int js, int ks, int i, int j, int k)
  {
    SArray<real,2,ndims,2> v;
    for(int d=0; d<ndims; d++)
    {
      if (d==0)
      {
      v(d,0) = vvar(d,k+ks,j+js,i+is);
      v(d,1) = vvar(d,k+ks,j+js,i+is+1);
      }
      if (d==1)
      {
      v(d,0) = vvar(d,k+ks,j+js,i+is);
      v(d,1) = vvar(d,k+ks,j+js+1,i+is);
      }
      if (d==2)
      {
      v(d,0) = vvar(d,k+ks,j+js,i+is);
      v(d,1) = vvar(d,k+ks+1,j+js,i+is);
      }
    }
    phiT(var, u, v);
  }
  
  real YAKL_INLINE phiTW(SArray<real,1,2> const &u, SArray<real,1,2> const &v) {
      return (v(0)*u(0) + v(1)*u(1))/2.0;
    }

template<ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_phiTW(real4d var, const real4d uwvar, const real4d wvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,2> u;
  SArray<real,1,2> v;
u(0) = uwvar(0,k+ks,j+js,i+is);
u(1) = uwvar(0,k+ks+1,j+js,i+is);
//Have to subtract 1 from k here since UW has an extra dof compared to w
v(0) = wvar(0,k+ks-1,j+js,i+is);
v(1) = wvar(0,k+ks,j+js,i+is);
if (addmode == ADD_MODE::REPLACE) {var(0, k+ks, j+js, i+is) = phiTW(u, v);}
if (addmode == ADD_MODE::ADD) {var(0, k+ks, j+js, i+is) += phiTW(u, v);}
}

template<ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_phiTW(SArray<real,1,1> &var, SArray<real,1,2> const &uw, const real4d wvar, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,2> v;
//SHOULD SUBTRACT HERE ALSO?
v(0) = wvar(0,k+ks-1,j+js,i+is);
v(1) = wvar(0,k+ks,j+js,i+is);
if (addmode == ADD_MODE::REPLACE) {var(0) = phiTW(uw, v);}
if (addmode == ADD_MODE::ADD) {var(0) += phiTW(uw, v);}
}

  #endif
