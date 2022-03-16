#ifndef _HAMILTONIAN_KE_H_
#define _HAMILTONIAN_KE_H_

#include "common.h"

// // This kinetic energy functional assumes that dens(0) holds the total density
// class Hamiltonian_Hk1D_rho {
// 
// public:
//   Geometry<ndims,1,1,1> *primal_geometry;
//   Geometry<ndims,1,1,1> *dual_geometry;
//   bool is_initialized;
// 
//    Hamiltonian_Hk1D_rho() {
//      this->is_initialized = false;
// }
// 
// void initialize(ModelParameters &params, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom)
// {
//   this->primal_geometry = &primal_geom;
//   this->dual_geometry = &dual_geom;
//   this->is_initialized = true;
// }
// 
// real YAKL_INLINE compute_KE(const realArr v, const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
// {
// 
//   real KE, PE, IE;
//   SArray<real,1> h0, h0im1;
//   SArray<real,1> U, he;
//   SArray<real,1,2> h0arr;
// 
//   //compute U = H v
//   compute_H<1,diff_ord> (U, v, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
// 
//   // Compute h0 = I h needed for phi calcs
//   compute_I<1,diff_ord> (h0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
//   compute_I<1,diff_ord> (h0im1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i-1, j, k);
// 
// //IS THIS CORRECT?
//   //compute he = phi h0
//   h0arr(0,0) = h0(0);
//   h0arr(0,1) = h0im1(0);
//   phi(he, h0arr);
// 
// //IS THIS CORRECT?
//   return 1./2. * (he(0) * ( U(0) * v(0,k+ks,j+js,i+is)));
// 
// }
// 
// 
//  void YAKL_INLINE compute_dKdv(realArr F, realArr K, realArr HE, const realArr v, const realArr U, const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k)
// {
//   SArray<real,1,2> D0;
//   SArray<real,1> he;
// 
// //IS THIS CORRECT?
//   //compute he = phi * h0
//   D0(0,0) = dens0(0, k+ks, j+js, i+is);
//   D0(0,1) = dens0(0, k+ks, j+js, i+is-1);
//   phi(he, D0);
//   HE(0, k+ks, j+js, i+is) = he(0);
// 
//   //compute F = he * U
//   F(0, k+ks, j+js, i+is) = U(0, k+ks, j+js, i+is) * he(0);
// 
//   //compute K = 1/2 * PhiT(U,V)
//   compute_phiT(K, U, v, is, js, ks, i, j, k);
//   K(0, k+ks, j+js, i+is) *= 0.5;
// 
// }
// 
// // Note that this ADDS to Bvar...
// void YAKL_INLINE compute_dKddens(realArr B, realArr Bfct, const realArr K, int is, int js, int ks, int i, int j, int k)
// {
//   SArray<real,1> K0;
//   compute_I<1, diff_ord>(K0, K, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
//   B(0, k+ks, j+js, i+is) += K0(0);
// }
// 
//  };
 
 
 
// This kinetic energy functional assumes that dens(0) holds the total density
// GENERALIZE WHEN VARIABLE SETS INTRODUCED
class Hamiltonian_Hk {
  
public:
  Geometry<1,1,1> *primal_geometry;
  Geometry<1,1,1> *dual_geometry;
  bool is_initialized;

   Hamiltonian_Hk() {
     this->is_initialized = false;
}

void initialize(ModelParameters &params, Geometry<1,1,1> &primal_geom, Geometry<1,1,1> &dual_geom)
{
  this->primal_geometry = &primal_geom;
  this->dual_geometry = &dual_geom;
  this->is_initialized = true;
}

real YAKL_INLINE compute_KE(const realArr v, const realArr dens, int is, int js, int ks, int i, int j, int k)
{
  
  real KE;
  SArray<real,1> h0, h0im1, h0jm1, h0km1;
  SArray<real,ndims> U, he;
  SArray<real,ndims,2> h0arr;

  //compute U = H v
  compute_H<1,diff_ord> (U, v, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

  // Compute h0 = I h needed for phi calcs
  compute_I<1,diff_ord> (h0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
  compute_I<1,diff_ord> (h0im1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i-1, j, k);
  if (ndims>=2) {compute_I<1,diff_ord> (h0jm1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j-1, k);}
  if (ndims>=3) {compute_I<1,diff_ord> (h0km1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k-1);}

  //compute he = phi h0
  for (int d=0; d<ndims; d++)
  {
    if (d == 0)
    {
      h0arr(d,0) = h0(0);
      h0arr(d,1) = h0im1(0);
    }
    if (d == 1)
    {
      h0arr(d,0) = h0(0);
      h0arr(d,1) = h0jm1(0);
    }
    if (d == 2)
    {
      h0arr(d,0) = h0(0);
      h0arr(d,1) = h0km1(0);
    }
  }
  phi(he, h0arr);
  
  KE = 0.;
  for (int d=0; d<ndims; d++)
    {
      KE = KE + he(d) * (U(d) * v(d,k+ks,j+js,i+is));
    }
  return 1./2. * KE;

}


 void YAKL_INLINE compute_dKdv(realArr F, realArr K, realArr HE, const realArr v, const realArr U, const realArr dens0, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,ndims,2> D0;
  SArray<real,ndims> he;

  //compute he = phi * h0
  for (int d=0; d<ndims; d++)
  {
    if (d == 0)
    {
      D0(d,0) = dens0(0, k+ks, j+js, i+is);
      D0(d,1) = dens0(0, k+ks, j+js, i+is-1);
    }
    if (d == 1)
    {
      D0(d,0) = dens0(0, k+ks, j+js, i+is);
      D0(d,1) = dens0(0, k+ks, j+js-1, i+is);
    }
    if (d == 2)
    {
      D0(d,0) = dens0(0, k+ks, j+js, i+is);
      D0(d,1) = dens0(0, k+ks-1, j+js, i+is);
    }
  }
  phi(he, D0);

  //compute F = he * U, set HE
  for (int d=0; d<ndims; d++)
  {
    HE(d, k+ks, j+js, i+is) = he(d);
    F(d, k+ks, j+js, i+is) = U(d, k+ks, j+js, i+is) * he(d);
  }
  
  //compute K = 1/2 * PhiT(U,V)
  compute_phiT(K, U, v, is, js, ks, i, j, k);
  K(0, k+ks, j+js, i+is) *= 0.5;
  
}

// Note that this ADDS to Bvar...
void YAKL_INLINE compute_dKddens(realArr B, const realArr K, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1> K0;
  compute_I<1, diff_ord>(K0, K, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
  B(0, k+ks, j+js, i+is) += K0(0);
}

 };
 
 
 
 
 
 // This kinetic energy functional assumes that dens(0) + densfct(0:2) make the total density
 class Hamiltonian_Hk_rhod {
   
 public:
   Geometry<1,1,1> *primal_geometry;
   Geometry<1,1,1> *dual_geometry;
   bool is_initialized;

    Hamiltonian_Hk_rhod() {
      this->is_initialized = false;
 }

 void initialize(ModelParameters &params, Geometry<1,1,1> &primal_geom, Geometry<1,1,1> &dual_geom)
 {
   this->primal_geometry = &primal_geom;
   this->dual_geometry = &dual_geom;
   this->is_initialized = true;
 }

 real YAKL_INLINE compute_KE(const realArr v, const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
 {
   

   real KE;
   SArray<real,1> rhod0, rhod0im1, rhod0jm1, rhod0km1;
   SArray<real,3> rhos0, rhos0im1, rhos0jm1, rhos0km1;
   SArray<real,ndims> U, he;
   SArray<real,ndims,2> h0arr;

   //compute U = H v
   compute_H<1,diff_ord> (U, v, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

   // Compute h0 = I h needed for phi calcs
   compute_I<1,diff_ord> (rhod0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
   compute_I<1,diff_ord> (rhod0im1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i-1, j, k);
   if (ndims>=2) {compute_I<1,diff_ord> (rhod0jm1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j-1, k);}
   if (ndims>=3) {compute_I<1,diff_ord> (rhod0km1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k-1);}
   compute_I<3,diff_ord> (rhos0im1, densfct, *this->primal_geometry, *this->dual_geometry, is, js, ks, i-1, j, k);
   compute_I<3,diff_ord> (rhos0im1, densfct, *this->primal_geometry, *this->dual_geometry, is, js, ks, i-1, j, k);
   if (ndims>=2) {compute_I<3,diff_ord> (rhos0jm1, densfct, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j-1, k);}
   if (ndims>=3) {compute_I<3,diff_ord> (rhos0km1, densfct, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k-1);}
   
   //compute he = phi h0
   // exploits linearity of I/Phi
   for (int d=0; d<ndims; d++)
   {
     if (d == 0)
     {
       h0arr(d,0) = rhod0(0) + rhos0(0) + rhos0(1) + rhos0(2);
       h0arr(d,1) = rhod0im1(0) + rhos0im1(0) + rhos0im1(1) + rhos0im1(2);
     }
     if (d == 1)
     {
       h0arr(d,0) = rhod0(0) + rhos0(0) + rhos0(1) + rhos0(2);
       h0arr(d,1) = rhod0jm1(0) + rhos0jm1(0) + rhos0jm1(1) + rhos0jm1(2);
     }
     if (d == 2)
     {
       h0arr(d,0) = rhod0(0) + rhos0(0) + rhos0(1) + rhos0(2);
       h0arr(d,1) = rhod0km1(0) + rhos0km1(0) + rhos0km1(1) + rhos0km1(2);
     }
   }
   phi(he, h0arr);
   
   KE = 0.;
   for (int d=0; d<ndims; d++)
     {
       KE = KE + he(d) * (U(d) * v(d,k+ks,j+js,i+is));
     }
   return 1./2. * KE;

 }


  void YAKL_INLINE compute_dKdv(realArr F, realArr K, realArr HE, const realArr v, const realArr U, const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k)
 {
   
   SArray<real,ndims,2> D0;
   SArray<real,ndims> he;

   //compute he = phi * h0
   // exploits linearity of I/phi
   for (int d=0; d<ndims; d++)
   {
     if (d == 0)
     {
       D0(d,0) = dens0(0, k+ks, j+js, i+is) + densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js, i+is);
       D0(d,1) = dens0(0, k+ks, j+js, i+is-1) + densfct0(0, k+ks, j+js, i+is-1) + densfct0(1, k+ks, j+js, i+is-1) + densfct0(2, k+ks, j+js, i+is-1);
     }
     if (d == 1)
     {
       D0(d,0) = dens0(0, k+ks, j+js, i+is) + densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js, i+is);
       D0(d,1) = dens0(0, k+ks, j+js-1, i+is) + densfct0(0, k+ks, j+js-1, i+is) + densfct0(1, k+ks, j+js-1, i+is) + densfct0(2, k+ks, j+js-1, i+is);
     }
     if (d == 2)
     {
       D0(d,0) = dens0(0, k+ks, j+js, i+is) + densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js, i+is);
       D0(d,1) = dens0(0, k+ks-1, j+js, i+is) + densfct0(0, k+ks-1, j+js, i+is) + densfct0(1, k+ks-1, j+js, i+is) + densfct0(2, k+ks-1, j+js, i+is);
     }
   }
   phi(he, D0);

   //compute F = he * U, set HE
   for (int d=0; d<ndims; d++)
   {
     HE(d, k+ks, j+js, i+is) = he(d);
     F(d, k+ks, j+js, i+is) = U(d, k+ks, j+js, i+is) * he(d);
   }
   
   //compute K = 1/2 * PhiT(U,V)
   compute_phiT(K, U, v, is, js, ks, i, j, k);
   K(0, k+ks, j+js, i+is) *= 0.5;
   
 }

 // Note that this ADDS to Bvar/Bfctvar...
 void YAKL_INLINE compute_dKddens(realArr B, realArr Bfct, const realArr K, int is, int js, int ks, int i, int j, int k)
 {
   SArray<real,1> K0;
   compute_I<1, diff_ord>(K0, K, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
   B(0, k+ks, j+js, i+is) += K0(0);
   Bfct(0, k+ks, j+js, i+is) += K0(0);
   Bfct(1, k+ks, j+js, i+is) += K0(0);
   Bfct(2, k+ks, j+js, i+is) += K0(0);
 }

  };









  // This kinetic energy functional assumes density is the sum of a compile time list of indices of DENS
  // EVENTUALLY, RIGHT NOW IT JUST USES DENS(0)...
  class Hamiltonian_Hk_extruded {
    
  public:
    Geometry<1,1,1> *primal_geometry;
    Geometry<1,1,1> *dual_geometry;
    bool is_initialized;

     Hamiltonian_Hk_extruded() {
       this->is_initialized = false;
  }

  void initialize(ModelParameters &params, Geometry<1,1,1> &primal_geom, Geometry<1,1,1> &dual_geom)
  {
    this->primal_geometry = &primal_geom;
    this->dual_geometry = &dual_geom;
    this->is_initialized = true;
  }

  real YAKL_INLINE compute_KE(const realArr v, const realArr w, const realArr dens, int is, int js, int ks, int i, int j, int k)
  {
    
//     SArray<real,1> K2;
//     SArray<real,ndims,1> U0, Uup, Uright;
//     SArray<real,1> UW0, UW1;
//     SArray<real,1> h0;
//     SArray<real,ndims,2> U;
//     SArray<real,2> UW;
// 
//     //compute U = H v, UW = Hv w
//     compute_Hext<1,diff_ord> (U0, v, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
//     compute_Hext<1,diff_ord> (Uright, v, *this->primal_geometry, *this->dual_geometry, is, js, ks, i+1, j, k);
//     //Have to subtract 1 from k here since UW has an extra dof compared to w!
// //THINK ABOUT THIS!
//     compute_Hv<1,vert_diff_ord> (UW0, w, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k-1);
//     compute_Hv<1,vert_diff_ord> (UW1, w, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
//     U(0,0) = U0(0);
//     U(0,1) = Uright(0);
//     UW(0) = UW0(0);
//     UW(1) = UW1(0);
//     if (ndims>=2) {
//       compute_Hext<1,diff_ord> (Uup, v, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j+1, k);
//     U(1,0) = U0(1);
//     U(1,1) = Uup(1);    
//   }
//     //compute K = 1/2 * PhiT(U,V) + 1/2 * PhiTW(UW,W)
//     compute_phiT(K2, U, v, is, js, ks, i, j, k);
//     compute_phiTW<ADD_MODE::ADD>(K2, UW, w, is, js, ks, i, j, k);
//     K2(0) *= 0.5;

      real K2 = 0.;

      // Have to subtract 1 from k here since UW has an extra dof compared to w!
      SArray<real,1> UW0, UW1;
      compute_Hv<1,vert_diff_ord> (UW0, w, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k-1);
      compute_Hv<1,vert_diff_ord> (UW1, w, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
      real w0, w1;
      //Have to subtract 1 from k here since UW has an extra dof compared to w
      w0 = w(0,k+ks-1,j+js,i+is);
      w1 = w(0,k+ks,j+js,i+is);
      K2 += 0.5 * (w0*UW0(0) + w1*UW1(0));

      real v0, v1;
      SArray<real,ndims> U0, U1;
      compute_Hext<1,diff_ord> (U0, v, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
      compute_Hext<1,diff_ord> (U1, v, *this->primal_geometry, *this->dual_geometry, is, js, ks, i+1, j, k);
      v0 = v(0,k+ks,j+js,i+is);
      v1 = v(0,k+ks,j+js,i+is+1);
      K2 += 0.5 * (v0*U0(0) + v1*U1(0));
      
      if (ndims==2) {
      real v0, v1;
      SArray<real,ndims> U0, U1;
      compute_Hext<1,diff_ord> (U0, v, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
      compute_Hext<1,diff_ord> (U1, v, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j+1, k);
      v0 = v(1,k+ks,j+js,i+is);
      v1 = v(1,k+ks,j+js+1,i+is);
      K2 += 0.5 * (v0*U0(1) + v1*U1(1));
      }
      
      K2 *= 0.5;
      
      // Compute h0 = I h needed for phi calcs
      SArray<real,1> h0;
      compute_Iext<1,diff_ord,vert_diff_ord> (h0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

      return h0(0) * K2;

  }


  void YAKL_INLINE compute_Fw(realArr FW, realArr HEw, const realArr UW, const realArr dens0, int is, int js, int ks, int i, int j, int k)
{
  //SArray<real,2> Dv;
  //compute hew = phiw * h0
   //Dv(0) = dens0(0, k+ks, j+js, i+is);
   //Dv(1) = dens0(0, k+ks-1, j+js, i+is);
   //real hew = phiW(Dv);
  //compute FW = hew * UW, set HEw
  real hew = (dens0(0, k+ks, j+js, i+is) + dens0(0, k+ks-1, j+js, i+is))/2.0;
  HEw(0, k+ks, j+js, i+is) = hew;
  FW(0, k+ks, j+js, i+is) = UW(0, k+ks, j+js, i+is) * hew;
  //std::cout << "HEw in Hk " << i << " " << j << " " << k << " " << HEw(0,k+ks,j+js,i+is) << "\n" << std::flush;
}

  void YAKL_INLINE compute_K(realArr K,  const realArr v, const realArr U, const realArr w, const realArr UW, int is, int js, int ks, int i, int j, int k)
{
  //compute K = 1/2 * PhiT(U,V) + 1/2 * PhiTW(UW,W)
  real K2 = 0.;
  
  real w0, w1, UW0, UW1;
  UW0 = UW(0,k+ks,j+js,i+is);
  UW1 = UW(0,k+ks+1,j+js,i+is);
  //Have to subtract 1 from k here since UW has an extra dof compared to w
  w0 = w(0,k+ks-1,j+js,i+is);
  w1 = w(0,k+ks,j+js,i+is);
  K2 += 0.5 * (w0*UW0 + w1*UW1);
  
  real v0, U0, v1, U1;
  U0 = U(0,k+ks,j+js,i+is);
  U1 = U(0,k+ks,j+js,i+is+1);
  v0 = v(0,k+ks,j+js,i+is);
  v1 = v(0,k+ks,j+js,i+is+1);
  K2 += 0.5 *(v0*U0 + v1*U1);
  
  if (ndims==2)
    {
  real v0, U0, v1, U1;
  U0 = U(1,k+ks,j+js,i+is);
  U1 = U(1,k+ks,j+js+1,i+is);
  v0 = v(1,k+ks,j+js,i+is);
  v1 = v(1,k+ks,j+js+1,i+is);
  K2 += 0.5 *(v0*U0 + v1*U1);
  }
  
  K(0, k+ks, j+js, i+is) = 0.5 * K2;
  //compute_phiT(K, U, v, is, js, ks, i, j, k);
  //compute_phiTW<ADD_MODE::ADD>(K, UW, w, is, js, ks, i, j, k);
  //K(0, k+ks, j+js, i+is) *= 0.5;
}

   void YAKL_INLINE compute_F(realArr F, realArr HE, const realArr U, const realArr dens0, int is, int js, int ks, int i, int j, int k)
  {
    SArray<real,ndims,2> D0;
    SArray<real,ndims> he;

    //compute he = phi * h0
    for (int d=0; d<ndims; d++)
    {
      if (d == 0)
      {
        D0(d,0) = dens0(0, k+ks, j+js, i+is);
        D0(d,1) = dens0(0, k+ks, j+js, i+is-1);
      }
      if (d == 1)
      {
        D0(d,0) = dens0(0, k+ks, j+js, i+is);
        D0(d,1) = dens0(0, k+ks, j+js-1, i+is);
      }
    }
    phi(he, D0);
    //compute F = he * U, set HE
    for (int d=0; d<ndims; d++)
    {
      HE(d, k+ks, j+js, i+is) = he(d);
      F(d, k+ks, j+js, i+is) = U(d, k+ks, j+js, i+is) * he(d);
    }
    //std::cout << "HE in Hk " << i << " " << j << " " << k << " " << HE(0,k+ks,j+js,i+is) << "\n" << std::flush;
  }

  template<ADD_MODE addmode=ADD_MODE::REPLACE> void YAKL_INLINE compute_dKddens(realArr B, const realArr K, int is, int js, int ks, int i, int j, int k)
  {
    SArray<real,1> K0;
    compute_Iext<1, diff_ord, vert_diff_ord>(K0, K, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
    if (addmode == ADD_MODE::REPLACE) {B(0, k+ks, j+js, i+is) = K0(0);}
    if (addmode == ADD_MODE::ADD) {B(0, k+ks, j+js, i+is) += K0(0);}
  }

   };
#endif
