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
class Hamiltonian_Hk_rho {
  
public:
  Geometry<ndims,1,1,1> *primal_geometry;
  Geometry<ndims,1,1,1> *dual_geometry;
  bool is_initialized;

   Hamiltonian_Hk_rho() {
     this->is_initialized = false;
}

void initialize(ModelParameters &params, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom)
{
  this->primal_geometry = &primal_geom;
  this->dual_geometry = &dual_geom;
  this->is_initialized = true;
}

real YAKL_INLINE compute_KE(const realArr v, const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
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


 void YAKL_INLINE compute_dKdv(realArr F, realArr K, realArr HE, const realArr v, const realArr U, const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k)
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
void YAKL_INLINE compute_dKddens(realArr B, realArr Bfct, const realArr K, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1> K0;
  compute_I<1, diff_ord>(K0, K, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
  B(0, k+ks, j+js, i+is) += K0(0);
}

 };
 
 
 
 
 
 // This kinetic energy functional assumes that dens(0) + densfct(0:2) make the total density
 class Hamiltonian_Hk_rhod {
   
 public:
   Geometry<ndims,1,1,1> *primal_geometry;
   Geometry<ndims,1,1,1> *dual_geometry;
   bool is_initialized;

    Hamiltonian_Hk_rhod() {
      this->is_initialized = false;
 }

 void initialize(ModelParameters &params, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom)
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
  
#endif
