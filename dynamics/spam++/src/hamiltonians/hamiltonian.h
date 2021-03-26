#ifndef _HAMILTONIAN_H_
#define _HAMILTONIAN_H_

#include "common.h"
#include "thermo.h"



struct pvpe {
  real pv=0., pe=0.;
};

class Functional_PVPE_rho {
public:
  bool is_initialized;

   Functional_PVPE_rho() {
     this->is_initialized = false;
}

void initialize(ModelParameters &params)
{
  this->is_initialized = true;
}

real YAKL_INLINE compute_hv(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1> hv;
  SArray<real,4> Dv;

  // compute hv = R h
  Dv(0) = dens(0, k+ks, j+js, i+is);
  Dv(1) = dens(0, k+ks, j+js, i+is-1);
  Dv(2) = dens(0, k+ks, j+js-1, i+is);
  Dv(3) = dens(0, k+ks, j+js-1, i+is-1);
  R(hv, Dv);

  return hv(0);
}

real YAKL_INLINE compute_zeta(const realArr v, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1> zeta;
  // compute zeta = D2 v
  compute_D2<1>(zeta, v, is, js, ks, i, j, k);
  return zeta(0);
}

real YAKL_INLINE compute_eta(const realArr v, const realArr coriolis, int is, int js, int ks, int i, int j, int k)
{
  real zeta = compute_zeta(v, is, js, ks, i, j, k);
  return zeta + coriolis(0, k+ks, j+js, i+is);
}


void YAKL_INLINE compute_q0f0(realArr q0, realArr f0, const realArr v, const realArr dens, const realArr densfct, const realArr coriolis, int is, int js, int ks, int i, int j, int k)
{

  real hv = compute_hv(dens, densfct, is, js, ks, i, j, k);
  real zeta = compute_zeta(v, is, js, ks, i, j, k);

  // compute q0 = zeta / hv and f0 = f / hv
    q0(0, k+ks, j+js, i+is) = zeta / hv;
    f0(0, k+ks, j+js, i+is) = coriolis(0, k+ks, j+js, i+is) / hv;

}

void YAKL_INLINE compute_q0(realArr q0, const realArr v, const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
{
real hv = compute_hv(dens, densfct, is, js, ks, i, j, k);
real zeta = compute_zeta(v, is, js, ks, i, j, k);
// compute q0 = zeta / hv and f0 = f / hv
  q0(0, k+ks, j+js, i+is) = zeta / hv;
}

pvpe YAKL_INLINE compute_PVPE(const realArr v, const realArr dens, const realArr densfct, const realArr coriolis, int is, int js, int ks, int i, int j, int k)
{
  pvpe vals;
  real eta = compute_eta(v, coriolis, is, js, ks, i, j, k);
  real hv = compute_hv(dens, densfct, is, js, ks, i, j, k);
  real q0 = eta / hv;

  vals.pv = eta;
  vals.pe = 0.5 * eta * q0;
  
  return vals;
}

};



class Functional_PVPE_rhod {
public:
  bool is_initialized;

   Functional_PVPE_rhod() {
     this->is_initialized = false;
}

void initialize(ModelParameters &params)
{
  this->is_initialized = true;
}

real YAKL_INLINE compute_hv(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1> hv;
  SArray<real,4> Dv;

  // compute hv = R h
  // Uses linearity of R
  Dv(0) = dens(0, k+ks, j+js, i+is) + densfct(0, k+ks, j+js, i+is) + densfct(1, k+ks, j+js, i+is) + densfct(2, k+ks, j+js, i+is);
  Dv(1) = dens(0, k+ks, j+js, i+is-1) + densfct(0, k+ks, j+js, i+is-1) + densfct(1, k+ks, j+js, i+is-1) + densfct(2, k+ks, j+js, i+is-1);
  Dv(2) = dens(0, k+ks, j+js-1, i+is) + densfct(0, k+ks, j+js-1, i+is) + densfct(1, k+ks, j+js-1, i+is) + densfct(2, k+ks, j+js-1, i+is);
  Dv(3) = dens(0, k+ks, j+js-1, i+is-1) + densfct(0, k+ks, j+js-1, i+is-1) + densfct(1, k+ks, j+js-1, i+is-1) + densfct(2, k+ks, j+js-1, i+is-1);
  R(hv, Dv);

  return hv(0);
}

real YAKL_INLINE compute_zeta(const realArr v, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1> zeta;
  // compute zeta = D2 v
  compute_D2<1>(zeta, v, is, js, ks, i, j, k);
  return zeta(0);
}

real YAKL_INLINE compute_eta(const realArr v, const realArr coriolis, int is, int js, int ks, int i, int j, int k)
{
  real zeta = compute_zeta(v, is, js, ks, i, j, k);
  return zeta + coriolis(0, k+ks, j+js, i+is);
}


void YAKL_INLINE compute_q0f0(realArr q0, realArr f0, const realArr v, const realArr dens, const realArr densfct, const realArr coriolis, int is, int js, int ks, int i, int j, int k)
{

  real hv = compute_hv(dens, densfct, is, js, ks, i, j, k);
  real zeta = compute_zeta(v, is, js, ks, i, j, k);

  // compute q0 = zeta / hv and f0 = f / hv
    q0(0, k+ks, j+js, i+is) = zeta / hv;
    f0(0, k+ks, j+js, i+is) = coriolis(0, k+ks, j+js, i+is) / hv;

}

void YAKL_INLINE compute_q0(realArr q0, const realArr v, const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
{
real hv = compute_hv(dens, densfct, is, js, ks, i, j, k);
real zeta = compute_zeta(v, is, js, ks, i, j, k);
// compute q0 = zeta / hv and f0 = f / hv
  q0(0, k+ks, j+js, i+is) = zeta / hv;
}

pvpe YAKL_INLINE compute_PVPE(const realArr v, const realArr dens, const realArr densfct, const realArr coriolis, int is, int js, int ks, int i, int j, int k)
{
  pvpe vals;
  real eta = compute_eta(v, coriolis, is, js, ks, i, j, k);
  real hv = compute_hv(dens, densfct, is, js, ks, i, j, k);
  real q0 = eta / hv;

  vals.pv = eta;
  vals.pe = 0.5 * eta * q0;
  
  return vals;
}

};






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
  
  real KE, PE, IE;
  SArray<real,1> h0, h0im1, h0jm1;
  SArray<real,2> U, he;
  SArray<real,2,2> h0arr;

  //compute U = H v
  compute_H<1,diff_ord> (U, v, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

  // Compute h0 = I h needed for phi calcs
  compute_I<1,diff_ord> (h0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
  compute_I<1,diff_ord> (h0im1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i-1, j, k);
  compute_I<1,diff_ord> (h0jm1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j-1, k);

  //compute he = phi h0
  h0arr(0,0) = h0(0);
  h0arr(1,0) = h0(0);
  h0arr(0,1) = h0im1(0);
  h0arr(1,1) = h0jm1(0);
  phi(he, h0arr);

  return 1./2. * (he(0) * ( U(0) * v(0,k+ks,j+js,i+is)) +
              + he(1) * ( U(1) * v(1,k+ks,j+js,i+is)));

}


 void YAKL_INLINE compute_dKdv(realArr F, realArr K, realArr HE, const realArr v, const realArr U, const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,2,2> D0;
  SArray<real,2> he;

  //compute he = phi * h0
  D0(0,0) = dens0(0, k+ks, j+js, i+is);
  D0(0,1) = dens0(0, k+ks, j+js, i+is-1);
  D0(1,0) = dens0(0, k+ks, j+js, i+is);
  D0(1,1) = dens0(0, k+ks, j+js-1, i+is);
  phi(he, D0);
  HE(0, k+ks, j+js, i+is) = he(0);
  HE(1, k+ks, j+js, i+is) = he(1);

  //compute F = he * U
  F(0, k+ks, j+js, i+is) = U(0, k+ks, j+js, i+is) * he(0);
  F(1, k+ks, j+js, i+is) = U(1, k+ks, j+js, i+is) * he(1);

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
   
   real KE, PE, IE;
   SArray<real,1> rhod0, rhod0im1, rhod0jm1;
   SArray<real,3> rhos0, rhos0im1, rhos0jm1;
   SArray<real,2> U, he;
   SArray<real,2,2> h0arr;

   //compute U = H v
   compute_H<1,diff_ord> (U, v, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

   // Compute h0 = I h needed for phi calcs
   compute_I<1,diff_ord> (rhod0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
   compute_I<1,diff_ord> (rhod0im1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i-1, j, k);
   compute_I<1,diff_ord> (rhod0jm1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j-1, k);
   compute_I<3,diff_ord> (rhos0, densfct, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
   compute_I<3,diff_ord> (rhos0im1, densfct, *this->primal_geometry, *this->dual_geometry, is, js, ks, i-1, j, k);
   compute_I<3,diff_ord> (rhos0jm1, densfct, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j-1, k);

   //compute he = phi h0
   // exploits linearity of I/Phi
   h0arr(0,0) = rhod0(0) + rhos0(0) + rhos0(1) + rhos0(2);
   h0arr(1,0) = rhod0(0) + rhos0(0) + rhos0(1) + rhos0(2);
   h0arr(0,1) = rhod0im1(0) + rhos0im1(0) + rhos0im1(1) + rhos0im1(2);
   h0arr(1,1) = rhod0jm1(0) + rhos0jm1(0) + rhos0jm1(1) + rhos0jm1(2);
   phi(he, h0arr);

   return 1./2. * (he(0) * ( U(0) * v(0,k+ks,j+js,i+is)) +
               + he(1) * ( U(1) * v(1,k+ks,j+js,i+is)));

 }


  void YAKL_INLINE compute_dKdv(realArr F, realArr K, realArr HE, const realArr v, const realArr U, const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k)
 {
   SArray<real,2,2> D0;
   SArray<real,2> he;

   //compute he = phi * h0
   // exploits linearity of I/phi
   D0(0,0) = dens0(0, k+ks, j+js, i+is) + densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js, i+is);
   D0(0,1) = dens0(0, k+ks, j+js, i+is-1) + densfct0(0, k+ks, j+js, i+is-1) + densfct0(1, k+ks, j+js, i+is-1) + densfct0(2, k+ks, j+js, i+is-1);
   D0(1,0) = dens0(0, k+ks, j+js, i+is) + densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js, i+is);
   D0(1,1) = dens0(0, k+ks, j+js-1, i+is) + densfct0(0, k+ks, j+js-1, i+is) + densfct0(1, k+ks, j+js-1, i+is) + densfct0(2, k+ks, j+js-1, i+is);
   phi(he, D0);
   HE(0, k+ks, j+js, i+is) = he(0);
   HE(1, k+ks, j+js, i+is) = he(1);

   //compute F = he * U
   F(0, k+ks, j+js, i+is) = U(0, k+ks, j+js, i+is) * he(0);
   F(1, k+ks, j+js, i+is) = U(1, k+ks, j+js, i+is) * he(1);

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
  
  
  
  
  
template<uint nt, uint ntfct> class Hamiltonian_TSWE_Hs {
 public:
   Geometry<ndims,1,1,1> *primal_geometry;
   Geometry<ndims,1,1,1> *dual_geometry;
   bool is_initialized;
   ThermoPotential *thermo;

    Hamiltonian_TSWE_Hs() {
      this->is_initialized = false;
 }
 
 void initialize(ModelParameters &params, ThermoPotential &thermodynamics, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom) 
 {
   this->primal_geometry = &primal_geom;
   this->dual_geometry = &dual_geom;
   this->is_initialized = true;
 }
 
 real YAKL_INLINE compute_PE(const realArr dens, const realArr densfct, const realArr hs, int is, int js, int ks, int i, int j, int k)
 {

  //P = S * hs + 1/2 S * h + sum_nt 1/2 h * t + sum_nt 1/2 h * tfct;
  SArray<real,2> dens0;
  compute_I<2,diff_ord> (dens0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

   real PE = hs(0, k+ks, j+js, i+is) * dens0(1) + 0.5 * dens0(0) * dens(1, k+ks, j+js, i+is);
   for (int l=2; l<2+nt; l++) { PE += 0.5 * dens0(0) * dens(l, k+ks, j+js, i+is);}
   for (int l=0; l<ntfct; l++) { PE += 0.5 * dens0(0) * densfct(l, k+ks, j+js, i+is);}
   return PE;
   
 }

 real YAKL_INLINE compute_IE(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
 {
   return 0;
 }

 //HOW SHOULD HS + OTHER CONSTANTS (GEOPOTENTIAL, ETC.) BE HANDLED IN A UNIFORM WAY
 //IDEALLY THEY ARE PART OF INITIALIZE FOR HAMILTONIAN 
 // MAYBE INITIALIZE TAKES CONSTANT VARS AS ARGUMENTS?
 // THEN WE CAN ALSO DO THINGS LIKE TAKE I HS ONCE, I GEOP ONCE, ETC.
// YES DO IT LIKE THIS...

//For now all Hs just take one extra argument: hs or geop!
//So it is ok for a bit

  void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens0, const realArr densfct0, const realArr hs, int is, int js, int ks, int i, int j, int k)
  {
    // dens0 is h,S,t; densfct is tfct


//ELIMINATE THIS EVENTUALLY...
    //compute I hs
    SArray<real,1> hs0;
    compute_I<1, diff_ord>(hs0, hs, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
    
    // Compute dHdh = 1/2 S + sum_nt 1/2 t + sum_nt 1/2 tfct
    B(0, k+ks, j+js, i+is) = dens0(1, k+ks, j+js, i+is)/2. ;
    for (int l=2; l<2+nt; l++) { B(0, k+ks, j+js, i+is) += dens0(l, k+ks, j+js, i+is)/2.;}
    for (int l=0; l<ntfct; l++) { B(0, k+ks, j+js, i+is) += densfct0(l, k+ks, j+js, i+is)/2.;}
          
    //Compute dHdS = 1/2 h + hs
    B(1, k+ks, j+js, i+is) = hs0(0) + dens0(0, k+ks, j+js, i+is)/2.;

    //Compute dHdt = 1/2 h
    for (int l=2; l<2+nt; l++) { B(l, k+ks, j+js, i+is) = dens0(0, k+ks, j+js, i+is)/2.;}
    
    //Compute dHdtfct = 1/2 h
    for (int l=0; l<ntfct; l++) { Bfct(0, k+ks, j+js, i+is) = dens0(0, k+ks, j+js, i+is)/2.;}
    
  }
};






template<uint nt, uint ntfct> class Hamiltonian_SWE_Hs {
 public:
   Geometry<ndims,1,1,1> *primal_geometry;
   Geometry<ndims,1,1,1> *dual_geometry;
   bool is_initialized;
   real g;
   ThermoPotential *thermo;
   
    Hamiltonian_SWE_Hs() {
      this->is_initialized = false;
 }
 
 void initialize(ModelParameters &params, ThermoPotential &thermodynamics, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom) 
 {
   this->primal_geometry = &primal_geom;
   this->dual_geometry = &dual_geom;
   this->is_initialized = true;
   g = params.g;
 }
 
 real YAKL_INLINE compute_PE(const realArr dens, const realArr densfct, const realArr hs, int is, int js, int ks, int i, int j, int k)
 {

  //P = g * h * hs + 1/2 g * h * h + sum_nt 1/2 h * t + sum_nt 1/2 h * tfct;
  SArray<real,1> h0;
  compute_I<1,diff_ord> (h0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

   real PE = g * hs(0, k+ks, j+js, i+is) * h0(0) + 0.5 * g * h0(0) * dens(0, k+ks, j+js, i+is);
  for (int l=1; l<1+nt; l++) { PE += 0.5 * h0(0) * dens(l, k+ks, j+js, i+is);}
   for (int l=0; l<ntfct; l++) { PE += 0.5 * h0(0) * densfct(l, k+ks, j+js, i+is);}
   return PE;
   
 }

 real YAKL_INLINE compute_IE(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
 {
   return 0;
 }

  void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens0, const realArr densfct0, const realArr hs, int is, int js, int ks, int i, int j, int k)
  {
    // dens0 is h,S,t; densfct is tfct


//ELIMINATE THIS EVENTUALLY...
    //compute I hs
    SArray<real,1> hs0;
    compute_I<1, diff_ord>(hs0, hs, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

    // Compute dHdh = g h + g hs + sum_nt 1/2 t + sum_nt 1/2 tfct
    B(0, k+ks, j+js, i+is) = g * hs(0) + g * dens0(0, k+ks, j+js, i+is);
    for (int l=1; l<1+nt; l++) { B(0, k+ks, j+js, i+is) += dens0(l, k+ks, j+js, i+is)/2.;}
    for (int l=0; l<ntfct; l++) { B(0, k+ks, j+js, i+is) += densfct0(l, k+ks, j+js, i+is)/2.;}
        

    //Compute dHdt = 1/2 h
    for (int l=1; l<1+nt; l++) { B(l, k+ks, j+js, i+is) = dens0(0, k+ks, j+js, i+is)/2.;}
    
    //Compute dHdtfct = 1/2 h
    for (int l=0; l<ntfct; l++) { Bfct(0, k+ks, j+js, i+is) = dens0(0, k+ks, j+js, i+is)/2.;}
    
  }
};




// ADD p-variants


class Hamiltonian_CE_Hs {
 public:
   Geometry<ndims,1,1,1> *primal_geometry;
   Geometry<ndims,1,1,1> *dual_geometry;
   bool is_initialized;
   ThermoPotential *thermo;
   
    Hamiltonian_CE_Hs() {
      this->is_initialized = false;
 }
 
 void initialize(ModelParameters &params, ThermoPotential &thermodynamics, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom)
 {
   this->thermo = &thermodynamics;
   this->primal_geometry = &primal_geom;
   this->dual_geometry = &dual_geom;
   this->is_initialized = true;
 }
 

 real YAKL_INLINE compute_PE(const realArr dens, const realArr densfct, const realArr geop, int is, int js, int ks, int i, int j, int k)
 {
   SArray<real,1> geop0;
   compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

   return dens(0, k+ks, j+js, i+is) * geop0(0);
   
 }
 
 //   alpha = 1./rho;
 //   entropic_var = entropic_density/rho;

     real YAKL_INLINE compute_alpha(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
       return 1./dens0(0, k+ks, j+js, i+is);
     }
     
     real YAKL_INLINE compute_entropic_var(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
       return dens0(1, k+ks, j+js, i+is)/dens0(0, k+ks, j+js, i+is);
     }

 real YAKL_INLINE compute_IE(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
 {
   SArray<real,2> dens0;
   compute_I<2, diff_ord>(dens0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
   real alpha = 1./dens0(0);
   real entropic_var = dens0(1)/dens0(0);
   return dens(0, k+ks, j+js, i+is) * thermo->compute_U(alpha, entropic_var, 0, 0, 0, 0);
 }

  void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens0, const realArr densfct0, const realArr geop, int is, int js, int ks, int i, int j, int k)
  {

    SArray<real,1> geop0;
    compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

    real alpha = compute_alpha(dens0, densfct0, is, js, ks, i, j, k);
    real entropic_var = compute_entropic_var(dens0, densfct0, is, js, ks, i, j, k);
    
    real U = thermo->compute_U(alpha, entropic_var, 0, 0, 0, 0);
    real p = -thermo->compute_dUdalpha(alpha, entropic_var, 0, 0, 0, 0);
    real generalized_Exner = thermo->compute_dUdentropic_var_density(alpha, entropic_var, 0, 0, 0, 0);

    B(0, k+ks, j+js, i+is) = geop0(0) + U + p * alpha - entropic_var * generalized_Exner;
    B(1, k+ks, j+js, i+is) = generalized_Exner;

  }
};







class Hamiltonian_MCE_rho_Hs {
 public:
   Geometry<ndims,1,1,1> *primal_geometry;
   Geometry<ndims,1,1,1> *dual_geometry;
   bool is_initialized;
   ThermoPotential *thermo;
   
    Hamiltonian_MCE_rho_Hs() {
      this->is_initialized = false;
 }
 
 void initialize(ModelParameters &params, ThermoPotential &thermodynamics, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom)
 {
   this->thermo = &thermodynamics;
   this->primal_geometry = &primal_geom;
   this->dual_geometry = &dual_geom;
   this->is_initialized = true;
 }
 

 real YAKL_INLINE compute_PE(const realArr dens, const realArr densfct, const realArr geop, int is, int js, int ks, int i, int j, int k)
 {
   SArray<real,1> geop0;
   compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
   return dens(0, k+ks, j+js, i+is) * geop0(0);
   
 }
 
 //   alpha = 1./rho;
 //   compute_entropic_var = entropic_density/rho;
 //   qd = (rho - \sum_n \rho_s)/rho;
 //   qv = rho_v/rho;
 //   ql = rho_l/rho;
 //   qi = rho_i/rho;
     real YAKL_INLINE compute_alpha(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
       return 1./dens0(0, k+ks, j+js, i+is);
     }
     
     real YAKL_INLINE compute_entropic_var(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
       return dens0(1, k+ks, j+js, i+is)/dens0(0, k+ks, j+js, i+is);
     }
     real YAKL_INLINE compute_qd(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
       return (dens0(0, k+ks, j+js, i+is) - densfct0(0, k+ks, j+js, i+is) - densfct0(1, k+ks, j+js, i+is) - densfct0(2, k+ks, j+js, i+is))/dens0(0, k+ks, j+js, i+is);
     }
     real YAKL_INLINE compute_qv(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
       return densfct0(0, k+ks, j+js, i+is)/dens0(0, k+ks, j+js, i+is);
     }
     real YAKL_INLINE compute_ql(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
       return densfct0(1, k+ks, j+js, i+is)/dens0(0, k+ks, j+js, i+is);
     }
     real YAKL_INLINE compute_qi(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
       return densfct0(2, k+ks, j+js, i+is)/dens0(0, k+ks, j+js, i+is);
     }


 real YAKL_INLINE compute_IE(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
 {

   SArray<real,2> dens0;
   SArray<real,3> densfct0;
   compute_I<2, diff_ord>(dens0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
   compute_I<3, diff_ord>(densfct0, densfct, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
   real alpha = 1./dens0(0);
   real entropic_var = dens0(1)/dens0(0);
   real qd = (dens0(0) - densfct(0) - densfct(1) - densfct(2))/dens0(0);
   real ql = densfct(0) / dens0(0);
   real qi = densfct(1) / dens0(0);
   real qv = densfct(2) / dens0(0);
   return dens(0, k+ks, j+js, i+is) * thermo->compute_U(alpha, entropic_var, qd, qv, ql, qi);
 }

  void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens0, const realArr densfct0, const realArr geop, int is, int js, int ks, int i, int j, int k)
  {

    SArray<real,1> geop0;
    compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
    
    real alpha = compute_alpha(dens0, densfct0, is, js, ks, i, j, k);
    real entropic_var = compute_entropic_var(dens0, densfct0, is, js, ks, i, j, k);
    real qd = compute_qd(dens0, densfct0, is, js, ks, i, j, k);
    real ql = compute_ql(dens0, densfct0, is, js, ks, i, j, k);
    real qi = compute_qi(dens0, densfct0, is, js, ks, i, j, k);
    real qv = compute_qv(dens0, densfct0, is, js, ks, i, j, k);
    
    real U = thermo->compute_U(alpha, entropic_var, qd, qv, ql, qi);
    real p = -thermo->compute_dUdalpha(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_Exner = thermo->compute_dUdentropic_var_density(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_d = thermo->compute_dUdqd(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_v = thermo->compute_dUdqv(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_l = thermo->compute_dUdql(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_i = thermo->compute_dUdqi(alpha, entropic_var, qd, qv, ql, qi);

    B(0, k+ks, j+js, i+is) = geop0(0) + U + p * alpha - entropic_var * generalized_Exner + 
    qv * (generalized_chemical_potential_d - generalized_chemical_potential_v) + 
    ql * (generalized_chemical_potential_d - generalized_chemical_potential_l) + 
    qi * (generalized_chemical_potential_d - generalized_chemical_potential_i);
    B(1, k+ks, j+js, i+is) = generalized_Exner;
    Bfct(0, k+ks, j+js, i+is) = generalized_chemical_potential_v - generalized_chemical_potential_d;
    Bfct(1, k+ks, j+js, i+is) = generalized_chemical_potential_l - generalized_chemical_potential_d;
    Bfct(2, k+ks, j+js, i+is) = generalized_chemical_potential_i - generalized_chemical_potential_d;
    
  }
};



class Hamiltonian_MCE_rhod_Hs {
 public:
   Geometry<ndims,1,1,1> *primal_geometry;
   Geometry<ndims,1,1,1> *dual_geometry;
   bool is_initialized;
   ThermoPotential *thermo;
   
    Hamiltonian_MCE_rhod_Hs() {
      this->is_initialized = false;
 }
 
 void initialize(ModelParameters &params, ThermoPotential &thermodynamics, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom)
 {
   this->thermo = &thermodynamics;
   this->primal_geometry = &primal_geom;
   this->dual_geometry = &dual_geom;
   this->is_initialized = true;
 }
 

 real YAKL_INLINE compute_PE(const realArr dens, const realArr densfct, const realArr geop, int is, int js, int ks, int i, int j, int k)
 {
   
   SArray<real,1> geop0;
   compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
   
   return (dens(0, k+ks, j+js, i+is) + densfct(0, k+ks, j+js, i+is) + densfct(1, k+ks, j+js, i+is) + densfct(2, k+ks, j+js, i+is) ) * geop0(0);
   
 }
 
 //   alpha = 1./rho;
 //   compute_entropic_var = entropic_density/rho;
 //   qd = rho_d/rho;
 //   qv = rho_v/rho;
 //   ql = rho_l/rho;
 //   qi = rho_i/rho;
 // with rho = rho_d + rho_v + rho_l + rho_i
     real YAKL_INLINE compute_alpha(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
       return 1./(dens0(0, k+ks, j+js, i+is) + densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js, i+is));
     }
     
     real YAKL_INLINE compute_entropic_var(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
       return dens0(1, k+ks, j+js, i+is)/(dens0(0, k+ks, j+js, i+is) + densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js, i+is));
     }
     real YAKL_INLINE compute_qd(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
       return (dens0(0, k+ks, j+js, i+is))/(dens0(0, k+ks, j+js, i+is) + densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js, i+is));
     }
     real YAKL_INLINE compute_qv(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
       return densfct0(0, k+ks, j+js, i+is)/(dens0(0, k+ks, j+js, i+is) + densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js, i+is));
     }
     real YAKL_INLINE compute_ql(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
       return densfct0(1, k+ks, j+js, i+is)/(dens0(0, k+ks, j+js, i+is) + densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js, i+is));
     }
     real YAKL_INLINE compute_qi(const realArr dens0, const realArr densfct0, int is, int js, int ks, int i, int j, int k) {
       return densfct0(2, k+ks, j+js, i+is)/(dens0(0, k+ks, j+js, i+is) + densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js, i+is));
     }


 real YAKL_INLINE compute_IE(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
 {

   SArray<real,2> dens0;
   SArray<real,3> densfct0;
   compute_I<2, diff_ord>(dens0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
   compute_I<3, diff_ord>(densfct0, densfct, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
   real rho = dens0(0) + densfct(0) + densfct(1) + densfct(2);
   real alpha = 1./rho;
   real entropic_var = dens0(1)/rho;
   real qd = dens0(0) / rho;
   real ql = densfct(0) / rho;
   real qi = densfct(1) / rho;
   real qv = densfct(2) / rho;
   return rho * thermo->compute_U(alpha, entropic_var, qd, qv, ql, qi);
 }

  void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens0, const realArr densfct0, const realArr geop, int is, int js, int ks, int i, int j, int k)
  {

    SArray<real,1> geop0;
    compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
    
    real alpha = compute_alpha(dens0, densfct0, is, js, ks, i, j, k);
    real entropic_var = compute_entropic_var(dens0, densfct0, is, js, ks, i, j, k);
    real qd = compute_qd(dens0, densfct0, is, js, ks, i, j, k);
    real ql = compute_ql(dens0, densfct0, is, js, ks, i, j, k);
    real qi = compute_qi(dens0, densfct0, is, js, ks, i, j, k);
    real qv = compute_qv(dens0, densfct0, is, js, ks, i, j, k);
    
    real U = thermo->compute_U(alpha, entropic_var, qd, qv, ql, qi);
    real p = -thermo->compute_dUdalpha(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_Exner = thermo->compute_dUdentropic_var_density(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_d = thermo->compute_dUdqd(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_v = thermo->compute_dUdqv(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_l = thermo->compute_dUdql(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_i = thermo->compute_dUdqi(alpha, entropic_var, qd, qv, ql, qi);

    real QNterm = qd * generalized_chemical_potential_d + qv * generalized_chemical_potential_v + qi * generalized_chemical_potential_i + qi * generalized_chemical_potential_i;
    B(0, k+ks, j+js, i+is) = geop0(0) + U + p * alpha - entropic_var * generalized_Exner - generalized_chemical_potential_d - QNterm; 
    Bfct(0, k+ks, j+js, i+is) = geop0(0) + U + p * alpha - entropic_var * generalized_Exner - generalized_chemical_potential_v - QNterm; 
    Bfct(1, k+ks, j+js, i+is) = geop0(0) + U + p * alpha - entropic_var * generalized_Exner - generalized_chemical_potential_l - QNterm; 
    Bfct(2, k+ks, j+js, i+is) = geop0(0) + U + p * alpha - entropic_var * generalized_Exner - generalized_chemical_potential_i - QNterm; 
    B(1, k+ks, j+js, i+is) = generalized_Exner;
  }
};




// p-variants
class Hamiltonian_CE_p_Hs : public Hamiltonian_CE_Hs
{
  real YAKL_INLINE compute_IE(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
  {
    SArray<real,2> dens0;
    compute_I<2, diff_ord>(dens0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
    real entropic_var = dens0(1)/dens0(0);
    real p = thermo->solve_p(dens0(0), entropic_var, 0, 0, 0, 0);
    return dens(0, k+ks, j+js, i+is) * thermo->compute_H(p, entropic_var, 0, 0, 0, 0) - p;
  }

   void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens0, const realArr densfct0, const realArr geop, int is, int js, int ks, int i, int j, int k)
   {

     SArray<real,1> geop0;
     compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

     real entropic_var = compute_entropic_var(dens0, densfct0, is, js, ks, i, j, k);
     real p = thermo->solve_p(dens0(0, k+ks, j+js, i+is), entropic_var, 0, 0, 0, 0);

     real H = thermo->compute_H(p, entropic_var, 0, 0, 0, 0);
     real generalized_Exner = thermo->compute_dHdentropic_var_density(p, entropic_var, 0, 0, 0, 0);

     B(0, k+ks, j+js, i+is) = geop0(0) + H - entropic_var * generalized_Exner;
     B(1, k+ks, j+js, i+is) = generalized_Exner;

   }
};

class Hamiltonian_MCE_rho_p_Hs : public Hamiltonian_MCE_rho_Hs
{
  real YAKL_INLINE compute_IE(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
  {

    SArray<real,2> dens0;
    SArray<real,3> densfct0;
    compute_I<2, diff_ord>(dens0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
    compute_I<3, diff_ord>(densfct0, densfct, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);

    real entropic_var = dens0(1)/dens0(0);
    real qd = (dens0(0) - densfct(0) - densfct(1) - densfct(2))/dens0(0);
    real ql = densfct(0) / dens0(0);
    real qi = densfct(1) / dens0(0);
    real qv = densfct(2) / dens0(0);
    real p = thermo->solve_p(dens0(0), entropic_var, qd, qv, ql, qi);
    return dens(0, k+ks, j+js, i+is) * thermo->compute_H(p, entropic_var, qd, qv, ql, qi) - p;
  }

   void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens0, const realArr densfct0, const realArr geop, int is, int js, int ks, int i, int j, int k)
   {

     SArray<real,1> geop0;
     compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
     
     real entropic_var = compute_entropic_var(dens0, densfct0, is, js, ks, i, j, k);
     real qd = compute_qd(dens0, densfct0, is, js, ks, i, j, k);
     real ql = compute_ql(dens0, densfct0, is, js, ks, i, j, k);
     real qi = compute_qi(dens0, densfct0, is, js, ks, i, j, k);
     real qv = compute_qv(dens0, densfct0, is, js, ks, i, j, k);

     real p = thermo->solve_p(dens0(0, k+ks, j+js, i+is), entropic_var, qd, qv, ql, qi);
     real H = thermo->compute_H(p, entropic_var, qd, qv, ql, qi);
     real generalized_Exner = thermo->compute_dHdentropic_var_density(p, entropic_var, qd, qv, ql, qi);
     real generalized_chemical_potential_d = thermo->compute_dHdqd(p, entropic_var, qd, qv, ql, qi);
     real generalized_chemical_potential_v = thermo->compute_dHdqv(p, entropic_var, qd, qv, ql, qi);
     real generalized_chemical_potential_l = thermo->compute_dHdql(p, entropic_var, qd, qv, ql, qi);
     real generalized_chemical_potential_i = thermo->compute_dHdqi(p, entropic_var, qd, qv, ql, qi);

     B(0, k+ks, j+js, i+is) = geop0(0) + H - entropic_var * generalized_Exner + 
     qv * (generalized_chemical_potential_d - generalized_chemical_potential_v) + 
     ql * (generalized_chemical_potential_d - generalized_chemical_potential_l) + 
     qi * (generalized_chemical_potential_d - generalized_chemical_potential_i);
     B(1, k+ks, j+js, i+is) = generalized_Exner;
     Bfct(0, k+ks, j+js, i+is) = generalized_chemical_potential_v - generalized_chemical_potential_d;
     Bfct(1, k+ks, j+js, i+is) = generalized_chemical_potential_l - generalized_chemical_potential_d;
     Bfct(2, k+ks, j+js, i+is) = generalized_chemical_potential_i - generalized_chemical_potential_d;
     
   }
};

class Hamiltonian_MCE_rhod_p_Hs : public Hamiltonian_MCE_rhod_Hs
{
  real YAKL_INLINE compute_IE(const realArr dens, const realArr densfct, int is, int js, int ks, int i, int j, int k)
  {

    SArray<real,2> dens0;
    SArray<real,3> densfct0;
    compute_I<2, diff_ord>(dens0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
    compute_I<3, diff_ord>(densfct0, densfct, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
    
    real rho = dens0(0) + densfct(0) + densfct(1) + densfct(2);
    real entropic_var = dens0(1)/rho;
    real qd = dens0(0) / rho;
    real ql = densfct(0) / rho;
    real qi = densfct(1) / rho;
    real qv = densfct(2) / rho;
    real p = thermo->solve_p(rho, entropic_var, qd, qv, ql, qi);
    return rho * thermo->compute_H(p, entropic_var, qd, qv, ql, qi) - p;
  }

   void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens0, const realArr densfct0, const realArr geop, int is, int js, int ks, int i, int j, int k)
   {

     SArray<real,1> geop0;
     compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
     
     real alpha = compute_alpha(dens0, densfct0, is, js, ks, i, j, k);
     real entropic_var = compute_entropic_var(dens0, densfct0, is, js, ks, i, j, k);
     real qd = compute_qd(dens0, densfct0, is, js, ks, i, j, k);
     real ql = compute_ql(dens0, densfct0, is, js, ks, i, j, k);
     real qi = compute_qi(dens0, densfct0, is, js, ks, i, j, k);
     real qv = compute_qv(dens0, densfct0, is, js, ks, i, j, k);
     real p = thermo->solve_p(1./alpha, entropic_var, qd, qv, ql, qi);

     real H = thermo->compute_H(p, entropic_var, qd, qv, ql, qi);
     real generalized_Exner = thermo->compute_dHdentropic_var_density(p, entropic_var, qd, qv, ql, qi);
     real generalized_chemical_potential_d = thermo->compute_dHdqd(p, entropic_var, qd, qv, ql, qi);
     real generalized_chemical_potential_v = thermo->compute_dHdqv(p, entropic_var, qd, qv, ql, qi);
     real generalized_chemical_potential_l = thermo->compute_dHdql(p, entropic_var, qd, qv, ql, qi);
     real generalized_chemical_potential_i = thermo->compute_dHdqi(p, entropic_var, qd, qv, ql, qi);

     real QNterm = qd * generalized_chemical_potential_d + qv * generalized_chemical_potential_v + qi * generalized_chemical_potential_i + qi * generalized_chemical_potential_i;
     B(0, k+ks, j+js, i+is) = geop0(0) + H - entropic_var * generalized_Exner - generalized_chemical_potential_d - QNterm; 
     Bfct(0, k+ks, j+js, i+is) = geop0(0) + H - entropic_var * generalized_Exner - generalized_chemical_potential_v - QNterm; 
     Bfct(1, k+ks, j+js, i+is) = geop0(0) + H - entropic_var * generalized_Exner - generalized_chemical_potential_l - QNterm; 
     Bfct(2, k+ks, j+js, i+is) = geop0(0) + H - entropic_var * generalized_Exner - generalized_chemical_potential_i - QNterm; 
     B(1, k+ks, j+js, i+is) = generalized_Exner;
   }
};



#endif
