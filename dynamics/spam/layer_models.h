#pragma once

#include "common.h"
#include "thermo.h"



class Hamiltonian_TSWE_Hs {
public:
 Geometry<1,1,1> *primal_geometry;
 Geometry<1,1,1> *dual_geometry;
 bool is_initialized;
 ThermoPotential *thermo;

  Hamiltonian_TSWE_Hs() {
    this->is_initialized = false;
}

void initialize(Parameters &params, ThermoPotential &thermodynamics, Geometry<1,1,1> &primal_geom, Geometry<1,1,1> &dual_geom) 
{
 this->primal_geometry = &primal_geom;
 this->dual_geometry = &dual_geom;
 this->is_initialized = true;
}

real YAKL_INLINE compute_PE(const real5d dens, const real5d hs, int is, int js, int ks, int i, int j, int k, int n)
{
// Assumes things are stored in dens as as 0=h, 1=S, tracers
//P = S * hs + 1/2 S * h + sum_nt 1/2 h * t + sum_nt 1/2 h * tfct;
SArray<real,1,2> dens0;
#ifdef _EXTRUDED   
compute_Iext<2,diff_ord,vert_diff_ord> (dens0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k, n);
#else
compute_I<2,diff_ord> (dens0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k, n);
#endif
 real PE = hs(0, k+ks, j+js, i+is, n) * dens0(1) + 0.5 * dens0(0) * dens(1, k+ks, j+js, i+is, n);
 for (int l=2; l<ndensity; l++) { PE += 0.5 * dens0(0) * dens(l, k+ks, j+js, i+is, n);}
 return PE;
 
}

real YAKL_INLINE compute_IE(const real5d dens, int is, int js, int ks, int i, int j, int k, int n)
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

void YAKL_INLINE compute_dHsdx(real5d B, const real5d dens0, const real5d hs, int is, int js, int ks, int i, int j, int k, int n)
{

  // Assumes things are stored in dens as as 0=h, 1=S, tracers

  //compute I hs
  SArray<real,1,1> hs0;
  #ifdef _EXTRUDED   
  compute_Iext<1, diff_ord, vert_diff_ord>(hs0, hs, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k, n);
  #else
  compute_I<1, diff_ord>(hs0, hs, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k, n);
  #endif
  //ELIMINATE THIS EVENTUALLY IE EITHER STORE AS A SEPARATE CONSTANT OR COMPUTE ONCE?...
  
  // Compute dHdh = 1/2 S + sum_nt 1/2 t
  B(0, k+ks, j+js, i+is, n) = dens0(1, k+ks, j+js, i+is, n)/2. ;
  for (int l=2; l<ndensity; l++) { B(0, k+ks, j+js, i+is, n) += dens0(l, k+ks, j+js, i+is, n)/2.;}
        
  //Compute dHdS = 1/2 h + hs
  B(1, k+ks, j+js, i+is, n) = hs0(0) + dens0(0, k+ks, j+js, i+is, n)/2.;

  //Compute dHdt = 1/2 h for tracers
  for (int l=2; l<ndensity; l++) { B(l, k+ks, j+js, i+is, n) = dens0(0, k+ks, j+js, i+is, n)/2.;}
  
}
};






class Hamiltonian_SWE_Hs {
public:
 Geometry<1,1,1> *primal_geometry;
 Geometry<1,1,1> *dual_geometry;
 bool is_initialized;
 real g;
 ThermoPotential *thermo;
 
  Hamiltonian_SWE_Hs() {
    this->is_initialized = false;
}

void initialize(Parameters &params, ThermoPotential &thermodynamics, Geometry<1,1,1> &primal_geom, Geometry<1,1,1> &dual_geom) 
{
 this->primal_geometry = &primal_geom;
 this->dual_geometry = &dual_geom;
 this->is_initialized = true;
 g = params.g;
}

real YAKL_INLINE compute_PE(const real5d dens, const real5d hs, int is, int js, int ks, int i, int j, int k, int n)
{
  // Assumes things are stored in dens as as 0=h, tracers

//P = g * h * hs + 1/2 g * h * h + sum_nt 1/2 h * t + sum_nt 1/2 h * tfct;
SArray<real,1,1> h0;
#ifdef _EXTRUDED   
compute_Iext<1,diff_ord,vert_diff_ord> (h0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k, n);
#else
compute_I<1,diff_ord> (h0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k, n);
#endif
//ELIMINATE THIS EVENTUALLY IE EITHER STORE AS A SEPARATE CONSTANT OR COMPUTE ONCE?...

 real PE = g * hs(0, k+ks, j+js, i+is, n) * h0(0) + 0.5 * g * h0(0) * dens(0, k+ks, j+js, i+is, n);
for (int l=1; l<ndensity; l++) { PE += 0.5 * h0(0) * dens(l, k+ks, j+js, i+is, n);}
 return PE;
 
}

real YAKL_INLINE compute_IE(const real5d dens, int is, int js, int ks, int i, int j, int k, int n)
{
 return 0;
}

void YAKL_INLINE compute_dHsdx(real5d B, const real5d dens0, const real5d hs, int is, int js, int ks, int i, int j, int k, int n)
{
  // Assumes things are stored in dens as as 0=h, tracers


//ELIMINATE THIS EVENTUALLY...
  //compute I hs
  SArray<real,1,1> hs0;
  #ifdef _EXTRUDED   
  compute_Iext<1, diff_ord,vert_diff_ord>(hs0, hs, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k, n);
  #else
  compute_I<1, diff_ord>(hs0, hs, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k, n);
  #endif

  // Compute dHdh = g h + g hs + sum_nt 1/2 t + sum_nt 1/2 tfct
  B(0, k+ks, j+js, i+is, n) = g * hs(0) + g * dens0(0, k+ks, j+js, i+is, n);
  for (int l=1; l<ndensity; l++) { B(0, k+ks, j+js, i+is) += dens0(l, k+ks, j+js, i+is, n)/2.;}

  //Compute dHdt = 1/2 h
  for (int l=1; l<ndensity; l++) { B(l, k+ks, j+js, i+is, n) = dens0(0, k+ks, j+js, i+is, n)/2.;}

}
};


