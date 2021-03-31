#ifndef _HAMILTONIAN_LAYERMODEL_H_
#define _HAMILTONIAN_LAYERMODEL_H_

#include "common.h"
#include "thermo.h"



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


#endif