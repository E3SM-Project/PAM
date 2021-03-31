#ifndef _HAMILTONIAN_EULER_H_
#define _HAMILTONIAN_EULER_H_

#include "common.h"
#include "thermo.h"



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
   real qv = densfct(0) / dens0(0);
   real ql = densfct(1) / dens0(0);
   real qi = densfct(2) / dens0(0);
   return dens(0, k+ks, j+js, i+is) * thermo->compute_U(alpha, entropic_var, qd, qv, ql, qi);
 }

  void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens0, const realArr densfct0, const realArr geop, int is, int js, int ks, int i, int j, int k)
  {

    SArray<real,1> geop0;
    compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
    
    real alpha = compute_alpha(dens0, densfct0, is, js, ks, i, j, k);
    real entropic_var = compute_entropic_var(dens0, densfct0, is, js, ks, i, j, k);
    real qd = compute_qd(dens0, densfct0, is, js, ks, i, j, k);
    real qv = compute_qv(dens0, densfct0, is, js, ks, i, j, k);
    real ql = compute_ql(dens0, densfct0, is, js, ks, i, j, k);
    real qi = compute_qi(dens0, densfct0, is, js, ks, i, j, k);
    
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
   real qv = densfct(0) / rho;
   real ql = densfct(1) / rho;
   real qi = densfct(2) / rho;
   return rho * thermo->compute_U(alpha, entropic_var, qd, qv, ql, qi);
 }

  void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens0, const realArr densfct0, const realArr geop, int is, int js, int ks, int i, int j, int k)
  {

    SArray<real,1> geop0;
    compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
    
    real alpha = compute_alpha(dens0, densfct0, is, js, ks, i, j, k);
    real entropic_var = compute_entropic_var(dens0, densfct0, is, js, ks, i, j, k);
    real qd = compute_qd(dens0, densfct0, is, js, ks, i, j, k);
    real qv = compute_qv(dens0, densfct0, is, js, ks, i, j, k);
    real ql = compute_ql(dens0, densfct0, is, js, ks, i, j, k);
    real qi = compute_qi(dens0, densfct0, is, js, ks, i, j, k);
    
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
    real qv = densfct(0) / dens0(0);
    real ql = densfct(1) / dens0(0);
    real qi = densfct(2) / dens0(0);
    real p = thermo->solve_p(dens0(0), entropic_var, qd, qv, ql, qi);
    return dens(0, k+ks, j+js, i+is) * thermo->compute_H(p, entropic_var, qd, qv, ql, qi) - p;
  }

   void YAKL_INLINE compute_dHsdx(realArr B, realArr Bfct, const realArr dens0, const realArr densfct0, const realArr geop, int is, int js, int ks, int i, int j, int k)
   {

     SArray<real,1> geop0;
     compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
     
     real entropic_var = compute_entropic_var(dens0, densfct0, is, js, ks, i, j, k);
     real qd = compute_qd(dens0, densfct0, is, js, ks, i, j, k);
     real qv = compute_qv(dens0, densfct0, is, js, ks, i, j, k);
     real ql = compute_ql(dens0, densfct0, is, js, ks, i, j, k);
     real qi = compute_qi(dens0, densfct0, is, js, ks, i, j, k);

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
    real qv = densfct(0) / rho;
    real ql = densfct(1) / rho;
    real qi = densfct(2) / rho;
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
     real qv = compute_qv(dens0, densfct0, is, js, ks, i, j, k);
     real ql = compute_ql(dens0, densfct0, is, js, ks, i, j, k);
     real qi = compute_qi(dens0, densfct0, is, js, ks, i, j, k);
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
