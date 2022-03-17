#ifndef _HAMILTONIAN_H_
#define _HAMILTONIAN_H_

#include "common.h"
#include "kinetic_energy.h"
#include "layer_models.h"
#include "compressible_euler.h"
#include "anelastic.h"
#include "functionals.h"

// 
// class VariableSet {
// 
// public: 
// 
//   virtual real YAKL_INLINE compute_alpha(const real4d dens0, int is, int js, int ks, int i, int j, int k) {};
//   virtual real YAKL_INLINE compute_specific_entropic_var(const real4d dens0, int is, int js, int ks, int i, int j, int k)  {};
//   virtual real YAKL_INLINE compute_qd(const real4d dens0, int is, int js, int ks, int i, int j, int k)  {};
//   virtual real YAKL_INLINE compute_qv(const real4d dens0, int is, int js, int ks, int i, int j, int k)  {};
//   virtual real YAKL_INLINE compute_ql(const real4d dens0, int is, int js, int ks, int i, int j, int k)  {};
//   virtual real YAKL_INLINE compute_qi(const real4d dens0, int is, int js, int ks, int i, int j, int k)  {};  
// 
//   virtual void YAKL_INLINE compute_D(const real4d dens0, int is, int js, int ks, int i, int j, int k) {};
//   virtual real YAKL_INLINE compute_entropic_var(const real4d dens0, int is, int js, int ks, int i, int j, int k)  {};
//   // Here we do assume that D is a linear combination of Dks. This is very reasonable.
//   virtual int YAKL_INLINE compute_partialDpartialDl(int l) {};

// };
// 
// 
// // dens is always [nofct, fct]
// 
// //This var set covers
// //swe/tswe with active tracers: [h, S, t_nofct], [t_fct,]
// //ce + mce_rho with inactive tracers: [rho, S], [rho_v, rho_l, rho_i]
// //cefct + mcefct_rho with inactive tracers: [], [rho, S, rho_v, rho_l, rho_i]
// // ie rho is the 1st variable
// // note that fct and nofct versions of ce/mce have same variable ordering! (but will have different ndensity_nofct and ndensity_fct...)
// class rhoVarSet {
//   real YAKL_INLINE compute_alpha(const real4d dens0, int is, int js, int ks, int i, int j, int k) {
//     return 1./dens0(0, k+ks, j+js, i+is);
//   };
//   real YAKL_INLINE compute_specific_entropic_var(const real4d dens0, int is, int js, int ks, int i, int j, int k) {
//     return dens0(1, k+ks, j+js, i+is)/dens0(0, k+ks, j+js, i+is);
//   };
//   real YAKL_INLINE compute_qd(const real4d dens0, int is, int js, int ks, int i, int j, int k) {
//     return (dens0(0, k+ks, j+js, i+is) - dens0(2, k+ks, j+js, i+is) - dens0(3, k+ks, j+js, i+is) - dens0(4, k+ks, j+js, i+is))/dens0(0, k+ks, j+js, i+is);
//   };
//   real YAKL_INLINE compute_qv(const real4d dens0, int is, int js, int ks, int i, int j, int k) {
//     return dens0(2, k+ks, j+js, i+is)/dens0(0, k+ks, j+js, i+is);
//   };
//   real YAKL_INLINE compute_ql(const real4d dens0, int is, int js, int ks, int i, int j, int k) {
//     return dens0(3, k+ks, j+js, i+is)/dens0(0, k+ks, j+js, i+is);
//   };
//   real YAKL_INLINE compute_qi(const real4d dens0, int is, int js, int ks, int i, int j, int k) {
//     return dens0(4, k+ks, j+js, i+is)/dens0(0, k+ks, j+js, i+is);
//   };
//   void YAKL_INLINE compute_D(const real4d dens, int is, int js, int ks, int i, int j, int k) {
//   return dens(0, k+ks, j+js, i+is);
//   };
//   real YAKL_INLINE compute_entropic_var(const real4d dens, int is, int js, int ks, int i, int j, int k)  {
//   return dens(1, k+ks, j+js, i+is);
//   };
//   real YAKL_INLINE compute_partialDpartialDl(int l) {
//    if (l==0) {return 1};
//    return 0;
//};

// 
// }
// 
// //This var set covers
// //mce_rhod with inactive tracers: [rhod, S], [rho_v, rho_l, rho_i]
// //mcefct_rhod with inactive tracers: [], [rhod, S, rho_v, rho_l, rho_i]
// // ie rho is linear combination of variables in dens
// // note that fct and nofct versions of ce/mce have same variable ordering! (but will have different ndensity_nofct and ndensity_fct...)
// class rhodVarSet {
// 
//   //   alpha = 1./rho;
//   //   compute_entropic_var = entropic_density/rho;
//   //   qd = rho_d/rho;
//   //   qv = rho_v/rho;
//   //   ql = rho_l/rho;
//   //   qi = rho_i/rho;
//   // with rho = rho_d + rho_v + rho_l + rho_i
//       real YAKL_INLINE compute_alpha(const real4d dens0, const real4d densfct0, int is, int js, int ks, int i, int j, int k) {
//         return 1./compute_D(dens0, is, js, ks, i, j, k);
//       }
// 
//       real YAKL_INLINE compute_specific_entropic_var(const real4d dens0, const real4d densfct0, int is, int js, int ks, int i, int j, int k) {
//         return dens0(1, k+ks, j+js, i+is)/compute_D(dens0, is, js, ks, i, j, k);
//       }
//       real YAKL_INLINE compute_qd(const real4d dens0, const real4d densfct0, int is, int js, int ks, int i, int j, int k) {
//         return dens0(0, k+ks, j+js, i+is))/compute_D(dens0, is, js, ks, i, j, k);
//       }
//       real YAKL_INLINE compute_qv(const real4d dens0, const real4d densfct0, int is, int js, int ks, int i, int j, int k) {
//         return dens0(2, k+ks, j+js, i+is)/compute_D(dens0, is, js, ks, i, j, k);
//       }
//       real YAKL_INLINE compute_ql(const real4d dens0, const real4d densfct0, int is, int js, int ks, int i, int j, int k) {
//         return dens0(3, k+ks, j+js, i+is)/compute_D(dens0, is, js, ks, i, j, k);
//       }
//       real YAKL_INLINE compute_qi(const real4d dens0, const real4d densfct0, int is, int js, int ks, int i, int j, int k) {
//         return dens0(4, k+ks, j+js, i+is)/compute_D(dens0, is, js, ks, i, j, k);
//       }
// 
//   void YAKL_INLINE compute_D(const real4d dens, int is, int js, int ks, int i, int j, int k) {
//   return dens0(0, k+ks, j+js, i+is) + dens0(2, k+ks, j+js, i+is) + dens0(3, k+ks, j+js, i+is) + dens0(4, k+ks, j+js, i+is);
//   };
//   real YAKL_INLINE compute_entropic_var(const real4d dens, int is, int js, int ks, int i, int j, int k)  {
//   return dens(1, k+ks, j+js, i+is);
//   };
//   real YAKL_INLINE compute_partialDpartialDl(int l) {
//    if (l==0 || l==2 || l ==3 || l==4) {return 1};
//    return 0;
//    };
// }
// 
// 
// // 
// class Hamiltonian {
// 
// public:
//   Geometry<1,1,1> *primal_geometry;
//   Geometry<1,1,1> *dual_geometry;
//   bool is_initialized;
// 
//    Hamiltonian() {
//      this->is_initialized = false;
// }
// 
// void initialize(ModelParameters &params, Geometry<1,1,1> &primal_geom, Geometry<1,1,1> &dual_geom)
// {
//   this->primal_geometry = &primal_geom;
//   this->dual_geometry = &dual_geom;
//   this->is_initialized = true;
// }
// };
// 
// 
// class Hk : public VariableSet, public Hamiltonian {
// public:
// real YAKL_INLINE compute_KE(const real4d v, const real4d dens, int is, int js, int ks, int i, int j, int k)
// {
// 
//   real KE;
//   SArray<real,1> h0, h0im1, h0jm1, h0km1;
//   SArray<real,ndims> U, he;
//   SArray<real,ndims,2> h0arr;
// 
//   //compute U = H v
//   compute_H<1,diff_ord> (U, v, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
// 
//   // Compute h0 = I h needed for phi calcs
//   compute_I<1,diff_ord> (h0, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
//   compute_I<1,diff_ord> (h0im1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i-1, j, k);
//   if (ndims>=2) {compute_I<1,diff_ord> (h0jm1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j-1, k);}
//   if (ndims>=3) {compute_I<1,diff_ord> (h0km1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k-1);}
// 
//   //compute he = phi h0
//   for (int d=0; d<ndims; d++)
//   {
//     if (d == 0)
//     {
//       h0arr(d,0) = h0(0);
//       h0arr(d,1) = h0im1(0);
//     }
//     if (d == 1)
//     {
//       h0arr(d,0) = h0(0);
//       h0arr(d,1) = h0jm1(0);
//     }
//     if (d == 2)
//     {
//       h0arr(d,0) = h0(0);
//       h0arr(d,1) = h0km1(0);
//     }
//   }
//   phi(he, h0arr);
// 
//   KE = 0.;
//   for (int d=0; d<ndims; d++)
//     {
//       KE = KE + he(d) * (U(d) * v(d,k+ks,j+js,i+is));
//     }
//   return 1./2. * KE;
// 
// }
// 
// 
//  void YAKL_INLINE compute_dKdv(real4d F, real4d K, real4d HE, const real4d v, const real4d U, const real4d dens0, int is, int js, int ks, int i, int j, int k)
// {
//   SArray<real,ndims,2> D0;
//   SArray<real,ndims> he;
// 
//   //compute he = phi * h0
//   for (int d=0; d<ndims; d++)
//   {
//     if (d == 0)
//     {
//       D0(d,0) = compute_D(dens0, is, js, ks, i, j, k);
//       D0(d,1) = compute_D(dens0, is, js, ks, i-1, j, k);
//     }
//     if (d == 1)
//     {
//       D0(d,0) = compute_D(dens0, is, js, ks, i, j, k);
//       D0(d,1) = compute_D(dens0, is, js, ks, i, j-1, k);
//     }
//   }
//   phi(he, D0);
// 
//   //compute F = he * U, set HE
//   for (int d=0; d<ndims; d++)
//   {
//     HE(d, k+ks, j+js, i+is) = he(d);
//     F(d, k+ks, j+js, i+is) = U(d, k+ks, j+js, i+is) * he(d);
//   }
// 
//   //compute K = 1/2 * PhiT(U,V)
//   compute_phiT(K, U, v, is, js, ks, i, j, k);
//   K(0, k+ks, j+js, i+is) *= 0.5;
// 
// }
// 
// // Note that this ADDS to Bvar...
// void YAKL_INLINE compute_dKddens(real4d B, const real4d K, int is, int js, int ks, int i, int j, int k)
// {
//   SArray<real,1> K0;
//   compute_I<1, diff_ord>(K0, K, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k);
// for(int l=0;l<ndensity;l++)
//{
  //   B(l, k+ks, j+js, i+is) += compute_partialDpartialDl(l) * K0(0);
//}
// }
// 
//  };

// class HamiltonianSWE: public SWE_HS, public Hk
// {
//   public:
// };
#endif
