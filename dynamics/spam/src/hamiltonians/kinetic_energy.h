#pragma once

#include "common.h"
#include "hodge_star.h"
#include "variableset.h"
#include "wedge.h"

// // This kinetic energy functional assumes that dens(0) holds the total
// density class Hamiltonian_Hk1D_rho {
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
// void initialize(Parameters &params, Geometry<ndims,1,1,1> &primal_geom,
// Geometry<ndims,1,1,1> &dual_geom)
// {
//   this->primal_geometry = &primal_geom;
//   this->dual_geometry = &dual_geom;
//   this->is_initialized = true;
// }
//
// real YAKL_INLINE compute_KE(const real5d v, const real5d dens, const real5d
// densfct, int is, int js, int ks, int i, int j, int k)
// {
//
//   real KE, PE, IE;
//   SArray<real,1> h0, h0im1;
//   SArray<real,1> U, he;
//   SArray<real,1,2> h0arr;
//
//   //compute U = H v
//   compute_H<1,diff_ord> (U, v, *this->primal_geometry, *this->dual_geometry,
//   is, js, ks, i, j, k);
//
//   // Compute h0 = I h needed for phi calcs
//   compute_I<1,diff_ord> (h0, dens, *this->primal_geometry,
//   *this->dual_geometry, is, js, ks, i, j, k); compute_I<1,diff_ord> (h0im1,
//   dens, *this->primal_geometry, *this->dual_geometry, is, js, ks, i-1, j, k);
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
//  void YAKL_INLINE compute_dKdv(real5d F, real5d K, real5d HE, const real5d v,
//  const real5d U, const real5d dens0, const real5d densfct0, int is, int js,
//  int ks, int i, int j, int k)
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
// void YAKL_INLINE compute_dKddens(real5d B, real5d Bfct, const real5d K, int
// is, int js, int ks, int i, int j, int k)
// {
//   SArray<real,1> K0;
//   compute_I<1, diff_ord>(K0, K, *this->primal_geometry, *this->dual_geometry,
//   is, js, ks, i, j, k); B(0, k+ks, j+js, i+is) += K0(0);
// }
//
//  };

class Hamiltonian_Hk {

public:
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;
  VariableSet varset;
  bool is_initialized;

  Hamiltonian_Hk() { this->is_initialized = false; }

  void initialize(VariableSet &variableset,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom) {
    this->primal_geometry = primal_geom;
    this->dual_geometry = dual_geom;
    this->is_initialized = true;
    this->varset = variableset;
  }

  real YAKL_INLINE compute_KE(const real5d v, const real5d dens, int is, int js,
                              int ks, int i, int j, int k, int n) const {

    real KE;
    SArray<real, 1, 1> h0, h0im1, h0jm1, h0km1;
    SArray<real, 1, ndims> U, he;
    SArray<real, 2, ndims, 2> h0arr;

    // compute U = H v
    compute_H<1, diff_ord>(U, v, this->primal_geometry, this->dual_geometry, is,
                           js, ks, i, j, k, n);

    // Compute h0 = I h needed for phi calcs
    compute_I<1, diff_ord>(h0, dens, this->primal_geometry, this->dual_geometry,
                           is, js, ks, i, j, k, n);
    compute_I<1, diff_ord>(h0im1, dens, this->primal_geometry,
                           this->dual_geometry, is, js, ks, i - 1, j, k, n);
    if (ndims >= 2) {
      compute_I<1, diff_ord>(h0jm1, dens, this->primal_geometry,
                             this->dual_geometry, is, js, ks, i, j - 1, k, n);
    }
    if (ndims >= 3) {
      compute_I<1, diff_ord>(h0km1, dens, this->primal_geometry,
                             this->dual_geometry, is, js, ks, i, j, k - 1, n);
    }

    // compute he = phi h0
    for (int d = 0; d < ndims; d++) {
      if (d == 0) {
        h0arr(d, 0) = h0(0);
        h0arr(d, 1) = h0im1(0);
      }
      if (d == 1) {
        h0arr(d, 0) = h0(0);
        h0arr(d, 1) = h0jm1(0);
      }
      if (d == 2) {
        h0arr(d, 0) = h0(0);
        h0arr(d, 1) = h0km1(0);
      }
    }
    phi(he, h0arr);

    KE = 0.;
    for (int d = 0; d < ndims; d++) {
      KE = KE + he(d) * (U(d) * v(d, k + ks, j + js, i + is, n));
    }
    return 0.5_fp * KE;
  }

  // FIX THIS TO GET TOTAL DENSITY FROM VARSET!

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void YAKL_INLINE compute_F_and_K(real5d F, real5d K, const real5d v,
                                   const real5d U, const real5d dens0, int is,
                                   int js, int ks, int i, int j, int k, int n,
                                   real fac = 1._fp) const {
    SArray<real, 2, ndims, 2> D0;
    SArray<real, 1, ndims> he;

    // compute he = phi * h0
    for (int d = 0; d < ndims; d++) {
      if (d == 0) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks, j + js, i + is - 1, n);
      }
      if (d == 1) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks, j + js - 1, i + is, n);
      }
      if (d == 2) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks - 1, j + js, i + is, n);
      }
    }
    phi(he, D0);

    // compute F = he * U
    for (int d = 0; d < ndims; d++) {
      if (addmode == ADD_MODE::REPLACE) {
        F(d, k + ks, j + js, i + is, n) =
            fac * U(d, k + ks, j + js, i + is, n) * he(d);
      } else if (addmode == ADD_MODE::ADD) {
        F(d, k + ks, j + js, i + is, n) +=
            fac * U(d, k + ks, j + js, i + is, n) * he(d);
      }
    }

    // compute K = 1/2 * PhiT(U,V)
    compute_phiT(K, U, v, is, js, ks, i, j, k, n);
    K(0, k + ks, j + js, i + is, n) *= 0.5_fp;
  }

  void YAKL_INLINE compute_F_and_he(real5d F, real5d HE, const real5d U,
                                    const real5d dens0, int is, int js, int ks,
                                    int i, int j, int k, int n) const {
    SArray<real, 2, ndims, 2> D0;
    SArray<real, 1, ndims> he;

    // compute he = phi * h0
    for (int d = 0; d < ndims; d++) {
      if (d == 0) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks, j + js, i + is - 1, n);
      }
      if (d == 1) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks, j + js - 1, i + is, n);
      }
      if (d == 2) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks - 1, j + js, i + is, n);
      }
    }
    phi(he, D0);

    // compute F = he * U, set HE
    for (int d = 0; d < ndims; d++) {
      F(d, k + ks, j + js, i + is, n) = U(d, k + ks, j + js, i + is, n) * he(d);
      HE(d, k + ks, j + js, i + is, n) = he(d);
    }
  }

  // FIX THIS TO GET TOTAL DENSITY FROM VARSET!
  //  Note that this ADDS to Bvar...
  void YAKL_INLINE compute_dKddens(real5d B, const real5d K, int is, int js,
                                   int ks, int i, int j, int k, int n,
                                   real fac = 1._fp) const {
    SArray<real, 1, 1> K0;
    compute_I<1, diff_ord>(K0, K, this->primal_geometry, this->dual_geometry,
                           is, js, ks, i, j, k, n);
    B(0, k + ks, j + js, i + is, n) += fac * K0(0);
  }
};

//  // This kinetic energy functional assumes that dens(0) + densfct(0:2) make
//  the total density class Hamiltonian_Hk_rhod {
//
//  public:
//    Geometry<1,1,1> *primal_geometry;
//    Geometry<1,1,1> *dual_geometry;
//    bool is_initialized;
//
//     Hamiltonian_Hk_rhod() {
//       this->is_initialized = false;
//  }
//
//  void initialize(Parameters &params, Geometry<1,1,1> &primal_geom,
//  Geometry<1,1,1> &dual_geom)
//  {
//    this->primal_geometry = &primal_geom;
//    this->dual_geometry = &dual_geom;
//    this->is_initialized = true;
//  }
//
//  real YAKL_INLINE compute_KE(const real5d v, const real5d dens, const real5d
//  densfct, int is, int js, int ks, int i, int j, int k, int n)
//  {
//
//
//    real KE;
//    SArray<real,1,1> rhod0, rhod0im1, rhod0jm1, rhod0km1;
//    SArray<real,1,3> rhos0, rhos0im1, rhos0jm1, rhos0km1;
//    SArray<real,1,ndims> U, he;
//    SArray<real,2,ndims,2> h0arr;
//
//    //compute U = H v
//    compute_H<1,diff_ord> (U, v, *this->primal_geometry, *this->dual_geometry,
//    is, js, ks, i, j, k, n);
//
//    // Compute h0 = I h needed for phi calcs
//    compute_I<1,diff_ord> (rhod0, dens, *this->primal_geometry,
//    *this->dual_geometry, is, js, ks, i, j, k, n); compute_I<1,diff_ord>
//    (rhod0im1, dens, *this->primal_geometry, *this->dual_geometry, is, js, ks,
//    i-1, j, k, n); if (ndims>=2) {compute_I<1,diff_ord> (rhod0jm1, dens,
//    *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j-1, k, n);}
//    if (ndims>=3) {compute_I<1,diff_ord> (rhod0km1, dens,
//    *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k-1, n);}
//    compute_I<3,diff_ord> (rhos0im1, densfct, *this->primal_geometry,
//    *this->dual_geometry, is, js, ks, i-1, j, k, n); compute_I<3,diff_ord>
//    (rhos0im1, densfct, *this->primal_geometry, *this->dual_geometry, is, js,
//    ks, i-1, j, k, n); if (ndims>=2) {compute_I<3,diff_ord> (rhos0jm1,
//    densfct, *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j-1,
//    k, n);} if (ndims>=3) {compute_I<3,diff_ord> (rhos0km1, densfct,
//    *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k-1, n);}
//
//    //compute he = phi h0
//    // exploits linearity of I/Phi
//    for (int d=0; d<ndims; d++)
//    {
//      if (d == 0)
//      {
//        h0arr(d,0) = rhod0(0) + rhos0(0) + rhos0(1) + rhos0(2);
//        h0arr(d,1) = rhod0im1(0) + rhos0im1(0) + rhos0im1(1) + rhos0im1(2);
//      }
//      if (d == 1)
//      {
//        h0arr(d,0) = rhod0(0) + rhos0(0) + rhos0(1) + rhos0(2);
//        h0arr(d,1) = rhod0jm1(0) + rhos0jm1(0) + rhos0jm1(1) + rhos0jm1(2);
//      }
//      if (d == 2)
//      {
//        h0arr(d,0) = rhod0(0) + rhos0(0) + rhos0(1) + rhos0(2);
//        h0arr(d,1) = rhod0km1(0) + rhos0km1(0) + rhos0km1(1) + rhos0km1(2);
//      }
//    }
//    phi(he, h0arr);
//
//    KE = 0.;
//    for (int d=0; d<ndims; d++)
//      {
//        KE = KE + he(d) * (U(d) * v(d,k+ks,j+js,i+is, n));
//      }
//    return 1./2. * KE;
//
//  }
//
//
//   void YAKL_INLINE compute_dKdv(real5d F, real5d K, real5d HE, const real5d
//   v, const real5d U, const real5d dens0, const real5d densfct0, int is, int
//   js, int ks, int i, int j, int k, int n)
//  {
//
//    SArray<real,2,ndims,2> D0;
//    SArray<real,1,ndims> he;
//
//    //compute he = phi * h0
//    // exploits linearity of I/phi
//    for (int d=0; d<ndims; d++)
//    {
//      if (d == 0)
//      {
//        D0(d,0) = dens0(0, k+ks, j+js, i+is, n) + densfct0(0, k+ks, j+js,
//        i+is, n) + densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js,
//        i+is); D0(d,1) = dens0(0, k+ks, j+js, i+is-1, n) + densfct0(0, k+ks,
//        j+js, i+is-1, n) + densfct0(1, k+ks, j+js, i+is-1) + densfct0(2, k+ks,
//        j+js, i+is-1);
//      }
//      if (d == 1)
//      {
//        D0(d,0) = dens0(0, k+ks, j+js, i+is) + densfct0(0, k+ks, j+js, i+is) +
//        densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js, i+is); D0(d,1)
//        = dens0(0, k+ks, j+js-1, i+is) + densfct0(0, k+ks, j+js-1, i+is) +
//        densfct0(1, k+ks, j+js-1, i+is) + densfct0(2, k+ks, j+js-1, i+is);
//      }
//      if (d == 2)
//      {
//        D0(d,0) = dens0(0, k+ks, j+js, i+is) + densfct0(0, k+ks, j+js, i+is) +
//        densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js, i+is); D0(d,1)
//        = dens0(0, k+ks-1, j+js, i+is) + densfct0(0, k+ks-1, j+js, i+is) +
//        densfct0(1, k+ks-1, j+js, i+is) + densfct0(2, k+ks-1, j+js, i+is);
//      }
//    }
//    phi(he, D0);
//
//    //compute F = he * U, set HE
//    for (int d=0; d<ndims; d++)
//    {
//      HE(d, k+ks, j+js, i+is) = he(d);
//      F(d, k+ks, j+js, i+is) = U(d, k+ks, j+js, i+is) * he(d);
//    }
//
//    //compute K = 1/2 * PhiT(U,V)
//    compute_phiT(K, U, v, is, js, ks, i, j, k);
//    K(0, k+ks, j+js, i+is) *= 0.5;
//
//  }
//
//  // Note that this ADDS to Bvar/Bfctvar...
//  void YAKL_INLINE compute_dKddens(real5d B, real5d Bfct, const real5d K, int
//  is, int js, int ks, int i, int j, int k, int n)
//  {
//    SArray<real,1,1> K0;
//    compute_I<1, diff_ord>(K0, K, *this->primal_geometry,
//    *this->dual_geometry, is, js, ks, i, j, k); B(0, k+ks, j+js, i+is) +=
//    K0(0); Bfct(0, k+ks, j+js, i+is) += K0(0); Bfct(1, k+ks, j+js, i+is) +=
//    K0(0); Bfct(2, k+ks, j+js, i+is) += K0(0);
//  }
//
//   };
//
//
//
//
//
//
//
//
//
class Hamiltonian_Hk_extruded {

public:
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;
  VariableSet varset;
  bool is_initialized;

  Hamiltonian_Hk_extruded() { this->is_initialized = false; }

  void initialize(VariableSet &variableset,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom) {
    this->primal_geometry = primal_geom;
    this->dual_geometry = dual_geom;
    this->is_initialized = true;
    this->varset = variableset;
  }

  real YAKL_INLINE compute_KE_top(const real5d v, const real5d w,
                                  const real5d dens, int is, int js, int ks,
                                  int i, int j, int k, int n) const {
    real K2 = 0.;
    // Have to subtract 1 from k here since UW has an extra dof compared to w!
    SArray<real, 1, 1> UW0;
    compute_Hv<1, vert_diff_ord>(UW0, w, this->primal_geometry,
                                 this->dual_geometry, is, js, ks, i, j, k - 1,
                                 n);
    real w0;
    // Have to subtract 1 from k here since UW has an extra dof compared to w
    w0 = w(0, k + ks - 1, j + js, i + is, n);
    K2 += 0.5 * w0 * UW0(0);
    return _compute_KE(K2, v, w, dens, is, js, ks, i, j, k, n);
  }

  real YAKL_INLINE compute_KE_bottom(const real5d v, const real5d w,
                                     const real5d dens, int is, int js, int ks,
                                     int i, int j, int k, int n) const {
    real K2 = 0.;
    // Have to subtract 1 from k here since UW has an extra dof compared to w!
    SArray<real, 1, 1> UW1;
    compute_Hv<1, vert_diff_ord>(UW1, w, this->primal_geometry,
                                 this->dual_geometry, is, js, ks, i, j, k, n);
    real w1;
    // Have to subtract 1 from k here since UW has an extra dof compared to w
    w1 = w(0, k + ks, j + js, i + is, n);
    K2 += 0.5 * w1 * UW1(0);
    return _compute_KE(K2, v, w, dens, is, js, ks, i, j, k, n);
  }

  real YAKL_INLINE compute_KE(const real5d v, const real5d w, const real5d dens,
                              int is, int js, int ks, int i, int j, int k,
                              int n) const {
    real K2 = 0.;
    // Have to subtract 1 from k here since UW has an extra dof compared to w!
    SArray<real, 1, 1> UW0, UW1;
    compute_Hv<1, vert_diff_ord>(UW0, w, this->primal_geometry,
                                 this->dual_geometry, is, js, ks, i, j, k - 1,
                                 n);
    compute_Hv<1, vert_diff_ord>(UW1, w, this->primal_geometry,
                                 this->dual_geometry, is, js, ks, i, j, k, n);
    real w0, w1;
    // Have to subtract 1 from k here since UW has an extra dof compared to w
    w0 = w(0, k + ks - 1, j + js, i + is, n);
    w1 = w(0, k + ks, j + js, i + is, n);
    K2 += 0.5 * (w0 * UW0(0) + w1 * UW1(0));
    return _compute_KE(K2, v, w, dens, is, js, ks, i, j, k, n);
  }

  real YAKL_INLINE _compute_KE(real K2, const real5d v, const real5d w,
                               const real5d dens, int is, int js, int ks, int i,
                               int j, int k, int n) const {
    real v0, v1;
    SArray<real, 1, ndims> U0, U1;
    compute_Hext<1, diff_ord>(U0, v, this->primal_geometry, this->dual_geometry,
                              is, js, ks, i, j, k, n);
    compute_Hext<1, diff_ord>(U1, v, this->primal_geometry, this->dual_geometry,
                              is, js, ks, i + 1, j, k, n);
    v0 = v(0, k + ks, j + js, i + is, n);
    v1 = v(0, k + ks, j + js, i + is + 1, n);
    K2 += 0.5 * (v0 * U0(0) + v1 * U1(0));

    if (ndims == 2) {
      real v0, v1;
      SArray<real, 1, ndims> U0, U1;
      compute_Hext<1, diff_ord>(U0, v, this->primal_geometry,
                                this->dual_geometry, is, js, ks, i, j, k, n);
      compute_Hext<1, diff_ord>(U1, v, this->primal_geometry,
                                this->dual_geometry, is, js, ks, i, j + 1, k,
                                n);
      v0 = v(1, k + ks, j + js, i + is, n);
      v1 = v(1, k + ks, j + js + 1, i + is, n);
      K2 += 0.5 * (v0 * U0(1) + v1 * U1(1));
    }

    K2 *= 0.5;

    // Compute h0 = I h needed for phi calcs
    SArray<real, 1, 1> h0;
    compute_Iext<1, diff_ord, vert_diff_ord>(h0, dens, this->primal_geometry,
                                             this->dual_geometry, is, js, ks, i,
                                             j, k, n);

    return h0(0) * K2;
  }

  // FIX THIS TO GET TOTAL DENSITY FROM VARSET!
  void YAKL_INLINE compute_Fw(real5d FW, real5d HEw, const real5d UW,
                              const real5d dens0, int is, int js, int ks, int i,
                              int j, int k, int n) const {
    // SArray<real,2> Dv;
    // compute hew = phiw * h0
    // Dv(0) = dens0(0, k+ks, j+js, i+is);
    // Dv(1) = dens0(0, k+ks-1, j+js, i+is);
    // real hew = phiW(Dv);
    // compute FW = hew * UW, set HEw
    real hew = (dens0(0, k + ks, j + js, i + is, n) +
                dens0(0, k + ks - 1, j + js, i + is, n)) /
               2.0;
    HEw(0, k + ks, j + js, i + is, n) = hew;
    FW(0, k + ks, j + js, i + is, n) = UW(0, k + ks, j + js, i + is, n) * hew;
    // std::cout << "HEw in Hk " << i << " " << j << " " << k << " " <<
    // HEw(0,k+ks,j+js,i+is) << "\n" << std::flush;
  }

  void YAKL_INLINE compute_K(real5d K, const real5d v, const real5d U,
                             const real5d w, const real5d UW, int is, int js,
                             int ks, int i, int j, int k, int n) const {
    // compute K = 1/2 * PhiT(U,V) + 1/2 * PhiTW(UW,W)
    real K2 = 0.;

    real w0, w1, UW0, UW1;
    UW0 = UW(0, k + ks, j + js, i + is, n);
    UW1 = UW(0, k + ks + 1, j + js, i + is, n);
    // Have to subtract 1 from k here since UW has an extra dof compared to w
    w0 = w(0, k + ks - 1, j + js, i + is, n);
    w1 = w(0, k + ks, j + js, i + is, n);
    K2 += 0.5 * (w0 * UW0 + w1 * UW1);

    real v0, U0, v1, U1;
    U0 = U(0, k + ks, j + js, i + is, n);
    U1 = U(0, k + ks, j + js, i + is + 1, n);
    v0 = v(0, k + ks, j + js, i + is, n);
    v1 = v(0, k + ks, j + js, i + is + 1, n);
    K2 += 0.5 * (v0 * U0 + v1 * U1);

    if (ndims == 2) {
      real v0, U0, v1, U1;
      U0 = U(1, k + ks, j + js, i + is, n);
      U1 = U(1, k + ks, j + js + 1, i + is, n);
      v0 = v(1, k + ks, j + js, i + is, n);
      v1 = v(1, k + ks, j + js + 1, i + is, n);
      K2 += 0.5 * (v0 * U0 + v1 * U1);
    }

    K(0, k + ks, j + js, i + is, n) = 0.5 * K2;
    // compute_phiT(K, U, v, is, js, ks, i, j, k);
    // compute_phiTW<ADD_MODE::ADD>(K, UW, w, is, js, ks, i, j, k);
    // K(0, k+ks, j+js, i+is) *= 0.5;
  }

  // FIX THIS TO GET TOTAL DENSITY FROM VARSET!
  void YAKL_INLINE compute_F(real5d F, real5d HE, const real5d U,
                             const real5d dens0, int is, int js, int ks, int i,
                             int j, int k, int n) const {
    SArray<real, 2, ndims, 2> D0;
    SArray<real, 1, ndims> he;

    // compute he = phi * h0
    for (int d = 0; d < ndims; d++) {
      if (d == 0) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks, j + js, i + is - 1, n);
      }
      if (d == 1) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks, j + js - 1, i + is, n);
      }
    }
    phi(he, D0);
    // compute F = he * U, set HE
    for (int d = 0; d < ndims; d++) {
      HE(d, k + ks, j + js, i + is, n) = he(d);
      F(d, k + ks, j + js, i + is, n) = U(d, k + ks, j + js, i + is, n) * he(d);
    }
    // std::cout << "HE in Hk " << i << " " << j << " " << k << " " <<
    // HE(0,k+ks,j+js,i+is) << "\n" << std::flush;
  }

  // FIX THIS TO GET TOTAL DENSITY FROM VARSET!
  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void YAKL_INLINE compute_dKddens(real5d B, const real5d K, int is, int js,
                                   int ks, int i, int j, int k, int n) const {
    SArray<real, 1, 1> K0;
    compute_Iext<1, diff_ord, vert_diff_ord>(K0, K, this->primal_geometry,
                                             this->dual_geometry, is, js, ks, i,
                                             j, k, n);
    if (addmode == ADD_MODE::REPLACE) {
      B(0, k + ks, j + js, i + is, n) = K0(0);
    }
    if (addmode == ADD_MODE::ADD) {
      B(0, k + ks, j + js, i + is, n) += K0(0);
    }
  }
};
