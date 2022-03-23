
#ifndef _FUNCTIONALS_H_
#define _FUNCTIONALS_H_

#include "common.h"

struct pvpe {
  real pv=0., pe=0.;
};

class Functional_PVPE {
public:
  bool is_initialized;

   Functional_PVPE() {
     this->is_initialized = false;
}

void initialize(ModelParameters &params)
{
  this->is_initialized = true;
}

real YAKL_INLINE compute_hv(const real4d dens, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,1> hv;
  SArray<real,1,4> Dv;

  // compute hv = R h
  Dv(0) = dens(0, k+ks, j+js, i+is);
  Dv(1) = dens(0, k+ks, j+js, i+is-1);
  Dv(2) = dens(0, k+ks, j+js-1, i+is);
  Dv(3) = dens(0, k+ks, j+js-1, i+is-1);
  R(hv, Dv);

  return hv(0);
}

real YAKL_INLINE compute_zeta(const real4d v, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,1> zeta;
  // compute zeta = D2 v
  compute_D2<1>(zeta, v, is, js, ks, i, j, k);
  return zeta(0);
}

real YAKL_INLINE compute_eta(const real4d v, const real4d coriolis, int is, int js, int ks, int i, int j, int k)
{
  real zeta = compute_zeta(v, is, js, ks, i, j, k);
  return zeta + coriolis(0, k+ks, j+js, i+is);
}

//This computes relative q0
void YAKL_INLINE compute_q0f0(real4d q0, real4d f0, const real4d v, const real4d dens, const real4d coriolis, int is, int js, int ks, int i, int j, int k)
{

  real hv = compute_hv(dens, is, js, ks, i, j, k);
  real zeta = compute_zeta(v, is, js, ks, i, j, k);

  // compute q0 = zeta / hv and f0 = f / hv
    q0(0, k+ks, j+js, i+is) = zeta / hv;
    f0(0, k+ks, j+js, i+is) = coriolis(0, k+ks, j+js, i+is) / hv;

}


//This computes TRUE q0
void YAKL_INLINE compute_q0(real4d q0, const real4d v, const real4d dens, const real4d coriolis, int is, int js, int ks, int i, int j, int k)
{
real hv = compute_hv(dens, is, js, ks, i, j, k);
real eta = compute_eta(v, coriolis, is, js, ks, i, j, k);
// compute q0 = zeta / hv and f0 = f / hv
  q0(0, k+ks, j+js, i+is) = eta / hv;
}

pvpe YAKL_INLINE compute_PVPE(const real4d v, const real4d dens, const real4d coriolis, int is, int js, int ks, int i, int j, int k)
{
  pvpe vals;
  real eta = compute_eta(v, coriolis, is, js, ks, i, j, k);
  real hv = compute_hv(dens, is, js, ks, i, j, k);
  real q0 = eta / hv;

  vals.pv = eta;
  vals.pe = 0.5 * eta * q0;
  
  return vals;
}

};


// 
// class Functional_PVPE_rhod {
// public:
//   bool is_initialized;
// 
//    Functional_PVPE_rhod() {
//      this->is_initialized = false;
// }
// 
// void initialize(ModelParameters &params)
// {
//   this->is_initialized = true;
// }
// 
// real YAKL_INLINE compute_hv(const real4d dens, const real4d densfct, int is, int js, int ks, int i, int j, int k)
// {
//   SArray<real,1> hv;
//   SArray<real,4> Dv;
// 
//   // compute hv = R h
//   // Uses linearity of R
//   Dv(0) = dens(0, k+ks, j+js, i+is) + densfct(0, k+ks, j+js, i+is) + densfct(1, k+ks, j+js, i+is) + densfct(2, k+ks, j+js, i+is);
//   Dv(1) = dens(0, k+ks, j+js, i+is-1) + densfct(0, k+ks, j+js, i+is-1) + densfct(1, k+ks, j+js, i+is-1) + densfct(2, k+ks, j+js, i+is-1);
//   Dv(2) = dens(0, k+ks, j+js-1, i+is) + densfct(0, k+ks, j+js-1, i+is) + densfct(1, k+ks, j+js-1, i+is) + densfct(2, k+ks, j+js-1, i+is);
//   Dv(3) = dens(0, k+ks, j+js-1, i+is-1) + densfct(0, k+ks, j+js-1, i+is-1) + densfct(1, k+ks, j+js-1, i+is-1) + densfct(2, k+ks, j+js-1, i+is-1);
//   R(hv, Dv);
// 
//   return hv(0);
// }
// 
// real YAKL_INLINE compute_zeta(const real4d v, int is, int js, int ks, int i, int j, int k)
// {
//   SArray<real,1> zeta;
//   // compute zeta = D2 v
//   compute_D2<1>(zeta, v, is, js, ks, i, j, k);
//   return zeta(0);
// }
// 
// real YAKL_INLINE compute_eta(const real4d v, const real4d coriolis, int is, int js, int ks, int i, int j, int k)
// {
//   real zeta = compute_zeta(v, is, js, ks, i, j, k);
//   return zeta + coriolis(0, k+ks, j+js, i+is);
// }
// 
// 
// void YAKL_INLINE compute_q0f0(real4d q0, real4d f0, const real4d v, const real4d dens, const real4d densfct, const real4d coriolis, int is, int js, int ks, int i, int j, int k)
// {
// 
//   real hv = compute_hv(dens, densfct, is, js, ks, i, j, k);
//   real zeta = compute_zeta(v, is, js, ks, i, j, k);
// 
//   // compute q0 = zeta / hv and f0 = f / hv
//     q0(0, k+ks, j+js, i+is) = zeta / hv;
//     f0(0, k+ks, j+js, i+is) = coriolis(0, k+ks, j+js, i+is) / hv;
// 
// }
// 
// void YAKL_INLINE compute_q0(real4d q0, const real4d v, const real4d dens, const real4d densfct, int is, int js, int ks, int i, int j, int k)
// {
// real hv = compute_hv(dens, densfct, is, js, ks, i, j, k);
// real zeta = compute_zeta(v, is, js, ks, i, j, k);
// // compute q0 = zeta / hv and f0 = f / hv
//   q0(0, k+ks, j+js, i+is) = zeta / hv;
// }
// 
// pvpe YAKL_INLINE compute_PVPE(const real4d v, const real4d dens, const real4d densfct, const real4d coriolis, int is, int js, int ks, int i, int j, int k)
// {
//   pvpe vals;
//   real eta = compute_eta(v, coriolis, is, js, ks, i, j, k);
//   real hv = compute_hv(dens, densfct, is, js, ks, i, j, k);
//   real q0 = eta / hv;
// 
//   vals.pv = eta;
//   vals.pe = 0.5 * eta * q0;
// 
//   return vals;
// }
// 
// };
// 
// 













class Functional_PVPE_extruded {

public:
  bool is_initialized;

   Functional_PVPE_extruded() {
     this->is_initialized = false;
}

void initialize(ModelParameters &params)
{
  this->is_initialized = true;
}

real YAKL_INLINE compute_hvxz(const real4d dens, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,1> hv;
  SArray<real,1,4> Dv;
  
  // compute hv = R h
  // Uses linearity of R
  Dv(0) = dens(0, k+ks, j+js, i+is);
  Dv(1) = dens(0, k+ks, j+js, i+is-1);
  Dv(2) = dens(0, k+ks+1, j+js, i+is);
  Dv(3) = dens(0, k+ks+1, j+js, i+is-1);
  R(hv, Dv);
  
  return hv(0);
}

real YAKL_INLINE compute_zetaxz(const real4d v, const real4d w, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,1,1> zeta;
  // compute zeta = Dxz "v"
  compute_Dxz<1> (zeta, v, w, is, js, ks, i, j, k);
  return zeta(0);
}

real YAKL_INLINE compute_etaxz(const real4d v, const real4d w, const real4d coriolisxz, int is, int js, int ks, int i, int j, int k)
{
  real zeta = compute_zetaxz(v, w, is, js, ks, i, j, k);
  return zeta + coriolisxz(0, k+ks, j+js, i+is);
}

//This computes true qxz
void YAKL_INLINE compute_qxz0(real4d qxz0, const real4d v, const real4d w, const real4d dens, const real4d coriolisxz, int is, int js, int ks, int i, int j, int k)
{
  //Need to subtract 1 here since d00(i,k) corresponds to p11(i,k)
 real hv = compute_hvxz(dens, is, js, ks, i, j, k-1);
 real eta = compute_etaxz(v, w, coriolisxz, is, js, ks, i, j, k-1);
 // compute q0 = zeta / hv and f0 = f / hv
qxz0(0, k+ks, j+js, i+is) = eta / hv;
}

//This computes relative qxz
void YAKL_INLINE compute_qxz0fxz0(real4d qxz0, real4d fxz0, const real4d v, const real4d w, const real4d dens, const real4d coriolisxz, int is, int js, int ks, int i, int j, int k)
{
  //Need to subtract 1 here since d00(i,k) corresponds to p11(i,k)
 real hv = compute_hvxz(dens, is, js, ks, i, j, k-1);
 real zeta = compute_zetaxz(v, w, is, js, ks, i, j, k-1);
 // compute q0 = zeta / hv and f0 = f / hv
qxz0(0, k+ks, j+js, i+is) = zeta / hv;
fxz0(0, k+ks, j+js, i+is) = coriolisxz(0, k+ks, j+js, i+is) / hv;
}

pvpe YAKL_INLINE compute_PVPE(const real4d v, const real4d w, const real4d dens, const real4d coriolisxz, int is, int js, int ks, int i, int j, int k)
{
  pvpe vals;
  //No subtraction here since this is called on primal cells p11
  real eta = compute_etaxz(v, w, coriolisxz, is, js, ks, i, j, k);
  real hv = compute_hvxz(dens, is, js, ks, i, j, k);
  real q0 = eta / hv;
  
  vals.pv = eta;
  vals.pe = 0.5 * eta * q0;
  
  return vals;
}

};

#endif
