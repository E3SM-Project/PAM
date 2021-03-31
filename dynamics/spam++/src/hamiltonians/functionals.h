
#ifndef _FUNCTIONALS_H_
#define _FUNCTIONALS_H_

#include "common.h"

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

#endif
