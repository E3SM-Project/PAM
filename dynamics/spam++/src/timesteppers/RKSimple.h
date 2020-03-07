

#ifndef _RKSIMPLE_H_
#define _RKSIMPLE_H_


#include "common.h"
#include "topology.h"
#include "variable_sets.h"
#include "exchange.h"

// EVENTUALLY THIS NEEDS TO BE SWAPPABLE...
#include "advection.h"


// THIS SHOULD LIVE IN A GENERAL timesteppers.h class
// Time scheme types
enum class TIME_TYPE { KGRK, ADER };


template<uint nstages> void set_stage_coefficients(TIME_TYPE rksimple_type, SArray<real, nstages> &stage_coeffs)
{
  if (rksimple_type == TIME_TYPE::KGRK && nstages == 4)
  {
  stage_coeffs(0) = 1./4.;
  stage_coeffs(1) = 1./3.;
  stage_coeffs(2) = 1./2.;
  stage_coeffs(3) = 1.;
  }
}

template<uint ndims, uint nprog, uint nconst, uint ndiag, uint nstages> class RKSimpleTimeIntegrator {

public:

  SArray<real, nstages> stage_coeffs;
  VariableSet<ndims, nprog> xtend;
  VariableSet<ndims, nprog> xtemp;
  VariableSet<ndims, nprog> *x;
  Tendencies<ndims, nprog, nconst, ndiag> *tendencies;
  const VariableSet<ndims, nconst> *const_vars;
  VariableSet<ndims, ndiag> *diagnostic_vars;
  ExchangeSet<ndims, nprog> *x_exchange;

  bool is_initialized;
  RKSimpleTimeIntegrator();
  RKSimpleTimeIntegrator( const RKSimpleTimeIntegrator<ndims,nprog,nconst,ndiag,nstages> &rksimple) = delete;
  RKSimpleTimeIntegrator& operator=( const RKSimpleTimeIntegrator<ndims,nprog,nconst,ndiag,nstages> &rksimple) = delete;
  void initialize(Tendencies<ndims, nprog, nconst, ndiag> &tend, VariableSet<ndims, nprog> &xvars, const VariableSet<ndims, nconst> &consts, VariableSet<ndims, ndiag> &diagnostics, ExchangeSet<ndims, nprog> &prog_exch);
  void stepForward(real dt);


};


    template<uint ndims, uint nprog, uint nconst, uint ndiag, uint nstages> RKSimpleTimeIntegrator<ndims,nprog,nconst,ndiag,nstages>::RKSimpleTimeIntegrator()
    {
      this->is_initialized = false;
      std::cout << "CREATED RKSIMPLE\n";
    }


  template<uint ndims, uint nprog, uint nconst, uint ndiag, uint nstages> void RKSimpleTimeIntegrator<ndims,nprog,nconst,ndiag,nstages>::initialize(Tendencies<ndims, nprog, nconst, ndiag> &tend, VariableSet<ndims, nprog> &xvars, const VariableSet<ndims, nconst> &consts, VariableSet<ndims, ndiag> &diagnostics, ExchangeSet<ndims, nprog> &prog_exch)
  {
    this->xtemp.initialize(xvars, "xtemp");
    this->xtend.initialize(xvars, "xtend");
    this->x = &xvars;
    this->tendencies = &tend;
    this->const_vars = &consts;
    this->diagnostic_vars = &diagnostics;
    this->x_exchange = &prog_exch;
    this->is_initialized = true;
  }

  template<uint ndims, uint nprog, uint nconst, uint ndiag, uint nstages> void RKSimpleTimeIntegrator<ndims,nprog,nconst,ndiag,nstages>::stepForward(real dt)
  {

    this->tendencies->compute_rhs(*this->const_vars, *this->x, *this->diagnostic_vars, this->xtend);
    this->xtemp.waxpy(dt * this->stage_coeffs(0), this->xtend, *this->x);

    for (int i=1; i<nstages; i++)
    {
      this->x_exchange->exchange_variable_set(this->xtemp);
      this->tendencies->compute_rhs(*this->const_vars, this->xtemp, *this->diagnostic_vars, this->xtend);
      this->xtemp.waxpy(dt * this->stage_coeffs(i), this->xtend, *this->x);
    }
    // THIS COPY CAN BE AVOIDED IF WE ARE CLEVER ie swap x and xtemp
    // Would require being careful with IO also?
    this->x->copy(this->xtemp);
    this->x_exchange->exchange_variable_set(*this->x);
  }

#endif
