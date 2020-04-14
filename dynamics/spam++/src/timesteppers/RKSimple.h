

#ifndef _RKSIMPLE_H_
#define _RKSIMPLE_H_


#include "common.h"
#include "topology.h"
#include "variable_sets.h"
#include "exchange.h"
#include "model.h"




template<uint ndims, uint nprog, uint nconst, uint naux, uint nstages> class RKSimpleTimeIntegrator {

public:

  SArray<real, nstages> stage_coeffs;
  VariableSet<ndims, nprog> xtend;
  VariableSet<ndims, nprog> xtemp;
  VariableSet<ndims, nprog> *x;
  Tendencies<ndims, nprog, nconst, naux> *tendencies;
  const VariableSet<ndims, nconst> *const_vars;
  VariableSet<ndims, naux> *auxiliary_vars;
  ExchangeSet<ndims, nprog> *x_exchange;

  bool is_initialized;
  RKSimpleTimeIntegrator();
  RKSimpleTimeIntegrator( const RKSimpleTimeIntegrator<ndims,nprog,nconst,naux,nstages> &rksimple) = delete;
  RKSimpleTimeIntegrator& operator=( const RKSimpleTimeIntegrator<ndims,nprog,nconst,naux,nstages> &rksimple) = delete;
  void initialize(Tendencies<ndims, nprog, nconst, naux> &tend, VariableSet<ndims, nprog> &xvars, const VariableSet<ndims, nconst> &consts, VariableSet<ndims, naux> &auxiliarys, ExchangeSet<ndims, nprog> &prog_exch);
  void stepForward(real dt);
  void set_stage_coefficients();

};


    template<uint ndims, uint nprog, uint nconst, uint naux, uint nstages> RKSimpleTimeIntegrator<ndims,nprog,nconst,naux,nstages>::RKSimpleTimeIntegrator()
    {
      this->is_initialized = false;
      std::cout << "CREATED RKSIMPLE\n";
    }

    // THIS IS IN FACT SPECIFIC TO KG RK
    // MUST GENERALIZE SOMEHOW...
    template<uint ndims, uint nprog, uint nconst, uint naux, uint nstages> void RKSimpleTimeIntegrator<ndims,nprog,nconst,naux,nstages>::set_stage_coefficients()
      {
        if (nstages == 4)
        {
        this->stage_coeffs(0) = 1./4.;
        this->stage_coeffs(1) = 1./3.;
        this->stage_coeffs(2) = 1./2.;
        this->stage_coeffs(3) = 1.;
        }
      }

  template<uint ndims, uint nprog, uint nconst, uint naux, uint nstages> void RKSimpleTimeIntegrator<ndims,nprog,nconst,naux,nstages>::initialize(Tendencies<ndims, nprog, nconst, naux> &tend, VariableSet<ndims, nprog> &xvars, const VariableSet<ndims, nconst> &consts, VariableSet<ndims, naux> &auxiliarys, ExchangeSet<ndims, nprog> &prog_exch)
  {
    this->xtemp.initialize(xvars, "xtemp");
    this->xtend.initialize(xvars, "xtend");
    this->x = &xvars;
    this->tendencies = &tend;
    this->const_vars = &consts;
    this->auxiliary_vars = &auxiliarys;
    this->x_exchange = &prog_exch;
    set_stage_coefficients();
    this->is_initialized = true;
  }




  template<uint ndims, uint nprog, uint nconst, uint naux, uint nstages> void RKSimpleTimeIntegrator<ndims,nprog,nconst,naux,nstages>::stepForward(real dt)
  {

    this->tendencies->compute_rhs(dt, *this->const_vars, *this->x, *this->auxiliary_vars, this->xtend);
    this->xtemp.waxpy(-1.*dt * this->stage_coeffs(0), this->xtend, *this->x);

    for (int i=1; i<nstages; i++)
    {
      this->x_exchange->exchange_variable_set(this->xtemp);
      this->tendencies->compute_rhs(dt, *this->const_vars, this->xtemp, *this->auxiliary_vars, this->xtend);
      this->xtemp.waxpy(-1.*dt * this->stage_coeffs(i), this->xtend, *this->x);
    }
    // THIS COPY CAN BE AVOIDED IF WE ARE CLEVER ie swap x and xtemp
    // Would require being careful with IO also?
    this->x->copy(this->xtemp);
    this->x_exchange->exchange_variable_set(*this->x);
  }

#endif
