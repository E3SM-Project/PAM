

#ifndef _RKSIMPLE_H_
#define _RKSIMPLE_H_


#include "common.h"
#include "topology.h"
#include "variable_sets.h"
#include "exchange.h"
#include "model.h"




template<uint nprog, uint nconst, uint naux, uint nstages> class RKSimpleTimeIntegrator {

public:

  SArray<real, nstages> stage_coeffs;
  VariableSet<nprog> xtend;
  VariableSet<nprog> xtemp;
  VariableSet<nprog> *x;
  Tendencies<nprog, nconst, naux> *tendencies;
  VariableSet<nconst> *const_vars;
  VariableSet<naux> *auxiliary_vars;
  ExchangeSet<nprog> *x_exchange;

  bool is_initialized;
  RKSimpleTimeIntegrator();
  RKSimpleTimeIntegrator( const RKSimpleTimeIntegrator<nprog,nconst,naux,nstages> &rksimple) = delete;
  RKSimpleTimeIntegrator& operator=( const RKSimpleTimeIntegrator<nprog,nconst,naux,nstages> &rksimple) = delete;
  void initialize(Tendencies<nprog, nconst, naux> &tend, VariableSet<nprog> &xvars, VariableSet<nconst> &consts, VariableSet<naux> &auxiliarys, ExchangeSet<nprog> &prog_exch);
  void stepForward(real dt);
  void set_stage_coefficients();

};


    template<uint nprog, uint nconst, uint naux, uint nstages> RKSimpleTimeIntegrator<nprog,nconst,naux,nstages>::RKSimpleTimeIntegrator()
    {
      this->is_initialized = false;
      std::cout << "CREATED RKSIMPLE\n";
    }

    // THIS IS IN FACT SPECIFIC TO KG RK
    // MUST GENERALIZE SOMEHOW...
    template<uint nprog, uint nconst, uint naux, uint nstages> void RKSimpleTimeIntegrator<nprog,nconst,naux,nstages>::set_stage_coefficients()
      {
        if (nstages == 4)
        {
        this->stage_coeffs(0) = 1./4.;
        this->stage_coeffs(1) = 1./3.;
        this->stage_coeffs(2) = 1./2.;
        this->stage_coeffs(3) = 1.;
        }
      }

  template<uint nprog, uint nconst, uint naux, uint nstages> void RKSimpleTimeIntegrator<nprog,nconst,naux,nstages>::initialize(Tendencies<nprog, nconst, naux> &tend, VariableSet<nprog> &xvars, VariableSet<nconst> &consts, VariableSet<naux> &auxiliarys, ExchangeSet<nprog> &prog_exch)
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




  template<uint nprog, uint nconst, uint naux, uint nstages> void RKSimpleTimeIntegrator<nprog,nconst,naux,nstages>::stepForward(real dt)
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
