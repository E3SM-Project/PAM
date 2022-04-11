
#pragma once


#include "common.h"
#include "topology.h"
#include "variable_sets.h"
#include "exchange.h"
#include "model.h"

#define NSTAGESMAX 6

template<uint nprog, uint nconst, uint naux> class RKSimpleTimeIntegrator {

public:

  SArray<real,1, NSTAGESMAX> stage_coeffs;
  VariableSet<nprog> xtend;
  VariableSet<nprog> xtemp;
  VariableSet<nprog> *x;
  Tendencies<nprog, nconst, naux> *tendencies;
  VariableSet<nconst> *const_vars;
  VariableSet<naux> *auxiliary_vars;
  ExchangeSet<nprog> *x_exchange;
  uint nstages;
  
  bool is_initialized;
  RKSimpleTimeIntegrator();
  RKSimpleTimeIntegrator( const RKSimpleTimeIntegrator<nprog,nconst,naux> &rksimple) = delete;
  RKSimpleTimeIntegrator& operator=( const RKSimpleTimeIntegrator<nprog,nconst,naux> &rksimple) = delete;
  void initialize(Parameters &params, Tendencies<nprog, nconst, naux> &tend, VariableSet<nprog> &xvars, VariableSet<nconst> &consts, VariableSet<naux> &auxiliarys, ExchangeSet<nprog> &prog_exch);
  void stepForward(real dt);
  void set_stage_coefficients(Parameters &params);

};


    template<uint nprog, uint nconst, uint naux> RKSimpleTimeIntegrator<nprog,nconst,naux>::RKSimpleTimeIntegrator()
    {
      this->is_initialized = false;
      std::cout << "CREATED RKSIMPLE\n";
    }

    template<uint nprog, uint nconst, uint naux> void RKSimpleTimeIntegrator<nprog,nconst,naux>::set_stage_coefficients(Parameters &params)
      {
        if (params.tstype == "kgrk4")
        {
        this->nstages = 4;
        this->stage_coeffs(0) = 1./4.;
        this->stage_coeffs(1) = 1./3.;
        this->stage_coeffs(2) = 1./2.;
        this->stage_coeffs(3) = 1.;
        }
      }

  template<uint nprog, uint nconst, uint naux> void RKSimpleTimeIntegrator<nprog,nconst,naux>::initialize(Parameters &params, Tendencies<nprog, nconst, naux> &tend, VariableSet<nprog> &xvars, VariableSet<nconst> &consts, VariableSet<naux> &auxiliarys, ExchangeSet<nprog> &prog_exch)
  {
    this->xtemp.initialize(xvars, "xtemp");
    this->xtend.initialize(xvars, "xtend");
    this->x = &xvars;
    this->tendencies = &tend;
    this->const_vars = &consts;
    this->auxiliary_vars = &auxiliarys;
    this->x_exchange = &prog_exch;
    set_stage_coefficients(params);
    this->is_initialized = true;
  }

  template<uint nprog, uint nconst, uint naux> void RKSimpleTimeIntegrator<nprog,nconst,naux>::stepForward(real dt)
  {

    this->tendencies->compute_rhs(dt, *this->const_vars, *this->x, *this->auxiliary_vars, this->xtend);
    this->xtemp.waxpy(-1.*dt * this->stage_coeffs(0), this->xtend, *this->x);

    for (int i=1; i<nstages; i++)
    {
      std::cout << "stage " << i << "\n";
      this->x_exchange->exchange_variable_set(this->xtemp);
      this->tendencies->compute_rhs(dt, *this->const_vars, this->xtemp, *this->auxiliary_vars, this->xtend);
      this->xtemp.waxpy(-1.*dt * this->stage_coeffs(i), this->xtend, *this->x);
    }
    // THIS COPY CAN BE AVOIDED IF WE ARE CLEVER ie swap x and xtemp
    // Would require being careful with IO also?
    this->x->copy(this->xtemp);
    this->x_exchange->exchange_variable_set(*this->x);
  }