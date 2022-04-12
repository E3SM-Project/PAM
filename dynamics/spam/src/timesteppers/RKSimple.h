
#pragma once


#include "common.h"
#include "topology.h"
#include "variable_sets.h"
#include "exchange.h"
#include "model.h"

#define NSTAGESMAX 6

class RKSimpleTimeIntegrator {

public:

  SArray<real,1, NSTAGESMAX> stage_coeffs;
  VariableSet<nprognostic> xtend;
  VariableSet<nprognostic> xtemp;
  VariableSet<nprognostic> *x;
  Tendencies *tendencies;
  VariableSet<nconstant> *const_vars;
  VariableSet<nauxiliary> *auxiliary_vars;
  ExchangeSet<nprognostic> *x_exchange;
  uint nstages;
  
  bool is_initialized;
  RKSimpleTimeIntegrator();
  RKSimpleTimeIntegrator( const RKSimpleTimeIntegrator &rksimple) = delete;
  RKSimpleTimeIntegrator& operator=( const RKSimpleTimeIntegrator &rksimple) = delete;
  void initialize(ModelParameters &params, Tendencies &tend, VariableSet<nprognostic> &xvars, VariableSet<nconstant> &consts, VariableSet<nauxiliary> &auxiliarys, ExchangeSet<nprognostic> &prog_exch);
  void stepForward(real dt);
  void set_stage_coefficients(ModelParameters &params);

};


    RKSimpleTimeIntegrator::RKSimpleTimeIntegrator()
    {
      this->is_initialized = false;
      std::cout << "CREATED RKSIMPLE\n";
    }

    void RKSimpleTimeIntegrator::set_stage_coefficients(ModelParameters &params)
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

  void RKSimpleTimeIntegrator::initialize(ModelParameters &params, Tendencies &tend, VariableSet<nprognostic> &xvars, VariableSet<nconstant> &consts, VariableSet<nauxiliary> &auxiliarys, ExchangeSet<nprognostic> &prog_exch)
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

  void RKSimpleTimeIntegrator::stepForward(real dt)
  {

    this->tendencies->compute_rhs(dt, *this->const_vars, *this->x, *this->auxiliary_vars, this->xtend);
    this->xtemp.waxpy(-1.*dt * this->stage_coeffs(0), this->xtend, *this->x);

    for (int i=1; i<nstages; i++)
    {
      //std::cout << "stage " << i << "\n";
      this->x_exchange->exchange_variable_set(this->xtemp);
      this->tendencies->compute_rhs(dt, *this->const_vars, this->xtemp, *this->auxiliary_vars, this->xtend);
      this->xtemp.waxpy(-1.*dt * this->stage_coeffs(i), this->xtend, *this->x);
    }
    // THIS COPY CAN BE AVOIDED IF WE ARE CLEVER ie swap x and xtemp
    // Would require being careful with IO also?
    this->x->copy(this->xtemp);
    this->x_exchange->exchange_variable_set(*this->x);
  }
