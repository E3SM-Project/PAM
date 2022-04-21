#pragma once


#include "common.h"
#include "topology.h"
#include "field_sets.h"
#include "exchange.h"
#include "model.h"

class SSPKKTimeIntegrator {

public:

    FieldSet<nprognostic> F1;
    FieldSet<nprognostic> x1;
    FieldSet<nprognostic> F2;
    FieldSet<nprognostic> x2;
    FieldSet<nprognostic> x3;
    FieldSet<nprognostic> F3;
    FieldSet<nprognostic> *x;
    Tendencies *tendencies;
    FieldSet<nconstant> *const_vars;
    FieldSet<nauxiliary> *auxiliary_vars;
    ExchangeSet<nprognostic> *x_exchange;
    std::string tstype;
    
  bool is_initialized;
  SSPKKTimeIntegrator();
  SSPKKTimeIntegrator( const SSPKKTimeIntegrator &ssprk) = delete;
  SSPKKTimeIntegrator& operator=( const SSPKKTimeIntegrator &ssprk) = delete;
  void initialize(ModelParameters &params, Tendencies &tend, FieldSet<nprognostic> &xvars, FieldSet<nconstant> &consts, FieldSet<nauxiliary> &auxiliarys, ExchangeSet<nprognostic> &prog_exch);
  void stepForward(real dt);


};


    SSPKKTimeIntegrator::SSPKKTimeIntegrator()
    {
      this->is_initialized = false;
    }


  void SSPKKTimeIntegrator::initialize(ModelParameters &params, Tendencies &tend, FieldSet<nprognostic> &xvars, FieldSet<nconstant> &consts, FieldSet<nauxiliary> &auxiliarys, ExchangeSet<nprognostic> &prog_exch)
  {
    tstype = params.tstype;
    
    this->x1.initialize(xvars, "x1");
    this->F1.initialize(xvars, "F1");
    this->x2.initialize(xvars, "x2");
    this->F2.initialize(xvars, "F2");
    if (tstype == "ssprk3")
    {
    this->x3.initialize(xvars, "x3");
    this->F3.initialize(xvars, "F3");
    }
    this->x = &xvars;
    this->tendencies = &tend;
    this->const_vars = &consts;
    this->auxiliary_vars = &auxiliarys;
    this->x_exchange = &prog_exch;
    this->is_initialized = true;
  }

  void SSPKKTimeIntegrator::stepForward(real dt)
  {

      this->tendencies->compute_rhs(dt, *this->const_vars, *this->x, *this->auxiliary_vars, this->F1);
      this->x1.waxpy(-1.*dt, this->F1, *this->x);
      this->x_exchange->exchange_variable_set(this->x1);
      this->tendencies->compute_rhs(dt, *this->const_vars, this->x1, *this->auxiliary_vars, this->F2);

     if (tstype == "ssprk2")
     {
    this->x2.waxpbypcz(0.5, 0.5, -0.5*dt, *this->x, this->x1, this->F2);
    this->x->copy(this->x2);
    }

    if (tstype == "ssprk3")
    {
   this->x2.waxpbypcz(0.75, 0.25, -0.25*dt, *this->x, this->x1, this->F2);
   this->x_exchange->exchange_variable_set(this->x2);
   this->tendencies->compute_rhs(dt, *this->const_vars, this->x2, *this->auxiliary_vars, this->F3);
   this->x3.waxpbypcz(1./3., 2./3., -2./3.*dt, *this->x, this->x2, this->F3);
   this->x->copy(this->x3);
   }

   this->x_exchange->exchange_variable_set(*this->x);

  }
