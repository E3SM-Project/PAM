

#ifndef _SSPRK_H_
#define _SSPRK_H_


#include "common.h"
#include "topology.h"
#include "variable_sets.h"
#include "exchange.h"
#include "model.h"



template<uint ndims, uint nprog, uint nconst, uint naux, uint nstages> class SSPKKTimeIntegrator {

public:

    VariableSet<ndims, nprog> F1;
    VariableSet<ndims, nprog> x1;
    VariableSet<ndims, nprog> F2;
    VariableSet<ndims, nprog> x2;
    VariableSet<ndims, nprog> x3;
    VariableSet<ndims, nprog> F3;
    VariableSet<ndims, nprog> *x;
    Tendencies<ndims, nprog, nconst, naux> *tendencies;
    const VariableSet<ndims, nconst> *const_vars;
    VariableSet<ndims, naux> *auxiliary_vars;
    ExchangeSet<ndims, nprog> *x_exchange;

  bool is_initialized;
  SSPKKTimeIntegrator();
  SSPKKTimeIntegrator( const SSPKKTimeIntegrator<ndims,nprog,nconst,naux,nstages> &ssprk) = delete;
  SSPKKTimeIntegrator& operator=( const SSPKKTimeIntegrator<ndims,nprog,nconst,naux,nstages> &ssprk) = delete;
  void initialize(Tendencies<ndims, nprog, nconst, naux> &tend, VariableSet<ndims, nprog> &xvars, const VariableSet<ndims, nconst> &consts, VariableSet<ndims, naux> &auxiliarys, ExchangeSet<ndims, nprog> &prog_exch);
  void stepForward(real dt);


};


    template<uint ndims, uint nprog, uint nconst, uint naux, uint nstages> SSPKKTimeIntegrator<ndims,nprog,nconst,naux,nstages>::SSPKKTimeIntegrator()
    {
      this->is_initialized = false;
      std::cout << "CREATED SSPKK\n";
    }


  template<uint ndims, uint nprog, uint nconst, uint naux, uint nstages> void SSPKKTimeIntegrator<ndims,nprog,nconst,naux,nstages>::initialize(Tendencies<ndims, nprog, nconst, naux> &tend, VariableSet<ndims, nprog> &xvars, const VariableSet<ndims, nconst> &consts, VariableSet<ndims, naux> &auxiliarys, ExchangeSet<ndims, nprog> &prog_exch)
  {
    this->x1.initialize(xvars, "x1");
    this->F1.initialize(xvars, "F1");
    this->x2.initialize(xvars, "x2");
    this->F2.initialize(xvars, "F2");
    if (nstages == 3)
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

  template<uint ndims, uint nprog, uint nconst, uint naux, uint nstages> void SSPKKTimeIntegrator<ndims,nprog,nconst,naux,nstages>::stepForward(real dt)
  {

      this->tendencies->compute_rhs(dt, *this->const_vars, *this->x, *this->auxiliary_vars, this->F1);
      this->x1.waxpy(-1.*dt, this->F1, *this->x);
      this->x_exchange->exchange_variable_set(this->x1);
      this->tendencies->compute_rhs(dt, *this->const_vars, this->x1, *this->auxiliary_vars, this->F2);

     if (nstages == 2)
     {
    this->x2.waxpbypcz(0.5, 0.5, -0.5*dt, *this->x, this->x1, this->F2);
    this->x->copy(this->x2);
    }

    if (nstages == 3)
    {

   this->x2.waxpbypcz(0.75, 0.25, -0.25*dt, *this->x, this->x1, this->F2);
   this->x_exchange->exchange_variable_set(this->x2);
   this->tendencies->compute_rhs(dt, *this->const_vars, this->x2, *this->auxiliary_vars, this->F3);
   this->x3.waxpbypcz(1./3., 2./3., -2./3.*dt, *this->x, this->x2, this->F3);
   this->x->copy(this->x3);
   }

   this->x_exchange->exchange_variable_set(*this->x);

  }

#endif
