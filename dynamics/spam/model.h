#pragma once

#include "common.h"
#include "topology.h"
#include "geometry.h"
#include "geometry.h"
#include "params.h"
#include "variable_sets.h"
#include "weno_func_recon.h" // needed to set TransformMatrices related stuff

//THIS NEEDS SLIGHT MODIFICATION FOR EXTRUDED- ADD MORE TRANSFORM RELATED STUFF?
template <uint nprog, uint nconst, uint ndiag> class Diagnostics {
public:

  const Topology *primal_topology;
  const Topology *dual_topology;
  Geometry<1,1,1> *primal_geometry;
  Geometry<1,1,1> *dual_geometry;
  ExchangeSet<ndiag> *diag_exchange;

  bool is_initialized;

   Diagnostics() {
     this->is_initialized = false;
     std::cout << "CREATED DIAGNOSTICS\n";
   }

   void initialize(const Topology &ptopo, const Topology &dtopo, Geometry<1,1,1> &pgeom, Geometry<1,1,1> &dgeom, ExchangeSet<ndiag> &diag_exchange)
   {
     this->primal_topology = &ptopo;
     this->dual_topology = &dtopo;
     this->primal_geometry = &pgeom;
     this->dual_geometry = &dgeom;
     this->diag_exchange = &diag_exchange;

     this->is_initialized = true;
   }
   
   virtual void compute_diag(const VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<ndiag> &diagnostic_vars) {};

};

//THIS NEEDS SLIGHT MODIFICATION FOR EXTRUDED- ADD MORE TRANSFORM RELATED STUFF?
template <uint nprog, uint nconst, uint naux> class Tendencies {
public:

  const Topology *primal_topology;
  const Topology *dual_topology;
  ExchangeSet<naux> *aux_exchange;
  ExchangeSet<nconst> *const_exchange;
  Geometry<1,1,1> *primal_geometry;
  Geometry<1,1,1> *dual_geometry;

  SArray<real,2,reconstruction_order,2> primal_to_gll;
  SArray<real,3,reconstruction_order,reconstruction_order,reconstruction_order> primal_wenoRecon;
  SArray<real,1,(reconstruction_order-1)/2+2> primal_wenoIdl;
  real primal_wenoSigma;

  SArray<real,2,dual_reconstruction_order,2> dual_to_gll;
  SArray<real,3,dual_reconstruction_order,dual_reconstruction_order,dual_reconstruction_order> dual_wenoRecon;
  SArray<real,1,(dual_reconstruction_order-1)/2+2> dual_wenoIdl;
  real dual_wenoSigma;

  SArray<real,2,coriolis_reconstruction_order,2> coriolis_to_gll;
  SArray<real,3,coriolis_reconstruction_order,coriolis_reconstruction_order,coriolis_reconstruction_order> coriolis_wenoRecon;
  SArray<real,1,(coriolis_reconstruction_order-1)/2+2> coriolis_wenoIdl;
  real coriolis_wenoSigma;

  bool is_initialized;
  
   Tendencies() {
     this->is_initialized = false;
     std::cout << "CREATED TENDENCIES\n";
   }

   void initialize(Parameters &params, const Topology &primal_topo, const Topology &dual_topo, Geometry<1,1,1> &primal_geom, Geometry<1,1,1> &dual_geom, ExchangeSet<naux> &aux_exchange, ExchangeSet<nconst> &const_exchange)
   {
     this->primal_topology = &primal_topo;
     this->dual_topology = &dual_topo;
     this->primal_geometry = &primal_geom;
     this->dual_geometry = &dual_geom;
     this->aux_exchange = &aux_exchange;
     this->const_exchange = &const_exchange;

    TransformMatrices::coefs_to_gll_lower( primal_to_gll );
    TransformMatrices::weno_sten_to_coefs(primal_wenoRecon);
    wenoSetIdealSigma<reconstruction_order>(primal_wenoIdl,primal_wenoSigma);

    TransformMatrices::coefs_to_gll_lower( dual_to_gll );
    TransformMatrices::weno_sten_to_coefs(dual_wenoRecon);
    wenoSetIdealSigma<dual_reconstruction_order>(dual_wenoIdl,dual_wenoSigma);

    TransformMatrices::coefs_to_gll_lower( coriolis_to_gll );
    TransformMatrices::weno_sten_to_coefs(coriolis_wenoRecon);
    wenoSetIdealSigma<coriolis_reconstruction_order>(coriolis_wenoIdl,coriolis_wenoSigma);
    
    this->is_initialized = true;
  }
  virtual void compute_constants(VariableSet<nconst> &const_vars, VariableSet<nprog> &x) {};
  virtual void YAKL_INLINE compute_rhs(real dt, VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<naux> &auxiliary_vars, VariableSet<nprog> &xtend) {};
};



