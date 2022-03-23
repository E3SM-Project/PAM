#ifndef _MODEL_H_
#define _MODEL_H_

#include "common.h"
#include "topology.h"
#include "variable_sets.h"
#include "ext_deriv.h"
#include "hodge_star.h"
#include "fct.h"
#include "recon.h"
#include "wedge.h"
#include "geometry.h"
#include "params.h"
//#include "string.h"


// Number of variables
// dens
uint constexpr nprognostic = 2;
#define DENSVAR 0
#define QXZVAR 1

uint constexpr nconstant = 7;
#define VVAR 0
#define WVAR 1
#define UVAR 2
#define UWVAR 3
#define VTVAR 4
#define WTVAR 5
#define D2VAR 6

//functional derivatives = F, FW, B, K, he, hew
//primal grid reconstruction stuff- U, W, dens0, edgerecon, recon, vertedgerecon, vertrecon
//fct stuff- Phi, Mf, edgeflux
//Q/W STUFF?

uint constexpr nauxiliary = 17;

#define DENS0VAR 0
#define DENSRECONVAR 1
#define DENSEDGERECONVAR 2
#define DENSVERTRECONVAR 3
#define DENSVERTEDGERECONVAR 4

#define PHIVAR 5
#define PHIVERTVAR 6
#define EDGEFLUXVAR 7
#define VERTEDGEFLUXVAR 8
#define MFVAR 9

#define QXZ0VAR 10
#define QXZRECONVAR 11
#define QXZVERTRECONVAR 12
#define QXZEDGERECONVAR 13
#define QXZVERTEDGERECONVAR 14
#define QXZFLUXVAR 15
#define QXZVERTFLUXVAR 16

// q, associated concentration 0-forms for den

uint constexpr ndiagnostic = 9;
#define DENSLDIAGVAR 0
#define QXZDIAGVAR 1

#define QXZRECONDIAGVAR 2
#define QXZVERTRECONDIAGVAR 3
#define QXZEDGERECONDIAGVAR 4
#define QXZVERTEDGERECONDIAGVAR 5
#define QXZFLUXDIAGVAR 6
#define QXZVERTFLUXDIAGVAR 7
#define D2DIAGVAR 8

//track total densities, dens min/max
uint constexpr nstats = 4;

#define DENSSTAT 0
#define DENSMINSTAT 1
#define DENSMAXSTAT 2
#define QXZSTAT 3


// *******   Model Specific Parameters   ***********//

void set_model_specific_params(std::string inFile, ModelParameters &params)
{
  params.etime = 0.0;
}

// ******* Diagnostics *************//

template <uint nprog, uint nconst, uint ndiag> class Diagnostics {
public:

  const Topology *primal_topology;
  const Topology *dual_topology;
  Geometry<1,1,1> *primal_geometry;
  Geometry<1,1,1> *dual_geometry;
  ExchangeSet<ndiag> *diag_exchange;

  //TransformMatrices<real> trans;

  SArray<real,2,reconstruction_order,2> primal_to_gll;
  SArray<real,3,reconstruction_order,reconstruction_order,reconstruction_order> primal_wenoRecon;
  SArray<real,1,(reconstruction_order-1)/2+2> primal_wenoIdl;
  real primal_wenoSigma;

  SArray<real,2,dual_reconstruction_order,2> dual_to_gll;
  SArray<real,3,dual_reconstruction_order,dual_reconstruction_order,dual_reconstruction_order> dual_wenoRecon;
  SArray<real,1,(dual_reconstruction_order-1)/2+2> dual_wenoIdl;
  real dual_wenoSigma;

  SArray<real,2,vert_reconstruction_order,2> primal_vert_to_gll;
  SArray<real,3,vert_reconstruction_order,vert_reconstruction_order,vert_reconstruction_order> primal_vert_wenoRecon;
  SArray<real,1,(vert_reconstruction_order-1)/2+2> primal_vert_wenoIdl;
  real primal_vert_wenoSigma;

  SArray<real,2,dual_vert_reconstruction_order,2> dual_vert_to_gll;
  SArray<real,3,dual_vert_reconstruction_order,dual_vert_reconstruction_order,dual_vert_reconstruction_order> dual_vert_wenoRecon;
  SArray<real,1,(dual_vert_reconstruction_order-1)/2+2> dual_vert_wenoIdl;
  real dual_vert_wenoSigma;
  
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

     TransformMatrices::coefs_to_gll_lower( primal_to_gll );
     TransformMatrices::weno_sten_to_coefs(primal_wenoRecon);
     wenoSetIdealSigma<reconstruction_order>(primal_wenoIdl,primal_wenoSigma);

     TransformMatrices::coefs_to_gll_lower( dual_to_gll );
     TransformMatrices::weno_sten_to_coefs(dual_wenoRecon);
     wenoSetIdealSigma<dual_reconstruction_order>(dual_wenoIdl,dual_wenoSigma);

     TransformMatrices::coefs_to_gll_lower( primal_vert_to_gll );
     TransformMatrices::weno_sten_to_coefs(primal_vert_wenoRecon);
     wenoSetIdealSigma<vert_reconstruction_order>(primal_vert_wenoIdl,primal_vert_wenoSigma);

     TransformMatrices::coefs_to_gll_lower( dual_vert_to_gll );
     TransformMatrices::weno_sten_to_coefs(dual_vert_wenoRecon);
     wenoSetIdealSigma<dual_vert_reconstruction_order>(dual_vert_wenoIdl,dual_vert_wenoSigma);
    
     this->is_initialized = true;

   }


   void compute_diag(const VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<ndiag> &diagnostic_vars)
   {

     // Compute T0
     int pis = primal_topology->is;
     int pjs = primal_topology->js;
     int pks = primal_topology->ks;
     
     // yakl::parallel_for("ComputeDens0", primal_topology->n_cells_interfaces, YAKL_LAMBDA (int iGlob) {
     //   int k, j, i;
     //   yakl::unpackIndices(iGlob, primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
       parallel_for( Bounds<3>( primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
     compute_Iext<ntdofs, diff_ord, vert_diff_ord>(diagnostic_vars.fields_arr[DENSLDIAGVAR].data, x.fields_arr[DENSVAR].data, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k);
     });


     // Compute q0
     int dis = dual_topology->is;
     int djs = dual_topology->js;
     int dks = dual_topology->ks;
     
     // yakl::parallel_for("ComputeQ0", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
     //   int k, j, i;
     //   yakl::unpackIndices(iGlob, dual_topology->ni -2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
       parallel_for( Bounds<3>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

     compute_Jext<nQdofs, diff_ord, vert_diff_ord>(diagnostic_vars.fields_arr[QXZDIAGVAR].data, x.fields_arr[QXZVAR].data, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k+1);
     });


//TEMPORARY        
     //compute qrecon
     // 
     // yakl::parallel_for("ComputeQEdgeRecon", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
     //   int k, j, i;
     //   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
     // 
     //    compute_straight_xz_edge_recon<nQdofs, reconstruction_type, reconstruction_order>(
     //      diagnostic_vars.fields_arr[QXZEDGERECONDIAGVAR].data, diagnostic_vars.fields_arr[QXZDIAGVAR].data, pis, pjs, pks, i, j, k,
     //      primal_wenoRecon, primal_to_gll, primal_wenoIdl, primal_wenoSigma);
     // 
     //    compute_straight_xz_vert_edge_recon<nQdofs, vert_reconstruction_type, vert_reconstruction_order>(
     //      diagnostic_vars.fields_arr[QXZVERTEDGERECONDIAGVAR].data, diagnostic_vars.fields_arr[QXZDIAGVAR].data, pis, pjs, pks, i, j, k,
     //      primal_vert_wenoRecon, primal_vert_to_gll, primal_vert_wenoIdl, primal_vert_wenoSigma);
     // 
     // });
     // 
     // this->diag_exchange->exchanges_arr[QXZEDGERECONDIAGVAR].exchange_field(diagnostic_vars.fields_arr[QXZEDGERECONDIAGVAR]);
     // this->diag_exchange->exchanges_arr[QXZVERTEDGERECONDIAGVAR].exchange_field(diagnostic_vars.fields_arr[QXZVERTEDGERECONDIAGVAR]);
     // 
     // yakl::parallel_for("ComputeQRECON", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
     //   int k, j, i;
     //   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
     // 
     //   compute_straight_xz_recon<nQdofs, reconstruction_type>(
     //     diagnostic_vars.fields_arr[QXZRECONDIAGVAR].data, diagnostic_vars.fields_arr[QXZEDGERECONDIAGVAR].data, const_vars.fields_arr[WTVAR].data, pis, pjs, pks, i, j, k);
     // 
     // });
     // yakl::parallel_for("ComputeQVERTRECON", primal_topology->n_cells_interfaces, YAKL_LAMBDA (int iGlob) {
     //   int k, j, i;
     //   yakl::unpackIndices(iGlob, primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
     // 
     //   compute_straight_xz_vert_recon<nQdofs, vert_reconstruction_type>(
     //     diagnostic_vars.fields_arr[QXZVERTRECONDIAGVAR].data, diagnostic_vars.fields_arr[QXZVERTEDGERECONDIAGVAR].data, const_vars.fields_arr[VTVAR].data, pis, pjs, pks, i, j, k);
     // 
     // });
     // this->diag_exchange->exchanges_arr[QXZRECONDIAGVAR].exchange_field(diagnostic_vars.fields_arr[QXZRECONDIAGVAR]);
     // this->diag_exchange->exchanges_arr[QXZVERTRECONDIAGVAR].exchange_field(diagnostic_vars.fields_arr[QXZVERTRECONDIAGVAR]);
     // 
     // // //compute Q flux
     //    yakl::parallel_for("ComputeQFLUX", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
     //      int k, j, i;
     //      yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
     // 
     //      if (qf_choice == QF_MODE::EC)
     //     { compute_Qxz_w_EC<nQdofs>(diagnostic_vars.fields_arr[QXZFLUXDIAGVAR].data, diagnostic_vars.fields_arr[QXZRECONDIAGVAR].data, diagnostic_vars.fields_arr[QXZVERTRECONDIAGVAR].data, const_vars.fields_arr[UVAR].data, pis, pjs, pks, i, j, k);}
     // 
     //      if (qf_choice == QF_MODE::NOEC)
     //     { compute_Qxz_w_nonEC<nQdofs>(diagnostic_vars.fields_arr[QXZFLUXDIAGVAR].data, diagnostic_vars.fields_arr[QXZRECONDIAGVAR].data, const_vars.fields_arr[UVAR].data, pis, pjs, pks, i, j, k);}
     // 
     // });
     //    yakl::parallel_for("ComputeQVERTFLUX", primal_topology->n_cells_interfaces, YAKL_LAMBDA (int iGlob) {
     //      int k, j, i;
     //      yakl::unpackIndices(iGlob, primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
     // 
     //      if (qf_choice == QF_MODE::EC)
     //      { compute_Qxz_u_EC<nQdofs>(diagnostic_vars.fields_arr[QXZVERTFLUXDIAGVAR].data, diagnostic_vars.fields_arr[QXZRECONDIAGVAR].data, diagnostic_vars.fields_arr[QXZVERTRECONDIAGVAR].data, const_vars.fields_arr[UWVAR].data, pis, pjs, pks, i, j, k);}
     // 
     //      if (qf_choice == QF_MODE::NOEC)
     //      { compute_Qxz_u_nonEC<nQdofs>(diagnostic_vars.fields_arr[QXZVERTFLUXDIAGVAR].data, diagnostic_vars.fields_arr[QXZVERTRECONDIAGVAR].data, const_vars.fields_arr[UWVAR].data, pis, pjs, pks, i, j, k);}
     // 
     // });
     //     this->diag_exchange->exchanges_arr[QXZFLUXDIAGVAR].exchange_field(diagnostic_vars.fields_arr[QXZFLUXDIAGVAR]);
     //     this->diag_exchange->exchanges_arr[QXZVERTFLUXDIAGVAR].exchange_field(diagnostic_vars.fields_arr[QXZVERTFLUXDIAGVAR]);        
     // 
     //     yakl::parallel_for("ComputeQXZTend", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
     //       int k, j, i;
     //       yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
     //     compute_Dxz<nQdofs> (diagnostic_vars.fields_arr[D2DIAGVAR].data, diagnostic_vars.fields_arr[QXZVERTFLUXDIAGVAR].data, diagnostic_vars.fields_arr[QXZFLUXDIAGVAR].data, pis, pjs, pks, i, j, k);
     //     //compute_Dxz<nQdofs> (diagnostic_vars.fields_arr[D2DIAGVAR].data, diagnostic_vars.fields_arr[QXZFLUXDIAGVAR].data, diagnostic_vars.fields_arr[QXZVERTFLUXDIAGVAR].data, pis, pjs, pks, i, j, k);
     //     });
        
}

};
// *******   Tendencies   ***********//

// THIS SHOULD BE GENERALIZABLE...
// ONLY COMPUTE FUNCTION NEEDS TO CHANGE
template <uint nprog, uint nconst, uint naux> class Tendencies {
public:

  const Topology *primal_topology;
  const Topology *dual_topology;
  ExchangeSet<naux> *aux_exchange;
  ExchangeSet<nconst> *const_exchange;
  Geometry<1,1,1> *primal_geometry;
  Geometry<1,1,1> *dual_geometry;

  //TransformMatrices<real> trans;

  SArray<real,2,reconstruction_order,2> primal_to_gll;
  SArray<real,3,reconstruction_order,reconstruction_order,reconstruction_order> primal_wenoRecon;
  SArray<real,1,(reconstruction_order-1)/2+2> primal_wenoIdl;
  real primal_wenoSigma;

  SArray<real,2,dual_reconstruction_order,2> dual_to_gll;
  SArray<real,3,dual_reconstruction_order,dual_reconstruction_order,dual_reconstruction_order> dual_wenoRecon;
  SArray<real,1,(dual_reconstruction_order-1)/2+2> dual_wenoIdl;
  real dual_wenoSigma;

  SArray<real,2,vert_reconstruction_order,2> primal_vert_to_gll;
  SArray<real,3,vert_reconstruction_order,vert_reconstruction_order,vert_reconstruction_order> primal_vert_wenoRecon;
  SArray<real,1,(vert_reconstruction_order-1)/2+2> primal_vert_wenoIdl;
  real primal_vert_wenoSigma;

  SArray<real,2,dual_vert_reconstruction_order,2> dual_vert_to_gll;
  SArray<real,3,dual_vert_reconstruction_order,dual_vert_reconstruction_order,dual_vert_reconstruction_order> dual_vert_wenoRecon;
  SArray<real,1,(dual_vert_reconstruction_order-1)/2+2> dual_vert_wenoIdl;
  real dual_vert_wenoSigma;
  
  bool is_initialized;
  
   Tendencies() {
     this->is_initialized = false;
     std::cout << "CREATED TENDENCIES\n";
   }

   void initialize(ModelParameters &params, const Topology &primal_topo, const Topology &dual_topo, Geometry<1,1,1> &primal_geom, Geometry<1,1,1> &dual_geom, ExchangeSet<naux> &aux_exchange, ExchangeSet<nconst> &const_exchange)
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

     TransformMatrices::coefs_to_gll_lower( primal_vert_to_gll );
     TransformMatrices::weno_sten_to_coefs(primal_vert_wenoRecon);
     wenoSetIdealSigma<vert_reconstruction_order>(primal_vert_wenoIdl,primal_vert_wenoSigma);

     TransformMatrices::coefs_to_gll_lower( dual_vert_to_gll );
     TransformMatrices::weno_sten_to_coefs(dual_vert_wenoRecon);
     wenoSetIdealSigma<dual_vert_reconstruction_order>(dual_vert_wenoIdl,dual_vert_wenoSigma);
    
    this->is_initialized = true;
  }


  void compute_constants(VariableSet<nconst> &const_vars, VariableSet<nprog> &x)
  {
    
    int pis = primal_topology->is;
    int pjs = primal_topology->js;
    int pks = primal_topology->ks;

    int dis = dual_topology->is;
    int djs = dual_topology->js;
    int dks = dual_topology->ks;

    // Compute U = H V and UW = Hv W

    //yakl::parallel_for("ComputeUVAR", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
      //int k, j, i;
      //yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
      parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

    compute_Hext<1, diff_ord>(const_vars.fields_arr[UVAR].data, const_vars.fields_arr[VVAR].data, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k);
    });
    this->const_exchange->exchanges_arr[UVAR].exchange_field(const_vars.fields_arr[UVAR]);
    
    //yakl::parallel_for("ComputeUWVAR", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
    //  int k, j, i;
    //  yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
    parallel_for( Bounds<3>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
    compute_Hv<1, vert_diff_ord>(const_vars.fields_arr[UWVAR].data, const_vars.fields_arr[WVAR].data, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k+1);
    });    
    this->const_exchange->exchanges_arr[UWVAR].exchange_field(const_vars.fields_arr[UWVAR]);

    //compute VT = "Wxz" Uw ie w/Uw at v
    //compute WT = "Wxz" U ie v/U at w 
  //  yakl::parallel_for("ComputeVTVAR", primal_topology->n_cells_interfaces, YAKL_LAMBDA (int iGlob) {
    //  int k, j, i;
    //  yakl::unpackIndices(iGlob, primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
      parallel_for( Bounds<3>( primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

compute_Wxz_u(const_vars.fields_arr[VTVAR].data,const_vars.fields_arr[UWVAR].data, pis, pjs, pks, i, j, k);
 });
    //yakl::parallel_for("ComputeWTVAR", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
    //  int k, j, i;
    //  yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
      parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

compute_Wxz_w(const_vars.fields_arr[WTVAR].data,const_vars.fields_arr[UVAR].data, pis, pjs, pks, i, j, k);
 });
    this->const_exchange->exchanges_arr[VTVAR].exchange_field(const_vars.fields_arr[VTVAR]);
    this->const_exchange->exchanges_arr[WTVAR].exchange_field(const_vars.fields_arr[WTVAR]);

    // yakl::parallel_for("ComputeQXZTend", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
    //   int k, j, i;
    //   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
    // compute_Dxz<1> (const_vars.fields_arr[D2VAR].data, const_vars.fields_arr[VVAR].data, const_vars.fields_arr[WVAR].data, pis, pjs, pks, i, j, k);
    // });
    
    
  }


  void YAKL_INLINE compute_rhs(real dt, VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<naux> &auxiliary_vars, VariableSet<nprog> &xtend)
  {

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

// Compute Q0

//yakl::parallel_for("ComputeQ0VAR", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
//  int k, j, i;
//  yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
parallel_for( Bounds<3>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
compute_Jext<nQdofs, diff_ord, vert_diff_ord>(auxiliary_vars.fields_arr[QXZ0VAR].data, x.fields_arr[QXZVAR].data, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k+1);
});

// Compute DENS0

//yakl::parallel_for("ComputeDens0VAR", primal_topology->n_cells_interfaces, YAKL_LAMBDA (int iGlob) {
//  int k, j, i;
//   yakl::unpackIndices(iGlob, primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
parallel_for( Bounds<3>( primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
compute_Iext<ntdofs, diff_ord, vert_diff_ord>(auxiliary_vars.fields_arr[DENS0VAR].data, x.fields_arr[DENSVAR].data, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k);
});

this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(auxiliary_vars.fields_arr[DENS0VAR]);
this->aux_exchange->exchanges_arr[QXZ0VAR].exchange_field(auxiliary_vars.fields_arr[QXZ0VAR]);



//compute qrecon

// yakl::parallel_for("ComputeQEdgeRecon", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

   compute_straight_xz_edge_recon<nQdofs, reconstruction_type, reconstruction_order>(
     auxiliary_vars.fields_arr[QXZEDGERECONVAR].data, auxiliary_vars.fields_arr[QXZ0VAR].data, pis, pjs, pks, i, j, k,
     primal_wenoRecon, primal_to_gll, primal_wenoIdl, primal_wenoSigma);

   compute_straight_xz_vert_edge_recon<nQdofs, vert_reconstruction_type, vert_reconstruction_order>(
     auxiliary_vars.fields_arr[QXZVERTEDGERECONVAR].data, auxiliary_vars.fields_arr[QXZ0VAR].data, pis, pjs, pks, i, j, k,
     primal_vert_wenoRecon, primal_vert_to_gll, primal_vert_wenoIdl, primal_vert_wenoSigma);
  
});

this->aux_exchange->exchanges_arr[QXZEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[QXZEDGERECONVAR]);
this->aux_exchange->exchanges_arr[QXZVERTEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[QXZVERTEDGERECONVAR]);

// yakl::parallel_for("ComputeQRECON", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
  compute_straight_xz_recon<nQdofs, reconstruction_type>(
    auxiliary_vars.fields_arr[QXZRECONVAR].data, auxiliary_vars.fields_arr[QXZEDGERECONVAR].data, const_vars.fields_arr[WTVAR].data, pis, pjs, pks, i, j, k);
});
// yakl::parallel_for("ComputeQVERTRECON", primal_topology->n_cells_interfaces, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
parallel_for( Bounds<3>( primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
  compute_straight_xz_vert_recon<nQdofs, vert_reconstruction_type>(
    auxiliary_vars.fields_arr[QXZVERTRECONVAR].data, auxiliary_vars.fields_arr[QXZVERTEDGERECONVAR].data, const_vars.fields_arr[VTVAR].data, pis, pjs, pks, i, j, k);
});
this->aux_exchange->exchanges_arr[QXZRECONVAR].exchange_field(auxiliary_vars.fields_arr[QXZRECONVAR]);
this->aux_exchange->exchanges_arr[QXZVERTRECONVAR].exchange_field(auxiliary_vars.fields_arr[QXZVERTRECONVAR]);

//compute Densrecon

// yakl::parallel_for("ComputeDensEdgeRecon", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

  compute_twisted_edge_recon<ntdofs, dual_reconstruction_type, dual_reconstruction_order>(
    auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data, dis, djs, dks, i, j, k,
    dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
    
    compute_twisted_vert_edge_recon<ntdofs, dual_vert_reconstruction_type, dual_vert_reconstruction_order>(
      auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data, dis, djs, dks, i, j, k, 
      dual_vert_wenoRecon, dual_vert_to_gll, dual_vert_wenoIdl, dual_vert_wenoSigma);
  });

this->aux_exchange->exchanges_arr[DENSEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSEDGERECONVAR]);
this->aux_exchange->exchanges_arr[DENSVERTEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR]);

// yakl::parallel_for("ComputeDensRECON", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

  compute_twisted_recon<ntdofs, dual_reconstruction_type>(
    auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, const_vars.fields_arr[UVAR].data, dis, djs, dks, i, j, k);
  });

    // yakl::parallel_for("ComputeDensVertRECON", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
    //   int k, j, i;
    //   yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
      parallel_for( Bounds<3>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
    compute_twisted_vert_recon<ntdofs, dual_vert_reconstruction_type>(
      auxiliary_vars.fields_arr[DENSVERTRECONVAR].data, auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data, const_vars.fields_arr[UWVAR].data, dis, djs, dks, i, j, k+1);
});
this->aux_exchange->exchanges_arr[DENSRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSRECONVAR]);
this->aux_exchange->exchanges_arr[DENSVERTRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSVERTRECONVAR]);


// compute FCT stuff
      // yakl::parallel_for("ComputeEdgeFlux", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
      //   int k, j, i;
      //   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
        parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

      compute_edgefluxes<ntdofs> (auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[DENSRECONVAR].data, const_vars.fields_arr[UVAR].data, dis, djs, dks, i, j, k);
      });
      
      // yakl::parallel_for("ComputeVertEdgeFlux", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
      //   int k, j, i;
      //   yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
        parallel_for( Bounds<3>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

      compute_vertedgefluxes<ntdofs> (auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data, auxiliary_vars.fields_arr[DENSVERTRECONVAR].data, const_vars.fields_arr[UWVAR].data, dis, djs, dks, i, j, k+1);
      });
      
this->aux_exchange->exchanges_arr[EDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[EDGEFLUXVAR]);
this->aux_exchange->exchanges_arr[VERTEDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[VERTEDGEFLUXVAR]);


// yakl::parallel_for("ComputeMf", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

  compute_Mfext<ntdofs> (auxiliary_vars.fields_arr[MFVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data, dt, dis, djs, dks, i, j, k);
});

    this->aux_exchange->exchanges_arr[MFVAR].exchange_field(auxiliary_vars.fields_arr[MFVAR]);

    // yakl::parallel_for("ComputePhiVert", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
    //   int k, j, i;
    //   yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
      parallel_for( Bounds<3>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

      compute_Phivert<ntdofs> (auxiliary_vars.fields_arr[PHIVERTVAR].data, auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k+1);
    });

    // yakl::parallel_for("ComputePhi", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
    //   int k, j, i;
    //   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
      parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
      compute_Phi<ntdofs> (auxiliary_vars.fields_arr[PHIVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k);
    });
    
//Don't do FCT for non-FCT vars
    for (int l=0; l<ntnofctdofs; l++)
    {
    auxiliary_vars.fields_arr[PHIVAR].set(l, 1.0);
    auxiliary_vars.fields_arr[PHIVERTVAR].set(l, 1.0);
    }
    
    this->aux_exchange->exchanges_arr[PHIVAR].exchange_field(auxiliary_vars.fields_arr[PHIVAR]);
    this->aux_exchange->exchanges_arr[PHIVERTVAR].exchange_field(auxiliary_vars.fields_arr[PHIVERTVAR]);



// //compute Q flux
   // yakl::parallel_for("ComputeQFLUX", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
   //   int k, j, i;
   //   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
     parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

     if (qf_choice == QF_MODE::EC)
    { compute_Qxz_w_EC<nQdofs>(auxiliary_vars.fields_arr[QXZFLUXVAR].data, auxiliary_vars.fields_arr[QXZRECONVAR].data, auxiliary_vars.fields_arr[QXZVERTRECONVAR].data, const_vars.fields_arr[UVAR].data, pis, pjs, pks, i, j, k);}

     if (qf_choice == QF_MODE::NOEC)
    { compute_Qxz_w_nonEC<nQdofs>(auxiliary_vars.fields_arr[QXZFLUXVAR].data, auxiliary_vars.fields_arr[QXZRECONVAR].data, const_vars.fields_arr[UVAR].data, pis, pjs, pks, i, j, k);}

});
   // yakl::parallel_for("ComputeQVERTFLUX", primal_topology->n_cells_interfaces, YAKL_LAMBDA (int iGlob) {
   //   int k, j, i;
   //   yakl::unpackIndices(iGlob, primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
   parallel_for( Bounds<3>( primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
     
     if (qf_choice == QF_MODE::EC)
     { compute_Qxz_u_EC<nQdofs>(auxiliary_vars.fields_arr[QXZVERTFLUXVAR].data, auxiliary_vars.fields_arr[QXZRECONVAR].data, auxiliary_vars.fields_arr[QXZVERTRECONVAR].data, const_vars.fields_arr[UWVAR].data, pis, pjs, pks, i, j, k);}

     if (qf_choice == QF_MODE::NOEC)
     { compute_Qxz_u_nonEC<nQdofs>(auxiliary_vars.fields_arr[QXZVERTFLUXVAR].data, auxiliary_vars.fields_arr[QXZVERTRECONVAR].data, const_vars.fields_arr[UWVAR].data, pis, pjs, pks, i, j, k);}

});
    this->aux_exchange->exchanges_arr[QXZFLUXVAR].exchange_field(auxiliary_vars.fields_arr[QXZFLUXVAR]);
    this->aux_exchange->exchanges_arr[QXZVERTFLUXVAR].exchange_field(auxiliary_vars.fields_arr[QXZVERTFLUXVAR]);


// // Compute Dxz "QXZFLUX"
// yakl::parallel_for("ComputeQXZTend", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
compute_Dxz<nQdofs> (xtend.fields_arr[QXZVAR].data, auxiliary_vars.fields_arr[QXZVERTFLUXVAR].data, auxiliary_vars.fields_arr[QXZFLUXVAR].data, pis, pjs, pks, i, j, k);
});

// compute D2bar "F"
    // yakl::parallel_for("ComputeDensTend", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
    //   int k, j, i;
    //   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
      parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
      compute_wDbar2_fct<ntdofs> (xtend.fields_arr[DENSVAR].data, auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[PHIVAR].data, const_vars.fields_arr[UVAR].data, dis, djs, dks, i, j, k);
      compute_wDvbar_fct<ntdofs, ADD_MODE::ADD> (xtend.fields_arr[DENSVAR].data, auxiliary_vars.fields_arr[DENSVERTRECONVAR].data, auxiliary_vars.fields_arr[PHIVERTVAR].data, const_vars.fields_arr[UWVAR].data, dis, djs, dks, i, j, k);
    });
    
}
};



// *******   Statistics Calculations   ***********//

// THIS STUFF SHOULD BE CLEANED UP AND GENERALIZED LIKE VARIABLE SETS IF POSSIBLE...
// ONLY COMPUTE FUNCTION NEEDS TO CHANGE!
class Stat
{
public:
  real2d data;
  std::string name;
  int ndofs;


void initialize(std::string statName, int ndof, ModelParameters &params, Parallel &par)
{
  name = statName;
  ndofs = ndof;

  if (par.masterproc)
  {
    data = real2d(name.c_str(), ndofs, params.Nsteps/params.Nstat + 1);
  }
}

};



template <uint nprog, uint nconst, uint nstats> class Stats
{
public:
  std::array<Stat,nstats> stats_arr;
  MPI_Request Req [nstats];
  MPI_Status  Status[nstats];
  int ierr;
  int statsize;
  int masterproc;
  const Topology *primal_topology;
  const Topology *dual_topology;
  Geometry<1,1,1> *primal_geometry;
  Geometry<1,1,1> *dual_geometry;
  real3d trimmed_dual_density, trimmed_primal_density;

  void initialize(ModelParameters &params, Parallel &par, const Topology &primal_topo, const Topology &dual_topo, Geometry<1,1,1> &primal_geom, Geometry<1,1,1> &dual_geom)
  {
    this->primal_topology = &primal_topo;
    this->dual_topology = &dual_topo;
    this->primal_geometry = &primal_geom;
    this->dual_geometry = &dual_geom;

    statsize = params.Nsteps/params.Nstat + 1;
    stats_arr[DENSSTAT].initialize("mass", ntdofs, params, par);
    stats_arr[DENSMAXSTAT].initialize("densmax", ntdofs, params, par);
    stats_arr[DENSMINSTAT].initialize("densmin", ntdofs, params, par);
    stats_arr[QXZSTAT].initialize("qmass", nQdofs, params, par);
    masterproc = par.masterproc;
  
    trimmed_dual_density = real3d("trimmed_dual_density", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
    trimmed_primal_density = real3d("trimmed_primal_density", this->primal_topology->nl, this->primal_topology->n_cells_y, this->primal_topology->n_cells_x);

  }



  void compute( VariableSet<nprog> &progvars,  VariableSet<nconst> &constvars, int tind)
  {


      SArray<real,1,ntdofs> masslocal, massglobal;
      SArray<real,1,ntdofs> densmaxlocal, densmaxglobal;
      SArray<real,1,ntdofs> densminlocal, densminglobal;
      SArray<real,1,nQdofs> Qlocal, Qglobal;

      for (int l=0;l<ntdofs;l++) {masslocal(l) = 0.; massglobal(l) = 0.;}
      for (int l=0;l<ntdofs;l++) {densmaxlocal(l) = 0.; densmaxglobal(l) = 0.;}
      for (int l=0;l<ntdofs;l++) {densminlocal(l) = 0.; densminglobal(l) = 0.;}
  
      int dis = dual_topology->is;
      int djs = dual_topology->js;
      int dks = dual_topology->ks;
      
    for (int l=0;l<ntdofs;l++)
    {
      parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
        trimmed_dual_density(k,j,i) = progvars.fields_arr[DENSVAR].data(l,k+dks,j+djs,i+dis);
      });

    masslocal(l) = yakl::intrinsics::sum(trimmed_dual_density);
    densmaxlocal(l) = yakl::intrinsics::maxval(trimmed_dual_density);
    densminlocal(l) = yakl::intrinsics::minval(trimmed_dual_density);
    }

    int pis = primal_topology->is;
    int pjs = primal_topology->js;
    int pks = primal_topology->ks;
    
    for (int l=0;l<nQdofs;l++)
    {
      parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
        trimmed_primal_density(k,j,i) = progvars.fields_arr[QXZVAR].data(l,k+pks,j+pjs,i+pis);
      });

    Qlocal(l) = yakl::intrinsics::sum(trimmed_primal_density);
    }
    
    
    //MPI sum/min/max
    this->ierr = MPI_Ireduce( &masslocal, &massglobal, ntdofs, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[DENSSTAT]);
    this->ierr = MPI_Ireduce( &densmaxlocal, &densmaxglobal, ntdofs, REAL_MPI, MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[DENSMAXSTAT]);
    this->ierr = MPI_Ireduce( &densminlocal, &densminglobal, ntdofs, REAL_MPI, MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[DENSMINSTAT]);
    this->ierr = MPI_Ireduce( &Qlocal, &Qglobal, nQdofs, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[QXZSTAT]);
    this->ierr = MPI_Waitall(nstats, this->Req, this->Status);


  if (masterproc)
  {
    for (int l=0;l<ntdofs;l++)
    {
  this->stats_arr[DENSSTAT].data(l,tind) = massglobal(l);
  this->stats_arr[DENSMAXSTAT].data(l,tind) = densmaxglobal(l);
  this->stats_arr[DENSMINSTAT].data(l,tind) = densminglobal(l);
}

for (int l=0;l<nQdofs;l++)
{
this->stats_arr[QXZSTAT].data(l,tind) = Qglobal(l);
}

  }
  }
};


template <uint nvars> void set_dofs_arr(SArray<int,2, nvars, 3> &dofs_arr, int var, int basedof, int extdof, int ndofs)
{
  dofs_arr(var, 0) = basedof;
  dofs_arr(var, 1) = extdof;
  dofs_arr(var, 2) = ndofs;
}

// *******   VariableSet Initialization   ***********//
template <uint nprog, uint nconst, uint naux, uint ndiag> void initialize_variables(const Topology &ptopo, const Topology &dtopo,
SArray<int,2, nprog, 3> &prog_ndofs_arr, SArray<int,2, nconst, 3> &const_ndofs_arr, SArray<int,2, naux, 3> &aux_ndofs_arr, SArray<int,2, ndiag, 3> &diag_ndofs_arr,
std::array<std::string, nprog> &prog_names_arr, std::array<std::string, nconst> &const_names_arr, std::array<std::string, naux> &aux_names_arr, std::array<std::string, ndiag> &diag_names_arr,
std::array<const Topology *, nprog> &prog_topo_arr, std::array<const Topology *, nconst> &const_topo_arr, std::array<const Topology *, naux> &aux_topo_arr, std::array<const Topology *, ndiag> &diag_topo_arr)
{

  //primal grid represents straight quantities, dual grid twisted quantities
// ndims is the BASEDIM size!

  // ADD Q RELATED stuff

  // dens, QXZ
  prog_topo_arr[DENSVAR] = &dtopo;
  prog_names_arr[DENSVAR] = "dens";
  set_dofs_arr(prog_ndofs_arr, DENSVAR, ndims, 1, ntdofs); //dens = twisted (n,1)-form
  prog_topo_arr[QXZVAR] = &ptopo;
  prog_names_arr[QXZVAR] = "QXZ";
  set_dofs_arr(prog_ndofs_arr, QXZVAR, ndims, 1, nQdofs); //Q = straight (n,1)-form
  
  // velocities
  const_topo_arr[VVAR] = &ptopo;
  const_topo_arr[WVAR] = &ptopo;
  const_names_arr[VVAR] = "v";
  const_names_arr[WVAR] = "w";
  set_dofs_arr(const_ndofs_arr, VVAR, 1, 0, 1); //v = straight (1,0)-form
  set_dofs_arr(const_ndofs_arr, WVAR, 0, 1, 1); //w = straight (0,1)-form
  const_topo_arr[UVAR] = &dtopo;
  const_topo_arr[UWVAR] = &dtopo;
  const_names_arr[UVAR] = "U";
  const_names_arr[UWVAR] = "Uw";
  set_dofs_arr(const_ndofs_arr, UVAR, ndims-1, 1, 1); //U = twisted (n-1,1)-form
  set_dofs_arr(const_ndofs_arr, UWVAR, ndims, 0, 1); //Uw = twisted (n,0)-form

  const_topo_arr[VTVAR] = &ptopo;
  const_topo_arr[WTVAR] = &ptopo;
  const_names_arr[VTVAR] = "vt";
  const_names_arr[WTVAR] = "wt";
  set_dofs_arr(const_ndofs_arr, VTVAR, 1, 0, 1); //VT = straight (1,0)-form ie w at v pts
  set_dofs_arr(const_ndofs_arr, WTVAR, 0, 1, 1); //WT = straight (0,1)-form ie v at w pts
  
  //dens primal grid reconstruction stuff- dens0, edgerecon, recon
  aux_topo_arr[DENSRECONVAR] = &dtopo;
  aux_topo_arr[DENSEDGERECONVAR] = &dtopo;
  aux_topo_arr[DENSVERTRECONVAR] = &dtopo;
  aux_topo_arr[DENSVERTEDGERECONVAR] = &dtopo;
  aux_names_arr[DENSVERTRECONVAR] = "densvertrecon";
  aux_names_arr[DENSVERTEDGERECONVAR] = "densvertedgerecon";
  aux_names_arr[DENSRECONVAR] = "densrecon";
  aux_names_arr[DENSEDGERECONVAR] = "densedgerecon";
  set_dofs_arr(aux_ndofs_arr, DENSRECONVAR, ndims-1, 1, ntdofs);  //densrecon lives on horiz dual edges, associated with F
  set_dofs_arr(aux_ndofs_arr, DENSEDGERECONVAR, ndims, 1, 2*ndims*ntdofs); //densedgerecon lives on dual cells, associated with F
  set_dofs_arr(aux_ndofs_arr, DENSVERTRECONVAR, ndims, 0, ntdofs);  //densvertrecon lives on vert dual edges, associated with Fw
  set_dofs_arr(aux_ndofs_arr, DENSVERTEDGERECONVAR, ndims, 1, 2*ntdofs); //densvertedgerecon lives on dual cells, associated with Fw

  aux_topo_arr[DENS0VAR] = &ptopo;
  aux_names_arr[DENS0VAR] = "dens0";
  set_dofs_arr(aux_ndofs_arr, DENS0VAR, 0, 0, ntdofs); //dens0 = straight (0,0)-form
  
  //Q stuff
  aux_topo_arr[QXZ0VAR] = &dtopo;
  aux_names_arr[QXZ0VAR] = "QXZ0";
  set_dofs_arr(aux_ndofs_arr, QXZ0VAR, 0, 0, nQdofs); //Q0 = twisted (0,0)-form

  aux_topo_arr[QXZRECONVAR] = &ptopo;
  aux_topo_arr[QXZEDGERECONVAR] = &ptopo;
  aux_topo_arr[QXZVERTRECONVAR] = &ptopo;
  aux_topo_arr[QXZVERTEDGERECONVAR] = &ptopo;
  aux_topo_arr[QXZFLUXVAR] = &ptopo;
  aux_topo_arr[QXZVERTFLUXVAR] = &ptopo;
  aux_names_arr[QXZRECONVAR] = "qxzrecon";
  aux_names_arr[QXZEDGERECONVAR] = "qxzedgerecon";
  aux_names_arr[QXZVERTRECONVAR] = "qxzvertrecon";
  aux_names_arr[QXZVERTEDGERECONVAR] = "qxzvertedgerecon";
  aux_names_arr[QXZFLUXVAR] = "qxzflux";
  aux_names_arr[QXZVERTFLUXVAR] = "qxzvertflux";
  set_dofs_arr(aux_ndofs_arr, QXZRECONVAR, 0, 1, nQdofs);  //qxzrecon lives on vert primal edges, associated with w
  set_dofs_arr(aux_ndofs_arr, QXZEDGERECONVAR, ndims, 1, 2*nQdofs); //qxzedgerecon lives on primal cells, associated with Fw/w
  set_dofs_arr(aux_ndofs_arr, QXZVERTRECONVAR, 1, 0, nQdofs);  //qxzsvertrecon lives on horiz primal edges, associated with v
  set_dofs_arr(aux_ndofs_arr, QXZVERTEDGERECONVAR, ndims, 1, 2*nQdofs); //qxzvertedgerecon lives on primal cells, associated with F/v
  set_dofs_arr(aux_ndofs_arr, QXZFLUXVAR, 0, 1, nQdofs); //qxzflux lives on vert primal edges, associated with w
  set_dofs_arr(aux_ndofs_arr, QXZVERTFLUXVAR, 1, 0, nQdofs); //qxzvertflux lives on horiz primal edges, associated with v

//TEMPORARY DIAGS...
  diag_topo_arr[QXZRECONDIAGVAR] = &ptopo;
  diag_topo_arr[QXZVERTRECONDIAGVAR] = &ptopo;
  diag_names_arr[QXZRECONDIAGVAR] = "qxzrecon";
  diag_names_arr[QXZVERTRECONDIAGVAR] = "qxzvertrecon";
  set_dofs_arr(diag_ndofs_arr, QXZRECONDIAGVAR, 0, 1, nQdofs);  //qxzrecon lives on vert primal edges, associated with w
  set_dofs_arr(diag_ndofs_arr, QXZVERTRECONDIAGVAR, 1, 0, nQdofs);  //qxzsvertrecon lives on horiz primal edges, associated with v
  
  diag_topo_arr[QXZEDGERECONDIAGVAR] = &ptopo;
  diag_topo_arr[QXZVERTEDGERECONDIAGVAR] = &ptopo;
  diag_names_arr[QXZEDGERECONDIAGVAR] = "qxzedgerecon";
  diag_names_arr[QXZVERTEDGERECONDIAGVAR] = "qxzvertedgerecon";
  set_dofs_arr(diag_ndofs_arr, QXZEDGERECONDIAGVAR, ndims, 1, 2*nQdofs); //qxzedgerecon lives on primal cells, associated with Fw/w
  set_dofs_arr(diag_ndofs_arr, QXZVERTEDGERECONDIAGVAR, ndims, 1, 2*nQdofs); //qxzvertedgerecon lives on primal cells, associated with F/v

  diag_topo_arr[QXZFLUXDIAGVAR] = &ptopo;
  diag_topo_arr[QXZVERTFLUXDIAGVAR] = &ptopo;
  diag_names_arr[QXZFLUXDIAGVAR] = "qxzflux";
  diag_names_arr[QXZVERTFLUXDIAGVAR] = "qxzvertflux";
  set_dofs_arr(diag_ndofs_arr, QXZFLUXDIAGVAR, 0, 1, nQdofs); //qxzflux lives on vert primal edges, associated with w
  set_dofs_arr(diag_ndofs_arr, QXZVERTFLUXDIAGVAR, 1, 0, nQdofs); //qxzvertflux lives on horiz primal edges, associated with v

  const_topo_arr[D2VAR] = &ptopo;
  const_names_arr[D2VAR] = "D2";
  set_dofs_arr(const_ndofs_arr, D2VAR, ndims, 1, 1); //Q = straight (n,1)-form
  diag_topo_arr[D2DIAGVAR] = &ptopo;
  diag_names_arr[D2DIAGVAR] = "D2Q";
  set_dofs_arr(diag_ndofs_arr, D2DIAGVAR, ndims, 1, nQdofs); //Q = straight (n,1)-form
  
  //fct stuff- Phi, Mf, edgeflux
  aux_topo_arr[PHIVAR] = &dtopo;
  aux_topo_arr[PHIVERTVAR] = &dtopo;
  aux_topo_arr[MFVAR] = &dtopo;
  aux_topo_arr[EDGEFLUXVAR] = &dtopo;
  aux_topo_arr[VERTEDGEFLUXVAR] = &dtopo;
  aux_names_arr[PHIVAR] = "Phi";
  aux_names_arr[PHIVERTVAR] = "PhiVert";
  aux_names_arr[MFVAR] = "Mf";
  aux_names_arr[EDGEFLUXVAR] = "edgeflux";
  aux_names_arr[VERTEDGEFLUXVAR] = "vertedgeflux";
  set_dofs_arr(aux_ndofs_arr, PHIVAR, ndims-1, 1, ntdofs); 
  set_dofs_arr(aux_ndofs_arr, PHIVERTVAR, ndims, 0, ntdofs); 
  set_dofs_arr(aux_ndofs_arr, MFVAR, ndims, 1, ntdofs);  
  set_dofs_arr(aux_ndofs_arr, EDGEFLUXVAR, ndims-1, 1, ntdofs); 
  set_dofs_arr(aux_ndofs_arr, VERTEDGEFLUXVAR, ndims, 0, ntdofs); 
  
  // concentration 0-forms for dens, Q
  diag_topo_arr[DENSLDIAGVAR] = &ptopo;
  diag_names_arr[DENSLDIAGVAR] = "densl";
  set_dofs_arr(diag_ndofs_arr, DENSLDIAGVAR, 0, 0, ntdofs); //densldiag = straight (0,0)-form
  
  diag_topo_arr[QXZDIAGVAR] = &dtopo;
  diag_names_arr[QXZDIAGVAR] = "QXZl";
  set_dofs_arr(diag_ndofs_arr, QXZDIAGVAR, 0, 0, nQdofs); //Qldiag = twisted (0,0)-form
  
}

#endif
