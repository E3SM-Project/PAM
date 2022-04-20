#pragma once

#include "common.h"
#include "stats.h"
#include "model.h"

#include "ext_deriv.h"
#include "hodge_star.h"
#include "fct.h"
#include "recon.h"
#include "wedge.h"


// *******   Diagnostics   ***********//

class ModelDiagnostics: public Diagnostics {
public:

 void compute_diag(ModelParameters &params, real time, const VariableSet<nconstant> &const_vars, VariableSet<nprognostic> &x, VariableSet<ndiagnostic> &diagnostic_vars)
 {
   // Compute T0
   int pis = primal_topology->is;
   int pjs = primal_topology->js;
   int pks = primal_topology->ks;
   
     parallel_for("Compute T0", Bounds<4>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
   compute_I<ntdofs, diff_ord>(diagnostic_vars.fields_arr[DENSLDIAGVAR].data, x.fields_arr[DENSVAR].data, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k, n);
   });


   // Compute q0
   int dis = dual_topology->is;
   int djs = dual_topology->js;
   int dks = dual_topology->ks;
   
     parallel_for("Compute Q0", Bounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
   compute_J<nQdofs, diff_ord>(diagnostic_vars.fields_arr[QXZDIAGVAR].data, x.fields_arr[QXZVAR].data, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k, n);
   });
}

};

// *******   Tendencies   ***********//

class ModelTendencies: public ExtrudedTendencies {
public:


     void compute_constants(VariableSet<nconstant> &const_vars, VariableSet<nprognostic> &x)
     {
      
     parallel_for("Compute Uvar", Bounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
     compute_H<1, diff_ord>(const_vars.fields_arr[UVAR].data, const_vars.fields_arr[VVAR].data, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k, n);
     });
     this->const_exchange->exchanges_arr[UVAR].exchange_field(const_vars.fields_arr[UVAR]);
       
       
     //compute UT = W U
     parallel_for("ComputeUTVAR", Bounds<4>(primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens), YAKL_LAMBDA(int k, int j, int i, int n) {
 compute_W(const_vars.fields_arr[UTVAR].data,const_vars.fields_arr[UVAR].data, pis, pjs, pks, i, j, k, n);
  });
     this->const_exchange->exchanges_arr[UTVAR].exchange_field(const_vars.fields_arr[UTVAR]);
     
     }

     void YAKL_INLINE compute_rhs(real dt, VariableSet<nconstant> &const_vars, VariableSet<nprognostic> &x, VariableSet<nauxiliary> &auxiliary_vars, VariableSet<nprognostic> &xtend)
     {
       
       int pis = primal_topology->is;
       int pjs = primal_topology->js;
       int pks = primal_topology->ks;

       int dis = dual_topology->is;
       int djs = dual_topology->js;
       int dks = dual_topology->ks;
       
       // Compute Q0
       parallel_for("Compute Q0", Bounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
       compute_J<nQdofs, diff_ord>(auxiliary_vars.fields_arr[Q0VAR].data, x.fields_arr[QVAR].data, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k, n);
       });

       // Compute T0
       parallel_for("Compute T0", Bounds<4>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
       compute_I<ntdofs, diff_ord>(auxiliary_vars.fields_arr[T0VAR].data, x.fields_arr[TVAR].data, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k, n);
       });

       this->aux_exchange->exchanges_arr[T0VAR].exchange_field(auxiliary_vars.fields_arr[T0VAR]);
       this->aux_exchange->exchanges_arr[Q0VAR].exchange_field(auxiliary_vars.fields_arr[Q0VAR]);

       //compute qrecon
       parallel_for("ComputeQEdgeRecon", Bounds<4>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
          compute_straight_edge_recon<nQdofs, reconstruction_type, reconstruction_order>(
            auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[Q0VAR].data, pis, pjs, pks, i, j, k, n,
            primal_wenoRecon, primal_to_gll, primal_wenoIdl, primal_wenoSigma);
       });

       this->aux_exchange->exchanges_arr[QEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[QEDGERECONVAR]);

       parallel_for("ComputeQRECON", Bounds<4>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
         compute_straight_recon<nQdofs, reconstruction_type>(
           auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[QEDGERECONVAR].data, const_vars.fields_arr[WTVAR].data, pis, pjs, pks, i, j, k, n);
       });
       this->aux_exchange->exchanges_arr[QRECONVAR].exchange_field(auxiliary_vars.fields_arr[QRECONVAR]);

       //compute Trecon
         parallel_for("ComputeTEdgeRecon", Bounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
         compute_twisted_edge_recon<ntdofs, dual_reconstruction_type, dual_reconstruction_order>(
           auxiliary_vars.fields_arr[TEDGERECONVAR].data, auxiliary_vars.fields_arr[T0VAR].data, dis, djs, dks, i, j, k, n,
           dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
         });
       this->aux_exchange->exchanges_arr[TEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[TEDGERECONVAR]);

         parallel_for("ComputeTRecon", Bounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
         compute_twisted_recon<ntdofs, dual_reconstruction_type>(
           auxiliary_vars.fields_arr[TRECONVAR].data, auxiliary_vars.fields_arr[TEDGERECONVAR].data, const_vars.fields_arr[UVAR].data, dis, djs, dks, i, j, k, n);
         });
       this->aux_exchange->exchanges_arr[TRECONVAR].exchange_field(auxiliary_vars.fields_arr[TRECONVAR]);

       // compute FCT stuff
               parallel_for("ComputeEdgeFlux", Bounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
             compute_edgefluxes<ntdofs> (auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[TRECONVAR].data, const_vars.fields_arr[UVAR].data, dis, djs, dks, i, j, k, n);
             });
       this->aux_exchange->exchanges_arr[EDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[EDGEFLUXVAR]);


         parallel_for("ComputeMf", Bounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
         compute_Mf<ntdofs> (auxiliary_vars.fields_arr[MFVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, dt, dis, djs, dks, i, j, k, n);
       });
           this->aux_exchange->exchanges_arr[MFVAR].exchange_field(auxiliary_vars.fields_arr[MFVAR]);

             parallel_for("ComputePhi", Bounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
             compute_Phi<ntdofs> (auxiliary_vars.fields_arr[PHIVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k, n);
           });
           
       //Don't do FCT for non-FCT vars
           for (int l=0; l<ntnofctdofs; l++)
           {
           auxiliary_vars.fields_arr[PHIVAR].set(l, 1.0);
           auxiliary_vars.fields_arr[PHIVERTVAR].set(l, 1.0);
           }
           this->aux_exchange->exchanges_arr[PHIVAR].exchange_field(auxiliary_vars.fields_arr[PHIVAR]);

           // //compute Q flux
                parallel_for("ComputeQFLUX", Bounds<4>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
                if (qf_choice == QF_MODE::EC)
               { compute_Q_EC<nQdofs>(auxiliary_vars.fields_arr[QFLUXVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, const_vars.fields_arr[UVAR].data, pis, pjs, pks, i, j, k, n);}
                if (qf_choice == QF_MODE::NOEC)
               { compute_Q_nonEC<nQdofs>(auxiliary_vars.fields_arr[QFLUXVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, const_vars.fields_arr[UVAR].data, pis, pjs, pks, i, j, k, n);}
           });
               this->aux_exchange->exchanges_arr[QFLUXVAR].exchange_field(auxiliary_vars.fields_arr[QFLUXVAR]);


           // // Compute D2 QFLUX
           parallel_for("ComputeQTend", Bounds<4>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
           compute_D2<nQdofs> (xtend.fields_arr[QVAR].data, auxiliary_vars.fields_arr[QFLUXVAR].data, pis, pjs, pks, i, j, k, n);
           });

           // compute D2bar "F"
                 parallel_for("ComputeTTend", Bounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
                 compute_wDbar2_fct<ntdofs> (xtend.fields_arr[DENSVAR].data, auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[PHIVAR].data, const_vars.fields_arr[UVAR].data, dis, djs, dks, i, j, k, n);
               });

     }
     
};



// *******   Statistics   ***********//

class ModelStats: public Stats {
public:
  
  real3d trimmed_dual_density, trimmed_primal_density;
   
  void initialize(Parameters &params, Parallel &par, const Topology &primal_topo, const Topology &dual_topo, Geometry &primal_geom, Geometry &dual_geom)
{
  Stats::initialize(params, par, primal_topo, dual_topo, primal_geom, dual_geom);
  this->stats_arr[TMASSSTAT].initialize("Tmass", ntdofs, this->statsize, this->nens, this->masterproc);
  this->stats_arr[TMAXSTAT].initialize("Tmax", ntdofs, this->statsize, this->nens, this->masterproc);
  this->stats_arr[TMINSTAT].initialize("Tmin", ntdofs, this->statsize, this->nens, this->masterproc);
  this->stats_arr[QMASSSTAT].initialize("Qmass", nQdofs, this->statsize, this->nens, this->masterproc);

  this->trimmed_dual_density = real3d("trimmed_dual_density", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
  this->trimmed_primal_density = real3d("trimmed_primal_density", this->primal_topology->nl, this->primal_topology->n_cells_y, this->primal_topology->n_cells_x);

}

   void compute( VariableSet<nprognostic> &progvars,  VariableSet<nconstant> &constvars, int tind)
   {

    for (int n=0;n<nens;n++)
    {
      SArray<real,1,ntdofs> masslocal, massglobal;
      SArray<real,1,ntdofs> Tmaxlocal, Tmaxglobal;
      SArray<real,1,ntdofs> Tminlocal, Tminglobal;
      SArray<real,1,nQdofs> Qlocal, Qglobal;

      for (int l=0;l<ntdofs;l++) {masslocal(l) = 0.; massglobal(l) = 0.;}
      for (int l=0;l<ntdofs;l++) {Tmaxlocal(l) = 0.; Tmaxglobal(l) = 0.;}
      for (int l=0;l<ntdofs;l++) {Tminlocal(l) = 0.; Tminglobal(l) = 0.;}
  
      int dis = dual_topology->is;
      int djs = dual_topology->js;
      int dks = dual_topology->ks;
      
    for (int l=0;l<ntdofs;l++)
    {
      parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
        trimmed_dual_density(k,j,i) = progvars.fields_arr[DENSVAR].data(l,k+dks,j+djs,i+dis,n);
      });

    masslocal(l) = yakl::intrinsics::sum(trimmed_dual_density);
    Tmaxlocal(l) = yakl::intrinsics::maxval(trimmed_dual_density);
    Tminlocal(l) = yakl::intrinsics::minval(trimmed_dual_density);
    }

    int pis = primal_topology->is;
    int pjs = primal_topology->js;
    int pks = primal_topology->ks;
    
    for (int l=0;l<nQdofs;l++)
    {
      parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
        trimmed_primal_density(k,j,i) = progvars.fields_arr[QXZVAR].data(l,k+pks,j+pjs,i+pis,n);
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
  this->stats_arr[DENSSTAT].data(l,tind,n) = massglobal(l);
  this->stats_arr[DENSMAXSTAT].data(l,tind,n) = Tmaxglobal(l);
  this->stats_arr[DENSMINSTAT].data(l,tind,n) = Tminglobal(l);
}

for (int l=0;l<nQdofs;l++)
{
this->stats_arr[QXZSTAT].data(l,tind) = Qglobal(l);
}
    }

}
};


// *******   VariableSet Initialization   ***********//
void initialize_variables(const Topology &ptopo, const Topology &dtopo,
SArray<int,2, nprognostic, 3> &prog_ndofs_arr, SArray<int,2, nconstant, 3> &const_ndofs_arr, SArray<int,2, nauxiliary, 3> &aux_ndofs_arr, SArray<int,2, ndiagnostic, 3> &diag_ndofs_arr,
std::array<std::string, nprognostic> &prog_names_arr, std::array<std::string, nconstant> &const_names_arr, std::array<std::string, nauxiliary> &aux_names_arr, std::array<std::string, ndiagnostic> &diag_names_arr,
std::array<const Topology *, nprognostic> &prog_topo_arr, std::array<const Topology *, nconstant> &const_topo_arr, std::array<const Topology *, nauxiliary> &aux_topo_arr, std::array<const Topology *, ndiagnostic> &diag_topo_arr)
{

  set_dofs_arr(prog_ndofs_arr, DENSVAR, ndims, 1, ndensity); //dens = twisted n-form

  prog_topo_arr[TVAR] = &dtopo;
  prog_topo_arr[QVAR] = &ptopo;

  const_topo_arr[VVAR] = &ptopo;
  const_topo_arr[UVAR] = &dtopo;
  const_topo_arr[UTVAR] = &ptopo;

  aux_topo_arr[T0VAR] = &ptopo;
  aux_topo_arr[TRECONVAR] = &dtopo;
  aux_topo_arr[TEDGERECONVAR] = &dtopo;
  aux_topo_arr[PHIVAR] = &dtopo;
  aux_topo_arr[MFVAR] = &dtopo;
  aux_topo_arr[EDGEFLUXVAR] = &dtopo;
  aux_topo_arr[Q0VAR] = &dtopo;
  aux_topo_arr[QRECONVAR] = &ptopo;
  aux_topo_arr[QEDGERECONVAR] = &ptopo;
  aux_topo_arr[QFLUXVAR] = &ptopo;

  diag_topo_arr[TDIAGVAR] = &ptopo;
  diag_topo_arr[QDIAGVAR] = &dtopo;

  prog_names_arr[TVAR] = "T";
  prog_names_arr[QVAR] = "Q";

  const_names_arr[VVAR] = "v";
  const_names_arr[UVAR] = "U";
  const_names_arr[UTVAR] = "UT";

  aux_names_arr[T0VAR] = "T0";
  aux_names_arr[TRECONVAR] = "Trecon";
  aux_names_arr[TEDGERECONVAR] = "Tedgerecon";
  aux_names_arr[PHIVAR] = "Phi";
  aux_names_arr[MFVAR] = "Mf";
  aux_names_arr[EDGEFLUXVAR] = "edgeflux";
  aux_names_arr[Q0VAR] = "q0";
  aux_names_arr[QRECONVAR] = "qrecon";
  aux_names_arr[QEDGERECONVAR] = "qedgerecon";
  aux_names_arr[QFLUXVAR] = "qflux";

  diag_names_arr[TDIAGVAR] = "T0";
  diag_names_arr[QDIAGVAR] = "Q0";

  set_dofs_arr(prog_ndofs_arr, TVAR, ndims, 1, ntdofs); //dens = twisted n-form
  set_dofs_arr(prog_ndofs_arr, QVAR, ndims, 1, nQdofs); //Q = straight n-form
  
  set_dofs_arr(const_ndofs_arr, VVAR, 1, 1, 1); //v = straight 1-form
  set_dofs_arr(const_ndofs_arr, UVAR, ndims-1, 1, 1); //U = twisted (n-1)-form
  set_dofs_arr(const_ndofs_arr, UTVAR, 1, 1, 1); //WT = straight 1-form
  

  set_dofs_arr(aux_ndofs_arr, TRECONVAR, ndims-1, 1, ntdofs);  //Trecon lives on dual edges, associated with F
  set_dofs_arr(aux_ndofs_arr, TEDGERECONVAR, ndims, 1, 2*ndims*ntdofs); //Tedgerecon lives on dual cells, associated with F
  set_dofs_arr(aux_ndofs_arr, T0VAR, 0, 1, ntdofs); //T0 = straight 0-form
  set_dofs_arr(aux_ndofs_arr, Q0VAR, 0, 1, nQdofs); //Q0 = twisted 0-form
  set_dofs_arr(aux_ndofs_arr, QRECONVAR, 1, 1, nQdofs);  //qrecon lives on primal edges, associated with UT
  set_dofs_arr(aux_ndofs_arr, QEDGERECONVAR, ndims, 1, 2*nQdofs); //qedgerecon lives on primal cells, associated with UT
  set_dofs_arr(aux_ndofs_arr, QFLUXVAR, 1, 1, nQdofs); //qxzflux lives on primal edges, associated with UT
  set_dofs_arr(aux_ndofs_arr, PHIVAR, ndims-1, 1, ntdofs); 
  set_dofs_arr(aux_ndofs_arr, MFVAR, ndims, 1, ntdofs);  
  set_dofs_arr(aux_ndofs_arr, EDGEFLUXVAR, ndims-1, 1, ntdofs); 
  

  set_dofs_arr(diag_ndofs_arr, TDIAGVAR, 0, 1, ntdofs); //Tdiag = straight 0-form
  set_dofs_arr(diag_ndofs_arr, QXZDIAGVAR, 0, 1, nQdofs); //Qdiag = twisted 0-form

}


//***************** Set Initial Conditions ***************************//

#define gaussian_2d(x,y)   (1._fp * exp(-100._fp * pow(x-0.5_fp ,2._fp)) * exp(-100._fp * pow(y-0.5_fp ,2._fp)))
real YAKL_INLINE gaussian(real x, real y)         { return gaussian_2d(x,y); }
#define vortex1_2d(x,y)   (1._fp *  exp(-100._fp * pow(x-0.75_fp ,2._fp)) * exp(-100._fp * pow(y-0.75_fp ,2._fp)))
#define vortex2_2d(x,y)   (0.5_fp * exp(-50._fp  * pow(x-0.25_fp ,2._fp)) * exp(-75._fp  * pow(y-0.25_fp , 2._fp)))
real YAKL_INLINE vortices(real x, real y)         { return vortex1_2d(x,y)   + vortex2_2d(x,y); }
real YAKL_INLINE square(real x, real y)         {return (x > 0.35_fp && x < 0.65_fp && y > 0.35_fp && y < 0.65_fp                        ) ? 1._fp : 0._fp;}
real YAKL_INLINE square_ur(real x, real y)         {return (x > 0.6_fp && x < 0.9_fp && y > 0.6_fp && y < 0.9_fp                        ) ? 1._fp : 0._fp;}
real YAKL_INLINE square_ll(real x, real y)         {return (x > 0.1_fp && x < 0.4_fp && y > 0.1_fp && y < 0.4_fp                        ) ? 1._fp : 0._fp;}
real YAKL_INLINE doublesquare(real x, real y)         {return square_ur(x,y) + square_ll(x,y);}



#define C_UNIFORM_WIND 1.

real YAKL_INLINE uniform_x_wind(real x) {
  return C_UNIFORM_WIND;
}
vec<2> YAKL_INLINE uniform_x_wind(real x, real y) {
  vec<2> vvec;
  vvec.u = C_UNIFORM_WIND;
  return vvec;
}



vec<2> YAKL_INLINE uniform_y_wind(real x, real y) {
  vec<2> vvec;
  vvec.v = -C_UNIFORM_WIND;
  return vvec;
}

vec<2> YAKL_INLINE uniform_xy_wind(real x, real y) {
  vec<2> vvec;
  vvec.u = -C_UNIFORM_WIND/sqrt(2._fp);
  vvec.v = C_UNIFORM_WIND/sqrt(2._fp);
  return vvec;
}

// FIX THIS
vec<2> YAKL_INLINE deformational_wind(real x, real y) {
  vec<2> vvec;
  vvec.v = C_UNIFORM_WIND;
  return vvec;
}

vec<2> YAKL_INLINE doublevortex_wind(real x, real y) {
vec<2> vvec;

real xprime1 = 1.0_fp / (pi * 3._fp / 40._fp ) * sin(pi / 1.0_fp * (x - 0.4_fp ));
real yprime1 = 1.0_fp / (pi * 3._fp / 40._fp ) * sin(pi / 1.0_fp * (y - 0.4_fp ));
real xprime2 = 1.0_fp / (pi * 3._fp / 40._fp ) * sin(pi / 1.0_fp * (x - 0.6_fp ));
real yprime2 = 1.0_fp / (pi * 3._fp / 40._fp ) * sin(pi / 1.0_fp * (y - 0.6_fp ));
real xprimeprime1 = 1.0_fp / (2.0_fp * pi * 3._fp / 40._fp ) * sin(2.0_fp * pi / 1.0_fp * (x - 0.4_fp ));
real yprimeprime1 = 1.0_fp / (2.0_fp * pi * 3._fp / 40._fp ) * sin(2.0_fp * pi / 1.0_fp * (y - 0.4_fp ));
real xprimeprime2 = 1.0_fp / (2.0_fp * pi * 3._fp / 40._fp ) * sin(2.0_fp * pi / 1.0_fp * (x - 0.6_fp ));
real yprimeprime2 = 1.0_fp / (2.0_fp * pi * 3._fp / 40._fp ) * sin(2.0_fp * pi / 1.0_fp * (y - 0.6_fp ));

vvec.u = -1.0_fp * (yprimeprime1 * exp(-0.5_fp *(xprime1 * xprime1 + yprime1 * yprime1)) + yprimeprime2 * exp(-0.5_fp *(xprime2 * xprime2 + yprime2 * yprime2)));
vvec.v = 1.0_fp * (xprimeprime1 * exp(-0.5_fp *(xprime1 * xprime1 + yprime1 * yprime1)) + xprimeprime2 * exp(-0.5_fp *(xprime2 * xprime2 + yprime2 * yprime2)));
return vvec;
}


void readModelParamsFile(std::string inFile, ModelParameters &params, Parallel &par, int nz)
{
  readParamsFile( inFile, params, par, nz);
  
  // Read the data initialization options
  params.windInitStr = config["initWind"].as<std::string>();
  for (int i=0;i<ntnofctdofs;i++)
  {params.TInitStr[i] = config["initT" + std::to_string(i+1)].as<std::string>();}
  for (int i=ntnofctdofs;i<ntdofs;i++)
  {params.TInitStr[i] = config["initTFCT" + std::to_string(i-ntnofctdofs+1)].as<std::string>();}
  for (int i=0;i<nQdofs;i++)
  {params.QInitStr[i] = config["initQ" + std::to_string(i+1)].as<std::string>();}  
}

void set_domain_sizes_ic (ModelParameters &params, std::string initData)
{
  params.xlen = 1.0_fp;
  params.ylen = 1.0_fp;
  params.xc = 0.5_fp;
  params.yc = 0.5_fp;
}

void set_initial_conditions (ModelParameters &params, VariableSet<nprognostic> &progvars, VariableSet<nconstant> &constvars, 
Geometry &primal_geom, Geometry &dual_geom)
{

  if (params.windInitStr == "uniform_x"    ) {pgeom.set_1form_values(uniform_x_wind,     constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);}
  if (params.windInitStr == "uniform_y"    ) {pgeom.set_1form_values(uniform_y_wind,     constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);}
  if (params.windInitStr == "uniform_xy"   ) {pgeom.set_1form_values(uniform_xy_wind,     constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);}
  if (params.windInitStr == "deformational") {pgeom.set_1form_values(deformational_wind, constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);}
  if (params.windInitStr == "doublevortex" ) {pgeom.set_1form_values(doublevortex_wind, constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);}
  
  for (int i=0; i<ntdofs; i++)
  {
  if (params.TInitStr[i] == "gaussian") {dgeom.set_2form_values(gaussian, progvars.fields_arr[TVAR], i);}
  if (params.TInitStr[i] == "vortices") {dgeom.set_2form_values(vortices, progvars.fields_arr[TVAR], i);}
  if (params.TInitStr[i] == "square")   {dgeom.set_2form_values(square,   progvars.fields_arr[TVAR], i);}
  if (params.TInitStr[i] == "doublesquare")   {dgeom.set_2form_values(doublesquare,   progvars.fields_arr[TVAR], i);}
  }
  for (int i=0; i<nQdofs; i++)
  {
    if (params.QInitStr[i] == "gaussian") {pgeom.set_2form_values(gaussian, progvars.fields_arr[QVAR], i);}
    if (params.QInitStr[i] == "vortices") {pgeom.set_2form_values(vortices, progvars.fields_arr[QVAR], i);}
    if (params.QInitStr[i] == "square")   {pgeom.set_2form_values(square,   progvars.fields_arr[QVAR], i);}
    if (params.QInitStr[i] == "doublesquare")   {pgeom.set_2form_values(doublesquare,   progvars.fields_arr[QVAR], i);}
  }
  
  
}
