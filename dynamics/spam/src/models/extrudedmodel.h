#pragma once

#include "common.h"
#include "stats.h"
#include "model.h"

#include "ext_deriv.h"
#include "hodge_star.h"
#include "fct.h"
#include "recon.h"
#include "wedge.h"
#include "hamiltonian.h"
#include "variableset.h"

// *******   Functionals/Hamiltonians   ***********//

Functional_PVPE_extruded PVPE;
Hamiltonian_Hk_extruded Hk;

#ifdef _SWE
Hamiltonian_SWE_Hs Hs;
VariableSet_SWE varset;
#elif _TSWE
Hamiltonian_TSWE_Hs Hs;
VariableSet_TSWE varset;
#elif _CE
Hamiltonian_CE_Hs Hs;
VariableSet_CE varset;
#elif _MCErho
Hamiltonian_MCE_Hs Hs;
VariableSet_MCE_rho varset;
#elif _MCErhod
Hamiltonian_MCE_Hs Hs;
VariableSet_MCE_rhod varset;
#elif _CEp
Hamiltonian_CE_p_Hs Hs;
VariableSet_CE_p varset;
#elif _MCErhop
Hamiltonian_MCE_p_Hs Hs;
VariableSet_MCE_rhop varset;
#elif _MCErhodp
Hamiltonian_MCE_p_Hs Hs;
VariableSet_MCE_rhodp varset;
#endif
//ADD ANELASTIC + MOIST ANELASTIC

#ifdef _THERMONONE
ThermoPotential thermo;
#elif _IDEAL_GAS_POTTEMP
IdealGas_Pottemp thermo;
#elif _IDEAL_GAS_ENTROPY
IdealGas_Entropy thermo;
#elif _CONST_KAPPA_VIRPOTTEMP
ConstantKappa_VirtualPottemp thermo;
#elif _UNAPPROX_POTTEMP
Unapprox_Pottemp thermo;
#elif _UNAPPROX_ENTROPY
Unapprox_Entropy thermo;
#endif

// *******   Diagnostics   ***********//

class ModelDiagnostics: public Diagnostics {
public:

 void compute_diag(const FieldSet<nconstant> &const_vars, FieldSet<nprognostic> &x, FieldSet<ndiagnostic> &diagnostic_vars)
 {

  int dis = dual_topology->is;
  int djs = dual_topology->js;
  int dks = dual_topology->ks;

  int pis = primal_topology->is;
  int pjs = primal_topology->js;
  int pks = primal_topology->ks;

  parallel_for("Compute DENS0 DIAG", SimpleBounds<4>( primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
compute_Iext<ndensity, diff_ord, vert_diff_ord>(diagnostic_vars.fields_arr[DENSLDIAGVAR].data, x.fields_arr[DENSVAR].data, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k, n);
});

  parallel_for("Compute Q0 DIAG", SimpleBounds<4>( dual_topology->ni-3, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
PVPE.compute_qxz0(diagnostic_vars.fields_arr[QXZDIAGVAR].data, x.fields_arr[VVAR].data, x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data, const_vars.fields_arr[CORIOLISXZVAR].data, dis, djs, dks, i, j, k+2, n);
});

//Bottom is k=1 and top is k=dual_topology->ni-2
parallel_for("Compute Q0 DIAG TOP/BOTTOM", SimpleBounds<3>(dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int j, int i, int n) { 
PVPE.compute_qxz0_bottom(diagnostic_vars.fields_arr[QXZDIAGVAR].data, x.fields_arr[VVAR].data, x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data, const_vars.fields_arr[CORIOLISXZVAR].data, dis, djs, dks, i, j, 1, n);
PVPE.compute_qxz0_top(diagnostic_vars.fields_arr[QXZDIAGVAR].data, x.fields_arr[VVAR].data, x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data, const_vars.fields_arr[CORIOLISXZVAR].data, dis, djs, dks, i, j, dual_topology->ni-2, n);
});

  diagnostic_vars.fields_arr[QXZDIAGVAR].set_bnd(0.0);
}

};

// *******   Tendencies   ***********//

class ModelTendencies: public ExtrudedTendencies {
public:

 void initialize(PamCoupler &coupler, ModelParameters &params, Topology &primal_topo, Topology &dual_topo, Geometry &primal_geom, Geometry &dual_geom, ExchangeSet<nauxiliary> &aux_exchange, ExchangeSet<nconstant> &const_exchange)
 {

ExtrudedTendencies::initialize(params, primal_topo, dual_topo, primal_geom, dual_geom, aux_exchange, const_exchange);
   varset.initialize(coupler, params, thermo,  *this->primal_topology, *this->dual_topology, *this->primal_geometry, *this->dual_geometry);
   PVPE.initialize(varset);
   Hk.initialize(varset, *this->primal_geometry, *this->dual_geometry);
   Hs.initialize(thermo, varset, *this->primal_geometry, *this->dual_geometry);
}

void convert_dynamics_to_coupler_state(PamCoupler &coupler, const FieldSet<nprognostic> &prog_vars, const FieldSet<nconstant> &const_vars)
{varset.convert_dynamics_to_coupler_state(coupler, prog_vars, const_vars);}
void convert_coupler_to_dynamics_state(PamCoupler &coupler, FieldSet<nprognostic> &prog_vars, const FieldSet<nconstant> &const_vars)
{varset.convert_coupler_to_dynamics_state(coupler, prog_vars, const_vars);}

     void compute_constants(FieldSet<nconstant> &const_vars, FieldSet<nprognostic> &x)
     {}

       void YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_I(
        real5d Uvar, real5d UWvar, real5d qxz0var, real5d fxz0var, real5d dens0var, 
        const real5d Vvar, const real5d Wvar, const real5d densvar, const real5d coriolisxzvar) {

    int pis = primal_topology->is;
    int pjs = primal_topology->js;
    int pks = primal_topology->ks;

    int dis = dual_topology->is;
    int djs = dual_topology->js;
    int dks = dual_topology->ks;

       parallel_for("Compute Dens0var", SimpleBounds<4>( primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
    compute_Iext<ndensity, diff_ord, vert_diff_ord>(dens0var, densvar, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k, n);
    });

      parallel_for("Compute Uvar", SimpleBounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
    compute_Hext<1, diff_ord>(Uvar, Vvar, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k, n);
    });
  
      parallel_for("Compute UWVar", SimpleBounds<4>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
    compute_Hv<1, vert_diff_ord>(UWvar, Wvar, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k+1, n);
    });    

      parallel_for("Compute Q0, F0", SimpleBounds<4>( dual_topology->ni-3, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
    PVPE.compute_qxz0fxz0(qxz0var, fxz0var, Vvar, Wvar, densvar, coriolisxzvar, dis, djs, dks, i, j, k+2, n);
    });
    parallel_for("Compute Q0, F0 bnd", SimpleBounds<3>( dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int j, int i, int n) { 
    PVPE.compute_qxz0fxz0_bottom(qxz0var, fxz0var, Vvar, Wvar, densvar, coriolisxzvar, dis, djs, dks, i, j, 1, n);
    PVPE.compute_qxz0fxz0_top(qxz0var, fxz0var, Vvar, Wvar, densvar, coriolisxzvar, dis, djs, dks, i, j,  dual_topology->ni-2, n);
    });
        }


        void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_II(
          real5d Fvar, real5d FWvar, real5d Kvar, real5d HEvar, real5d HEWvar, const real5d Vvar, const real5d Uvar, const real5d Wvar, const real5d UWvar, const real5d dens0var) {

            int dis = dual_topology->is;
            int djs = dual_topology->js;
            int dks = dual_topology->ks;

            int pis = primal_topology->is;
            int pjs = primal_topology->js;
            int pks = primal_topology->ks;
            
      //THIS WILL NEED SOME SLIGHT MODIFICATIONS FOR CASE OF NON-ZERO UWVAR_B IE BOUNDARY FLUXES
      //BUT FOR NOW IT IS FINE SINCE UWVAR=0 on BND AND THEREFORE K COMPUTATIONS IGNORE IT
              parallel_for("Compute Fvar, Kvar", SimpleBounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
            Hk.compute_F(Fvar, HEvar, Uvar, dens0var, dis, djs, dks, i, j, k, n);
            Hk.compute_K(Kvar, Vvar, Uvar, Wvar, UWvar, dis, djs, dks, i, j, k, n);
          });
            parallel_for("Compute FWvar", SimpleBounds<4>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
          Hk.compute_Fw(FWvar, HEWvar, UWvar, dens0var, dis, djs, dks, i, j, k+1, n);
        });
        
          }

      void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_III(
        real5d Bvar, real5d FTvar, real5d FTWvar,
        const real5d Fvar, const real5d Uvar, const real5d FWvar, const real5d UWvar,
        const real5d Kvar, const real5d densvar, const real5d HSvar) {

    int pis = primal_topology->is;
    int pjs = primal_topology->js;
    int pks = primal_topology->ks;

    parallel_for("Compute FTvar", SimpleBounds<4>( primal_topology->ni-2, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
    compute_Wxz_u(FTvar, FWvar, pis, pjs, pks, i, j, k+1, n);
    });
    parallel_for("Compute FTvar bnd", SimpleBounds<3>( primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int j, int i, int n) { 
    compute_Wxz_u_bottom(FTvar, FWvar, pis, pjs, pks, i, j, 0, n);
    compute_Wxz_u_top(FTvar, FWvar, pis, pjs, pks, i, j, primal_topology->ni-1, n);
    });
      parallel_for("Compute FTWvar", SimpleBounds<4>( primal_topology->nl-2, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
    compute_Wxz_w(FTWvar, Fvar, pis, pjs, pks, i, j, k+1, n);
    });
    parallel_for("Compute FTWvar bnd", SimpleBounds<3>( primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int j, int i, int n) { 
    compute_Wxz_w_bottom(FTWvar, Fvar, pis, pjs, pks, i, j, 0, n);
    compute_Wxz_w_top(FTWvar, Fvar, pis, pjs, pks, i, j, primal_topology->nl-1, n);
    });
            parallel_for("Compute Bvar", SimpleBounds<4>( primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
    Hs.compute_dHsdx(Bvar, densvar, HSvar, pis, pjs, pks, i, j, k, n);
    Hk.compute_dKddens<ADD_MODE::ADD>(Bvar, Kvar, pis, pjs, pks, i, j, k, n);
        });

        }

      void YAKL_INLINE compute_edge_reconstructions(real5d densedgereconvar, real5d densvertedgereconvar, 
        real5d qxzedgereconvar, real5d qxzvertedgereconvar, 
        real5d coriolisxzedgereconvar, real5d coriolisxzvertedgereconvar,
        const real5d dens0var, const real5d qxz0var, const real5d fxz0var) {

          int dis = dual_topology->is;
          int djs = dual_topology->js;
          int dks = dual_topology->ks;

          int pis = primal_topology->is;
          int pjs = primal_topology->js;
          int pks = primal_topology->ks;
          
           parallel_for("ComputeDensEdgeRecon", SimpleBounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
            compute_twisted_edge_recon<ndensity, dual_reconstruction_type, dual_reconstruction_order>(
              densedgereconvar, dens0var, dis, djs, dks, i, j, k, n,
              dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
              compute_twisted_vert_edge_recon<ndensity, dual_vert_reconstruction_type, dual_vert_reconstruction_order>(
                densvertedgereconvar, dens0var, dis, djs, dks, i, j, k, n, 
                dual_vert_wenoRecon, dual_vert_to_gll, dual_vert_wenoIdl, dual_vert_wenoSigma);
            });

              parallel_for("ComputeQEdgeRecon",SimpleBounds<4>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
               compute_straight_xz_edge_recon<1, reconstruction_type, reconstruction_order>(
                 qxzedgereconvar, qxz0var, pis, pjs, pks, i, j, k, n,
                 primal_wenoRecon, primal_to_gll, primal_wenoIdl, primal_wenoSigma);
               compute_straight_xz_vert_edge_recon<1, vert_reconstruction_type, vert_reconstruction_order>(
                 qxzvertedgereconvar, qxz0var, pis, pjs, pks, i, j, k, n,
                 primal_vert_wenoRecon, primal_vert_to_gll, primal_vert_wenoIdl, primal_vert_wenoSigma);
                 compute_straight_xz_edge_recon<1, coriolis_reconstruction_type, coriolis_reconstruction_order>(
                   coriolisxzedgereconvar, fxz0var, pis, pjs, pks, i, j, k, n,
                   coriolis_wenoRecon, coriolis_to_gll, coriolis_wenoIdl, coriolis_wenoSigma);
                 compute_straight_xz_vert_edge_recon<1, coriolis_vert_reconstruction_type, coriolis_vert_reconstruction_order>(
                   coriolisxzvertedgereconvar, fxz0var, pis, pjs, pks, i, j, k, n,
                   coriolis_vert_wenoRecon, coriolis_vert_to_gll, coriolis_vert_wenoIdl, coriolis_vert_wenoSigma);
            });
            
      }
      
      void YAKL_INLINE compute_recons(
      real5d densreconvar,  real5d densvertreconvar,
      real5d qxzreconvar,  real5d qxzvertreconvar,
      real5d coriolisxzreconvar,  real5d coriolisxzvertreconvar,
      const real5d densedgereconvar, const real5d densvertedgereconvar, 
      const real5d qxzedgereconvar, const real5d qxzvertedgereconvar,
      const real5d coriolisxzedgereconvar, const real5d coriolisxzvertedgereconvar,
      const real5d HEvar, const real5d HEWvar,
      const real5d Uvar, const real5d UWvar,
      const real5d FTvar, const real5d FTWvar) {

    int dis = dual_topology->is;
    int djs = dual_topology->js;
    int dks = dual_topology->ks;

    int pis = primal_topology->is;
    int pjs = primal_topology->js;
    int pks = primal_topology->ks;

      parallel_for("ComputeDensRECON", SimpleBounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
      compute_twisted_recon<ndensity, dual_reconstruction_type>(
        densreconvar, densedgereconvar, Uvar, dis, djs, dks, i, j, k, n);
        //scale twisted recons
        for (int d=0;d<ndims;d++) {
        for (int l=0;l<ndensity;l++) {
        densreconvar(l+d*ndensity,k+dks,j+djs,i+dis,n) = densreconvar(l+d*ndensity,k+dks,j+djs,i+dis,n) / HEvar(d,k+dks,j+djs,i+dis,n);
      }}
      });


          parallel_for("ComputeDensVertRECON", SimpleBounds<4>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
        compute_twisted_vert_recon<ndensity, dual_vert_reconstruction_type>(
          densvertreconvar, densvertedgereconvar, UWvar, dis, djs, dks, i, j, k+1, n);
          //scale twisted recons
          for (int l=0;l<ndensity;l++) {
          densvertreconvar(l,k+dks+1,j+djs,i+dis,n) = densvertreconvar(l,k+dks+1,j+djs,i+dis,n) / HEWvar(0,k+dks+1,j+djs,i+dis,n);
        }
    });

      parallel_for("ComputeQRECON", SimpleBounds<4>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
      compute_straight_xz_recon<1, reconstruction_type>(
        qxzreconvar, qxzedgereconvar, FTWvar, pis, pjs, pks, i, j, k, n);
        compute_straight_xz_recon<1, coriolis_reconstruction_type>(
          coriolisxzreconvar, coriolisxzedgereconvar, FTWvar, pis, pjs, pks, i, j, k, n);
    });
      parallel_for("ComputeQVERTRECON",  SimpleBounds<4>( primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
      compute_straight_xz_vert_recon<1, vert_reconstruction_type>(
        qxzvertreconvar, qxzvertedgereconvar, FTvar, pis, pjs, pks, i, j, k, n);
        compute_straight_xz_vert_recon<1, coriolis_vert_reconstruction_type>(
          coriolisxzvertreconvar, coriolisxzvertedgereconvar, FTvar, pis, pjs, pks, i, j, k, n);
    });

    }

      void YAKL_INLINE compute_tendencies(
      real5d denstendvar, real5d Vtendvar, real5d Wtendvar,
      const real5d densreconvar, const real5d densvertreconvar,
      const real5d qxzreconvar, const real5d qxzvertreconvar,
      const real5d coriolisxzreconvar, const real5d coriolisxzvertreconvar,
      const real5d Bvar, const real5d Fvar, const real5d FWvar, const real5d Phivar, const real5d Phivertvar) {

    int pis = primal_topology->is;
    int pjs = primal_topology->js;
    int pks = primal_topology->ks;

    int dis = dual_topology->is;
    int djs = dual_topology->js;
    int dks = dual_topology->ks;


      parallel_for("Compute Wtend", SimpleBounds<4>( primal_topology->nl-2, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
      compute_wDv_fct<ndensity> (Wtendvar, densvertreconvar, Phivertvar, Bvar, pis, pjs, pks, i, j, k+1, n);
       if (qf_choice == QF_MODE::EC)
      { compute_Qxz_w_EC<1, ADD_MODE::ADD>(Wtendvar, qxzreconvar, qxzvertreconvar, Fvar, pis, pjs, pks, i, j, k+1, n);}
       if (qf_choice == QF_MODE::NOEC)
      { compute_Qxz_w_nonEC<1, ADD_MODE::ADD>(Wtendvar, qxzreconvar, Fvar, pis, pjs, pks, i, j, k+1, n);}
      compute_Qxz_w_EC<1, ADD_MODE::ADD>(Wtendvar, coriolisxzreconvar, coriolisxzvertreconvar, Fvar, pis, pjs, pks, i, j, k+1, n);
    });
    parallel_for("Compute Wtend Bnd", SimpleBounds<3>(primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int j, int i, int n) { 
    compute_wDv_fct<ndensity> (Wtendvar, densvertreconvar, Phivertvar, Bvar, pis, pjs, pks, i, j, 0, n);
    compute_wDv_fct<ndensity> (Wtendvar, densvertreconvar, Phivertvar, Bvar, pis, pjs, pks, i, j, primal_topology->nl-1, n);
     if (qf_choice == QF_MODE::EC)
    { compute_Qxz_w_EC_bottom<1, ADD_MODE::ADD>(Wtendvar, qxzreconvar, qxzvertreconvar, Fvar, pis, pjs, pks, i, j, 0, n);
      compute_Qxz_w_EC_top<1, ADD_MODE::ADD>(Wtendvar, qxzreconvar, qxzvertreconvar, Fvar, pis, pjs, pks, i, j, primal_topology->nl-1, n);}
     if (qf_choice == QF_MODE::NOEC)
    { compute_Qxz_w_nonEC_bottom<1, ADD_MODE::ADD>(Wtendvar, qxzreconvar, Fvar, pis, pjs, pks, i, j, 0, n);
      compute_Qxz_w_nonEC_top<1, ADD_MODE::ADD>(Wtendvar, qxzreconvar, Fvar, pis, pjs, pks, i, j, primal_topology->nl-1, n);}
    compute_Qxz_w_EC_bottom<1, ADD_MODE::ADD>(Wtendvar, coriolisxzreconvar, coriolisxzvertreconvar, Fvar, pis, pjs, pks, i, j, 0, n);
    compute_Qxz_w_EC_top<1, ADD_MODE::ADD>(Wtendvar, coriolisxzreconvar, coriolisxzvertreconvar, Fvar, pis, pjs, pks, i, j, primal_topology->nl-1, n);
    });

    parallel_for("Compute Vtend",SimpleBounds<4>( primal_topology->ni-2, primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
      compute_wD1_fct<ndensity> (Vtendvar, densreconvar, Phivar, Bvar, pis, pjs, pks, i, j, k+1, n);
      if (qf_choice == QF_MODE::EC)
      { compute_Qxz_u_EC<1, ADD_MODE::ADD>(Vtendvar,qxzreconvar, qxzvertreconvar,FWvar, pis, pjs, pks, i, j, k+1, n);}
      if (qf_choice == QF_MODE::NOEC)
      { compute_Qxz_u_nonEC<1, ADD_MODE::ADD>(Vtendvar, qxzvertreconvar, FWvar, pis, pjs, pks, i, j, k+1, n);}
      compute_Qxz_u_EC<1, ADD_MODE::ADD>(Vtendvar,coriolisxzreconvar, coriolisxzvertreconvar,FWvar, pis, pjs, pks, i, j, k+1, n);
    });
    parallel_for("Compute Vtend Bnd",SimpleBounds<3>( primal_topology->n_cells_y, primal_topology->n_cells_x, primal_topology->nens) , YAKL_LAMBDA(int j, int i, int n) { 
    compute_wD1_fct<ndensity> (Vtendvar, densreconvar, Phivar, Bvar, pis, pjs, pks, i, j, 0, n);
    compute_wD1_fct<ndensity> (Vtendvar, densreconvar, Phivar, Bvar, pis, pjs, pks, i, j, primal_topology->ni-1, n);
    if (qf_choice == QF_MODE::EC)
    { compute_Qxz_u_EC_bottom<1, ADD_MODE::ADD>(Vtendvar,qxzreconvar, qxzvertreconvar,FWvar, pis, pjs, pks, i, j, 0, n);
      compute_Qxz_u_EC_top<1, ADD_MODE::ADD>(Vtendvar,qxzreconvar, qxzvertreconvar,FWvar, pis, pjs, pks, i, j, primal_topology->ni-1, n);}
    if (qf_choice == QF_MODE::NOEC)
    { compute_Qxz_u_nonEC_bottom<1, ADD_MODE::ADD>(Vtendvar, qxzvertreconvar, FWvar, pis, pjs, pks, i, j, 0, n);
      compute_Qxz_u_nonEC_top<1, ADD_MODE::ADD>(Vtendvar, qxzvertreconvar, FWvar, pis, pjs, pks, i, j, primal_topology->ni-1, n);}
    compute_Qxz_u_EC_bottom<1, ADD_MODE::ADD>(Vtendvar,coriolisxzreconvar, coriolisxzvertreconvar,FWvar, pis, pjs, pks, i, j, 0, n);
    compute_Qxz_u_EC_top<1, ADD_MODE::ADD>(Vtendvar,coriolisxzreconvar, coriolisxzvertreconvar,FWvar, pis, pjs, pks, i, j, primal_topology->ni-1, n);
    });

      parallel_for("Compute Dens Tend", SimpleBounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
      compute_wDbar2_fct<ndensity> (denstendvar, densreconvar, Phivar, Fvar, dis, djs, dks, i, j, k, n);
      compute_wDvbar_fct<ndensity, ADD_MODE::ADD> (denstendvar, densvertreconvar, Phivertvar, FWvar, dis, djs, dks, i, j, k, n);
      });

      }



      
      void YAKL_INLINE compute_rhs(real dt, FieldSet<nconstant> &const_vars, FieldSet<nprognostic> &x, FieldSet<nauxiliary> &auxiliary_vars, FieldSet<nprognostic> &xtend)
      {
           
          //Compute U, W, q0, dens0
          compute_functional_derivatives_and_diagnostic_quantities_I(
          auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[UWVAR].data, auxiliary_vars.fields_arr[QXZ0VAR].data, auxiliary_vars.fields_arr[FXZ0VAR].data,
          auxiliary_vars.fields_arr[DENS0VAR].data,
          x.fields_arr[VVAR].data, x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data, const_vars.fields_arr[CORIOLISXZVAR].data);
      
          auxiliary_vars.fields_arr[QXZ0VAR].set_bnd(0.0);
          auxiliary_vars.fields_arr[FXZ0VAR].set_bnd(0.0);
          auxiliary_vars.fields_arr[UWVAR].set_bnd(0.0);
          this->aux_exchange->exchanges_arr[UVAR].exchange_field(auxiliary_vars.fields_arr[UVAR]);
          this->aux_exchange->exchanges_arr[UWVAR].exchange_field(auxiliary_vars.fields_arr[UWVAR]);
          this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(auxiliary_vars.fields_arr[DENS0VAR]);
          this->aux_exchange->exchanges_arr[QXZ0VAR].exchange_field(auxiliary_vars.fields_arr[QXZ0VAR]);
          this->aux_exchange->exchanges_arr[FXZ0VAR].exchange_field(auxiliary_vars.fields_arr[FXZ0VAR]);

          //Compute K, F, FW, he, hew
          compute_functional_derivatives_and_diagnostic_quantities_II(
          auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[FWVAR].data, auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[HEVAR].data, auxiliary_vars.fields_arr[HEWVAR].data,
          x.fields_arr[VVAR].data, auxiliary_vars.fields_arr[UVAR].data, x.fields_arr[WVAR].data, auxiliary_vars.fields_arr[UWVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data);

          auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);
          this->aux_exchange->exchanges_arr[FVAR].exchange_field(auxiliary_vars.fields_arr[FVAR]);
          this->aux_exchange->exchanges_arr[FWVAR].exchange_field(auxiliary_vars.fields_arr[FWVAR]);
          this->aux_exchange->exchanges_arr[KVAR].exchange_field(auxiliary_vars.fields_arr[KVAR]);
          this->aux_exchange->exchanges_arr[HEVAR].exchange_field(auxiliary_vars.fields_arr[HEVAR]);
          this->aux_exchange->exchanges_arr[HEWVAR].exchange_field(auxiliary_vars.fields_arr[HEWVAR]);
              
          //Compute FT, B
          compute_functional_derivatives_and_diagnostic_quantities_III(
          auxiliary_vars.fields_arr[BVAR].data,  auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[FTWVAR].data,
          auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[UVAR].data,
          auxiliary_vars.fields_arr[FWVAR].data, auxiliary_vars.fields_arr[UWVAR].data,
          auxiliary_vars.fields_arr[KVAR].data, x.fields_arr[DENSVAR].data, const_vars.fields_arr[HSVAR].data);
          //auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data, const_vars.fields_arr[HSVAR].data);

          this->aux_exchange->exchanges_arr[BVAR].exchange_field(auxiliary_vars.fields_arr[BVAR]);
          this->aux_exchange->exchanges_arr[FTVAR].exchange_field(auxiliary_vars.fields_arr[FTVAR]);
          this->aux_exchange->exchanges_arr[FTWVAR].exchange_field(auxiliary_vars.fields_arr[FTWVAR]);
              
          // Compute densrecon, densvertrecon, qrecon and frecon
          compute_edge_reconstructions(
          auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data, 
          auxiliary_vars.fields_arr[QXZEDGERECONVAR].data, auxiliary_vars.fields_arr[QXZVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISXZEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[QXZ0VAR].data, auxiliary_vars.fields_arr[FXZ0VAR].data);

          this->aux_exchange->exchanges_arr[DENSEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSEDGERECONVAR]);
          this->aux_exchange->exchanges_arr[DENSVERTEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR]);
          this->aux_exchange->exchanges_arr[QXZEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[QXZEDGERECONVAR]);
          this->aux_exchange->exchanges_arr[QXZEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[QXZEDGERECONVAR]);
          this->aux_exchange->exchanges_arr[CORIOLISXZVERTEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR]);
          this->aux_exchange->exchanges_arr[CORIOLISXZVERTEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR]);

          compute_recons(
          auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSVERTRECONVAR].data, 
          auxiliary_vars.fields_arr[QXZRECONVAR].data, auxiliary_vars.fields_arr[QXZVERTRECONVAR].data, 
          auxiliary_vars.fields_arr[CORIOLISXZRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISXZVERTRECONVAR].data, 
          auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[QXZEDGERECONVAR].data, auxiliary_vars.fields_arr[QXZVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISXZEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[HEVAR].data, auxiliary_vars.fields_arr[HEWVAR].data,
          auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[UWVAR].data,
          auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[FTWVAR].data);
          
          this->aux_exchange->exchanges_arr[DENSRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSRECONVAR]);
          this->aux_exchange->exchanges_arr[DENSVERTRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSVERTRECONVAR]);
          this->aux_exchange->exchanges_arr[QXZRECONVAR].exchange_field(auxiliary_vars.fields_arr[QXZRECONVAR]);
          this->aux_exchange->exchanges_arr[QXZVERTRECONVAR].exchange_field(auxiliary_vars.fields_arr[QXZVERTRECONVAR]);
          this->aux_exchange->exchanges_arr[CORIOLISXZRECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISXZRECONVAR]);
          this->aux_exchange->exchanges_arr[CORIOLISXZVERTRECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISXZVERTRECONVAR]);


    //Compute fct
    int dis = dual_topology->is;
    int djs = dual_topology->js;
    int dks = dual_topology->ks;

      parallel_for("ComputeEdgeFlux", SimpleBounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
    compute_edgefluxes<ndensity> (auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, dis, djs, dks, i, j, k, n);
    });
      parallel_for("ComputeVertEdgeFlux", SimpleBounds<4>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
      compute_vertedgefluxes<ndensity> (auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data, auxiliary_vars.fields_arr[DENSVERTRECONVAR].data, auxiliary_vars.fields_arr[FWVAR].data, dis, djs, dks, i, j, k, n);
    });
    this->aux_exchange->exchanges_arr[EDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[EDGEFLUXVAR]);
    this->aux_exchange->exchanges_arr[VERTEDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[VERTEDGEFLUXVAR]);

    parallel_for("ComputeMf", SimpleBounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
    compute_Mfext<ndensity> (auxiliary_vars.fields_arr[MFVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data, dt, dis, djs, dks, i, j, k, n);
    });

    this->aux_exchange->exchanges_arr[MFVAR].exchange_field(auxiliary_vars.fields_arr[MFVAR]);
    
      parallel_for("ComputePhiVert", SimpleBounds<4>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
      compute_Phivert<ndensity> (auxiliary_vars.fields_arr[PHIVERTVAR].data, auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k+1, n);
    });
      parallel_for("ComputePhi", SimpleBounds<4>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA(int k, int j, int i, int n) { 
      compute_Phi<ndensity> (auxiliary_vars.fields_arr[PHIVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k, n);
    });

    //Don't do FCT for non-FCT vars
        for (int l=0; l<ndensity; l++)
        {
          //if (not varset.dens_pos(l))
          if (!varset.dens_pos[l])
          {
          auxiliary_vars.fields_arr[PHIVAR].set(l, 1.0);
          auxiliary_vars.fields_arr[PHIVERTVAR].set(l, 1.0);
          }
        }
        
    this->aux_exchange->exchanges_arr[PHIVAR].exchange_field(auxiliary_vars.fields_arr[PHIVAR]);
    this->aux_exchange->exchanges_arr[PHIVERTVAR].exchange_field(auxiliary_vars.fields_arr[PHIVERTVAR]);

          // Compute tendencies
          compute_tendencies(
          xtend.fields_arr[DENSVAR].data, xtend.fields_arr[VVAR].data, xtend.fields_arr[WVAR].data,
          auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSVERTRECONVAR].data,
          auxiliary_vars.fields_arr[QXZRECONVAR].data, auxiliary_vars.fields_arr[QXZVERTRECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISXZRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISXZVERTRECONVAR].data,
          auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[FWVAR].data, 
          auxiliary_vars.fields_arr[PHIVAR].data, auxiliary_vars.fields_arr[PHIVERTVAR].data);

  }

};



// *******   Statistics   ***********//

class ModelStats: public Stats {
public:
  
   real3d TEarr, KEarr, PEarr, IEarr, PVarr, PENSarr, trimmed_density;
   
  void initialize(ModelParameters &params, Parallel &par, const Topology &primal_topo, const Topology &dual_topo, Geometry &primal_geom, Geometry &dual_geom)
{
  Stats::initialize(params, par, primal_topo, dual_topo, primal_geom, dual_geom);
  this->stats_arr[DENSSTAT].initialize("mass", ndensity, this->statsize, this->nens, this->masterproc);
  this->stats_arr[DENSMAXSTAT].initialize("densmax", ndensity, this->statsize, this->nens, this->masterproc);
  this->stats_arr[DENSMINSTAT].initialize("densmin", ndensity, this->statsize, this->nens, this->masterproc);
  this->stats_arr[ESTAT].initialize("energy", 4, this->statsize, this->nens, this->masterproc);
  this->stats_arr[PVSTAT].initialize("pv", 1, this->statsize, this->nens, this->masterproc);
  this->stats_arr[PESTAT].initialize("pens", 1, this->statsize, this->nens, this->masterproc);
  
  this->TEarr = real3d("TE", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
  this->KEarr = real3d("KE", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
  this->IEarr = real3d("IE", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
  this->PEarr = real3d("PE", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
  this->PVarr = real3d("PV", this->primal_topology->nl, this->primal_topology->n_cells_y, this->primal_topology->n_cells_x);
  this->PENSarr = real3d("PENS", this->primal_topology->nl, this->primal_topology->n_cells_y, this->primal_topology->n_cells_x);
  this->trimmed_density = real3d("trimmed_density", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
  
}

   void compute( FieldSet<nprognostic> &progvars,  FieldSet<nconstant> &constvars, int tind)
   {

    for (int n=0;n<nens;n++)
    {

      SArray<real,1,ndensity> masslocal, massglobal;
      SArray<real,1,ndensity> densmaxlocal, densmaxglobal;
      SArray<real,1,ndensity> densminlocal, densminglobal;
      SArray<real,1,1> pvlocal, pvglobal;
      SArray<real,1,4> elocal, eglobal;
      SArray<real,1,1> pelocal, peglobal;


      pvlocal(0) = 0.;
      pvglobal(0) = 0.;
      pelocal(0) = 0.;
      peglobal(0) = 0.;
      elocal(0) = 0.;
      elocal(1) = 0.;
      elocal(2) = 0.;
      elocal(3) = 0.;
      eglobal(0) = 0.;
      eglobal(1) = 0.;
      eglobal(2) = 0.;
      eglobal(3) = 0.;
      for (int l=0;l<ndensity;l++) {masslocal(l) = 0.; massglobal(l) = 0.;}
      for (int l=0;l<ndensity;l++) {densmaxlocal(l) = 0.; densmaxglobal(l) = 0.;}
      for (int l=0;l<ndensity;l++) {densminlocal(l) = 0.; densminglobal(l) = 0.;}

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

parallel_for("Compute energetics stats", SimpleBounds<3>( dual_topology->nl-2, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
   real KE, PE, IE;
KE = Hk.compute_KE(progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data, progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k+1, n);
PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data, constvars.fields_arr[HSVAR].data, dis, djs, dks, i, j, k+1, n);
IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k+1, n);
TEarr(k+1, j, i) = KE + PE + IE;
KEarr(k+1, j, i) = KE;
PEarr(k+1, j, i) = PE;
IEarr(k+1, j, i) = IE;
});
parallel_for("Compute energetics stats bnd", SimpleBounds<2>( dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int j, int i) { 
real KE, PE, IE;
KE = Hk.compute_KE_bottom(progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data, progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, 0, n);
PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data, constvars.fields_arr[HSVAR].data, dis, djs, dks, i, j, 0, n);
IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, 0, n);
TEarr(0, j, i) = KE + PE + IE;
KEarr(0, j, i) = KE;
PEarr(0, j, i) = PE;
IEarr(0, j, i) = IE;
KE = Hk.compute_KE_top(progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data, progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, dual_topology->nl-1, n);
PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data, constvars.fields_arr[HSVAR].data, dis, djs, dks, i, j, dual_topology->nl-1, n);
IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, dual_topology->nl-1, n);
TEarr(dual_topology->nl-1, j, i) = KE + PE + IE;
KEarr(dual_topology->nl-1, j, i) = KE;
PEarr(dual_topology->nl-1, j, i) = PE;
IEarr(dual_topology->nl-1, j, i) = IE;
});

elocal(0) = yakl::intrinsics::sum(TEarr);
elocal(1) = yakl::intrinsics::sum(KEarr);
elocal(2) = yakl::intrinsics::sum(PEarr);
elocal(3) = yakl::intrinsics::sum(IEarr);

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

parallel_for("Compute PV/PE stats", SimpleBounds<3>( primal_topology->nl-2, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
 pvpe vals_pvpe;
 vals_pvpe = PVPE.compute_PVPE(progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data, progvars.fields_arr[DENSVAR].data, constvars.fields_arr[CORIOLISXZVAR].data, pis, pjs, pks, i, j, k+1, n);
 PVarr(k+1, j, i) = vals_pvpe.pv;
 PENSarr(k+1, j, i) = vals_pvpe.pe;
  });
  parallel_for("Compute PV/PE stats bnd", SimpleBounds<2>( primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int j, int i) { 
   pvpe vals_pvpe;
   vals_pvpe = PVPE.compute_PVPE_bottom(progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data, progvars.fields_arr[DENSVAR].data, constvars.fields_arr[CORIOLISXZVAR].data, pis, pjs, pks, i, j, 0, n);
   PVarr(0, j, i) = vals_pvpe.pv;
   PENSarr(0, j, i) = vals_pvpe.pe;
   vals_pvpe = PVPE.compute_PVPE_top(progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data, progvars.fields_arr[DENSVAR].data, constvars.fields_arr[CORIOLISXZVAR].data, pis, pjs, pks, i, j, primal_topology->nl-1, n);
   PVarr(primal_topology->nl-1, j, i) = vals_pvpe.pv;
   PENSarr(primal_topology->nl-1, j, i) = vals_pvpe.pe;
    });

    pvlocal(0) = yakl::intrinsics::sum(PVarr);
    pelocal(0) = yakl::intrinsics::sum(PENSarr);

    for (int l=0;l<ndensity;l++)
    {
      parallel_for("Compute trimmed density", SimpleBounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
        trimmed_density(k,j,i) = progvars.fields_arr[DENSVAR].data(l,k+dks,j+djs,i+dis,n);
      });


    masslocal(l) = yakl::intrinsics::sum(trimmed_density);
    densmaxlocal(l) = yakl::intrinsics::maxval(trimmed_density);
    densminlocal(l) = yakl::intrinsics::minval(trimmed_density);
  }

    //MPI sum/min/max
    this->ierr = MPI_Ireduce( &masslocal, &massglobal, ndensity, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[DENSSTAT]);
    this->ierr = MPI_Ireduce( &densmaxlocal, &densmaxglobal, ndensity, REAL_MPI, MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[DENSMAXSTAT]);
    this->ierr = MPI_Ireduce( &densminlocal, &densminglobal, ndensity, REAL_MPI, MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[DENSMINSTAT]);
    this->ierr = MPI_Ireduce( &pvlocal, &pvglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PVSTAT]);
    this->ierr = MPI_Ireduce( &pelocal, &peglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PESTAT]);
    this->ierr = MPI_Ireduce( &elocal, &eglobal, 4, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[ESTAT]);

    this->ierr = MPI_Waitall(nstats, this->Req, this->Status);


  if (masterproc)
  {
    for (int l=0;l<ndensity;l++)
    {
  this->stats_arr[DENSSTAT].data(l,tind,n) = massglobal(l);
  this->stats_arr[DENSMAXSTAT].data(l,tind,n) = densmaxglobal(l);
  this->stats_arr[DENSMINSTAT].data(l,tind,n) = densminglobal(l);
}

  this->stats_arr[ESTAT].data(0,tind,n) = eglobal(0);
  this->stats_arr[ESTAT].data(1,tind,n) = eglobal(1);
  this->stats_arr[ESTAT].data(2,tind,n) = eglobal(2);
  this->stats_arr[ESTAT].data(3,tind,n) = eglobal(3);
  this->stats_arr[PVSTAT].data(0,tind,n) = pvglobal(0);
  this->stats_arr[PESTAT].data(0,tind,n) = peglobal(0);

}

}
}
};


// *******   FieldSet Initialization   ***********//
void initialize_variables(const Topology &ptopo, const Topology &dtopo,
SArray<int,2, nprognostic, 3> &prog_ndofs_arr, SArray<int,2, nconstant, 3> &const_ndofs_arr, SArray<int,2, nauxiliary, 3> &aux_ndofs_arr, SArray<int,2, ndiagnostic, 3> &diag_ndofs_arr,
std::array<std::string, nprognostic> &prog_names_arr, std::array<std::string, nconstant> &const_names_arr, std::array<std::string, nauxiliary> &aux_names_arr, std::array<std::string, ndiagnostic> &diag_names_arr,
std::array<const Topology *, nprognostic> &prog_topo_arr, std::array<const Topology *, nconstant> &const_topo_arr, std::array<const Topology *, nauxiliary> &aux_topo_arr, std::array<const Topology *, ndiagnostic> &diag_topo_arr)
{

  //primal grid represents straight quantities, dual grid twisted quantities
// ndims is the BASEDIM size!

  // v, w, dens
  prog_topo_arr[VVAR] = &ptopo;
  prog_topo_arr[WVAR] = &ptopo;
  prog_topo_arr[DENSVAR] = &dtopo;
  prog_names_arr[VVAR] = "v";
  prog_names_arr[WVAR] = "w";
  prog_names_arr[DENSVAR] = "dens";
  set_dofs_arr(prog_ndofs_arr, VVAR, 1, 0, 1); //v = straight (1,0)-form
  set_dofs_arr(prog_ndofs_arr, WVAR, 0, 1, 1); //w = straight (0,1)-form
  set_dofs_arr(prog_ndofs_arr, DENSVAR, ndims, 1, ndensity); //dens = twisted (n,1)-form

  // hs  
  const_topo_arr[HSVAR] = &dtopo;
  const_names_arr[HSVAR] = "hs";
  set_dofs_arr(const_ndofs_arr, HSVAR, ndims, 1, 1); //hs = twisted (n,1)-form
  const_topo_arr[CORIOLISXZVAR] = &ptopo;
  const_names_arr[CORIOLISXZVAR] = "coriolisxz";
  set_dofs_arr(const_ndofs_arr, CORIOLISXZVAR, 1, 1, 1); //f = straight (1,1)-form
  
  //functional derivatives = F, B, K, he, U
  aux_topo_arr[BVAR] = &ptopo;
  aux_topo_arr[FVAR] = &dtopo;
  aux_topo_arr[UVAR] = &dtopo;
  aux_topo_arr[HEVAR] = &dtopo;
  aux_topo_arr[FWVAR] = &dtopo;
  aux_topo_arr[UWVAR] = &dtopo;
  aux_topo_arr[HEWVAR] = &dtopo;
  aux_topo_arr[KVAR] = &dtopo;
  aux_names_arr[KVAR] = "K";
  aux_names_arr[BVAR] = "B";
  aux_names_arr[FVAR] = "F";
  aux_names_arr[UVAR] = "U";
  aux_names_arr[HEVAR] = "he";
  aux_names_arr[FWVAR] = "Fw";
  aux_names_arr[UWVAR] = "Uw";
  aux_names_arr[HEWVAR] = "hew";
  set_dofs_arr(aux_ndofs_arr, BVAR, 0, 0, ndensity); //B = straight (0,0)-form
  set_dofs_arr(aux_ndofs_arr, KVAR, ndims, 1, 1);   //K = twisted (n,1)-form
  set_dofs_arr(aux_ndofs_arr, FVAR, ndims-1, 1, 1);  //F = twisted (n-1,1)-form
  set_dofs_arr(aux_ndofs_arr, UVAR, ndims-1, 1, 1); //U = twisted (n-1,1)-form
  set_dofs_arr(aux_ndofs_arr, HEVAR, ndims-1, 1, 1); //he lives on horiz dual edges, associated with F
  set_dofs_arr(aux_ndofs_arr, FWVAR, ndims, 0, 1);  //Fw = twisted (n,0)-form
  set_dofs_arr(aux_ndofs_arr, UWVAR, ndims, 0, 1); //Uw = twisted (n,0)-form
  set_dofs_arr(aux_ndofs_arr, HEWVAR, ndims, 0, 1); //hew lives on vert dual edges, associated with Fw
  
  //dens primal grid reconstruction stuff- dens0, edgerecon, recon
  aux_topo_arr[DENSRECONVAR] = &dtopo;
  aux_topo_arr[DENSEDGERECONVAR] = &dtopo;
  aux_topo_arr[DENSVERTRECONVAR] = &dtopo;
  aux_topo_arr[DENSVERTEDGERECONVAR] = &dtopo;
  aux_topo_arr[DENS0VAR] = &ptopo;
  aux_names_arr[DENS0VAR] = "dens0";
  aux_names_arr[DENSVERTRECONVAR] = "densvertrecon";
  aux_names_arr[DENSVERTEDGERECONVAR] = "densvertedgerecon";
  aux_names_arr[DENSRECONVAR] = "densrecon";
  aux_names_arr[DENSEDGERECONVAR] = "densedgerecon";
  set_dofs_arr(aux_ndofs_arr, DENSRECONVAR, ndims-1, 1, ndensity);  //densrecon lives on horiz dual edges, associated with F
  set_dofs_arr(aux_ndofs_arr, DENSEDGERECONVAR, ndims, 1, 2*ndims*ndensity); //densedgerecon lives on dual cells, associated with F
  set_dofs_arr(aux_ndofs_arr, DENSVERTRECONVAR, ndims, 0, ndensity);  //densvertrecon lives on vert dual edges, associated with Fw
  set_dofs_arr(aux_ndofs_arr, DENSVERTEDGERECONVAR, ndims, 1, 2*ndensity); //densedgerecon lives on dual cells, associated with Fw
  set_dofs_arr(aux_ndofs_arr, DENS0VAR, 0, 0, ndensity); //dens0 = straight (0,0)-form

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
  set_dofs_arr(aux_ndofs_arr, PHIVAR, ndims-1, 1, ndensity); 
  set_dofs_arr(aux_ndofs_arr, PHIVERTVAR, ndims, 0, ndensity); 
  set_dofs_arr(aux_ndofs_arr, MFVAR, ndims, 1, ndensity);  
  set_dofs_arr(aux_ndofs_arr, EDGEFLUXVAR, ndims-1, 1, ndensity); 
  set_dofs_arr(aux_ndofs_arr, VERTEDGEFLUXVAR, ndims, 0, ndensity); 
  
  
  //Q stuff
  aux_topo_arr[QXZ0VAR] = &dtopo;
  aux_names_arr[QXZ0VAR] = "QXZ0";
  set_dofs_arr(aux_ndofs_arr, QXZ0VAR, 0, 0, 1); //Q0 = twisted (0,0)-form
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
  set_dofs_arr(aux_ndofs_arr, QXZRECONVAR, 0, 1, 1);  //qxzrecon lives on vert primal edges, associated with w
  set_dofs_arr(aux_ndofs_arr, QXZEDGERECONVAR, ndims, 1, 2*1); //qxzedgerecon lives on primal cells, associated with Fw/w
  set_dofs_arr(aux_ndofs_arr, QXZVERTRECONVAR, 1, 0, 1);  //qxzsvertrecon lives on horiz primal edges, associated with v
  set_dofs_arr(aux_ndofs_arr, QXZVERTEDGERECONVAR, ndims, 1, 2*1); //qxzvertedgerecon lives on primal cells, associated with F/v
  set_dofs_arr(aux_ndofs_arr, QXZFLUXVAR, 0, 1, 1); //qxzflux lives on vert primal edges, associated with w
  set_dofs_arr(aux_ndofs_arr, QXZVERTFLUXVAR, 1, 0, 1); //qxzvertflux lives on horiz primal edges, associated with v

  aux_topo_arr[FTVAR] = &ptopo;
  aux_topo_arr[FTWVAR] = &ptopo;
  aux_names_arr[FTVAR] = "FT";
  aux_names_arr[FTWVAR] = "FTW";
  set_dofs_arr(aux_ndofs_arr, FTVAR, 1, 0, 1); //FT = straight (1,0)-form ie Fw at v pts
  set_dofs_arr(aux_ndofs_arr, FTWVAR, 0, 1, 1); //FTW = straight (0,1)-form ie F at w pts
  
  aux_topo_arr[FXZ0VAR] = &dtopo;
  aux_topo_arr[CORIOLISXZRECONVAR] = &ptopo;
  aux_topo_arr[CORIOLISXZEDGERECONVAR] = &ptopo; 
  aux_topo_arr[CORIOLISXZVERTRECONVAR] = &ptopo;
  aux_topo_arr[CORIOLISXZVERTEDGERECONVAR] = &ptopo; 
  aux_names_arr[FXZ0VAR] = "fxz0";
  aux_names_arr[CORIOLISXZRECONVAR] = "coriolisxzrecon";
  aux_names_arr[CORIOLISXZEDGERECONVAR] = "coriolisxzedgerecon";
  aux_names_arr[CORIOLISXZVERTRECONVAR] = "coriolisxzvertrecon";
  aux_names_arr[CORIOLISXZVERTEDGERECONVAR] = "coriolisxzvertedgerecon";
  set_dofs_arr(aux_ndofs_arr, FXZ0VAR, 0, 0, 1);  //fxz0 is a twisted 0,0 form
  set_dofs_arr(aux_ndofs_arr, CORIOLISXZRECONVAR, 0, 1, 1);  //coriolisxzrecon lives on vert primal edges, associated with w
  set_dofs_arr(aux_ndofs_arr, CORIOLISXZEDGERECONVAR, ndims, 1, 2*1); //coriolisxzedgerecon lives on primal cells, associated with Fw/w
  set_dofs_arr(aux_ndofs_arr, CORIOLISXZVERTRECONVAR, 1, 0, 1);  //coriolisxzsvertrecon lives on horiz primal edges, associated with v
  set_dofs_arr(aux_ndofs_arr, CORIOLISXZVERTEDGERECONVAR, ndims, 1, 2*1); //coriolisxzvertedgerecon lives on primal cells, associated with F/v  
  
  // q, concentration 0-forms for dens
  diag_topo_arr[DENSLDIAGVAR] = &ptopo;
  diag_names_arr[DENSLDIAGVAR] = "densl";
  set_dofs_arr(diag_ndofs_arr, DENSLDIAGVAR, 0, 0, ndensity); //densldiag = straight (0,0)-form
  diag_topo_arr[QXZDIAGVAR] = &dtopo;
  diag_names_arr[QXZDIAGVAR] = "QXZl";
  set_dofs_arr(diag_ndofs_arr, QXZDIAGVAR, 0, 0, 1); //Qldiag = twisted (0,0)-form

  // #if defined _AN || defined _MAN
  // aux_topo_arr[PVAR] = &ptopo; //p = straight 0-form
  // aux_names_arr[PVAR] = "p";
  // set_dofs_arr(aux_ndofs_arr, PVAR, 0, 1, 1);  //p = straight 0-form
  // #endif
  
}


//***************** Set Initial Conditions ***************************//

struct dbv_constants {
real const g = 9.80616_fp;
real const Lx = 5000000._fp;
real const Ly = 5000000._fp;
real const coriolis = 0.00006147_fp;
real const H0 = 750.0_fp;
real const ox = 0.1_fp;
real const oy = 0.1_fp;
real const sigmax = 3._fp/40._fp*Lx;
real const sigmay = 3._fp/40._fp*Ly;
real const dh = 75.0_fp;
real const xc1 = (0.5_fp - ox) * Lx;
real const yc1 = (0.5_fp - oy) * Ly;
real const xc2 = (0.5_fp + ox) * Lx;
real const yc2 = (0.5_fp + oy) * Ly;
real const xc = 0.5_fp * Lx;
real const yc = 0.5_fp * Ly;
real const c = 0.05_fp;
real const a = 1.0_fp/3.0_fp;
real const D = 0.5_fp * Lx;
};
dbv_constants dbl_vortex_constants;


real YAKL_INLINE double_vortex_coriolis(real x, real y)
{
  return dbl_vortex_constants.coriolis;
}

real YAKL_INLINE double_vortex_h(real x, real y)
{
  real xprime1 = dbl_vortex_constants.Lx / (pi * dbl_vortex_constants.sigmax) * sin(pi / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc1));
  real yprime1 = dbl_vortex_constants.Ly / (pi * dbl_vortex_constants.sigmay) * sin(pi / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc1));
  real xprime2 = dbl_vortex_constants.Lx / (pi * dbl_vortex_constants.sigmax) * sin(pi / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc2));
  real yprime2 = dbl_vortex_constants.Ly / (pi * dbl_vortex_constants.sigmay) * sin(pi / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc2));
  real xprimeprime1 = dbl_vortex_constants.Lx / (2.0_fp * pi * dbl_vortex_constants.sigmax) * sin(2.0_fp * pi / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc1));
  real yprimeprime1 = dbl_vortex_constants.Ly / (2.0_fp * pi * dbl_vortex_constants.sigmay) * sin(2.0_fp * pi / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc1));
  real xprimeprime2 = dbl_vortex_constants.Lx / (2.0_fp * pi * dbl_vortex_constants.sigmax) * sin(2.0_fp * pi / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc2));
  real yprimeprime2 = dbl_vortex_constants.Ly / (2.0_fp * pi * dbl_vortex_constants.sigmay) * sin(2.0_fp * pi / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc2));

  return dbl_vortex_constants.H0 - dbl_vortex_constants.dh * (exp(-0.5_fp * (xprime1 * xprime1 + yprime1 * yprime1)) + exp(-0.5_fp * (xprime2 * xprime2 + yprime2 * yprime2)) - 4._fp * pi * dbl_vortex_constants.sigmax * dbl_vortex_constants.sigmay / dbl_vortex_constants.Lx / dbl_vortex_constants.Ly);
}

vecext<2> YAKL_INLINE double_vortex_v(real x, real y) {
vecext<2> vvec;

real xprime1 = dbl_vortex_constants.Lx / (pi * dbl_vortex_constants.sigmax) * sin(pi / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc1));
real yprime1 = dbl_vortex_constants.Ly / (pi * dbl_vortex_constants.sigmay) * sin(pi / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc1));
real xprime2 = dbl_vortex_constants.Lx / (pi * dbl_vortex_constants.sigmax) * sin(pi / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc2));
real yprime2 = dbl_vortex_constants.Ly / (pi * dbl_vortex_constants.sigmay) * sin(pi / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc2));
real xprimeprime1 = dbl_vortex_constants.Lx / (2.0_fp * pi * dbl_vortex_constants.sigmax) * sin(2.0_fp * pi / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc1));
real yprimeprime1 = dbl_vortex_constants.Ly / (2.0_fp * pi * dbl_vortex_constants.sigmay) * sin(2.0_fp * pi / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc1));
real xprimeprime2 = dbl_vortex_constants.Lx / (2.0_fp * pi * dbl_vortex_constants.sigmax) * sin(2.0_fp * pi / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc2));
real yprimeprime2 = dbl_vortex_constants.Ly / (2.0_fp * pi * dbl_vortex_constants.sigmay) * sin(2.0_fp * pi / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc2));

vvec.u = - dbl_vortex_constants.g * dbl_vortex_constants.dh / dbl_vortex_constants.coriolis / dbl_vortex_constants.sigmay * (yprimeprime1 * exp(-0.5_fp*(xprime1 * xprime1 + yprime1 * yprime1)) + yprimeprime2 * exp(-0.5_fp*(xprime2 * xprime2 + yprime2 * yprime2)));
vvec.w = dbl_vortex_constants.g * dbl_vortex_constants.dh / dbl_vortex_constants.coriolis / dbl_vortex_constants.sigmax * (xprimeprime1 * exp(-0.5_fp*(xprime1 * xprime1 + yprime1 * yprime1)) + xprimeprime2 * exp(-0.5_fp*(xprime2 * xprime2 + yprime2 * yprime2)));
return vvec;
}

real YAKL_INLINE double_vortex_S(real x, real y)
{
  //real sval = g * (1. + c * sin(2. * M_PI / Lx * (x - xc)) * sin(2. * M_PI / Ly * (y - yc)) * exp(-((x-xc)*(x-xc) + (y-yc)*(y-yc))/(a*a*D*D)));
  real sval = dbl_vortex_constants.g * (1._fp + dbl_vortex_constants.c * exp(-((x-dbl_vortex_constants.xc)*(x-dbl_vortex_constants.xc) + (y-dbl_vortex_constants.yc)*(y-dbl_vortex_constants.yc))/(dbl_vortex_constants.a*dbl_vortex_constants.a*dbl_vortex_constants.D*dbl_vortex_constants.D)));
  //real sval = g * (1. + c * sin(2. * M_PI / Lx * (x- xc)));
  //real sval = g;
  //real sval = g * (1. + c * ((x > 0.35 * Lx && x < 0.65 * Lx && y > 0.35 * Ly && y < 0.65 * Ly ) ? 1. : 0.));
  return sval * double_vortex_h(x,y);
}


// CAN WE GENERALIZE THESE? HOW? NEED TO CALL OUT TO THE CORRECT HEIGHT FUNCTION...
// ALSO NEED TO SET VARIOUS LX/LY SIZES
// MAYBE THIS LATTER IS DOABLE VIA PARAMS?
// SIMILAR WITH GAUSSIANS
// MAYBE WHAT WE DO IS HAVE A TRACER STRUCT, AND THEN SET THE VALUES FOR THE TRACER STRUCTURE ACCORDINGLY?
// AND THE TRACER STRUCT CAN HAVE A LINK TO DENSITY/HEIGHT FUNCTION!
// IN FACT, I THINK ALL TEST CASES STRUCTURES SHOULD INHERIT FROM A MASTER STRUCTURE THAT HAS LX/LY/ETC ALREADY SET UP?
// SOMETHING LIKE THIS...

real YAKL_INLINE tracer_square_cent(real x, real y)         {return (x > 0.35_fp*dbl_vortex_constants.Lx && x < 0.65_fp*dbl_vortex_constants.Lx && y > 0.35_fp*dbl_vortex_constants.Ly && y < 0.65_fp*dbl_vortex_constants.Ly                        ) ? 0.005_fp : 0.;}
real YAKL_INLINE tracer_square_ur(real x, real y)         {return (x > 0.6_fp*dbl_vortex_constants.Lx && x < 0.9_fp*dbl_vortex_constants.Lx && y > 0.6_fp*dbl_vortex_constants.Ly && y < 0.9_fp*dbl_vortex_constants.Ly                        ) ? 0.005_fp : 0.;}
real YAKL_INLINE tracer_square_ll(real x, real y)         {return (x > 0.1_fp*dbl_vortex_constants.Lx && x < 0.4_fp*dbl_vortex_constants.Lx && y > 0.1_fp*dbl_vortex_constants.Ly && y < 0.4_fp*dbl_vortex_constants.Ly                        ) ? 0.005_fp : 0.;}
real YAKL_INLINE tracer_square_urpll(real x, real y)         {return tracer_square_ur(x,y) + tracer_square_ll(x,y);}

real YAKL_INLINE double_vortex_tracer_square_cent (real x, real y) {return double_vortex_h(x,y) * tracer_square_cent(x, y);}
real YAKL_INLINE double_vortex_tracer_square_urpll (real x, real y) {return double_vortex_h(x,y) * tracer_square_urpll(x, y);}

real YAKL_INLINE double_vortex_tracer_gaussian(real x, real y)     {return double_vortex_h(x,y) * 0.005_fp * exp(-((x-dbl_vortex_constants.xc)*(x-dbl_vortex_constants.xc) + (y-dbl_vortex_constants.yc)*(y-dbl_vortex_constants.yc))/(dbl_vortex_constants.a*dbl_vortex_constants.a*dbl_vortex_constants.D*dbl_vortex_constants.D));}
//{ return 0.005 *double_vortex_h(x,y) * exp(-100. * pow((x-xc)/Lx,2.)) * exp(-100. * pow((y-yc)/Ly,2.)); }



struct smallbubble_constants {
real const g = 9.80616_fp;
real const Lx = 1000._fp;
real const Ly = 1500._fp;
real const yc = 0.5_fp * Ly;
real const xc = 0.5_fp * Lx;
real const theta0 = 300.0_fp;
real const zc = 350._fp;
real const dss = 0.5_fp;
real const rc = 250._fp;
real const rh0 = 0.8_fp;
};
smallbubble_constants rb_constants;

struct largebubble_constants  {
  real const g = 9.80616_fp;
  real const Lx = 20000._fp;
  real const Lz = 20000._fp;
  real const xc = 0.5_fp * Lx;
  real const zc = 0.5_fp * Lz;
  real const theta0 = 300.0_fp;
  real const bzc = 2000._fp;
  real const xrad = 2000._fp;
  real const zrad = 2000._fp;
  real const amp_theta = 2.0_fp;
  real const amp_vapor = 0.8_fp;
  real const Cpv = 1859._fp;
  real const Cpd = 1003._fp;
};
largebubble_constants lrb_constants;

// Universal

real YAKL_INLINE isentropic_T(real x, real z, real theta0, real g)
{
  return theta0 - z * g / thermo.cst.Cpd;
}

real YAKL_INLINE isentropic_p(real x, real z, real theta0, real g)
{
  return thermo.cst.pr * pow(isentropic_T(x, z, theta0, g) / theta0, 1./thermo.cst.kappa_d);
}

real YAKL_INLINE isentropic_rho(real x, real z, real theta0, real g) {
  real p = isentropic_p(x, z, theta0, g);
  real T = isentropic_T(x, z, theta0, g);
  real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
  return 1._fp/alpha;
}


real YAKL_INLINE linear_ellipsoid(real x, real z, real x0, real z0, real xrad, real zrad, real amp)
{
  real xn = (x-x0)/xrad;
  real zn = (z-z0)/zrad;
  real dist = sqrt( xn*xn + zn*zn );
  return amp * std::max( 1._fp - dist , 0._fp );  
}

real YAKL_INLINE flat_geop(real x, real z, real g)
{
  return g * z;
}

// Returns saturation vapor pressure
 real YAKL_INLINE saturation_vapor_pressure(real temp) {
   real tc = temp - 273.15_fp;
   return 610.94_fp * exp( 17.625_fp*tc / (243.04_fp + tc) );
 }






 // (Moist) Rising Bubble (small-scale)

 real YAKL_INLINE rb_entropicvar(real x, real z) {
   real p = isentropic_p(x, z, rb_constants.theta0, rb_constants.g);
   real T = isentropic_T(x, z, rb_constants.theta0, rb_constants.g);
   real r = sqrt((x-rb_constants.xc)*(x-rb_constants.xc) + (z-rb_constants.zc)*(z-rb_constants.zc));
   real dtheta = (r<rb_constants.rc) ? rb_constants.dss * 0.5_fp * (1._fp + cos(pi * r/rb_constants.rc)) : 0._fp;
   real dT = dtheta * pow(p/thermo.cst.pr, thermo.cst.kappa_d);
   return thermo.compute_entropic_var(p, T+dT, 0, 0, 0, 0);
 }

 real YAKL_INLINE rb_rho(real x, real z)
 {
   return isentropic_rho(x, z, rb_constants.theta0, rb_constants.g);
 }

 real YAKL_INLINE rb_entropicdensity(real x, real z) {
   return rb_entropicvar(x,z) * rb_rho(x,z);
 }

 real YAKL_INLINE rb_rho_acousticbalance(real x, real z) {
   real rho_b = rb_rho(x,z);
   real theta = rb_entropicvar(x,z);
   return rho_b * rb_constants.theta0 / theta;

 }

 real YAKL_INLINE rb_entropicdensity_acousticbalance(real x, real z) {
   return rb_entropicvar(x,z) * rb_rho_acousticbalance(x,z);
 }

 real YAKL_INLINE rb_geop(real x, real z)
 {
   return flat_geop(x,z,rb_constants.g);
 }

 // The moist bubble is just the dry bubble with an added moisture perturbation
 // There is no attempt at balancing anything, just rho_d and theta_h are in balance!
 // This is not quite hydrostatic balance even...


 // We assume a formula here for SVP that might not be consistent with the thermodynamics
 real YAKL_INLINE mrb_rho_v(real x, real z) {
   real r = sqrt((x-rb_constants.xc)*(x-rb_constants.xc) + (z-rb_constants.zc)*(z-rb_constants.zc));
   real rh = (r<rb_constants.rc) ? rb_constants.rh0 * (1._fp + cos(pi * r/rb_constants.rc)) : 0._fp;
   real Th = isentropic_T(x, z, rb_constants.theta0, rb_constants.g);
   real svp = saturation_vapor_pressure(Th);
   real pv = svp * rh;
   return pv / (thermo.cst.Rv * Th);
 }

 real YAKL_INLINE mrb_rho_d(real x, real z) {
   real p = isentropic_p(x, z, rb_constants.theta0, rb_constants.g);
   real T = isentropic_T(x, z, rb_constants.theta0, rb_constants.g);
   real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
   return 1._fp/alpha;
 }

 real YAKL_INLINE mrb_rho(real x, real z) {
   real rhod = mrb_rho_d(x,z);
   real rhov = mrb_rho_v(x,z);
   return rhod + rhov;
 }

 real YAKL_INLINE mrb_entropicdensity(real x, real z) {
   real p = isentropic_p(x, z, rb_constants.theta0, rb_constants.g);
   real T = isentropic_T(x, z, rb_constants.theta0, rb_constants.g);
   real r = sqrt((x-rb_constants.xc)*(x-rb_constants.xc) + (z-rb_constants.zc)*(z-rb_constants.zc));
   real dtheta = (r<rb_constants.rc) ? rb_constants.dss * 0.5_fp * (1. + cos(pi * r/rb_constants.rc)) : 0._fp;
   real dT = dtheta * pow(p/thermo.cst.pr, thermo.cst.kappa_d);
   real theta = thermo.compute_entropic_var(p, T+dT, 1, 0, 0, 0);
   return theta * mrb_rho(x,z);
 }








 // (Moist) Rising Bubble (large-scale)

 real YAKL_INLINE lrb_entropicvar(real x, real z) {
   
   real p = isentropic_p(x, z, lrb_constants.theta0, lrb_constants.g);
   real T0 = isentropic_T(x, z, lrb_constants.theta0, lrb_constants.g);
   real dtheta = linear_ellipsoid(x, z, lrb_constants.xc, lrb_constants.bzc, lrb_constants.xrad, lrb_constants.zrad, lrb_constants.amp_theta);
   real dT = dtheta * pow(p/thermo.cst.pr, thermo.cst.kappa_d);
   return thermo.compute_entropic_var(p, T0+dT, 0, 0, 0, 0);

 }

 real YAKL_INLINE lrb_rho(real x, real z)
 {
   return isentropic_rho(x, z, lrb_constants.theta0, lrb_constants.g);
 }

 real YAKL_INLINE lrb_entropicdensity(real x, real z) {
   return lrb_entropicvar(x,z) * lrb_rho(x,z);
 }

 // real YAKL_INLINE lrb_rho_acousticbalance(real x, real z) {
 //   real rho_b = lrb_rho(x,z);
 //   real theta = lrb_entropicvar(x,z);
 //   return rho_b * lrb_constants.theta0 / theta;
 // 
 // }
 // 
 // real YAKL_INLINE lrb_entropicdensity_acousticbalance(real x, real z) {
 //   return lrb_entropicvar(x,z) * lrb_rho_acousticbalance(x,z);
 // }

 real YAKL_INLINE lrb_geop(real x, real z)
 {
   return flat_geop(x,z,lrb_constants.g);
 }

 real YAKL_INLINE mlrb_rho_d(real x, real z) {
   real p = isentropic_p(x, z, lrb_constants.theta0, lrb_constants.g);
   real T = isentropic_T(x, z, lrb_constants.theta0, lrb_constants.g);
   real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
   return 1._fp/alpha;
 }

 real YAKL_INLINE mlrb_rho_v(real x, real z) {
  
   real pert = linear_ellipsoid(x, z, lrb_constants.xc, lrb_constants.bzc, lrb_constants.xrad, lrb_constants.zrad, lrb_constants.amp_vapor);
   real Th = isentropic_T(x, z, lrb_constants.theta0, rb_constants.g);
   real svp = saturation_vapor_pressure(Th);
   real pv = svp * pert;
   return pv / (thermo.cst.Rv * Th);
 }


 real YAKL_INLINE mlrb_rho(real x, real z) {
   real rhod = mlrb_rho_d(x,z);
   real rhov = mlrb_rho_v(x,z);
   return rhod + rhov;
 }

 real YAKL_INLINE mlrb_entropicvar(real x, real z) {
   
   real p = isentropic_p(x, z, lrb_constants.theta0, lrb_constants.g);
   real T0 = isentropic_T(x, z, lrb_constants.theta0, lrb_constants.g);
   real dtheta = linear_ellipsoid(x, z, lrb_constants.xc, lrb_constants.bzc, lrb_constants.xrad, lrb_constants.zrad, lrb_constants.amp_theta);
   real dT = dtheta * pow(p/thermo.cst.pr, thermo.cst.kappa_d);
   return thermo.compute_entropic_var(p, T0+dT, 1, 0, 0, 0);
 }
 
 real YAKL_INLINE mlrb_entropicdensity(real x, real z) {
   return mlrb_entropicvar(x,z) * mlrb_rho(x,z);
 }




 void readModelParamsFile(std::string inFile, ModelParameters &params, Parallel &par, int nz)
 {
   readParamsFile( inFile, params, par, nz);

   //Read config file
   YAML::Node config = YAML::LoadFile(inFile);
      
   params.acoustic_balance = config["balance_initial_density"].as<bool>(false);

   // Read the data initialization options
   params.initdataStr = config["initData"].as<std::string>();

   serial_print("IC: " + params.initdataStr, par.masterproc);
   serial_print("acoustically balanced: " + std::to_string(params.acoustic_balance), par.masterproc);
   
   for (int i=0;i<ntracers_dycore;i++)
   {
     params.tracerdataStr[i] = config["initTracer" + std::to_string(i)].as<std::string>();
     params.dycore_tracerpos[i] = config["initTracerPos" + std::to_string(i)].as<bool>();
     serial_print("Dycore Tracer" + std::to_string(i) + " IC: " + params.tracerdataStr[i], par.masterproc);
   }
   

 }


void set_domain_sizes_ic (ModelParameters &params, std::string initData)
{

  params.ylen = 1.0;
  params.yc = 0.5;
  
  if (initData == "doublevortex")
{
  
  params.xlen = dbl_vortex_constants.Lx;
  params.zlen = dbl_vortex_constants.Ly;
  params.xc = dbl_vortex_constants.xc;
  params.zc = dbl_vortex_constants.yc;
}
if (initData == "risingbubble" or initData == "moistrisingbubble")
{
params.xlen = rb_constants.Lx;
params.zlen = rb_constants.Ly;
params.xc = rb_constants.Lx * 0.5_fp;
params.zc = rb_constants.Ly  * 0.5_fp;
}

if (initData == "largerisingbubble" or initData == "moistlargerisingbubble")
{
params.xlen = lrb_constants.Lx;
params.zlen = lrb_constants.Lz;
params.xc = lrb_constants.Lx * 0.5_fp;
params.zc = lrb_constants.Lz  * 0.5_fp;
}

}
  
void set_initial_conditions (ModelParameters &params, FieldSet<nprognostic> &progvars, FieldSet<nconstant> &constvars, 
Geometry &primal_geom, Geometry &dual_geom)
{

  if (params.initdataStr == "doublevortex")
  {
      dual_geom.set_11form_values(double_vortex_h, progvars.fields_arr[DENSVAR], 0);
#ifdef _TSWE
      dual_geom.set_11form_values(double_vortex_S, progvars.fields_arr[DENSVAR], 1);
#endif
primal_geom.set_10form_values(double_vortex_v, progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
primal_geom.set_01form_values(double_vortex_v,progvars.fields_arr[WVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
  
  //RESTORE ONCE PRIMAL GRID CFV RECON IS FIXED IN PARALLEL
  //primal_geom.set_11form_values(double_vortex_coriolis, constvars.fields_arr[CORIOLISVAR], 0);

// HOW DO GENERALIZE THESE?
// WANT TO SCALE TRACER FIELDS BY ACTUAL HEIGHT FIELDS...
// SHOULD BE USABLE FOR ANY IC!
for (int i=0; i<ntracers_dycore; i++)
{
if (params.tracerdataStr[i] == "gaussian") {dual_geom.set_11form_values(double_vortex_tracer_gaussian, progvars.fields_arr[DENSVAR], i+ndensity_dycore);}
if (params.tracerdataStr[i] == "square") {dual_geom.set_11form_values(double_vortex_tracer_square_cent, progvars.fields_arr[DENSVAR], i+ndensity_dycore);}
if (params.tracerdataStr[i] == "doublesquare") {dual_geom.set_11form_values(double_vortex_tracer_square_urpll, progvars.fields_arr[DENSVAR], i+ndensity_dycore);}
}
Hs.set_parameters(dbl_vortex_constants.g);
  }
  
  if (params.initdataStr == "risingbubble")
  {
    if (params.acoustic_balance)
    {
    dual_geom.set_11form_values(rb_rho_acousticbalance, progvars.fields_arr[DENSVAR], 0);
    dual_geom.set_11form_values(rb_entropicdensity_acousticbalance, progvars.fields_arr[DENSVAR], 1);
    }
    else
    {
      dual_geom.set_11form_values(rb_rho, progvars.fields_arr[DENSVAR], 0);
      dual_geom.set_11form_values(rb_entropicdensity, progvars.fields_arr[DENSVAR], 1);      
    }
    dual_geom.set_11form_values(rb_geop, constvars.fields_arr[HSVAR], 0);
  }

  
  // ADD ACOUSTIC BALANCING HERE
  if (params.initdataStr == "moistrisingbubble")
  {
    // ADD ANELASTIC/MOIST ANELASTIC?
    #if defined _MCErhod || defined _MCErhodp
    dual_geom.set_11form_values(mrb_rho_d, progvars.fields_arr[DENSVAR], 0);
    #elif defined _MCErho || defined _MCErhop
    dual_geom.set_11form_values(mrb_rho, progvars.fields_arr[DENSVAR], 0);    
    #endif
    dual_geom.set_11form_values(mrb_entropicdensity, progvars.fields_arr[DENSVAR], 1);
    dual_geom.set_11form_values(mrb_rho_v, progvars.fields_arr[DENSVAR], varset.dm_id_vap + ndensity_nophysics);
    dual_geom.set_11form_values(rb_geop, constvars.fields_arr[HSVAR], 0);
  }
  
  if (params.initdataStr == "largerisingbubble")
  {
    // if (params.acoustic_balance)
    // {
    // dual_geom.set_11form_values(lrb_rho_acousticbalance, progvars.fields_arr[DENSVAR], 0);
    // dual_geom.set_11form_values(lrb_entropicdensity_acousticbalance, progvars.fields_arr[DENSVAR], 1);
    // }
    // else
    // {
      dual_geom.set_11form_values(lrb_rho, progvars.fields_arr[DENSVAR], 0);
      dual_geom.set_11form_values(lrb_entropicdensity, progvars.fields_arr[DENSVAR], 1);      
    //}
    dual_geom.set_11form_values(lrb_geop, constvars.fields_arr[HSVAR], 0);
  }

  
  // ADD ACOUSTIC BALANCING HERE
  if (params.initdataStr == "moistlargerisingbubble")
  {
    // ADD ANELASTIC/MOIST ANELASTIC?
    #if defined _MCErhod || defined _MCErhodp
    dual_geom.set_11form_values(mlrb_rho_d, progvars.fields_arr[DENSVAR], 0);
    #elif defined _MCErho || defined _MCErhop
    dual_geom.set_11form_values(mlrb_rho, progvars.fields_arr[DENSVAR], 0);    
    #endif
    dual_geom.set_11form_values(mlrb_entropicdensity, progvars.fields_arr[DENSVAR], 1);
    dual_geom.set_11form_values(mlrb_rho_v, progvars.fields_arr[DENSVAR], varset.dm_id_vap + ndensity_nophysics);
    dual_geom.set_11form_values(lrb_geop, constvars.fields_arr[HSVAR], 0);
  }
  
}