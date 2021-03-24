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
#include "string.h"

//DO/SHOULD WE HAVE OUT OF SLICE VELOCITY AND ASSOCIATED CORIOLIS?
//NOT FOR INITIAL VERSION....

// rho + entropic variable density
// AN/MAN don't actually evolve rho though, it stays equal to initial constant value
#ifdef _AN
uint constexpr ndensity = 2 + NTRACERS;
#endif
#ifdef _MAN
uint constexpr ndensity = 2 + NTRACERS;
#endif
#ifdef _FC
uint constexpr ndensity = 2 + NTRACERS;
#endif
#ifdef _MFC
uint constexpr ndensity = 2 + NTRACERS;
#endif

//WHAT IS THE BEST WAY TO HANDLE ADDITIONAL ACTIVE TRACER IE TRACE SPECIES DENSITIES?
//IN RUN SCRIPT I THINK, BUT THIS IS REALLY A FUNCTION OF EQN SET...
//do we advect around non-active tracers?
//for the forseeable future I think we can just worry about vapor, liquid, (ice)

#if NTRACERS > 0
uint constexpr ntracers = NTRACERS;
#endif

#if NTRACERS_FCT > 0
uint constexpr ndensityfct = NTRACERS_FCT;
#endif

// Number of variables
// v, w, dens, densfct
#if NTRACERS_FCT > 0
uint constexpr nprognostic = 4;
#define DENSFCTVAR 3
#else
uint constexpr nprognostic = 3;
#endif

#define VVAR 0
#define WVAR 1
#define DENSVAR 2


// MORE NEEDED HERE? EVENTUALLY OUT OF SLICE CORIOLIS, ANYTHING ELSE? GRAVITY/GEOPOTENIAL?
// topography
uint constexpr nconstant = 1;
#define TOPOGRAPHYVAR 0

// I THINK THIS STUFF IS BASICALLY THE SAME...

//functional derivatives + aux quantities= F, Fw, B, BFCT, K
//primal grid reconstruction stuff- he, U, Uw, dens0, densfct0, edgerecon, recon, edgereconfct, reconfct + vert stuff
//dual grid reconstruction stuff- q0, f0, FT, qedgerecon, qrecon + vert stuff

//EVENTUALLY NEED CORIOLIS STUFF IF WE DO OUT OF SLICE VELOCITY...

//fct stuff- Phi, Mf, edgeflux

#if NTRACERS_FCT > 0
uint constexpr nauxiliary = 27;
#else
uint constexpr nauxiliary = 18;
#endif


#define FVAR 0
#define FWVAR 1
#define BVAR 2
#define KVAR 3

#define HEVAR 4
#define UVAR 5
#define HEWVAR 6
#define UWVAR 7

#define DENS0VAR 8
#define DENSRECONVAR 9
#define DENSEDGERECONVAR 10
#define DENSVERTRECONVAR 11
#define DENSVERTEDGERECONVAR 12

#define FTVAR 13
#define Q0VAR 14
#define F0VAR 15
#define QRECONVAR 16
#define QEDGERECONVAR 17

#if NTRACERS_FCT > 0
#define BFCTVAR 18
#define DENSFCT0VAR 19
#define DENSFCTRECONVAR 20
#define DENSFCTEDGERECONVAR 21
#define DENSFCTVERTRECONVAR 22
#define DENSFCTVERTEDGERECONVAR 23
#define PHIVAR 24
#define EDGEFLUXVAR 25
#define MFVAR 26
#endif

// q, associated concentration 0-forms for den

#if NTRACERS_FCT > 0
uint constexpr ndiagnostic = 3;
#define DENSFCTLDIAGVAR 2
#else
uint constexpr ndiagnostic = 2;
#endif
#define QDIAGVAR 0
#define DENSLDIAGVAR 1


//track total densities, dens min/max, densfct min/max, energy (P+K+I+total), PV
#if NTRACERS_FCT > 0
uint constexpr nstats = 8;
#else
uint constexpr nstats = 5;
#endif

#define DENSSTAT 0
#define DENSMINSTAT 1
#define DENSMAXSTAT 2
#define ESTAT 3
#define PVSTAT 4
#if NTRACERS_FCT > 0
#define DENSFCTSTAT 5
#define DENSFCTMAXSTAT 6
#define DENSFCTMINSTAT 7
#endif





// *******   Model Specific Parameters   ***********//

void set_model_specific_params(std::string inFile, ModelParameters &params)
{


  std::string strDataInit = "";

  // Read in equals-separated key = value file line by line
  std::ifstream fInStream(inFile);
  std::string line;
  while (std::getline(fInStream, line)) {
    // Remove spaces and tabs from the line
    line.erase (std::remove(line.begin(), line.end(), ' '), line.end());
    line.erase (std::remove(line.begin(), line.end(), '\t'), line.end());

    // If the line isn't empty and doesn't begin with a comment specifier, split it based on the colon
    if (!line.empty() && line.find("//",0) != 0) {
      // Find the colon
      uint splitloc = line.find('=',0);
      // Store the key and value strings
      std::string key   = line.substr(0,splitloc);
      std::string value = line.substr(splitloc+1,line.length()-splitloc);

      // Transform the value into a string stream for convenience
      std::stringstream ssVal(value);

      // Match the key, and store the value
           if ( !strcmp( "dataInit"   , key.c_str() ) ) { ssVal >> strDataInit       ;}
      //else {
      //  std::cout << "Error: key " << key << " not understood in file " << inFile << "\n";
      //}
    }
  }

  if (!strcmp("",strDataInit.c_str())) { std::cout << "Error: key " << "dataInit" << " not set.\n"; exit(-1); }

    size_t splitloc = strDataInit.find("//",0);
    std::string sub_str;
    if (splitloc != std::string::npos){
      sub_str = strDataInit.substr(0,splitloc);
    } else {
      sub_str = strDataInit;
    }
    if      ( !strcmp(sub_str.c_str(),"rb" ) ) { params.data_init_cond = DATA_INIT::RB  ; }
    else if ( !strcmp(sub_str.c_str(),"nhgw" ) ) { params.data_init_cond = DATA_INIT::NHGW  ; }
    else if ( !strcmp(sub_str.c_str(),"hgw"   ) ) { params.data_init_cond = DATA_INIT::HGW    ; }
    else  {
      std::cout << "Error: unrecognized dataInit " << strDataInit << "\n";
      exit(-1);
    }

  params.etime = 0.0;

//FIX THIS!!!
  if (params.data_init_cond == DATA_INIT::RB)
  {
  params.xlen = Lx;
  params.xc = Lx/2.;
  params.ylen = Ly;
  params.yc = Ly/2.;
  }

}

// ******* Diagnostics *************//

// THIS SHOULD BE GENERALIZABLE...
// ONLY COMPUTE FUNCTION NEEDS TO CHANGE

template <uint nprog, uint nconst, uint ndiag> class Diagnostics {
public:

  const Topology *primal_topology;
  const Topology *dual_topology;
  Geometry<ndims,1,1,1> *primal_geometry;
  Geometry<ndims,1,1,1> *dual_geometry;

  bool is_initialized;

   Diagnostics() {
     this->is_initialized = false;
     std::cout << "CREATED DIAGNOSTICS\n";
   }

   void initialize(const Topology &ptopo, const Topology &dtopo, Geometry<ndims,1,1,1> &pgeom, Geometry<ndims,1,1,1> &dgeom)
   {
     this->primal_topology = &ptopo;
     this->dual_topology = &dtopo;
     this->primal_geometry = &pgeom;
     this->dual_geometry = &dgeom;
     this->is_initialized = true;
   }


#if NTRACERS_FCT > 0
   void YAKL_INLINE compute_diagnostic_quantities(
     realArr Q0var, realArr dens0var, realArr densfct0var,
     const realArr Vvar, const realArr densvar, const realArr densfctvar) {
#else
void YAKL_INLINE compute_diagnostic_quantities(
  realArr Q0var, realArr dens0var,
  const realArr Vvar, const realArr densvar) {
#endif

//PRIMAL VS. DUAL HERE?
// IE WHAT DO WE ACTUALLY LOOP OVER?
       int is = primal_topology->is;
       int js = primal_topology->js;
       int ks = primal_topology->ks;

       yakl::parallel_for("ComputeDiagI", topology->n_cells, YAKL_LAMBDA (int iGlob) {
         SArray<real,1> zeta;
         real hv;

         int k, j, i;
         yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

// DO THESE COMPUTATIONS NEED TO BE ADJUSTED FOR 2D SPLIT-TYPE STUFF?
 // compute zeta = D2 v
 compute_D2<1>(zeta, Vvar, is, js, ks, i, j, k);

 // compute q0 = zeta / R h
   compute_R<0> (hv, densvar, is, js, ks, i, j, k);
   Q0var(0, k+ks, j+js, i+is) = (zeta(0)) / hv;

// compute dens0var
for(int l=0; l<ndensity; l++)
{dens0var(l, k+ks, j+js, i+is) = densvar(l, k+ks, j+js, i+is) / densvar(0, k+ks, j+js, i+is);}

// compute densfct0var
#if NTRACERS_FCT > 0
for(int l=0; l<ndensityfct; l++)
{densfct0var(l, k+ks, j+js, i+is) = densfctvar(l, k+ks, j+js, i+is) / densvar(0, k+ks, j+js, i+is);}
#endif
       });

     }



   void compute_diag(const VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<ndiag> &diagnostic_vars)
   {

#if NTRACERS_FCT > 0
   compute_diagnostic_quantities(
   diagnostic_vars.fields_arr[QDIAGVAR].data, diagnostic_vars.fields_arr[DENSLDIAGVAR].data, diagnostic_vars.fields_arr[DENSFCTLDIAGVAR].data,
   x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data, x.fields_arr[DENSFCTVAR].data);
#else
compute_diagnostic_quantities(
diagnostic_vars.fields_arr[QDIAGVAR].data, diagnostic_vars.fields_arr[DENSLDIAGVAR].data,
x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data,);
#endif
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
  Geometry<ndims,1,1,1> *primal_geometry;
  Geometry<ndims,1,1,1> *dual_geometry;

  TransformMatrices<real> trans;

  SArray<real,reconstruction_order,2> primal_to_gll;
  SArray<real,reconstruction_order,reconstruction_order,reconstruction_order> primal_wenoRecon;
  SArray<real,(reconstruction_order-1)/2+2> primal_wenoIdl;
  real primal_wenoSigma;

  SArray<real,dual_reconstruction_order,2> dual_to_gll;
  SArray<real,dual_reconstruction_order,dual_reconstruction_order,dual_reconstruction_order> dual_wenoRecon;
  SArray<real,(dual_reconstruction_order-1)/2+2> dual_wenoIdl;
  real dual_wenoSigma;

  SArray<real,vert_reconstruction_order,2> primal_vert_to_gll;
  SArray<real,vert_reconstruction_order,vert_reconstruction_order,vert_reconstruction_order> primal_vert_wenoRecon;
  SArray<real,(vert_reconstruction_order-1)/2+2> primal_vert_wenoIdl;
  real primal_vert_wenoSigma;

  SArray<real,dual_vert_reconstruction_order,2> dual_vert_to_gll;
  SArray<real,dual_vert_reconstruction_order,dual_vert_reconstruction_order,dual_vert_reconstruction_order> dual_vert_wenoRecon;
  SArray<real,(dual_vert_reconstruction_order-1)/2+2> dual_vert_wenoIdl;
  real dual_vert_wenoSigma;

  bool is_initialized;

   Tendencies() {
     this->is_initialized = false;
     std::cout << "CREATED TENDENCIES\n";
   }

  void initialize(const Topology &primal_topo, const Topology &dual_topo, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom, ExchangeSet<naux> &aux_exchange, ExchangeSet<nconst> &const_exchange)
  {
    this->primal_topology = &primal_topo;
    this->dual_topology = &dual_topo;
    this->primal_geometry = &primal_geom;
    this->dual_geometry = &dual_geom;
    this->aux_exchange = &aux_exchange;
    this->const_exchange = &const_exchange;

    trans.coefs_to_gll_lower( primal_to_gll );
    trans.weno_sten_to_coefs(primal_wenoRecon);
    wenoSetIdealSigma<reconstruction_order>(primal_wenoIdl,primal_wenoSigma);

    trans.coefs_to_gll_lower( dual_to_gll );
    trans.weno_sten_to_coefs(dual_wenoRecon);
    wenoSetIdealSigma<dual_reconstruction_order>(dual_wenoIdl,dual_wenoSigma);

    trans.coefs_to_gll_lower( primal_vert_to_gll );
    trans.weno_sten_to_coefs(primal_vert_wenoRecon);
    wenoSetIdealSigma<vert_reconstruction_order>(primal_vert_wenoIdl,primal_vert_wenoSigma);

    trans.coefs_to_gll_lower( dual_vert_to_gll );
    trans.weno_sten_to_coefs(dual_vert_wenoRecon);
    wenoSetIdealSigma<dual_vert_reconstruction_order>(dual_vert_wenoIdl,dual_vert_wenoSigma);

    this->is_initialized = true;
  }


//HOW DO WE HANDLE ANELASTIC HERE
//DO WE EVOLVE RHO AND HOPE CONSTRAINT EQN GIVES RHO=CONST?
//WOULD BASICALLY REQUIRE EXACT SOLN OF CONSTRAINT- ROUND-OFF WILL BREAK THIS EVENTUALLY...
//OR NOT EVOLVE IT?
//UNDERLYING HAMILTONIAN STRUCTURE HAS RHO AS ACTIVE VARIABLE

#if NTRACERS_FCT > 0
  void YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_I(
    realArr Uvar, realArr Uwvar, realArr Q0var, realArr dens0var, realArr densfct0var,
    const realArr Vvar,  const realArr Wvar, const realArr densvar, const realArr densfctvar) {
#else
void YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_I(
  realArr Uvar, realArr Uwvar, realArr Q0var, realArr dens0var,
  const realArr Vvar, const realArr Wvar, const realArr densvar) {
#endif

//PRIMAL VS. DUAL HERE?
// IE WHAT DO WE ACTUALLY LOOP OVER?
      int is = topology->is;
      int js = topology->js;
      int ks = topology->ks;

      yakl::parallel_for("ComputeDiagI", topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        real hv;
        yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

//DO WE NEED A HORIZ/VERT HODGE STAR SPLITTING?
//MAYBE FOR HEVI, BUT IT REALLY DOESN'T MAKE A LOT OF SENSE...

        // compute dens0var = I densvar, U = H v, densfct0var = I densfctvar
//FIX THESE COMPUTATIONS...
        compute_H<1, diff_ord>(Uvar, Vvar, *this->geom, is, js, ks, i, j, k);
        compute_H<1, diff_ord>(Uwvar, Wvar, *this->geom, is, js, ks, i, j, k);


//PRIMAL VS. DUAL GEOM HERE?
        compute_I<ndensity, diff_ord>(dens0var, densvar, *this->geom, is, js, ks, i, j, k);
#if NTRACERS_FCT > 0
        compute_I<ndensityfct, diff_ord>(densfct0var, densfctvar, *this->geom, is, js, ks, i, j, k);
#endif

// DOES THIS STUFF NEED HORIZ-VERT SPLITTING?

// compute zeta = D2 v
compute_D2<1>(Q0var, Vvar, is, js, ks, i, j, k);

// compute q0 = zeta / R h
  compute_R<0> (hv, densvar, is, js, ks, i, j, k);
  Q0var(0, k+ks, j+js, i+is) = Q0var(0, k+ks, j+js, i+is) / hv;

      });

    }

    void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_II(
      realArr Fvar, realArr Fwvar, realArr Kvar, realArr HEvar, realArr HEwvar, const realArr Vvar, const realArr Wvar, const realArr Uvar, const realArr Uwvar, const realArr dens0var) {

//PRIMAL VS. DUAL HERE?
// IE WHAT DO WE ACTUALLY LOOP OVER?
        int is = topology->is;
        int js = topology->js;
        int ks = topology->ks;

        yakl::parallel_for("ComputeDiagII", topology->n_cells, YAKL_LAMBDA (int iGlob) {
          int k, j, i;
          yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

//FIX THIS FOR SPLIT
        // compute he = phi * h0
        compute_phi<0>(HEvar, dens0var, is, js, ks, i, j, k);
        compute_phi<0>(HEwvar, dens0var, is, js, ks, i, j, k);

        //compute F = he * U, Fw = hew * Uw
        Fvar(0, k+ks, j+js, i+is) = Uvar(0, k+ks, j+js, i+is) * HEvar(0, k+ks, j+js, i+is);
        Fwvar(0, k+ks, j+js, i+is) = Uwvar(0, k+ks, j+js, i+is) * HEwvar(0, k+ks, j+js, i+is);

//FIX THIS FOR SPLIT
        //compute KE = 0.5 * phiT(u,v)
        compute_phiT(Kvar, Uvar, Vvar, is, js, ks, i, j, k);
        Kvar(0, k+ks, j+js, i+is) *= 0.5;

      });

      }


//HOW SHOULD FT BE HANLDED FOR SPLIT?
#if NTRACERS_FCT > 0
  void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_III(
    realArr FTvar, realArr Bvar, realArr Bfctvar,
    const realArr Fvar, const realArr Fwvar, const realArr Uvar, const realArr Uwvar,
    const realArr Kvar, const realArr dens0var, const realArr densfct0var) {
#else
void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_III(
  realArr FTvar, realArr Bvar,
  const realArr Fvar, const realArr Fwvar, const realArr Uvar, const realArr Uwvar,
  const realArr Kvar, const realArr dens0var) {
#endif

//PRIMAL VS. DUAL HERE?
// IE WHAT DO WE ACTUALLY LOOP OVER?
      int is = topology->is;
      int js = topology->js;
      int ks = topology->ks;

      yakl::parallel_for("ComputeDiagIII", topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        SArray<real,1> hs0;
        yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

//HOW SHOULD FT BE HANLDED FOR SPLIT?
      compute_W(FTvar, Fvar, is, js, ks, i, j, k);

//FIX THESE....
//HERE WE SHOULD CALL OUT TO THERMODYNAMICS MODULE...
#ifdef _AN
//WHICH GEOM?
      // Compute B = IK + S0/2.
      compute_I<1, diff_ord>(Bvar, Kvar, *this->geom, is, js, ks, i, j, k);
      Bvar(0, k+ks, j+js, i+is) += dens0var(1, k+ks, j+js, i+is)/2.;
      //Compute T = Ihs + h0/2;
      compute_I<1, diff_ord>(hs0, HSvar, *this->geom, is, js, ks, i, j, k);
      Bvar(1, k+ks, j+js, i+is) = dens0var(0, k+ks, j+js, i+is)/2. + hs0(0);
#endif

#ifdef _FC
      // Compute B = IK + gh0 + ghs0
      compute_I<1, diff_ord>(Bvar, Kvar, *this->geom, is, js, ks, i, j, k);
      compute_I<1, diff_ord>(hs0, HSvar, *this->geom, is, js, ks, i, j, k);
      Bvar(0, k+ks, j+js, i+is) += g*dens0var(0, k+ks, j+js, i+is) + g*hs0(0);
#endif

#ifdef _MAN
      // Compute B = IK + gh0 + ghs0
      compute_I<1, diff_ord>(Bvar, Kvar, *this->geom, is, js, ks, i, j, k);
      compute_I<1, diff_ord>(hs0, HSvar, *this->geom, is, js, ks, i, j, k);
      Bvar(0, k+ks, j+js, i+is) += g*dens0var(0, k+ks, j+js, i+is) + g*hs0(0);
#endif

#ifdef _MFC
      // Compute B = IK + gh0 + ghs0
      compute_I<1, diff_ord>(Bvar, Kvar, *this->geom, is, js, ks, i, j, k);
      compute_I<1, diff_ord>(hs0, HSvar, *this->geom, is, js, ks, i, j, k);
      Bvar(0, k+ks, j+js, i+is) += g*dens0var(0, k+ks, j+js, i+is) + g*hs0(0);
#endif


    });

    }




#if NTRACERS_FCT > 0
  void YAKL_INLINE compute_edge_reconstructions(
realArr densedgereconvar, realArr densvertedgereconvar, realArr densfctedgereconvar, realArr densfctvertedgereconvar, 
realArr Qedgereconvar,  realArr Qvertedgereconvar,
const realArr dens0var, const realArr densfct0var, const realArr Q0var) {
#else
  void YAKL_INLINE compute_edge_reconstructions(
realArr densedgereconvar, realArr densvertedgereconvar, 
realArr Qedgereconvar,  realArr Qvertedgereconvar,
const realArr dens0var, const realArr Q0var) {
#endif

//PRIMAL VS. DUAL HERE?
// IE WHAT DO WE ACTUALLY LOOP OVER?
    int is = topology->is;
    int js = topology->js;
    int ks = topology->ks;

    yakl::parallel_for("ComputeEdgeRecon", topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);


//FIX THESE UP...
      compute_primal_edge_recon<ndensity, reconstruction_type, reconstruction_order>(densedgereconvar, dens0var, is, js, ks, i, j, k, primal_wenoRecon, primal_to_gll, primal_wenoIdl, primal_wenoSigma);
      compute_primal_edge_recon<ndensity, vert_reconstruction_type, vert_reconstruction_order>(densvertedgereconvar, dens0var, is, js, ks, i, j, k, primal_vert_wenoRecon, primal_vert_to_gll, primal_vert_wenoIdl, primal_vert_wenoSigma);
#if NTRACERS_FCT > 0
      compute_primal_edge_recon<ndensityfct, reconstruction_type, reconstruction_order>(densfctedgereconvar, densfct0var, is, js, ks, i, j, k, primal_wenoRecon, primal_to_gll, primal_wenoIdl, primal_wenoSigma);
      compute_primal_edge_recon<ndensityfct, vert_reconstruction_type, vert_reconstruction_order>(densfctvertedgereconvar, densfct0var, is, js, ks, i, j, k, primal_vert_wenoRecon, primal_vert_to_gll, primal_vert_wenoIdl, primal_vert_wenoSigma);
#endif
      compute_dual_edge_recon<1, dual_reconstruction_type, dual_reconstruction_order>(Qedgereconvar, Q0var, is, js, ks, i, j, k, dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
      compute_dual_edge_recon<1, dual_vert_reconstruction_type, dual_vert_reconstruction_order>(Qvertedgereconvar, Q0var, is, js, ks, i, j, k, dual_vert_wenoRecon, dual_vert_to_gll, dual_vert_wenoIdl, dual_vert_wenoSigma);
    });

  }

//HOW SHOULD FT STUFF BE HANDLED?
#if NTRACERS_FCT > 0
  void YAKL_INLINE compute_recons(
realArr densreconvar, realArr densvertreconvar, realArr densfctreconvar, realArr densfctvertreconvar, 
realArr Qreconvar, realArr Qvertreconvar,
const realArr densedgereconvar,   const realArr densvertedgereconvar, const realArr densfctedgereconvar, const realArr densfctvertedgereconvar, 
const realArr Qedgereconvar, const realArr Qvertedgereconvar, const realArr HEvar, const realArr HEwvar,
const realArr FTvar, const realArr Uvar, const realArr Uwvar) {
#else
void YAKL_INLINE compute_recons(
realArr densreconvar, realArr densvertreconvar,
realArr Qreconvar, realArr Qvertreconvar,
const realArr densedgereconvar,   const realArr densvertedgereconvar,
const realArr Qedgereconvar, const realArr Qvertedgereconvar, const realArr HEvar, const realArr HEwvar,
const realArr FTvar, const realArr Uvar, const realArr Uwvar) {
#endif


//PRIMAL VS. DUAL HERE?
// IE WHAT DO WE ACTUALLY LOOP OVER?
    int is = topology->is;
    int js = topology->js;
    int ks = topology->ks;

    yakl::parallel_for("ComputeRecon", topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

//FIX ALL THIS UP...
      compute_primal_recon<ndensity, reconstruction_type>(densreconvar, densedgereconvar, Uvar, is, js, ks, i, j, k);
      compute_primal_recon<ndensity, vert_reconstruction_type>(densvertreconvar, densvertedgereconvar, Uwvar, is, js, ks, i, j, k);
#if NTRACERS_FCT > 0
      compute_primal_recon<ndensityfct, reconstruction_type>(densfctreconvar, densfctedgereconvar, Uvar, is, js, ks, i, j, k);
      compute_primal_recon<ndensityfct, vert_reconstruction_type>(densfctvertreconvar, densfctvertedgereconvar, Uwvar, is, js, ks, i, j, k);
#endif
//HOW SHOULD FT STUFF BE HANDLED?
      compute_dual_recon<1, dual_reconstruction_type>(Qreconvar, Qedgereconvar, FTvar, is, js, ks, i, j, k);
      compute_dual_recon<1, dual_vert_reconstruction_type>(Qvertreconvar, Qvertedgereconvar, FTvar, is, js, ks, i, j, k);


    //scale primal recons

    for (int l=0;l<ndensity;l++) {
    for (int d=0;d<1;d++) {
    densreconvar(l+d*ndensity,k+ks,j+js,i+is) = densreconvar(l+d*ndensity,k+ks,j+js,i+is) / HEvar(d,k+ks,j+js,i+is);
     }
    densvertreconvar(l+ndensity,k+ks,j+js,i+is) = densvertreconvar(l+ndensity,k+ks,j+js,i+is) / HEwvar(0,k+ks,j+js,i+is);
  }
#if NTRACERS_FCT > 0
  for (int l=0;l<ndensityfct;l++) {
  for (int d=0;d<1;d++) {
  densfctreconvar(l+d*ndensityfct,k+ks,j+js,i+is) = densfctreconvar(l+d*ndensityfct,k+ks,j+js,i+is) / HEvar(d,k+ks,j+js,i+is);
  }
  densfctvertreconvar(l+ndensityfct,k+ks,j+js,i+is) = densfctvertreconvar(l+ndensityfct,k+ks,j+js,i+is) / HEwvar(0,k+ks,j+js,i+is);
}
#endif

    });

}


#if NTRACERS_FCT > 0
  void YAKL_INLINE compute_tendencies(
  realArr denstendvar, realArr densfcttendvar, realArr Vtendvar, realArr Wtendvar,
  const realArr densreconvar, const realArr densvertreconvar, const realArr densfctreconvar, const realArr densfctvertreconvar, const realArr Qreconvar,
  const realArr Bvar, const realArr Bfctvar, const realArr Fvar, const realArr Fwvar, const realArr Phivar) {
#else
  void YAKL_INLINE compute_tendencies(
  realArr denstendvar, realArr Vtendvar, realArr Wtendvar,
  const realArr densreconvar, const realArr densvertreconvar, const realArr Qreconvar,
  const realArr Bvar, const realArr Bfctvar, const realArr Fvar, const realArr Fwvar, const realArr Phivar) {
#endif

//PRIMAL VS. DUAL HERE?
// IE WHAT DO WE ACTUALLY LOOP OVER?
    int is = topology->is;
    int js = topology->js;
    int ks = topology->ks;

      yakl::parallel_for("ComputeTendencies", topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

//WHAT NEEDS TO HAPPEN FOR HORIZ-VERT SPLITTING HERE?
    compute_wDbar2<ndensity> (denstendvar, densreconvar, Fvar, is, js, ks, i, j, k);
    compute_wD1<ndensity> (Vtendvar, densreconvar, Bvar, is, js, ks, i, j, k);
    compute_wDbar2<ndensity, ADD_MODE::ADD> (denstendvar, densreconvar, Fwvar, is, js, ks, i, j, k);
    compute_wD1<ndensity, ADD_MODE::ADD> (Wtendvar, densreconvar, Bvar, is, js, ks, i, j, k);

#if NTRACERS_FCT > 0
//WHAT NEEDS TO HAPPEN FOR HORIZ-VERT SPLITTING HERE?
    compute_wDbar2_fct<ndensityfct> (densfcttendvar, densfctreconvar, Phivar, Fvar, is, js, ks, i, j, k);
    compute_wD1_fct<ndensityfct, ADD_MODE::ADD> (Vtendvar, densfctreconvar, Phivar, Bfctvar, is, js, ks, i, j, k);
    compute_wDbar2_fct<ndensityfct, ADD_MODE::ADD> (densfcttendvar, densfctreconvar, Phivar, Fwvar, is, js, ks, i, j, k);
    compute_wD1_fct<ndensityfct, ADD_MODE::ADD> (Wtendvar, densfctreconvar, Phivar, Bfctvar, is, js, ks, i, j, k);
#endif

//HOW DOES THIS GET HANDLED?
    if (qf_choice == QF_MODE::EC)
    { compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, Qreconvar, Fvar, is, js, ks, i, j, k);}
    if (qf_choice == QF_MODE::NOEC)
    { compute_Q_nonEC<1, ADD_MODE::ADD>(Vtendvar, Qreconvar, Fvar, is, js, ks, i, j, k);}


  });

  }




  void YAKL_INLINE compute_rhs(real dt, VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<naux> &auxiliary_vars, VariableSet<nprog> &xtend)
  {

      //Compute U, Uw, q0, he, hew, dens0, densfct0
#if NTRACERS_FCT > 0
      compute_functional_derivatives_and_diagnostic_quantities_I(
      auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[UWVAR].data, auxiliary_vars.fields_arr[Q0VAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data,
      x.fields_arr[VVAR].data, x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data, x.fields_arr[DENSFCTVAR].data);

      this->aux_exchange->exchanges_arr[DENSFCT0VAR].exchange_field(auxiliary_vars.fields_arr[DENSFCT0VAR]);
#else
      compute_functional_derivatives_and_diagnostic_quantities_I(
      auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[UWVAR].data, auxiliary_vars.fields_arr[Q0VAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data,
      x.fields_arr[VVAR].data, x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data);
#endif

      this->aux_exchange->exchanges_arr[UVAR].exchange_field(auxiliary_vars.fields_arr[UVAR]);
      this->aux_exchange->exchanges_arr[UWVAR].exchange_field(auxiliary_vars.fields_arr[UWVAR]);
      this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(auxiliary_vars.fields_arr[DENS0VAR]);
      this->aux_exchange->exchanges_arr[Q0VAR].exchange_field(auxiliary_vars.fields_arr[Q0VAR]);


      //Compute K, F, Fw, he, hew
      compute_functional_derivatives_and_diagnostic_quantities_II(
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[FWVAR].data,auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[HEVAR].data, auxiliary_vars.fields_arr[HEWVAR].data,
      x.fields_arr[VVAR].data, x.fields_arr[WVAR].data, auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[UWVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data);

      this->aux_exchange->exchanges_arr[FVAR].exchange_field(auxiliary_vars.fields_arr[FVAR]);
      this->aux_exchange->exchanges_arr[FWVAR].exchange_field(auxiliary_vars.fields_arr[FWVAR]);
      this->aux_exchange->exchanges_arr[KVAR].exchange_field(auxiliary_vars.fields_arr[KVAR]);
      this->aux_exchange->exchanges_arr[HEVAR].exchange_field(auxiliary_vars.fields_arr[HEVAR]);
      this->aux_exchange->exchanges_arr[HEWVAR].exchange_field(auxiliary_vars.fields_arr[HEWVAR]);

//HOW THE HELL SHOULD FT STUFF BE HANDLED?
      //Compute FT, B, Bfct
#if NTRACERS_FCT > 0
      compute_functional_derivatives_and_diagnostic_quantities_III(
      auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[BFCTVAR].data,
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[FWVAR].data, auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[UWVAR].data,
      auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data);

      this->aux_exchange->exchanges_arr[BFCTVAR].exchange_field(auxiliary_vars.fields_arr[BFCTVAR]);
#else
      compute_functional_derivatives_and_diagnostic_quantities_III(
      auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[BVAR].data,
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[FWVAR].data, auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[UWVAR].data,
      auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data);
#endif

      this->aux_exchange->exchanges_arr[FTVAR].exchange_field(auxiliary_vars.fields_arr[FTVAR]);
      this->aux_exchange->exchanges_arr[BVAR].exchange_field(auxiliary_vars.fields_arr[BVAR]);

//FIX HERE DOWN
      // Compute densrecon, densfctrecon, qrecon and frecon
#if NTRACERS_FCT > 0
      compute_edge_reconstructions(
      auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR].data, auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data, auxiliary_vars.fields_arr[Q0VAR].data, auxiliary_vars.fields_arr[F0VAR].data);

      this->aux_exchange->exchanges_arr[DENSFCTEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR]);
#else
compute_edge_reconstructions(
auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data,
auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[Q0VAR].data, auxiliary_vars.fields_arr[F0VAR].data);
#endif

      this->aux_exchange->exchanges_arr[DENSEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[QEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[QEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[CORIOLISEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR]);

#if NTRACERS_FCT > 0
      compute_recons(
      auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSFCTRECONVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
      auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR].data,
      auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data, auxiliary_vars.fields_arr[HEVAR].data,
      auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[UVAR].data);
      this->aux_exchange->exchanges_arr[DENSFCTRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSFCTRECONVAR]);
#else
    compute_recons(
    auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
    auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
    auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data, auxiliary_vars.fields_arr[HEVAR].data,
    auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[UVAR].data);
#endif
      this->aux_exchange->exchanges_arr[DENSRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSRECONVAR]);
      this->aux_exchange->exchanges_arr[QRECONVAR].exchange_field(auxiliary_vars.fields_arr[QRECONVAR]);
      this->aux_exchange->exchanges_arr[CORIOLISRECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISRECONVAR]);


//Compute fct
#if NTRACERS_FCT > 0

//PRIMAL VS. DUAL HERE?
// IE WHAT DO WE ACTUALLY LOOP OVER?

int is = topology->is;
int js = topology->js;
int ks = topology->ks;

yakl::parallel_for("ComputeEdgeFlux", topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
compute_edgefluxes<ndensityfct> (auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[DENSFCTRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, is, js, ks, i, j, k);
});
this->aux_exchange->exchanges_arr[EDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[EDGEFLUXVAR]);


yakl::parallel_for("ComputeMf", topology->n_cells, YAKL_LAMBDA (int iGlob) {
int k, j, i;
yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
compute_Mf<ndensityfct> (auxiliary_vars.fields_arr[MFVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, dt, is, js, ks, i, j, k);
});

this->aux_exchange->exchanges_arr[MFVAR].exchange_field(auxiliary_vars.fields_arr[MFVAR]);

yakl::parallel_for("ComputePhi", topology->n_cells, YAKL_LAMBDA (int iGlob) {
int k, j, i;
yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
compute_Phi<ndensityfct> (auxiliary_vars.fields_arr[PHIVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[DENSFCTVAR].data, is, js, ks, i, j, k);
});


this->aux_exchange->exchanges_arr[PHIVAR].exchange_field(auxiliary_vars.fields_arr[PHIVAR]);
#endif

      // Compute tendencies
#if NTRACERS_FCT > 0
      compute_tendencies(
      xtend.fields_arr[DENSVAR].data, xtend.fields_arr[DENSFCTVAR].data, xtend.fields_arr[VVAR].data,
      auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSFCTRECONVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
      auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[BFCTVAR].data, auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[PHIVAR].data);
#else
  compute_tendencies(
  xtend.fields_arr[DENSVAR].data, xtend.fields_arr[VVAR].data,
  auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
  auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[FVAR].data);
#endif
};

};
// *******   Statistics Calculations   ***********//







// THIS STUFF SHOULD BE CLEANED UP AND GENERALIZED LIKE VARIABLE SETS IF POSSIBLE...
// ONLY COMPUTE FUNCTION NEEDS TO CHANGE!
class Stat
{
public:
  realArr data;
  std::string name;
  int ndofs;


void initialize(std::string statName, int ndof, ModelParameters &params, Parallel &par)
{
  name = statName;
  ndofs = ndof;

  if (par.masterproc)
  {
    data = realArr(name.c_str(), ndofs, params.Nsteps/params.Nstat + 1);
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
  Geometry<ndims,1,1,1> *primal_geometry;
  Geometry<ndims,1,1,1> *dual_geometry;

  void initialize(ModelParameters &params, Parallel &par, const Topology &primal_topo, const Topology &dual_topo, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom)
  {
      this->primal_topology = &primal_topo;
      this->dual_topology = &dual_topo;
      this->primal_geometry = &primal_geom;
      this->dual_geometry = &dual_geom;

    statsize = params.Nsteps/params.Nstat + 1;
    stats_arr[DENSSTAT].initialize("mass", ndensity, params, par);
    stats_arr[DENSMAXSTAT].initialize("densmax", ndensity, params, par);
    stats_arr[DENSMINSTAT].initialize("densmin", ndensity, params, par);
#if NTRACERS_FCT > 0
    stats_arr[DENSFCTSTAT].initialize("massfct", ndensityfct, params, par);
    stats_arr[DENSFCTMAXSTAT].initialize("densfctmax", ndensityfct, params, par);
    stats_arr[DENSFCTMINSTAT].initialize("densfctmin", ndensityfct, params, par);
#endif
    stats_arr[ESTAT].initialize("energy", 4, params, par);
    stats_arr[PVSTAT].initialize("pv", 1, params, par);
    masterproc = par.masterproc;
  }



  void compute( VariableSet<nprog> &progvars,  VariableSet<nconst> &constvars, int i)
  {


      SArray<real,ndensity> masslocal, massglobal;
      SArray<real,ndensity> densmaxlocal, densmaxglobal;
      SArray<real,ndensity> densminlocal, densminglobal;
#if NTRACERS_FCT > 0
      SArray<real,ndensityfct> massfctlocal, massfctglobal;
      SArray<real,ndensityfct> densfctmaxlocal, densfctmaxglobal;
      SArray<real,ndensityfct> densfctminlocal, densfctminglobal;
#endif
      SArray<real,4> elocal, eglobal;
      SArray<real,1> pvlocal, pvglobal;
      pvlocal(0) = 0.;
      pvglobal(0) = 0.;
      for (int l=0;l<4;l++) {elocal(l) = 0.; eglobal(l) = 0.;}
      for (int l=0;l<ndensity;l++) {masslocal(l) = 0.; massglobal(l) = 0.;}
      for (int l=0;l<ndensity;l++) {densmaxlocal(l) = 0.; densmaxglobal(l) = 0.;}
      for (int l=0;l<ndensity;l++) {densminlocal(l) = 0.; densminglobal(l) = 0.;}
#if NTRACERS_FCT > 0
      for (int l=0;l<ndensityfct;l++) {massfctlocal(l) = 0.; massfctglobal(l) = 0.;}
      for (int l=0;l<ndensityfct;l++) {densfctmaxlocal(l) = 0.; densfctmaxglobal(l) = 0.;}
      for (int l=0;l<ndensityfct;l++) {densfctminlocal(l) = 0.; densfctminglobal(l) = 0.;}
#endif


//PRIMAL VS. DUAL HERE?
// IE WHAT DO WE ACTUALLY LOOP OVER?
      int is = topology->is;
      int js = topology->js;
      int ks = topology->ks;

      yakl::parallel_for("ComputeStats", topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

//ALL OF THIS NEEDS TO BE UPDATED...
        real eta, hv, q0, KE, PE;
        SArray<real,ndensity> dens0;
        SArray<real,1> zeta, h0im1, h0jm1;
        SArray<real,2> U, he;
        SArray<real,2,2> h0arr;

//WHICH GEOM HERE?
        //compute stats locally
        compute_I<ndensity,diff_ord> (dens0, progvars.fields_arr[DENSVAR].data, *this->geom, is, js, ks, i, j, k);
        compute_H<1,diff_ord> (U, progvars.fields_arr[VVAR].data, *this->geom, is, js, ks, i, j, k);

        //Slight hack that only computes I for first dof
        compute_I<1,diff_ord> (h0im1, progvars.fields_arr[DENSVAR].data, *this->geom, is, js, ks, i-1, j, k);
        compute_I<1,diff_ord> (h0jm1, progvars.fields_arr[DENSVAR].data, *this->geom, is, js, ks, i, j-1, k);

h0arr(0,0) = dens0(0);
h0arr(1,0) = dens0(0);
h0arr(0,1) = h0im1(0);
h0arr(1,1) = h0jm1(0);
phi(he, h0arr);

        KE = 1./2. * (he(0) * ( U(0) * progvars.fields_arr[VVAR].data(0,k+ks,j+js,i+is)) +
                    + he(1) * ( U(1) * progvars.fields_arr[VVAR].data(1,k+ks,j+js,i+is)));

#ifdef _AN
IE = 0.5*progvars.fields_arr[DENSVAR].data(0,k+ks,j+js,i+is)*dens0(1) + dens0(1)*constvars.fields_arr[HSVAR].data(0,k+ks,j+js,i+is);
#endif
#ifdef _FC
IE = 0.5*g*progvars.fields_arr[DENSVAR].data(0,k+ks,j+js,i+is)*dens0(0) + g*dens0(0)*constvars.fields_arr[HSVAR].data(0,k+ks,j+js,i+is);
#endif
#ifdef _MAN
IE = 0.5*g*progvars.fields_arr[DENSVAR].data(0,k+ks,j+js,i+is)*dens0(0) + g*dens0(0)*constvars.fields_arr[HSVAR].data(0,k+ks,j+js,i+is);
#endif
#ifdef _MFC
IE = 0.5*g*progvars.fields_arr[DENSVAR].data(0,k+ks,j+js,i+is)*dens0(0) + g*dens0(0)*constvars.fields_arr[HSVAR].data(0,k+ks,j+js,i+is);
#endif

        compute_D2<1>(zeta, progvars.fields_arr[VVAR].data, is, js, ks, i, j, k);
        eta = zeta(0) + constvars.fields_arr[CORIOLISVAR].data(0,k+ks,j+js,i+is);
        compute_R<0> (hv, progvars.fields_arr[DENSVAR].data, is, js, ks, i, j, k);
        q0 = eta / hv;

       elocal(0) += KE;
       elocal(1) += IE;
       elocal(2) += PE;
       elocal(3) += KE + IE + PE;
       pvlocal(0) += eta;

    });

    for (int l=0;l<ndensity;l++)
    {
    masslocal(l) = progvars.fields_arr[DENSVAR].sum(l);
    densmaxlocal(l) = progvars.fields_arr[DENSVAR].max(l);
    densminlocal(l) = progvars.fields_arr[DENSVAR].min(l);
    }

#if NTRACERS_FCT > 0
    for (int l=0;l<ndensityfct;l++)
    {
      massfctlocal(l) = progvars.fields_arr[DENSFCTVAR].sum(l);
      densfctmaxlocal(l) = progvars.fields_arr[DENSFCTVAR].max(l);
      densfctminlocal(l) = progvars.fields_arr[DENSFCTVAR].min(l);
    }
#endif

    //MPI sum/min/max
    this->ierr = MPI_Ireduce( &masslocal, &massglobal, ndensity, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[DENSSTAT]);
    this->ierr = MPI_Ireduce( &densmaxlocal, &densmaxglobal, ndensity, REAL_MPI, MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[DENSMAXSTAT]);
    this->ierr = MPI_Ireduce( &densminlocal, &densminglobal, ndensity, REAL_MPI, MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[DENSMINSTAT]);

#if NTRACERS_FCT > 0
    this->ierr = MPI_Ireduce( &massfctlocal, &massfctglobal, ndensityfct, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[DENSFCTSTAT]);
    this->ierr = MPI_Ireduce( &densfctmaxlocal, &densfctmaxglobal, ndensityfct, REAL_MPI, MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[DENSFCTMAXSTAT]);
    this->ierr = MPI_Ireduce( &densfctminlocal, &densfctminglobal, ndensityfct, REAL_MPI, MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[DENSFCTMINSTAT]);
#endif
    this->ierr = MPI_Ireduce( &pvlocal, &pvglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PVSTAT]);
    this->ierr = MPI_Ireduce( &elocal, &eglobal, 4, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[ESTAT]);

    this->ierr = MPI_Waitall(nstats, this->Req, this->Status);


  if (masterproc)
  {
    for (int l=0;l<ndensity;l++)
    {
  this->stats_arr[DENSSTAT].data(l,i) = massglobal(l);
  this->stats_arr[DENSMAXSTAT].data(l,i) = densmaxglobal(l);
  this->stats_arr[DENSMINSTAT].data(l,i) = densminglobal(l);
}

#if NTRACERS_FCT > 0
  for (int l=0;l<ndensityfct;l++)
  {
  this->stats_arr[DENSFCTSTAT].data(l,i) = massfctglobal(l);
  this->stats_arr[DENSFCTMAXSTAT].data(l,i) = densfctmaxglobal(l);
  this->stats_arr[DENSFCTMINSTAT].data(l,i) = densfctminglobal(l);
}
#endif
  this->stats_arr[ESTAT].data(0,i) = eglobal(0);
  this->stats_arr[ESTAT].data(1,i) = eglobal(1);
  this->stats_arr[ESTAT].data(2,i) = eglobal(2);
  this->stats_arr[ESTAT].data(3,i) = eglobal(3);
  this->stats_arr[PVSTAT].data(0,i) = pvglobal(0);
  }
  }
};




// *******   VariableSet Initialization   ***********//
template <uint nprog, uint nconst, uint naux, uint ndiag> void initialize_variables(const Topology &ptopo, const Topology &dtopo,
SArray<int, nprog, 6> &prog_ndofs_arr, SArray<int, nconst, 6> &const_ndofs_arr, SArray<int, naux, 6> &aux_ndofs_arr, SArray<int, ndiag, 6> &diag_ndofs_arr,
std::array<std::string, nprog> &prog_names_arr, std::array<std::string, nconst> &const_names_arr, std::array<std::string, naux> &aux_names_arr, std::array<std::string, ndiag> &diag_names_arr,
std::array<const Topology *, nprog> &prog_topo_arr, std::array<const Topology *, nconst> &const_topo_arr, std::array<const Topology *, naux> &aux_topo_arr, std::array<const Topology *, ndiag> &diag_topo_arr)
{

//primal grid represents straight quantities, dual grid twisted quantities
// 0,1,2,3 ndofs for basemesh, 0,1 ndofs for 1D vert mesh (stored at 4,5)

//CAREFULLY CHECK/THINK ABOUT ALL THIS...
  // v, w, dens
  prog_topo_arr[VVAR] = &ptopo;
  prog_topo_arr[WVAR] = &ptopo;
  prog_topo_arr[DENSVAR] = &dtopo;
  prog_names_arr[VVAR] = "v";
  prog_names_arr[VVAR] = "w";
  prog_names_arr[DENSVAR] = "dens";
//ALL FUCKED
  prog_ndofs_arr(VVAR,1) = 1; //v = straight 1-form/1-form
  prog_ndofs_arr(VVAR,5) = 1; //v = straight 1-form/1-form
  prog_ndofs_arr(WVAR,1) = 1; //w = straight n-form/0-form
  prog_ndofs_arr(WVAR,1) = 1; //w = straight n-form/0-form
  prog_ndofs_arr(DENSVAR,2) = ndensity; //dens = twisted 2-form
  prog_ndofs_arr(DENSVAR,2) = ndensity; //dens = twisted 2-form

  // hs, coriolis
  const_topo_arr[TOPOGRAPHYVAR] = &topo;
  const_topo_arr[CORIOLISVAR] = &topo;
  const_names_arr[HSVAR] = "hs";
  const_names_arr[CORIOLISVAR] = "coriolis";
  const_ndofs_arr(HSVAR,2) = 1; //hs = twisted 2-form
  const_ndofs_arr(CORIOLISVAR,0) = 1; //f = straight 2-form

  //functional derivatives = F, B, K, he, U
  aux_topo_arr[BVAR] = &topo;
  aux_topo_arr[FVAR] = &topo;
  aux_topo_arr[UVAR] = &topo;
  aux_topo_arr[HEVAR] = &topo;
  aux_topo_arr[KVAR] = &topo;
  aux_names_arr[KVAR] = "K";
  aux_names_arr[BVAR] = "B";
  aux_names_arr[FVAR] = "F";
  aux_names_arr[UVAR] = "U";
  aux_names_arr[HEVAR] = "he";
  aux_ndofs_arr(BVAR,2) = ndensity; //B = straight 0-form
  aux_ndofs_arr(KVAR,2) = 1; //K = twisted 2-form
  aux_ndofs_arr(FVAR,1) = 1; //F = twisted 1-form
  aux_ndofs_arr(UVAR,1) = 1; //U = twisted 1-form
  aux_ndofs_arr(HEVAR,1) = 1; //he lives on edges

  //dens primal grid reconstruction stuff- dens0, edgerecon, recon
  aux_topo_arr[DENSRECONVAR] = &topo;
  aux_topo_arr[DENSEDGERECONVAR] = &topo;
  aux_topo_arr[DENS0VAR] = &topo;
  aux_names_arr[DENS0VAR] = "dens0";
  aux_names_arr[DENSRECONVAR] = "densrecon";
  aux_names_arr[DENSEDGERECONVAR] = "densedgerecon";
  aux_ndofs_arr(DENSRECONVAR,1) = ndensity; //densrecon lives on edges
  aux_ndofs_arr(DENSEDGERECONVAR,2) = 4*ndensity; //densedgerecon lives on cells
  aux_ndofs_arr(DENS0VAR,2) = ndensity; //dens0 = straight 0-form

  //dual grid reconstruction stuff- q0, f0, FT, qedgerecon, qrecon, coriolisedgercon, coriolisrecon
  aux_topo_arr[FTVAR] = &topo;
  aux_topo_arr[CORIOLISRECONVAR] = &topo;
  aux_topo_arr[CORIOLISEDGERECONVAR] = &topo;
  aux_topo_arr[Q0VAR] = &topo;
  aux_topo_arr[F0VAR] = &topo;
  aux_topo_arr[QRECONVAR] = &topo;
  aux_topo_arr[QEDGERECONVAR] = &topo;
  aux_names_arr[FTVAR] = "FT";
  aux_names_arr[CORIOLISRECONVAR] = "coriolisrecon";
  aux_names_arr[CORIOLISEDGERECONVAR] = "coriolisedgerecon";
  aux_names_arr[Q0VAR] = "q";
  aux_names_arr[F0VAR] = "f";
  aux_names_arr[QRECONVAR] = "qrecon";
  aux_names_arr[QEDGERECONVAR] = "qedgerecon";
  aux_ndofs_arr(FTVAR,1) = 1; //FT = straight 1-form
  aux_ndofs_arr(Q0VAR,0) = 1; //q0 = twisted 0-form
  aux_ndofs_arr(F0VAR,0) = 1; //f0 = twisted 0-form
  aux_ndofs_arr(QRECONVAR,1) = 1; //qrecon lives on edges
  aux_ndofs_arr(QEDGERECONVAR,0) = 4; //qedgerecon lives on dual cells
  aux_ndofs_arr(CORIOLISRECONVAR,1) = 1; //coriolisrecon lives on edges
  aux_ndofs_arr(CORIOLISEDGERECONVAR,0) = 4; //coriolisedgerecon lives on dual cells

  // q, concentration 0-forms for dens
  diag_topo_arr[QDIAGVAR] = &topo;
  diag_topo_arr[DENSLDIAGVAR] = &topo;
  diag_names_arr[QDIAGVAR] = "q";
  diag_names_arr[DENSLDIAGVAR] = "densl";
  diag_ndofs_arr(QDIAGVAR,0) = 1; //qdiag = twisted 0-form
  diag_ndofs_arr(DENSLDIAGVAR,2) = ndensity; //densldiag = straight 0-form

  //densfct stuff- densfct, BFCT, densfct0, edgereconfct, reconfct, Phi, Mf, edgeflux, concentration 0-forms for densfct
#if NTRACERS_FCT > 0
  prog_topo_arr[DENSFCTVAR] = &topo;
  prog_names_arr[DENSFCTVAR] = "densfct";
  prog_ndofs_arr(DENSFCTVAR,2) = ndensityfct; //densfct = twisted 2-form

  aux_topo_arr[BFCTVAR] = &topo;
  aux_topo_arr[DENSFCTRECONVAR] = &topo;
  aux_topo_arr[DENSFCTEDGERECONVAR] = &topo;
  aux_topo_arr[DENSFCT0VAR] = &topo;
  aux_topo_arr[PHIVAR] = &topo;
  aux_topo_arr[MFVAR] = &topo;
  aux_topo_arr[EDGEFLUXVAR] = &topo;
  aux_names_arr[BFCTVAR] = "Bfct";
  aux_names_arr[DENSFCT0VAR] = "densfct0";
  aux_names_arr[DENSFCTRECONVAR] = "densfctrecon";
  aux_names_arr[DENSFCTEDGERECONVAR] = "densfctedgerecon";
  aux_names_arr[PHIVAR] = "Phi";
  aux_names_arr[MFVAR] = "Mf";
  aux_names_arr[EDGEFLUXVAR] = "edgeflux";
  aux_ndofs_arr(BFCTVAR,2) = ndensityfct; //Bfct = straight 0-form
  aux_ndofs_arr(DENSFCTRECONVAR,1) = ndensityfct; //densfctrecon lives on edges
  aux_ndofs_arr(DENSFCTEDGERECONVAR,2) = 4*ndensityfct; //densfctedgerecon lives on cells
  aux_ndofs_arr(DENSFCT0VAR,2) = ndensityfct; //densfct0 = straight 0-form
  aux_ndofs_arr(PHIVAR,1) = ndensityfct;
  aux_ndofs_arr(MFVAR,2) = ndensityfct;
  aux_ndofs_arr(EDGEFLUXVAR,1) = ndensityfct;

  diag_topo_arr[DENSFCTLDIAGVAR] = &topo;
  diag_names_arr[DENSFCTLDIAGVAR] = "densfctl";
  diag_ndofs_arr(DENSFCTLDIAGVAR,2) = ndensityfct; //densfctldiag = straight 0-form
#endif


}

  // *******   Initial Conditions   ***********//

//ALL BROKEN
real YAKL_INLINE rb_p(real x, real y, real z)
{
    real xprime1 = Lx / (M_PI * sigmax) * sin(M_PI / Lx * (x - xc1));
    real yprime1 = Ly / (M_PI * sigmay) * sin(M_PI / Ly * (y - yc1));
    real xprime2 = Lx / (M_PI * sigmax) * sin(M_PI / Lx * (x - xc2));
    real yprime2 = Ly / (M_PI * sigmay) * sin(M_PI / Ly * (y - yc2));
    real xprimeprime1 = Lx / (2.0 * M_PI * sigmax) * sin(2 * M_PI / Lx * (x - xc1));
    real yprimeprime1 = Ly / (2.0 * M_PI * sigmay) * sin(2 * M_PI / Ly * (y - yc1));
    real xprimeprime2 = Lx / (2.0 * M_PI * sigmax) * sin(2 * M_PI / Lx * (x - xc2));
    real yprimeprime2 = Ly / (2.0 * M_PI * sigmay) * sin(2 * M_PI / Ly * (y - yc2));

    return H0 - dh * (exp(-0.5 * (xprime1 * xprime1 + yprime1 * yprime1)) + exp(-0.5 * (xprime2 * xprime2 + yprime2 * yprime2)) - 4. * M_PI * sigmax * sigmay / Lx / Ly);
}

vec<2> YAKL_INLINE rb_v(real x, real y) {
  vec<2> vvec;

  real xprime1 = Lx / (M_PI * sigmax) * sin(M_PI / Lx * (x - xc1));
  real yprime1 = Ly / (M_PI * sigmay) * sin(M_PI / Ly * (y - yc1));
  real xprime2 = Lx / (M_PI * sigmax) * sin(M_PI / Lx * (x - xc2));
  real yprime2 = Ly / (M_PI * sigmay) * sin(M_PI / Ly * (y - yc2));
  real xprimeprime1 = Lx / (2.0 * M_PI * sigmax) * sin(2 * M_PI / Lx * (x - xc1));
  real yprimeprime1 = Ly / (2.0 * M_PI * sigmay) * sin(2 * M_PI / Ly * (y - yc1));
  real xprimeprime2 = Lx / (2.0 * M_PI * sigmax) * sin(2 * M_PI / Lx * (x - xc2));
  real yprimeprime2 = Ly / (2.0 * M_PI * sigmay) * sin(2 * M_PI / Ly * (y - yc2));

  vvec.u = - g * dh / coriolis / sigmay * (yprimeprime1 * exp(-0.5*(xprime1 * xprime1 + yprime1 * yprime1)) + yprimeprime2 * exp(-0.5*(xprime2 * xprime2 + yprime2 * yprime2)));
  vvec.v = g * dh / coriolis / sigmax * (xprimeprime1 * exp(-0.5*(xprime1 * xprime1 + yprime1 * yprime1)) + xprimeprime2 * exp(-0.5*(xprime2 * xprime2 + yprime2 * yprime2)));
  return vvec;
}


template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<nprog> &progvars, VariableSet<nconst> &constvars, 
  Geometry<2, nquadx, nquady, nquadz> &primal_geom, Geometry<2, nquadx, nquady, nquadz> &dual_geom)
{

//THIS IS ALL BROKEN
//HOW SHOULD IT REALLY WORK

//Each initial condition gives a p,T,v,moisture (q? rho_s?) fields as a function of position; and also sets topography + coriolis
//From this, use THERMO to compute whatever our entropic variable is
//Then can either:
//1) get rho from p,T using THERMO, or 
//2) Compute a hydrostatically balanced rho given entropic variable + moisture
//Then get rho_s from rho and q

//What happens for AN? We set p=pref in initial conditions, I think?
//Generally we have IC's in hydrostatic balance + some perturbation
//So use p0 as pref in AN? etc.

    if (params.data_init_cond == DATA_INIT::RB)
    {
        dual_geom.set_primal_2form_values(rb_rho, progvars.fields_arr[DENSVAR], 0);
        dual_geom.set_primal_2form_values(rb_S, progvars.fields_arr[DENSVAR], 0);
        primal_geom.set_dual_1form_values(double_vortex_v, progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
        primal_geom.set_dual_1form_values(double_vortex_v, progvars.fields_arr[WVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);

//here we set density...
#ifdef _FC
#endif
#ifdef _MFC
#endif

    }




}



#endif
