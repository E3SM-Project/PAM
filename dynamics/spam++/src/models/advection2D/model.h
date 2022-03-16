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


// We assume that NTRACERS, NTRACERS_FCT and NTRACERS_Q are all >0
// THIS BREAKS 1D ADVECTION
// IS THIS OK?
// Can do a hacky fix that ignores all the Q var stuff for 1D?
// Really the interesting 1D advection would be a slice-type variant?
// Maybe merge with this?

// Number of variables

uint constexpr nprognostic = 3;
uint constexpr nconstant = 3;
uint constexpr nauxiliary = 13;
uint constexpr ndiagnostic = 3;
uint constexpr nstats = 7;

#define TVAR 0
#define TFCTVAR 1
#define QVAR 2

#define VVAR 0
#define UVAR 1
#define UTVAR 2

#define T0VAR 0
#define TRECONVAR 1
#define TEDGERECONVAR 2
#define TFCT0VAR 3
#define TFCTRECONVAR 4
#define TFCTEDGERECONVAR 5
#define PHIVAR 6
#define EDGEFLUXVAR 7
#define MFVAR 8
#define Q0VAR 9
#define QRECONVAR 10
#define QEDGERECONVAR 11
#define QFLUXVAR 12

#define TDIAGVAR 0
#define TFCTDIAGVAR 1
#define QDIAGVAR 2

#define TMASSSTAT 0
#define TMINSTAT 1
#define TMAXSTAT 2
#define TFCTMASSSTAT 3
#define TFCTMINSTAT 4
#define TFCTMAXSTAT 5
#define QMASSSTAT 6




// Initial conditions related variables and functions
#define gaussian_1d(x)     (1. * exp(-100. * pow(x-0.5,2.)))
#define gaussian_2d(x,y)   (1. * exp(-100. * pow(x-0.5,2.)) * exp(-100. * pow(y-0.5,2.)))
real YAKL_INLINE gaussian(real x)                 { return gaussian_1d(x); }
real YAKL_INLINE gaussian(real x, real y)         { return gaussian_2d(x,y); }

#define vortex1_1d(x)     (1. *  exp(-100. * pow(x-0.75,2.)))
#define vortex1_2d(x,y)   (1. *  exp(-100. * pow(x-0.75,2.)) * exp(-100. * pow(y-0.75,2.)))
#define vortex2_1d(x)     (0.5 * exp(-50.  * pow(x-0.25,2.)))
#define vortex2_2d(x,y)   (0.5 * exp(-50.  * pow(x-0.25,2.)) * exp(-75.  * pow(y-0.25,2.)))
real YAKL_INLINE vortices(real x)                 { return vortex1_1d(x)     + vortex2_1d(x); }
real YAKL_INLINE vortices(real x, real y)         { return vortex1_2d(x,y)   + vortex2_2d(x,y); }

real YAKL_INLINE square(real x)                 {return (x > 0.35 && x < 0.65                                                ) ? 1. : 0.;}
real YAKL_INLINE square(real x, real y)         {return (x > 0.35 && x < 0.65 && y > 0.35 && y < 0.65                        ) ? 1. : 0.;}





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
  vvec.u = -C_UNIFORM_WIND/sqrt(2.);
  vvec.v = C_UNIFORM_WIND/sqrt(2.);
  return vvec;
}




// FIX THIS
vec<2> YAKL_INLINE deformational_wind(real x, real y) {
  vec<2> vvec;
  vvec.v = C_UNIFORM_WIND;
  return vvec;
}

// *******   Model Specific Parameters   ***********//


void set_model_specific_params(std::string inFile, ModelParameters &params)
{


  std::string strDataInit1 = "";
  std::string strDataInit2 = "";
  std::string strDataInit3 = "";
  std::string strDataFCTInit1 = "";
  std::string strDataFCTInit2 = "";
  std::string strDataFCTInit3 = "";
  std::string strDataQInit1 = "";
  std::string strDataQInit2 = "";
  std::string strDataQInit3 = "";
  std::string strWindInit = "";

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
           if ( !strcmp( "dataInit1"   , key.c_str() ) ) { ssVal >> strDataInit1       ;}
      else if ( !strcmp( "dataInit2"   , key.c_str() ) ) { ssVal >> strDataInit2       ;}
      else if ( !strcmp( "dataInit3"   , key.c_str() ) ) { ssVal >> strDataInit3       ;}
      else if ( !strcmp( "dataFCTInit1"   , key.c_str() ) ) { ssVal >> strDataFCTInit1       ;}
      else if ( !strcmp( "dataFCTInit2"   , key.c_str() ) ) { ssVal >> strDataFCTInit2       ;}
      else if ( !strcmp( "dataFCTInit3"   , key.c_str() ) ) { ssVal >> strDataFCTInit3       ;}
      else if ( !strcmp( "dataQInit1"   , key.c_str() ) ) { ssVal >> strDataQInit1       ;}
      else if ( !strcmp( "dataQInit2"   , key.c_str() ) ) { ssVal >> strDataQInit2       ;}
      else if ( !strcmp( "dataQInit3"   , key.c_str() ) ) { ssVal >> strDataQInit3       ;}
      else if ( !strcmp( "windInit"    , key.c_str() ) ) { ssVal >> strWindInit        ;}


      //else {
      //  std::cout << "Error: key " << key << " not understood in file " << inFile << "\n";
      //}
    }
  }

  if (!strcmp("",strDataInit1.c_str())) { std::cout << "Error: key " << "dataInit1" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataInit2.c_str())) { std::cout << "Error: key " << "dataInit2" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataInit3.c_str())) { std::cout << "Error: key " << "dataInit3" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataFCTInit1.c_str())) { std::cout << "Error: key " << "dataFCTInit1" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataFCTInit2.c_str())) { std::cout << "Error: key " << "dataFCTInit2" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataFCTInit3.c_str())) { std::cout << "Error: key " << "dataFCTInit3" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataQInit1.c_str())) { std::cout << "Error: key " << "dataQInit1" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataQInit2.c_str())) { std::cout << "Error: key " << "dataQInit2" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataQInit3.c_str())) { std::cout << "Error: key " << "dataQInit3" << " not set.\n"; exit(-1); }
  if (!strcmp("",strWindInit.c_str()))  { std::cout << "Error: key " << "windInit"  << " not set.\n"; exit(-1); }

    size_t splitloc = strDataInit1.find("//",0);
    std::string sub_str;
    if (splitloc != std::string::npos){
      sub_str = strDataInit1.substr(0,splitloc);
    } else {
      sub_str = strDataInit1;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.data_init_cond[0] = DATA_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.data_init_cond[0] = DATA_INIT::VORTICES  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.data_init_cond[0] = DATA_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized dataInit " << strDataInit1 << "\n";
      exit(-1);
    }

    splitloc = strDataInit2.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strDataInit2.substr(0,splitloc);
    } else {
      sub_str = strDataInit2;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.data_init_cond[1] = DATA_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.data_init_cond[1] = DATA_INIT::VORTICES  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.data_init_cond[1] = DATA_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized dataInit " << strDataInit2 << "\n";
      exit(-1);
    }

    splitloc = strDataInit3.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strDataInit3.substr(0,splitloc);
    } else {
      sub_str = strDataInit3;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.data_init_cond[2] = DATA_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.data_init_cond[2] = DATA_INIT::VORTICES  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.data_init_cond[2] = DATA_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized dataInit " << strDataInit3 << "\n";
      exit(-1);
    }

    splitloc = strDataFCTInit1.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strDataFCTInit1.substr(0,splitloc);
    } else {
      sub_str = strDataFCTInit1;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.dataFCT_init_cond[0] = DATA_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.dataFCT_init_cond[0] = DATA_INIT::VORTICES  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.dataFCT_init_cond[0] = DATA_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized dataInit " << strDataFCTInit1 << "\n";
      exit(-1);
    }

    splitloc = strDataFCTInit2.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strDataFCTInit2.substr(0,splitloc);
    } else {
      sub_str = strDataFCTInit2;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.dataFCT_init_cond[1] = DATA_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.dataFCT_init_cond[1] = DATA_INIT::VORTICES  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.dataFCT_init_cond[1] = DATA_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized dataInit " << strDataFCTInit2 << "\n";
      exit(-1);
    }

    splitloc = strDataFCTInit3.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strDataFCTInit3.substr(0,splitloc);
    } else {
      sub_str = strDataFCTInit3;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.dataFCT_init_cond[2] = DATA_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.dataFCT_init_cond[2] = DATA_INIT::VORTICES  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.dataFCT_init_cond[2] = DATA_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized dataInit " << strDataFCTInit3 << "\n";
      exit(-1);
    }

    splitloc = strDataQInit1.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strDataQInit1.substr(0,splitloc);
    } else {
      sub_str = strDataQInit1;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.dataQ_init_cond[0] = DATA_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.dataQ_init_cond[0] = DATA_INIT::VORTICES  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.dataQ_init_cond[0] = DATA_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized dataInit " << strDataQInit1 << "\n";
      exit(-1);
    }

    splitloc = strDataQInit2.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strDataQInit2.substr(0,splitloc);
    } else {
      sub_str = strDataQInit2;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.dataQ_init_cond[1] = DATA_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.dataQ_init_cond[1] = DATA_INIT::VORTICES  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.dataQ_init_cond[1] = DATA_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized dataInit " << strDataQInit2 << "\n";
      exit(-1);
    }

    splitloc = strDataQInit3.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strDataQInit3.substr(0,splitloc);
    } else {
      sub_str = strDataQInit3;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.dataQ_init_cond[2] = DATA_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.dataQ_init_cond[2] = DATA_INIT::VORTICES  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.dataQ_init_cond[2] = DATA_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized dataInit " << strDataQInit3 << "\n";
      exit(-1);
    }

    splitloc = strWindInit.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strWindInit.substr(0,splitloc);
    } else {
      sub_str = strWindInit;
    }
    if      ( !strcmp(sub_str.c_str(),"uniform_x"     ) ) { params.wind_init_cond = WIND_INIT::UNIFORM_X     ; }
    else if ( !strcmp(sub_str.c_str(),"uniform_y"     ) ) { params.wind_init_cond = WIND_INIT::UNIFORM_Y     ; }
    else if ( !strcmp(sub_str.c_str(),"uniform_z"     ) ) { params.wind_init_cond = WIND_INIT::UNIFORM_Z     ; }
    else if ( !strcmp(sub_str.c_str(),"uniform_xy"    ) ) { params.wind_init_cond = WIND_INIT::UNIFORM_XY    ; }
    else if ( !strcmp(sub_str.c_str(),"deformational" ) ) { params.wind_init_cond = WIND_INIT::DEFORMATIONAL ; }
    else  {
      std::cout << "Error: unrecognized windInit " << strWindInit << "\n";
      exit(-1);
    }

  params.etime = 0.0;


  params.xlen = 1.0;
  params.xc = 0.5;
  if (ndims>=2)
  {

    params.ylen = 1.0;
    params.yc = 0.5;
  }
  if (ndims ==3)
  {
  params.zlen = 1.0;
  params.zc = 0.5;
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

      void compute_diag(const VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<ndiag> &diagnostic_vars)
      {

        // Compute T0 and TFCT0
        int pis = primal_topology->is;
        int pjs = primal_topology->js;
        int pks = primal_topology->ks;
        
        yakl::parallel_for("ComputeT0", primal_topology->n_cells, YAKL_LAMBDA (int iGlob) {
          int k, j, i;
          yakl::unpackIndices(iGlob, primal_topology->n_cells_z, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);

        compute_I<ntdofs, diff_ord>(diagnostic_vars.fields_arr[TDIAGVAR].data, x.fields_arr[TVAR].data, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k);
        compute_I<ntfctdofs, diff_ord>(diagnostic_vars.fields_arr[TFCTDIAGVAR].data, x.fields_arr[TFCTVAR].data, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k);
        });


        // Compute q0
        int dis = dual_topology->is;
        int djs = dual_topology->js;
        int dks = dual_topology->ks;

        yakl::parallel_for("ComputeQ0", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
          int k, j, i;
          yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

        compute_J<nQdofs, diff_ord>(diagnostic_vars.fields_arr[QDIAGVAR].data, x.fields_arr[QVAR].data, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k);

        });


      }
};

// *******   Tendencies   ***********//

template <uint nprog, uint nconst, uint naux> class Tendencies {
public:

  const Topology *primal_topology;
  const Topology *dual_topology;
  ExchangeSet<naux> *aux_exchange;
  ExchangeSet<nconst> *const_exchange;
  Geometry<ndims,1,1,1> *primal_geometry;
  Geometry<ndims,1,1,1> *dual_geometry;

  TransformMatrices<real> trans;

  SArray<real,reconstruction_order,2> to_gll;
  SArray<real,reconstruction_order,reconstruction_order,reconstruction_order> wenoRecon;
  SArray<real,(reconstruction_order-1)/2+2> wenoIdl;
  real wenoSigma;

  SArray<real,dual_reconstruction_order,2> dual_to_gll;
  SArray<real,dual_reconstruction_order,dual_reconstruction_order,dual_reconstruction_order> dualwenoRecon;
  SArray<real,(dual_reconstruction_order-1)/2+2> dualwenoIdl;
  real dualwenoSigma;


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

    // Setup the matrix to transform a stencil of ord cell averages into tord GLL points
    trans.coefs_to_gll_lower( to_gll );
    trans.weno_sten_to_coefs(wenoRecon);
    wenoSetIdealSigma<reconstruction_order>(wenoIdl,wenoSigma);

    // Setup the matrix to transform a stencil of ord cell averages into tord GLL points
    trans.coefs_to_gll_lower( dual_to_gll );
    trans.weno_sten_to_coefs(dualwenoRecon);
    wenoSetIdealSigma<dual_reconstruction_order>(dualwenoIdl,dualwenoSigma);
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

    // Compute U = H V

    yakl::parallel_for("ComputeUVAR", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

    compute_H<1, diff_ord>(const_vars.fields_arr[UVAR].data, const_vars.fields_arr[VVAR].data, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k);

    });
    this->const_exchange->exchanges_arr[UVAR].exchange_field(const_vars.fields_arr[UVAR]);

    //compute UT = W U

    yakl::parallel_for("ComputeUTVAR", primal_topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, primal_topology->n_cells_z, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);

compute_W(const_vars.fields_arr[UTVAR].data,const_vars.fields_arr[UVAR].data, pis, pjs, pks, i, j, k);

 });
    this->const_exchange->exchanges_arr[UTVAR].exchange_field(const_vars.fields_arr[UTVAR]);


  }
  
  
  void compute_rhs(real dt, VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<naux> &auxiliary_vars, VariableSet<nprog> &xtend)
  {

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

// Compute Q0

yakl::parallel_for("ComputeQ0VAR", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

compute_J<nQdofs, diff_ord>(auxiliary_vars.fields_arr[Q0VAR].data, x.fields_arr[QVAR].data, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k);

});

// Compute T0, TFCT0

yakl::parallel_for("ComputeT0TFCT0VAR", primal_topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, primal_topology->n_cells_z, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);

compute_I<ntdofs, diff_ord>(auxiliary_vars.fields_arr[T0VAR].data, x.fields_arr[TVAR].data, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k);
compute_I<ntfctdofs, diff_ord>(auxiliary_vars.fields_arr[TFCT0VAR].data, x.fields_arr[TFCTVAR].data, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k);

});

this->aux_exchange->exchanges_arr[T0VAR].exchange_field(auxiliary_vars.fields_arr[T0VAR]);
this->aux_exchange->exchanges_arr[TFCT0VAR].exchange_field(auxiliary_vars.fields_arr[TFCT0VAR]);
this->aux_exchange->exchanges_arr[Q0VAR].exchange_field(auxiliary_vars.fields_arr[Q0VAR]);



//compute qrecon

yakl::parallel_for("ComputeQEdgeRecon", primal_topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, primal_topology->n_cells_z, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);

  compute_straight_edge_recon<nQdofs, reconstruction_type, reconstruction_order>(
    auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[Q0VAR].data, pis, pjs, pks, i, j, k,
    wenoRecon, to_gll, wenoIdl, wenoSigma);

});

this->aux_exchange->exchanges_arr[QEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[QEDGERECONVAR]);

yakl::parallel_for("ComputeQRECON", primal_topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, primal_topology->n_cells_z, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);

  compute_straight_recon<nQdofs, reconstruction_type>(
    auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[QEDGERECONVAR].data, const_vars.fields_arr[UTVAR].data, pis, pjs, pks, i, j, k);

});
this->aux_exchange->exchanges_arr[QRECONVAR].exchange_field(auxiliary_vars.fields_arr[QRECONVAR]);


//compute Trecon/Tfctrecon

yakl::parallel_for("ComputeTEdgeRecon", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

  compute_twisted_edge_recon<ntdofs, dual_reconstruction_type, dual_reconstruction_order>(
    auxiliary_vars.fields_arr[TEDGERECONVAR].data, auxiliary_vars.fields_arr[T0VAR].data, dis, djs, dks, i, j, k,
    dualwenoRecon, dual_to_gll, dualwenoIdl, dualwenoSigma);

  compute_twisted_edge_recon<ntfctdofs, dual_reconstruction_type, dual_reconstruction_order>(
    auxiliary_vars.fields_arr[TFCTEDGERECONVAR].data, auxiliary_vars.fields_arr[TFCT0VAR].data, dis, djs, dks, i, j, k,
    dualwenoRecon, dual_to_gll, dualwenoIdl, dualwenoSigma);

});

this->aux_exchange->exchanges_arr[TEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[TEDGERECONVAR]);
this->aux_exchange->exchanges_arr[TFCTEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[TFCTEDGERECONVAR]);

yakl::parallel_for("ComputeTRECON", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

  compute_twisted_recon<ntdofs, dual_reconstruction_type>(
    auxiliary_vars.fields_arr[TRECONVAR].data, auxiliary_vars.fields_arr[TEDGERECONVAR].data, const_vars.fields_arr[UVAR].data, dis, djs, dks, i, j, k);

  compute_twisted_recon<ntfctdofs, dual_reconstruction_type>(
    auxiliary_vars.fields_arr[TFCTRECONVAR].data, auxiliary_vars.fields_arr[TFCTEDGERECONVAR].data, const_vars.fields_arr[UVAR].data, dis, djs, dks, i, j, k);

});
this->aux_exchange->exchanges_arr[TRECONVAR].exchange_field(auxiliary_vars.fields_arr[TRECONVAR]);
this->aux_exchange->exchanges_arr[TFCTRECONVAR].exchange_field(auxiliary_vars.fields_arr[TFCTRECONVAR]);


// compute FCT stuff
      yakl::parallel_for("ComputeEdgeFlux", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
      compute_edgefluxes<ntfctdofs> (auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[TFCTRECONVAR].data, const_vars.fields_arr[UVAR].data, dis, djs, dks, i, j, k);
      });
this->aux_exchange->exchanges_arr[EDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[EDGEFLUXVAR]);


yakl::parallel_for("ComputeMf", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  compute_Mf<ntfctdofs> (auxiliary_vars.fields_arr[MFVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, dt, dis, djs, dks, i, j, k);
});

    this->aux_exchange->exchanges_arr[MFVAR].exchange_field(auxiliary_vars.fields_arr[MFVAR]);

    yakl::parallel_for("ComputePhi", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
      compute_Phi<ntfctdofs> (auxiliary_vars.fields_arr[PHIVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[TFCTVAR].data, dis, djs, dks, i, j, k);
    });


    this->aux_exchange->exchanges_arr[PHIVAR].exchange_field(auxiliary_vars.fields_arr[PHIVAR]);


//compute Q flux
   yakl::parallel_for("ComputeQFLUX", primal_topology->n_cells, YAKL_LAMBDA (int iGlob) {
     int k, j, i;
     yakl::unpackIndices(iGlob, primal_topology->n_cells_z, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);

     if (qf_choice == QF_MODE::EC)
     { compute_Q_EC<nQdofs>(auxiliary_vars.fields_arr[QFLUXVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, const_vars.fields_arr[UVAR].data, pis, pjs, pks, i, j, k);}

     if (qf_choice == QF_MODE::NOEC)
     { compute_Q_nonEC<nQdofs>(auxiliary_vars.fields_arr[QFLUXVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, const_vars.fields_arr[UVAR].data, pis, pjs, pks, i, j, k);}

});
   this->aux_exchange->exchanges_arr[QFLUXVAR].exchange_field(auxiliary_vars.fields_arr[QFLUXVAR]);


// Compute D2 Q
yakl::parallel_for("ComputeQTend", primal_topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, primal_topology->n_cells_z, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
compute_D2<nQdofs> (xtend.fields_arr[QVAR].data, auxiliary_vars.fields_arr[QFLUXVAR].data, pis, pjs, pks, i, j, k);
});

// compute D2bar "F"
    yakl::parallel_for("ComputeTTend", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

      compute_wDbar2<ntdofs> (xtend.fields_arr[TVAR].data, auxiliary_vars.fields_arr[TRECONVAR].data, const_vars.fields_arr[UVAR].data, dis, djs, dks, i, j, k);
      compute_wDbar2_fct<ntfctdofs> (xtend.fields_arr[TFCTVAR].data, auxiliary_vars.fields_arr[TFCTRECONVAR].data, auxiliary_vars.fields_arr[PHIVAR].data, const_vars.fields_arr[UVAR].data, dis, djs, dks, i, j, k);
    });
    
}
};

// *******   Statistics Calculations   ***********//
class Stat
{
public:
  realArr data;
  int ndofs;

  std::string name;

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

  void initialize(ModelParameters &params, Parallel &par, const Topology &primal_topo, const Topology &dual_topo, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom)
  {
    statsize = params.Nsteps/params.Nstat + 1;
    stats_arr[TMASSSTAT].initialize("Tmass", ntdofs, params, par);
    stats_arr[TMINSTAT].initialize("Tmin", ntdofs, params, par);
    stats_arr[TMAXSTAT].initialize("Tmax", ntdofs, params, par);
    stats_arr[TFCTMASSSTAT].initialize("Tfctmass", ntfctdofs, params, par);
    stats_arr[TFCTMINSTAT].initialize("Tfctmin", ntfctdofs, params, par);
    stats_arr[TFCTMAXSTAT].initialize("Tfctmax", ntfctdofs, params, par);
    stats_arr[QMASSSTAT].initialize("Qmass", nQdofs, params, par);
    masterproc = par.masterproc;
  }



  void compute( VariableSet<nprog> &progvars,  VariableSet<nconst> &constvars, int i)
  {

      SArray<real,ntdofs> masslocal, massglobal;
      SArray<real,ntdofs> maxlocal, maxglobal;
      SArray<real,ntdofs> minlocal, minglobal;
      SArray<real,ntfctdofs> massfctlocal, massfctglobal;
      SArray<real,ntfctdofs> maxfctlocal, maxfctglobal;
      SArray<real,ntfctdofs> minfctlocal, minfctglobal;
      SArray<real,nQdofs> massQlocal, massQglobal;

      for (int l=0;l<ntdofs;l++)
      {
        masslocal(l) = progvars.fields_arr[TVAR].sum(l);
        maxlocal(l) = progvars.fields_arr[TVAR].max(l);
        minlocal(l) = progvars.fields_arr[TVAR].min(l);
      }
      for (int l=0;l<ntfctdofs;l++)
      {
        massfctlocal(l) = progvars.fields_arr[TFCTVAR].sum(l);
        maxfctlocal(l) = progvars.fields_arr[TFCTVAR].max(l);
        minfctlocal(l) = progvars.fields_arr[TFCTVAR].min(l);
      }
      for (int l=0;l<nQdofs;l++)
      {
        massQlocal(l) = progvars.fields_arr[QVAR].sum(l);
      }
    //MPI sum/min/max
    this->ierr = MPI_Ireduce( &masslocal, &massglobal, ntdofs, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[TMASSSTAT]);
    this->ierr = MPI_Ireduce( &maxlocal, &maxglobal, ntdofs, REAL_MPI, MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[TMAXSTAT]);
    this->ierr = MPI_Ireduce( &minlocal, &minglobal, ntdofs, REAL_MPI, MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[TMINSTAT]);
    this->ierr = MPI_Ireduce( &massfctlocal, &massfctglobal, ntfctdofs, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[TFCTMASSSTAT]);
    this->ierr = MPI_Ireduce( &maxfctlocal, &maxfctglobal, ntfctdofs, REAL_MPI, MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[TFCTMAXSTAT]);
    this->ierr = MPI_Ireduce( &minfctlocal, &minfctglobal, ntfctdofs, REAL_MPI, MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[TFCTMINSTAT]);
    this->ierr = MPI_Ireduce( &massQlocal, &massQglobal, nQdofs, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[QMASSSTAT]);
    this->ierr = MPI_Waitall(nstats, this->Req, this->Status);

  if (masterproc)
  {
    for (int l=0;l<ntdofs;l++)
    {
    this->stats_arr[TMASSSTAT].data(l,i) = massglobal(l);
    this->stats_arr[TMAXSTAT].data(l,i) = maxglobal(l);
    this->stats_arr[TMINSTAT].data(l,i) = minglobal(l);
  }
    for (int l=0;l<ntfctdofs;l++)
    {
    this->stats_arr[TFCTMASSSTAT].data(l,i) = massfctglobal(l);
    this->stats_arr[TFCTMAXSTAT].data(l,i) = maxfctglobal(l);
    this->stats_arr[TFCTMINSTAT].data(l,i) = minfctglobal(l);
  }
  for (int l=0;l<nQdofs;l++)
  {
    this->stats_arr[QMASSSTAT].data(l,i) = massQglobal(l);
  }
  }
}
};



// *******   VariableSet Initialization   ***********//

template <uint nprog, uint nconst, uint naux, uint ndiag> void initialize_variables(const Topology &ptopo, const Topology &dtopo,
SArray<int, nprog, 4> &prog_ndofs_arr, SArray<int, nconst, 4> &const_ndofs_arr, SArray<int, naux, 4> &aux_ndofs_arr, SArray<int, ndiag, 4> &diag_ndofs_arr,
std::array<std::string, nprog> &prog_names_arr, std::array<std::string, nconst> &const_names_arr, std::array<std::string, naux> &aux_names_arr, std::array<std::string, ndiag> &diag_names_arr,
std::array<const Topology *, nprog> &prog_topo_arr, std::array<const Topology *, nconst> &const_topo_arr, std::array<const Topology *, naux> &aux_topo_arr, std::array<const Topology *, ndiag> &diag_topo_arr)
{
  prog_topo_arr[TVAR] = &dtopo;
  prog_topo_arr[TFCTVAR] = &dtopo;
  prog_topo_arr[QVAR] = &ptopo;

  const_topo_arr[VVAR] = &ptopo;
  const_topo_arr[UVAR] = &dtopo;
  const_topo_arr[UTVAR] = &ptopo;

  aux_topo_arr[T0VAR] = &ptopo;
  aux_topo_arr[TRECONVAR] = &dtopo;
  aux_topo_arr[TEDGERECONVAR] = &dtopo;
  aux_topo_arr[TFCT0VAR] = &ptopo;
  aux_topo_arr[TFCTRECONVAR] = &dtopo;
  aux_topo_arr[TFCTEDGERECONVAR] = &dtopo;
  aux_topo_arr[PHIVAR] = &dtopo;
  aux_topo_arr[MFVAR] = &dtopo;
  aux_topo_arr[EDGEFLUXVAR] = &dtopo;
  aux_topo_arr[Q0VAR] = &dtopo;
  aux_topo_arr[QRECONVAR] = &ptopo;
  aux_topo_arr[QEDGERECONVAR] = &ptopo;
  aux_topo_arr[QFLUXVAR] = &ptopo;

  diag_topo_arr[TDIAGVAR] = &ptopo;
  diag_topo_arr[TFCTDIAGVAR] = &ptopo;
  diag_topo_arr[QDIAGVAR] = &dtopo;

  prog_names_arr[TVAR] = "T";
  prog_names_arr[TFCTVAR] = "Tfct";
  prog_names_arr[QVAR] = "Q";

  const_names_arr[VVAR] = "v";
  const_names_arr[UVAR] = "U";
  const_names_arr[UTVAR] = "UT";

  aux_names_arr[T0VAR] = "T0";
  aux_names_arr[TRECONVAR] = "Trecon";
  aux_names_arr[TEDGERECONVAR] = "Tedgerecon";
  aux_names_arr[TFCT0VAR] = "Tfct0";
  aux_names_arr[TFCTRECONVAR] = "Tfctrecon";
  aux_names_arr[TFCTEDGERECONVAR] = "Tfctedgerecon";
  aux_names_arr[PHIVAR] = "Phi";
  aux_names_arr[MFVAR] = "Mf";
  aux_names_arr[EDGEFLUXVAR] = "edgeflux";
  aux_names_arr[Q0VAR] = "q0";
  aux_names_arr[QRECONVAR] = "qrecon";
  aux_names_arr[QEDGERECONVAR] = "qedgerecon";
  aux_names_arr[QFLUXVAR] = "qflux";

  diag_names_arr[TDIAGVAR] = "T0";
  diag_names_arr[TFCTDIAGVAR] = "Tfct0";
  diag_names_arr[QDIAGVAR] = "Q0";

    prog_ndofs_arr(TVAR,ndims) = ntdofs;
    prog_ndofs_arr(TFCTVAR,ndims) = ntfctdofs;
    prog_ndofs_arr(QVAR,ndims) = nQdofs;

    const_ndofs_arr(VVAR,1) = 1;
    const_ndofs_arr(UVAR,ndims-1) = 1;
    const_ndofs_arr(UTVAR,1) = 1;

    aux_ndofs_arr(T0VAR,0) = ntdofs;
    aux_ndofs_arr(TRECONVAR,ndims-1) = ntdofs;
    aux_ndofs_arr(TEDGERECONVAR,ndims) = 2*ndims*ntdofs;
    aux_ndofs_arr(TFCT0VAR,0) = ntfctdofs;
    aux_ndofs_arr(TFCTRECONVAR,ndims-1) = ntfctdofs;
    aux_ndofs_arr(TFCTEDGERECONVAR,ndims) = 2*ndims*ntfctdofs;
    aux_ndofs_arr(PHIVAR,ndims-1) = ntfctdofs;
    aux_ndofs_arr(MFVAR,ndims) = ntfctdofs;
    aux_ndofs_arr(EDGEFLUXVAR,ndims-1) = ntfctdofs;
    aux_ndofs_arr(Q0VAR,0) = nQdofs;
    aux_ndofs_arr(QRECONVAR,ndims-1) = nQdofs;
    aux_ndofs_arr(QEDGERECONVAR,ndims) = 2*ndims*nQdofs;
    aux_ndofs_arr(QFLUXVAR,ndims-1) = nQdofs;

    diag_ndofs_arr(TDIAGVAR,0) = ntdofs;
    diag_ndofs_arr(TFCTDIAGVAR,0) = ntfctdofs;
    diag_ndofs_arr(QDIAGVAR,0) = nQdofs;
}

  // *******   Initial Conditions   ***********//

template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<nprog> &progvars, VariableSet<nconst> &constvars, 
Geometry<1, nquadx, nquady, nquadz> &pgeom, Geometry<1, nquadx, nquady, nquadz> &dgeom)
{

    for (int i=0; i<ntdofs; i++)
    {
    if (params.data_init_cond[i] == DATA_INIT::GAUSSIAN) {dgeom.set_1form_values(gaussian, progvars.fields_arr[TVAR], i);}
    if (params.data_init_cond[i] == DATA_INIT::SQUARE)   {dgeom.set_1form_values(square,   progvars.fields_arr[TVAR], i);}
    }
    for (int i=0; i<ntfctdofs; i++)
    {
    if (params.dataFCT_init_cond[i] == DATA_INIT::GAUSSIAN) {dgeom.set_1form_values(gaussian, progvars.fields_arr[TFCTVAR], i);}
    if (params.dataFCT_init_cond[i] == DATA_INIT::SQUARE)   {dgeom.set_1form_values(square,   progvars.fields_arr[TFCTVAR], i);}
    }
    if (params.wind_init_cond == WIND_INIT::UNIFORM_X ) {pgeom.set_1form_values(uniform_x_wind,     constvars.fields_arr[VVAR], 0);}
}

template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<nprog> &progvars, VariableSet<nconst> &constvars, 
Geometry<2, nquadx, nquady, nquadz> &pgeom, Geometry<2, nquadx, nquady, nquadz> &dgeom)
{
    for (int i=0; i<ntdofs; i++)
    {
    if (params.data_init_cond[i] == DATA_INIT::GAUSSIAN) {dgeom.set_2form_values(gaussian, progvars.fields_arr[TVAR], i);}
    if (params.data_init_cond[i] == DATA_INIT::VORTICES) {dgeom.set_2form_values(vortices, progvars.fields_arr[TVAR], i);}
    if (params.data_init_cond[i] == DATA_INIT::SQUARE)   {dgeom.set_2form_values(square,   progvars.fields_arr[TVAR], i);}
    }
    for (int i=0; i<ntfctdofs; i++)
    {
    if (params.dataFCT_init_cond[i] == DATA_INIT::GAUSSIAN) {dgeom.set_2form_values(gaussian, progvars.fields_arr[TFCTVAR], i);}
    if (params.dataFCT_init_cond[i] == DATA_INIT::VORTICES) {dgeom.set_2form_values(vortices, progvars.fields_arr[TFCTVAR], i);}
    if (params.dataFCT_init_cond[i] == DATA_INIT::SQUARE)   {dgeom.set_2form_values(square,   progvars.fields_arr[TFCTVAR], i);}
    }
    for (int i=0; i<nQdofs; i++)
    {
    if (params.dataQ_init_cond[i] == DATA_INIT::GAUSSIAN) {pgeom.set_2form_values(gaussian, progvars.fields_arr[QVAR], i);}
    if (params.dataQ_init_cond[i] == DATA_INIT::VORTICES) {pgeom.set_2form_values(vortices, progvars.fields_arr[QVAR], i);}
    if (params.dataQ_init_cond[i] == DATA_INIT::SQUARE)   {pgeom.set_2form_values(square,   progvars.fields_arr[QVAR], i);}
    }
    if (params.wind_init_cond == WIND_INIT::UNIFORM_X    ) {pgeom.set_1form_values(uniform_x_wind,     constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);}
    if (params.wind_init_cond == WIND_INIT::UNIFORM_Y    ) {pgeom.set_1form_values(uniform_y_wind,     constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);}
    if (params.wind_init_cond == WIND_INIT::UNIFORM_XY   ) {pgeom.set_1form_values(uniform_xy_wind,    constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);}
    if (params.wind_init_cond == WIND_INIT::DEFORMATIONAL) {pgeom.set_1form_values(deformational_wind, constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);}

}


#endif
