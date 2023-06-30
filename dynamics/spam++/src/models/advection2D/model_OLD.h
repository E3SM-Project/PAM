#ifndef _MODEL_H_
#define _MODEL_H_


#include "common.h"
#include "topology.h"
#include "variable_sets.h"
#include "ext_deriv.h"
#include "hodge_star.h"
#include "fct.h"
#include "recon.h"
#include "geometry.h"
#include "params.h"
#include "string.h"



// Number of variables
uint constexpr nprognostic = 1;
uint constexpr nconstant = 2;
uint constexpr nauxiliary = 6;
uint constexpr ndiagnostic = 1;
uint constexpr nstats = 3;

#define QVAR 0

#define VVAR 0
#define UVAR 1

#define Q0VAR 0
#define QRECONVAR 1
#define QEDGERECONVAR 2
#define PHIVAR 3
#define EDGEFLUXVAR 4
#define MFVAR 5

#define QDIAGVAR 0

#define MSTAT 0
#define MINSTAT 1
#define MAXSTAT 2


// Initial conditions related variables and functions
#define gaussian_1d(x)     (1. * exp(-100. * pow(x-0.5,2.)))
#define gaussian_2d(x,y)   (1. * exp(-100. * pow(x-0.5,2.)) * exp(-100. * pow(y-0.5,2.)))
#define gaussian_3d(x,y,z) (1. * exp(-100. * pow(x-0.5,2.)) * exp(-100. * pow(y-0.5,2.)) * exp(-100. * pow(z-0.5,2.)))
real YAKL_INLINE gaussian(real x)                 { return gaussian_1d(x); }
real YAKL_INLINE gaussian(real x, real y)         { return gaussian_2d(x,y); }
real YAKL_INLINE gaussian(real x, real y, real z) { return gaussian_3d(x,y,z); }

#define vortex1_1d(x)     (1. *  exp(-100. * pow(x-0.75,2.)))
#define vortex1_2d(x,y)   (1. *  exp(-100. * pow(x-0.75,2.)) * exp(-100. * pow(y-0.75,2.)))
#define vortex1_3d(x,y,z) (1. *  exp(-100. * pow(x-0.75,2.)) * exp(-100. * pow(y-0.75,2.)) * exp(-100. * pow(z-0.75,2.)))
#define vortex2_1d(x)     (0.5 * exp(-50.  * pow(x-0.25,2.)))
#define vortex2_2d(x,y)   (0.5 * exp(-50.  * pow(x-0.25,2.)) * exp(-75.  * pow(y-0.25,2.)))
#define vortex2_3d(x,y,z) (0.5 * exp(-50.  * pow(x-0.25,2.)) * exp(-75.  * pow(y-0.25,2.)) * exp(-100. * pow(z-0.25,2.)))
real YAKL_INLINE vortices(real x)                 { return vortex1_1d(x)     + vortex2_1d(x); }
real YAKL_INLINE vortices(real x, real y)         { return vortex1_2d(x,y)   + vortex2_2d(x,y); }
real YAKL_INLINE vortices(real x, real y, real z) { return vortex1_3d(x,y,z) + vortex2_3d(x,y,z); }

real YAKL_INLINE square(real x)                 {return (x > 0.35 && x < 0.65                                                ) ? 1. : 0.;}
real YAKL_INLINE square(real x, real y)         {return (x > 0.35 && x < 0.65 && y > 0.35 && y < 0.65                        ) ? 1. : 0.;}
real YAKL_INLINE square(real x, real y, real z) {return (x > 0.35 && x < 0.65 && y > 0.35 && y < 0.65 && z > 0.35 && z < 0.65) ? 1. : 0.;}





#define C_UNIFORM_WIND 1.

real YAKL_INLINE uniform_x_wind(real x) {
  return C_UNIFORM_WIND;
}
vec<2> YAKL_INLINE uniform_x_wind(real x, real y) {
  vec<2> vvec;
  vvec.u = C_UNIFORM_WIND;
  return vvec;
}
vec<3> YAKL_INLINE uniform_x_wind(real x, real y, real z) {
  vec<3> vvec;
  vvec.u = C_UNIFORM_WIND;
  return vvec;
}



vec<2> YAKL_INLINE uniform_y_wind(real x, real y) {
  vec<2> vvec;
  vvec.v = C_UNIFORM_WIND;
  return vvec;
}
vec<3> YAKL_INLINE uniform_y_wind(real x, real y, real z) {
  vec<3> vvec;
  vvec.v = C_UNIFORM_WIND;
  return vvec;
}

vec<2> YAKL_INLINE uniform_xy_wind(real x, real y) {
  vec<2> vvec;
  vvec.u = sqrt(C_UNIFORM_WIND/2.);
  vvec.v = sqrt(C_UNIFORM_WIND/2.);
  return vvec;
}


vec<3> YAKL_INLINE uniform_z_wind(real x, real y, real z) {
  vec<3> vvec;
  vvec.w = C_UNIFORM_WIND;
  return vvec;
}

// FIX THIS
vec<2> YAKL_INLINE deformational_wind(real x, real y) {
  vec<2> vvec;
  vvec.v = C_UNIFORM_WIND;
  return vvec;
}
vec<3> YAKL_INLINE deformational_wind(real x, real y, real z) {
  vec<3> vvec;
  vvec.v = C_UNIFORM_WIND;
  return vvec;
}

// *******   Model Specific Parameters   ***********//


void set_model_specific_params(std::string inFile, ModelParameters &params)
{


  std::string strDataInit1 = "";
  std::string strDataInit2 = "";
  std::string strDataInit3 = "";
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
      else if ( !strcmp( "windInit"    , key.c_str() ) ) { ssVal >> strWindInit        ;}
      //else {
      //  std::cout << "Error: key " << key << " not understood in file " << inFile << "\n";
      //}
    }
  }

  if (!strcmp("",strDataInit1.c_str())) { std::cout << "Error: key " << "dataInit1" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataInit2.c_str())) { std::cout << "Error: key " << "dataInit2" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataInit3.c_str())) { std::cout << "Error: key " << "dataInit3" << " not set.\n"; exit(-1); }
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

  const Topology *topology;
  Geometry<ndims,1,1,1> *geom;

  bool is_initialized;

   Diagnostics() {
     this->is_initialized = false;
     std::cout << "CREATED DIAGNOSTICS\n";
   }

   void initialize(const Topology &topo, Geometry<ndims,1,1,1> &geom)
   {
     this->topology = &topo;
     this->geom = &geom;

     this->is_initialized = true;
   }

      void compute_diag(const VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<ndiag> &diagnostic_vars)
      {

        int is = topology->is;
        int js = topology->js;
        int ks = topology->ks;

        // Compute q0
        yakl::parallel_for("ComputeQ0", topology->n_cells, YAKL_LAMBDA (int iGlob) {
          int k, j, i;
          yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

        compute_I<nqdofs, diff_ord>(diagnostic_vars.fields_arr[Q0VAR].data, x.fields_arr[QVAR].data, *this->geom, is, js, ks, i, j, k);

        });


      }
};

// *******   Tendencies   ***********//

template <uint nprog, uint nconst, uint naux> class Tendencies {
public:

  const Topology *topology;
  ExchangeSet<naux> *aux_exchange;
  ExchangeSet<nconst> *const_exchange;
  Geometry<ndims,1,1,1> *geom;

  TransformMatrices<real> trans;
  SArray<real,reconstruction_order,2> to_gll;
  SArray<real,reconstruction_order,reconstruction_order,reconstruction_order> wenoRecon;
  SArray<real,(reconstruction_order-1)/2+2> wenoIdl;
  real wenoSigma;




  bool is_initialized;

   Tendencies() {
     this->is_initialized = false;
     std::cout << "CREATED TENDENCIES\n";
   }

  void initialize(const Topology &topo, Geometry<ndims,1,1,1> &geom, ExchangeSet<naux> &aux_exchange, ExchangeSet<nconst> &const_exchange)
  {
    this->topology = &topo;
    this->geom = &geom;
    this->aux_exchange = &aux_exchange;
    this->const_exchange = &const_exchange;

    // Setup the matrix to transform a stencil of ord cell averages into tord GLL points
    trans.coefs_to_gll_lower( to_gll );
    trans.weno_sten_to_coefs(wenoRecon);
    wenoSetIdealSigma<reconstruction_order>(wenoIdl,wenoSigma);

    this->is_initialized = true;
  }



  void compute_rhs(real dt, VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<naux> &auxiliary_vars, VariableSet<nprog> &xtend)
  {


int is = topology->is;
int js = topology->js;
int ks = topology->ks;

// Compute q0 and U
yakl::parallel_for("ComputeQ0UVAR", topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

compute_H<nqdofs, diff_ord>(const_vars.fields_arr[UVAR].data, const_vars.fields_arr[VVAR].data, *this->geom, is, js, ks, i, j, k);
compute_I<nqdofs, diff_ord>(auxiliary_vars.fields_arr[Q0VAR].data, x.fields_arr[QVAR].data, *this->geom, is, js, ks, i, j, k);

});


this->aux_exchange->exchanges_arr[Q0VAR].exchange_field(auxiliary_vars.fields_arr[Q0VAR]);
this->const_exchange->exchanges_arr[UVAR].exchange_field(const_vars.fields_arr[UVAR]);

//compute qrecon

yakl::parallel_for("ComputeEdgeRecon", topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

  compute_primal_edge_recon<nqdofs, reconstruction_type, reconstruction_order>(
    auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[Q0VAR].data, is, js, ks, i, j, k,
    wenoRecon, to_gll, wenoIdl, wenoSigma);

});

this->aux_exchange->exchanges_arr[QEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[QEDGERECONVAR]);

yakl::parallel_for("ComputeQRECON", topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

  compute_primal_recon<nqdofs, reconstruction_type>(
    auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[QEDGERECONVAR].data, const_vars.fields_arr[UVAR].data, is, js, ks, i, j, k);

});
this->aux_exchange->exchanges_arr[QRECONVAR].exchange_field(auxiliary_vars.fields_arr[QRECONVAR]);




    if (fct)
    {

      yakl::parallel_for("ComputeEdgeFlux", topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
      compute_edgefluxes<nqdofs> (auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, const_vars.fields_arr[UVAR].data, is, js, ks, i, j, k);
      });
this->aux_exchange->exchanges_arr[EDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[EDGEFLUXVAR]);


yakl::parallel_for("ComputeMf", topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
  compute_Mf<nqdofs> (auxiliary_vars.fields_arr[MFVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, dt, is, js, ks, i, j, k);
});

    this->aux_exchange->exchanges_arr[MFVAR].exchange_field(auxiliary_vars.fields_arr[MFVAR]);

    yakl::parallel_for("ComputePhi", topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
      compute_Phi<nqdofs> (auxiliary_vars.fields_arr[PHIVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[QVAR].data, is, js, ks, i, j, k);
    });


    this->aux_exchange->exchanges_arr[PHIVAR].exchange_field(auxiliary_vars.fields_arr[PHIVAR]);
  }

    yakl::parallel_for("ComputeQTend", topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

      if (fct)
      {compute_wDbar2_fct<nqdofs> (xtend.fields_arr[QVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[PHIVAR].data, const_vars.fields_arr[UVAR].data, is, js, ks, i, j, k);}
      else
      {compute_wDbar2<nqdofs> (xtend.fields_arr[QVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, const_vars.fields_arr[UVAR].data, is, js, ks, i, j, k);}

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

  void initialize(ModelParameters &params, Parallel &par, const Topology &topo, Geometry<ndims,1,1,1> &geom)
  {
    statsize = params.Nsteps/params.Nstat + 1;
    stats_arr[MSTAT].initialize("qmass", nqdofs, params, par);
    stats_arr[MINSTAT].initialize("qmin", nqdofs, params, par);
    stats_arr[MAXSTAT].initialize("qmax", nqdofs, params, par);
    masterproc = par.masterproc;
  }




  void compute( VariableSet<nprog> &progvars,  VariableSet<nconst> &constvars, int i)
  {

      SArray<real,nqdofs> masslocal, massglobal;
      SArray<real,nqdofs> maxlocal, maxglobal;
      SArray<real,nqdofs> minlocal, minglobal;

      for (int l=0;l<nqdofs;l++)
      {
        masslocal(l) = progvars.fields_arr[QVAR].sum(l);
        maxlocal(l) = progvars.fields_arr[QVAR].max(l);
        minlocal(l) = progvars.fields_arr[QVAR].min(l);
      }

    //MPI sum/min/max
    this->ierr = MPI_Ireduce( &masslocal, &massglobal, nqdofs, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[MSTAT]);
    this->ierr = MPI_Ireduce( &maxlocal, &maxglobal, nqdofs, REAL_MPI, MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[MAXSTAT]);
    this->ierr = MPI_Ireduce( &minlocal, &minglobal, nqdofs, REAL_MPI, MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[MINSTAT]);
    this->ierr = MPI_Waitall(nstats, this->Req, this->Status);

  if (masterproc)
  {
    for (int l=0;l<nqdofs;l++)
    {
    this->stats_arr[MSTAT].data(l,i) = massglobal(l);
    this->stats_arr[MAXSTAT].data(l,i) = maxglobal(l);
    this->stats_arr[MINSTAT].data(l,i) = minglobal(l);
  }
  }
  }
};



// *******   VariableSet Initialization   ***********//

template <uint nprog, uint nconst, uint naux, uint ndiag> void initialize_variables(const Topology &topo,
SArray<int, nprog, 4> &prog_ndofs_arr, SArray<int, nconst, 4> &const_ndofs_arr, SArray<int, naux, 4> &aux_ndofs_arr, SArray<int, ndiag, 4> &diag_ndofs_arr,
std::array<std::string, nprog> &prog_names_arr, std::array<std::string, nconst> &const_names_arr, std::array<std::string, naux> &aux_names_arr, std::array<std::string, ndiag> &diag_names_arr,
std::array<const Topology *, nprog> &prog_topo_arr, std::array<const Topology *, nconst> &const_topo_arr, std::array<const Topology *, naux> &aux_topo_arr, std::array<const Topology *, ndiag> &diag_topo_arr)
{
  prog_topo_arr[QVAR] = &topo;

  const_topo_arr[UVAR] = &topo;
  const_topo_arr[VVAR] = &topo;

  aux_topo_arr[Q0VAR] = &topo;
  aux_topo_arr[QRECONVAR] = &topo;
  aux_topo_arr[QEDGERECONVAR] = &topo;
  aux_topo_arr[PHIVAR] = &topo;
  aux_topo_arr[MFVAR] = &topo;
  aux_topo_arr[EDGEFLUXVAR] = &topo;

  diag_topo_arr[QDIAGVAR] = &topo;


  prog_names_arr[QVAR] = "q";

  const_names_arr[UVAR] = "u";
  const_names_arr[VVAR] = "v";

  aux_names_arr[Q0VAR] = "q0";
  aux_names_arr[QRECONVAR] = "qrecon";
  aux_names_arr[QEDGERECONVAR] = "qedgerecon";
  aux_names_arr[PHIVAR] = "Phi";
  aux_names_arr[MFVAR] = "Mf";
  aux_names_arr[EDGEFLUXVAR] = "edgeflux";

  diag_names_arr[QDIAGVAR] = "q0";


    prog_ndofs_arr(QVAR,ndims) = nqdofs;

    const_ndofs_arr(VVAR,1) = 1;
    const_ndofs_arr(UVAR,ndims-1) = 1;

    aux_ndofs_arr(Q0VAR,ndims) = nqdofs;
    aux_ndofs_arr(QRECONVAR,ndims-1) = nqdofs;
    aux_ndofs_arr(QEDGERECONVAR,ndims) = 2*ndims*nqdofs;
    aux_ndofs_arr(PHIVAR,ndims-1) = nqdofs;
    aux_ndofs_arr(MFVAR,ndims) = nqdofs;
    aux_ndofs_arr(EDGEFLUXVAR,ndims-1) = nqdofs;

    diag_ndofs_arr(QDIAGVAR,ndims) = nqdofs;


}

  // *******   Initial Conditions   ***********//

template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<nprog> &progvars, VariableSet<nconst> &constvars, Geometry<1, nquadx, nquady, nquadz> &geom)
{

    for (int i=0; i<nqdofs; i++)
    {
    if (params.data_init_cond[i] == DATA_INIT::GAUSSIAN) {geom.set_primal_1form_values(gaussian, progvars.fields_arr[QVAR], i);}
    if (params.data_init_cond[i] == DATA_INIT::SQUARE)   {geom.set_primal_1form_values(square,   progvars.fields_arr[QVAR], i);}
    }
    if (params.wind_init_cond == WIND_INIT::UNIFORM_X ) {geom.set_dual_1form_values(uniform_x_wind,     constvars.fields_arr[VVAR], 0);}
}

template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<nprog> &progvars, VariableSet<nconst> &constvars, Geometry<2, nquadx, nquady, nquadz> &geom)
{
    for (int i=0; i<nqdofs; i++)
    {
    if (params.data_init_cond[i] == DATA_INIT::GAUSSIAN) {geom.set_primal_2form_values(gaussian, progvars.fields_arr[QVAR], i);}
    if (params.data_init_cond[i] == DATA_INIT::VORTICES) {geom.set_primal_2form_values(vortices, progvars.fields_arr[QVAR], i);}
    if (params.data_init_cond[i] == DATA_INIT::SQUARE)   {geom.set_primal_2form_values(square,   progvars.fields_arr[QVAR], i);}
    }
    if (params.wind_init_cond == WIND_INIT::UNIFORM_X    ) {geom.set_dual_1form_values(uniform_x_wind,     constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);}
    if (params.wind_init_cond == WIND_INIT::UNIFORM_Y    ) {geom.set_dual_1form_values(uniform_y_wind,     constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);}
    if (params.wind_init_cond == WIND_INIT::UNIFORM_XY   ) {geom.set_dual_1form_values(uniform_xy_wind,    constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);}
    if (params.wind_init_cond == WIND_INIT::DEFORMATIONAL) {geom.set_dual_1form_values(deformational_wind, constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);}

}


template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<nprog> &progvars, VariableSet<nconst> &constvars, Geometry<3, nquadx, nquady, nquadz> &geom)
{
    for (int i=0; i<nqdofs; i++)
    {
    if (params.data_init_cond[i] == DATA_INIT::GAUSSIAN) {geom.set_primal_3form_values(gaussian, progvars.fields_arr[QVAR], i);}
    if (params.data_init_cond[i] == DATA_INIT::SQUARE)   {geom.set_primal_3form_values(square,   progvars.fields_arr[QVAR], i);}
    }
    if (params.wind_init_cond == WIND_INIT::UNIFORM_X    ) {geom.set_dual_1form_values(uniform_x_wind,     constvars.fields_arr[VVAR], 0);}
    if (params.wind_init_cond == WIND_INIT::UNIFORM_Y    ) {geom.set_dual_1form_values(uniform_y_wind,     constvars.fields_arr[VVAR], 0);}
    if (params.wind_init_cond == WIND_INIT::UNIFORM_Z    ) {geom.set_dual_1form_values(uniform_z_wind,     constvars.fields_arr[VVAR], 0);}
    if (params.wind_init_cond == WIND_INIT::DEFORMATIONAL) {geom.set_dual_1form_values(deformational_wind, constvars.fields_arr[VVAR], 0);}




}

#endif