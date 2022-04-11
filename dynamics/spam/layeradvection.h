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

 void compute_diag(const VariableSet<nconstant> &const_vars, VariableSet<nprognostic> &x, VariableSet<ndiagnostic> &diagnostic_vars)
 {

}

};

// *******   Tendencies   ***********//

class ModelTendencies: public Tendencies {
public:


     void compute_constants(VariableSet<nconstant> &const_vars, VariableSet<nprognostic> &x)
     {}

     void YAKL_INLINE compute_rhs(real dt, VariableSet<nconstant> &const_vars, VariableSet<nprognostic> &x, VariableSet<nauxiliary> &auxiliary_vars, VariableSet<nprognostic> &xtend)
     {
};



// *******   Statistics   ***********//

class ModelStats: public Stats {
public:
  
   real3d trimmed_density;
   
  void initialize(Parameters &params, Parallel &par, const Topology &primal_topo, const Topology &dual_topo, Geometry &primal_geom, Geometry &dual_geom)
{
  Stats::initialize(params, par, primal_topo, dual_topo, primal_geom, dual_geom);
  this->stats_arr[DENSSTAT].initialize("mass", ndensity, this->statsize, this->nens, this->masterproc);
  this->stats_arr[DENSMAXSTAT].initialize("densmax", ndensity, this->statsize, this->nens, this->masterproc);
  this->stats_arr[DENSMINSTAT].initialize("densmin", ndensity, this->statsize, this->nens, this->masterproc);

//BROKEN
  stats_arr[TMASSSTAT].initialize("Tmass", ntdofs, params, par);
  stats_arr[TMINSTAT].initialize("Tmin", ntdofs, params, par);
  stats_arr[TMAXSTAT].initialize("Tmax", ntdofs, params, par);
  stats_arr[TFCTMASSSTAT].initialize("Tfctmass", ntfctdofs, params, par);
  stats_arr[TFCTMINSTAT].initialize("Tfctmin", ntfctdofs, params, par);
  stats_arr[TFCTMAXSTAT].initialize("Tfctmax", ntfctdofs, params, par);
  stats_arr[QMASSSTAT].initialize("Qmass", nQdofs, params, par);
  
  this->PVarr = real3d("PV", this->primal_topology->nl, this->primal_topology->n_cells_y, this->primal_topology->n_cells_x);
  this->trimmed_density = real3d("trimmed_density", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
  
}

   void compute( VariableSet<nprognostic> &progvars,  VariableSet<nconstant> &constvars, int tind)
   {

    for (int n=0;n<nens;n++)
    {
      //ADD THIS- CURRENTLY BROKEN
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

//THIS STUFF IS BROKEN
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


//***************** Set Initial Conditions ***************************//

//BROKEN
void set_initial_conditions (Parameters &params, VariableSet<nprognostic> &progvars, VariableSet<nconstant> &constvars, 
Geometry &primal_geom, Geometry &dual_geom)
{

  if (params.initdataStr == "doublevortex")
  {
    std::cout << "IC: double vortex " << "\n";

      dual_geom.set_2form_values(double_vortex_h, progvars.fields_arr[DENSVAR], 0);
#ifdef _TSWE
      dual_geom.set_2form_values(double_vortex_S, progvars.fields_arr[DENSVAR], 1);
#endif
  primal_geom.set_1form_values(double_vortex_v, progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
  primal_geom.set_2form_values(double_vortex_coriolis, constvars.fields_arr[CORIOLISVAR], 0);

// HOW DO GENERALIZE THESE?
// WANT TO SCALE TRACER FIELDS BY ACTUAL HEIGHT FIELDS...
// SHOULD BE USABLE FOR ANY IC!
for (int i=0; i<ntracers; i++)
{
if (params.tracerdataStr[i] == "gaussian") {dual_geom.set_2form_values(double_vortex_tracer_gaussian, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
if (params.tracerdataStr[i] == "square") {dual_geom.set_2form_values(double_vortex_tracer_square_cent, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
if (params.tracerdataStr[i] == "doublesquare") {dual_geom.set_2form_values(double_vortex_tracer_square_urpll, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
}


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
  
}