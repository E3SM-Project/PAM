#pragma once

#include "DataManager.h"
#include "MultipleFields.h"
#include "pam_coupler.h" //Has DataManager and pam_const
#include "common.h"
#include "field_sets.h"
#include "fileio.h"
#include "topology.h"
#include "geometry.h"
#include "params.h"
#include "parallel.h"
#if defined _HAMILTONIAN && defined _LAYER
#include "layermodel.h"
#elif defined _HAMILTONIAN && defined _EXTRUDED
#include "extrudedmodel.h"
#elif defined _ADVECTION && defined _LAYER 
#include "layeradvection.h"
#elif defined _ADVECTION && defined _EXTRUDED 
#include "extrudedadvection.h"
#endif
#include "SSPRK.h"
#include "RKSimple.h"

using pam::PamCoupler;

class Dycore {
  public:

#ifdef _LAYER
    UniformRectangularStraightGeometry primal_geometry;
    UniformRectangularTwistedGeometry dual_geometry;
#elif _EXTRUDED    
    UniformRectangularStraightExtrudedGeometry primal_geometry;
    UniformRectangularTwistedExtrudedGeometry dual_geometry;
#endif

    ModelStats stats;
    FieldSet<nprognostic> prognostic_vars;
    FieldSet<nconstant> constant_vars;
    FieldSet<ndiagnostic> diagnostic_vars;
    FieldSet<nauxiliary> auxiliary_vars;
    ExchangeSet<nprognostic> prog_exchange;
    ExchangeSet<nconstant> const_exchange;
    ExchangeSet<nauxiliary> aux_exchange;
    ExchangeSet<ndiagnostic> diag_exchange;
    FileIO io;
    ModelTendencies tendencies;
    ModelDiagnostics diagnostics;
    Topology primal_topology;
    Topology dual_topology;
    ModelParameters params;
    Parallel par;    
    #if _TIME_TYPE==0
          RKSimpleTimeIntegrator tint;
    #endif
    #if _TIME_TYPE==1
          SSPKKTimeIntegrator tint;
    #endif
    
    int ierr;
    real etime = 0.0;
    uint prevstep = 0;

  void init(PamCoupler &coupler) {

    serial_print("setting up dycore", par.masterproc);

    //Set parameters
    debug_print("reading parameters", par.masterproc);
    std::string inFile = coupler.get_option<std::string>( "standalone_input_file" );
    readModelParamsFile(inFile, params, par, coupler.get_nz());
    debug_print("read parameters", par.masterproc);

// HOW DO WE HANDLE VERTICAL LEVELS STUFF?
// BASICALLY A NEW GEOMETRY, I THINK?

    // Initialize the grid
    debug_print("start init topo/geom", par.masterproc);
    primal_topology.initialize(par,true);
    primal_geometry.initialize(primal_topology, params);
    dual_topology.initialize(par,false);
    dual_geometry.initialize(dual_topology, params);
    debug_print("finish init topo/geom", par.masterproc);

    // Allocate the variables
    debug_print("start init field/exchange sets", par.masterproc);
    //this gives basedof, extdof and ndofs
    SArray<int,2, nprognostic, 3> prog_dofs_arr;
    SArray<int,2, nconstant, 3> const_dofs_arr;
    SArray<int,2, ndiagnostic, 3> diag_dofs_arr;
    SArray<int,2, nauxiliary, 3> aux_dofs_arr;
    std::array<std::string, nprognostic> prog_names_arr;
    std::array<std::string, nconstant> const_names_arr;
    std::array<std::string, ndiagnostic> diag_names_arr;
    std::array<std::string, nauxiliary> aux_names_arr;
    std::array<const Topology *, nprognostic> prog_topo_arr;
    std::array<const Topology *, nconstant> const_topo_arr;
    std::array<const Topology *, ndiagnostic> diag_topo_arr;
    std::array<const Topology *, nauxiliary> aux_topo_arr;
    initialize_variables(primal_topology, dual_topology,
        prog_dofs_arr, const_dofs_arr, aux_dofs_arr, diag_dofs_arr,
        prog_names_arr, const_names_arr, aux_names_arr, diag_names_arr,
        prog_topo_arr, const_topo_arr, aux_topo_arr, diag_topo_arr);
    
    prognostic_vars.initialize("x", prog_names_arr, prog_topo_arr, prog_dofs_arr);
    constant_vars.initialize("cons", const_names_arr, const_topo_arr, const_dofs_arr);
    diagnostic_vars.initialize("diag", diag_names_arr, diag_topo_arr, diag_dofs_arr);
    auxiliary_vars.initialize("aux", aux_names_arr, aux_topo_arr, aux_dofs_arr);
    prog_exchange.initialize(prog_topo_arr, prog_dofs_arr);
    const_exchange.initialize(const_topo_arr, const_dofs_arr);
    aux_exchange.initialize(aux_topo_arr, aux_dofs_arr);
    diag_exchange.initialize(diag_topo_arr, diag_dofs_arr);
    debug_print("finish init field/exchange sets", par.masterproc);

    // Initialize the statistics
    debug_print("start stats init", par.masterproc);
    stats.initialize(params, par, primal_topology, dual_topology, primal_geometry, dual_geometry);
    debug_print("end stats init", par.masterproc);
    
    // Create the outputter
    debug_print("start io init", par.masterproc);
    io.initialize(params.outputName, primal_topology, dual_topology, par, prognostic_vars, constant_vars, diagnostic_vars, stats);
    debug_print("finish io init", par.masterproc);

    // // Initialize the tendencies and diagnostics
    debug_print("start tendencies/diagnostic init", par.masterproc);
    tendencies.initialize(coupler, params, primal_topology, dual_topology, primal_geometry, dual_geometry, aux_exchange, const_exchange);
    diagnostics.initialize(primal_topology, dual_topology, primal_geometry, dual_geometry, diag_exchange);
    debug_print("end tendencies/diagnostic init", par.masterproc);

    //EVENTUALLY THIS NEEDS TO BE MORE CLEVER IE POSSIBLY DO NOTHING BASED ON THE IC STRING?
    // set the initial conditions and compute initial stats
    debug_print("start ic setting", par.masterproc);
    set_initial_conditions(params, prognostic_vars, constant_vars, primal_geometry, dual_geometry);
    prog_exchange.exchange_variable_set(prognostic_vars);
    const_exchange.exchange_variable_set(constant_vars);
    tendencies.compute_constants(constant_vars, prognostic_vars);
    const_exchange.exchange_variable_set(constant_vars);
    stats.compute(prognostic_vars, constant_vars, 0);
    debug_print("end ic setting", par.masterproc);
        
    // // Initialize the time stepper 
    debug_print("start ts init", par.masterproc);
    tint.initialize(params, tendencies, prognostic_vars, constant_vars, auxiliary_vars, prog_exchange);
    debug_print("end ts init", par.masterproc);

    // convert dynamics state to Coupler state
    tendencies.convert_dynamics_to_coupler_state(coupler, prognostic_vars, constant_vars);
    
    // Output the initial model state
    debug_print("start initial io", par.masterproc);
    diagnostics.compute_diag(constant_vars, prognostic_vars, diagnostic_vars);
    io.outputInit(etime);
    debug_print("end initial io", par.masterproc);
    
    prevstep = 1;
  };
    
    
  // Given the model data and CFL value, compute the maximum stable time step
  real compute_time_step(PamCoupler const &coupler, real cfl_in = -1) {return 0._fp;};
        
  void timeStep( PamCoupler &coupler , real dtphys ) {

    serial_print("taking a dycore dtphys step", par.masterproc);

    // convert Coupler state to dynamics state
    tendencies.convert_coupler_to_dynamics_state(coupler, prognostic_vars, constant_vars);
    
    // Time stepping loop
    debug_print("start time stepping loop", par.masterproc);
    for (uint nstep = 0; nstep<params.crm_per_phys; nstep++) {
       yakl::fence();
       tint.stepForward(params.dtcrm);
       yakl::fence();
       
       etime += params.dtcrm;
      if ((nstep+prevstep)%params.Nout == 0) {
        serial_print("dycore step " + std::to_string((nstep+prevstep)) + " time " + std::to_string(etime), par.masterproc);
        diagnostics.compute_diag(constant_vars, prognostic_vars, diagnostic_vars);
        io.output(etime);
        io.outputStats(stats);
      }
      
       if (nstep%params.Nstat == 0)
       {
         stats.compute(prognostic_vars, constant_vars, (nstep+prevstep)/params.Nstat);
       }
       
    }
    prevstep += params.crm_per_phys;

    // convert dynamics state to Coupler state
    tendencies.convert_dynamics_to_coupler_state(coupler, prognostic_vars, constant_vars);
    
//ADD COUPLER DYCORE FUNCTIONS HERE AS WELL

debug_print("end time stepping loop", par.masterproc);
  };



  void finalize(PamCoupler &coupler) { 
    io.outputStats(stats);
  }
  
  const char * dycore_name() const { return "SPAM++"; }
  
  void set_domain_sizes(std::string initData, real &xlen, real &ylen, real &zlen)
  {
    set_domain_sizes_ic(params, initData);
    xlen = params.xlen;
    ylen = params.ylen;
    zlen = params.zlen;
  }
  
  void partition_domain(std::string infile, int &crm_nx, int &crm_ny)
  {
    _partition_domain(infile, params, par);
    crm_nx = par.nx;
    if (ndims>=2) {crm_ny = par.ny;}
  }
};


