#pragma once

#include "MultipleFields.h"
#include "pam_coupler.h" //Has DataManager and pam_const
#include "common.h"
#include "variable_sets.h"
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
    VariableSet<nprognostic> prognostic_vars;
    VariableSet<nconstant> constant_vars;
    VariableSet<ndiagnostic> diagnostic_vars;
    VariableSet<nauxiliary> auxiliary_vars;
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

  void init(PamCoupler const &coupler) {

    // Get MPI Info
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&par.nranks);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&par.myrank);
    
    //Set parameters
    std::string inFile = coupler.get_option<std::string>( "standalone_input_file" );
    std::cout << "reading parameters\n" << std::flush;
    readModelParamsFile(inFile, params, par, coupler.get_nz());
    std::cout << "read parameters\n" << std::flush;
    
// HOW DO WE HANDLE LX/LY/XC/YC;
// ALSO VERTICAL LEVELS STUFF?

    // Initialize the grid
    std::cout << "start init topo/geom\n" << std::flush;
    primal_topology.initialize(par,true);
    primal_geometry.initialize(primal_topology, params);
    dual_topology.initialize(par,false);
    dual_geometry.initialize(dual_topology, params);
    std::cout << "finish init topo/geom\n" << std::flush;

    // Allocate the variables
    std::cout << "start init variable sets and exchange sets\n" << std::flush;
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
    std::cout << "finish init variable sets\n" << std::flush;    

    // Initialize the statistics
    stats.initialize(params, par, primal_topology, dual_topology, primal_geometry, dual_geometry);
    
    // Create the outputter
    std::cout << "start io init\n" << std::flush;
    io.initialize(params.outputName, primal_topology, dual_topology, par, prognostic_vars, constant_vars, diagnostic_vars, stats);
    std::cout << "end io init\n" << std::flush;
    // 
    // set the initial conditions
    std::cout << "start set initial conditions\n" << std::flush;
    set_initial_conditions(params, prognostic_vars, constant_vars, primal_geometry, dual_geometry);
    // Do a boundary exchange
    prog_exchange.exchange_variable_set(prognostic_vars);
    const_exchange.exchange_variable_set(constant_vars);
    std::cout << "end set initial conditions\n" << std::flush;
    // 
    // // Initialize the time stepper, tendencies, diagnostics; compute initial stats
    std::cout << "start ts init\n" << std::flush;
    tendencies.initialize(params, primal_topology, dual_topology, primal_geometry, dual_geometry, aux_exchange, const_exchange);
    tendencies.compute_constants(constant_vars, prognostic_vars);
    diagnostics.initialize(primal_topology, dual_topology, primal_geometry, dual_geometry, diag_exchange);
    tint.initialize(params, tendencies, prognostic_vars, constant_vars, auxiliary_vars, prog_exchange);
    std::cout << "end ts init\n" << std::flush;
    stats.compute(prognostic_vars, constant_vars, 0);

    // Output the initial model state
    std::cout << "start initial output\n" << std::flush;
    diagnostics.compute_diag(constant_vars, prognostic_vars, diagnostic_vars);
    io.outputInit(etime);
    std::cout << "end initial output\n" << std::flush;
    
    prevstep = 1;
  };
    
    
  //void convert_dynamics_to_coupler_state( PamCoupler &coupler , real5d state , real5d tracers ) const {};
  //void convert_coupler_state_to_dynamics( PamCoupler const &coupler , real5d const &state , real5d const &tracers ) {};
    
  // Given the model data and CFL value, compute the maximum stable time step
  real compute_time_step(PamCoupler const &coupler, real cfl_in = -1) {return 0._fp;};
        
  void timeStep( PamCoupler &coupler , real dtphys ) {
        
    // Time stepping loop
    std::cout << "start timestepping loop\n" << std::flush;
    for (uint nstep = 0; nstep<params.crm_per_phys; nstep++) {
       yakl::fence();
       tint.stepForward(params.dtcrm);
       yakl::fence();
       
       etime += params.dtcrm;
      if ((nstep+prevstep)%params.Nout == 0) {
        std::cout << "step " << (nstep+prevstep) << " time " << etime << "\n";
        diagnostics.compute_diag(constant_vars, prognostic_vars, diagnostic_vars);
        io.output(etime);
      }
      
       if (nstep%params.Nstat == 0)
       {
         stats.compute(prognostic_vars, constant_vars, (nstep+prevstep)/params.Nstat);
       }
       
    }
    prevstep += params.crm_per_phys;

    std::cout << "end timestepping loop\n" << std::flush;
  };
        
  void finalize(PamCoupler &coupler) { 
    io.outputStats(stats);
  }
  
  const char * dycore_name() const { return "SPAM++"; }
  
  void set_domain_sizes(std::string initData, real &xlen, real &ylen)
  {
    //std::cout << "setting domain sizes\n";
    set_domain_sizes_ic(params, initData);
    xlen = params.xlen;
    ylen = params.ylen;
    //std::cout << params.xlen << " " << params.ylen << "\n";
  }
};


