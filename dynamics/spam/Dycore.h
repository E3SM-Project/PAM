#pragma once

//INCLUDE PAM CONST HERE???

#include "MultipleFields.h"
#include "DataManager.h"
#include "pam_coupler.h"

using pam::PamCoupler;

class Dycore {
  
  void convert_dynamics_to_coupler_state( PamCoupler &coupler , real5d state , real5d tracers ) const {};
  void convert_coupler_state_to_dynamics( PamCoupler const &coupler , real5d const &state , real5d const &tracers ) {};
    
  // Given the model data and CFL value, compute the maximum stable time step
  real compute_time_step(PamCoupler const &coupler, real cfl_in = -1) {return 0.;};
    
  void init(PamCoupler const &coupler) {};
    
  void init_idealized_state_and_tracers( PamCoupler &coupler ) {};
    
  void output(PamCoupler const &coupler, real etime) {};
      
  void timeStep( PamCoupler &coupler , real dtphys ) {};
        
  void finalize(PamCoupler &coupler) { }
  
  const char * dycore_name() const { return "SPAM++"; }
};


