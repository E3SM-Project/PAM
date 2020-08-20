
#pragma once

#include "const.h"

/* REQUIRED:
int  static constexpr nTimeDerivs = [# tendency time derivatives needed by the time stepping scheme];
bool static constexpr timeAvg     = [whether the spatial operator should return a time-averaged tendency];

static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

void init(std::string inFile):
    - Process input file with filename "inFile"
    - Allocate and initialize internal stuff

void timeStep( StateArr &state , real dt ): 
    - Perform a single time step

const char * getTemporalName():
    - Return the name and info about this temporal operator
*/
class Temporal : public Spatial {
public:

  #ifndef TIME_LSRK3
    static_assert(false,"ERROR: You included the wrong Temporal_*_defines.h file")
  #endif
  static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

  StateArr stateTmp;
  TendArr tendArr;

  
  void init(std::string inFile) {
    Spatial::init(inFile);
    stateTmp = createStateArr();
    tendArr = createTendArr();
  }


  void timeStep( StateArr &stateArr , real dt ) {
    // TODO: pass an MPI communicator to computeTendencies

    /////////////////////////////////////////////////////////////////
    // STAGE 1
    /////////////////////////////////////////////////////////////////
    // Loop over different items in the spatial splitting
    for (int spl = 0 ; spl < numSplit() ; spl++) {
      computeTendencies( stateArr , tendArr , dt , spl );

      auto &tendArr  = this->tendArr ;
      auto &stateTmp = this->stateTmp;
      auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
        real &state  = get(stateTmp,loc  ,spl);
        real &state0 = get(stateArr,loc  ,spl);
        real &tend   = get(tendArr ,loc,0,spl);
        state = state0 + dt / 3 * tend;
      };

      applyTendencies( applySingleTendency , spl );
    }

    /////////////////////////////////////////////////////////////////
    // STAGE 2
    /////////////////////////////////////////////////////////////////
    // Loop over different items in the spatial splitting
    for (int spl = 0 ; spl < numSplit() ; spl++) {
      computeTendencies( stateTmp , tendArr , dt , spl );

      auto &tendArr = this->tendArr;
      auto &stateTmp = this->stateTmp;
      auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
        real &state  = get(stateTmp,loc  ,spl);
        real &state0 = get(stateArr,loc  ,spl);
        real &tend   = get(tendArr ,loc,0,spl);
        state = state0 + dt / 2 * tend;
      };

      applyTendencies( applySingleTendency , spl );
    }

    /////////////////////////////////////////////////////////////////
    // STAGE 3
    /////////////////////////////////////////////////////////////////
    // Loop over different items in the spatial splitting
    for (int spl = 0 ; spl < numSplit() ; spl++) {
      computeTendencies( stateTmp , tendArr , dt , spl );

      auto &tendArr = this->tendArr;
      auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
        real &state = get(stateArr,loc  ,spl);
        real &tend  = get(tendArr ,loc,0,spl);
        state += dt / 1 * tend;
      };

      applyTendencies( applySingleTendency , spl );
    }
  }


  const char * getTemporalName() { return "ADER-DT"; }

};

