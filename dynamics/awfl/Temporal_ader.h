
#pragma once

#include "const.h"

int  constexpr nTimeDerivs = 1;
bool constexpr timeAvg     = true;
int  constexpr nAder       = ngll;

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
template <class Spatial> class Temporal_ader {
public:
  static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

  typedef typename Spatial::TendArr  TendArr ;
  typedef typename Spatial::StateArr StateArr;
  typedef typename Spatial::Location Location;

  TendArr tendArr;

  Spatial spaceOp;
  
  void init(std::string inFile) {
    spaceOp.init(inFile);
    tendArr = spaceOp.createTendArr();
  }


  void timeStep( StateArr &stateArr , real dt ) {
    // TODO: pass an MPI communicator to computeTendencies

    // Loop over different items in the spatial splitting
    for (int spl = 0 ; spl < spaceOp.numSplit() ; spl++) {
      real dtloc = dt;
      spaceOp.computeTendencies( stateArr , tendArr , dtloc , spl );

      auto &tendArr = this->tendArr;
      auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
        real &state = spaceOp.get(stateArr,loc  ,spl);
        real &tend  = spaceOp.get(tendArr ,loc,0,spl);
        state += dtloc * tend;
      };

      spaceOp.applyTendencies( applySingleTendency , spl );
    }
  }


  void finalize(StateArr &state) {
    spaceOp.finalize(state);
  }


  const char * getTemporalName() { return "ADER-DT"; }

};

