
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

  #ifndef TIME_SSPRK3
    static_assert(false,"ERROR: You included the wrong Temporal_*_defines.h file");
  #endif
  static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

  StateArr stateTmp1;
  StateArr stateTmp2;
  TendArr tendArr;

  
  void init(std::string inFile) {
    Spatial::init(inFile);
    stateTmp1 = createStateArr();
    stateTmp2 = createStateArr();
    tendArr   = createTendArr();
  }


  void timeStep( StateArr &stateArr , real dt ) {
    // TODO: pass an MPI communicator to computeTendencies

    for (int spl = 0 ; spl < numSplit() ; spl++) {
      /////////////////////////////////////////////////////////////////
      // STAGE 1
      /////////////////////////////////////////////////////////////////
      {
        computeTendencies( stateArr , tendArr , dt , spl );

        auto &tendArr   = this->tendArr  ;
        auto &stateTmp1 = this->stateTmp1;
        auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
          real &state0 = get(stateArr ,loc  ,spl);
          real &state1 = get(stateTmp1,loc  ,spl);
          real &tend   = get(tendArr  ,loc,0,spl);
          state1 = state0 + dt * tend;
        };

        applyTendencies( applySingleTendency , spl );
      }

      /////////////////////////////////////////////////////////////////
      // STAGE 2
      /////////////////////////////////////////////////////////////////
      {
        computeTendencies( stateTmp1 , tendArr , dt , spl );

        auto &tendArr   = this->tendArr ;
        auto &stateTmp1 = this->stateTmp1;
        auto &stateTmp2 = this->stateTmp2;
        auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
          real &state0 = get(stateArr ,loc  ,spl);
          real &state1 = get(stateTmp1,loc  ,spl);
          real &state2 = get(stateTmp2,loc  ,spl);
          real &tend   = get(tendArr  ,loc,0,spl);
          state2 = 3.*state0/4. + state1/4. + dt*tend/4.;
        };

        applyTendencies( applySingleTendency , spl );
      }

      /////////////////////////////////////////////////////////////////
      // STAGE 3
      /////////////////////////////////////////////////////////////////
      {
        computeTendencies( stateTmp2 , tendArr , dt , spl );

        auto &tendArr   = this->tendArr ;
        auto &stateTmp1 = this->stateTmp1;
        auto &stateTmp2 = this->stateTmp2;
        auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
          real &state0 = get(stateArr ,loc  ,spl);
          real &state1 = get(stateTmp1,loc  ,spl);
          real &state2 = get(stateTmp2,loc  ,spl);
          real &tend   = get(tendArr  ,loc,0,spl);
          state0 = state0/3. + 2.*state2/3. + 2.*dt*tend/3.;
        };

        applyTendencies( applySingleTendency , spl );
      }
    }
  }


  const char * getTemporalName() { return "ADER-DT"; }

};

