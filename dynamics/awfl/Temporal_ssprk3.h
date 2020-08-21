
#pragma once

#include "const.h"

int  constexpr nTimeDerivs = 1;
bool constexpr timeAvg     = false;
int  constexpr nAder       = 1;

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
template <class Spatial> class Temporal_ssprk3 {
public:
  static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

  typedef typename Spatial::StateTendArr  StateTendArr;
  typedef typename Spatial::StateArr      StateArr;
  typedef typename Spatial::TracerTendArr TracerTendArr;
  typedef typename Spatial::TracerArr     TracerArr;
  typedef typename Spatial::Location      Location;

  StateTendArr  stateTendArr;
  TracerTendArr tracerTendArr;
  StateTendArr  stateTendArrTmp;
  TracerTendArr tracerTendArrTmp;

  Spatial spaceOp;

  StateArr stateTmp1;
  StateArr stateTmp2;

  TracerArr tracerTmp1;
  TracerArr tracerTmp2;

  
  void init(std::string inFile) {
    spaceOp.init(inFile);
    stateTmp1        = spaceOp.createStateArr     ();
    stateTmp2        = spaceOp.createStateArr     ();
    stateTendArr     = spaceOp.createStateTendArr ();
    stateTendArrTmp  = spaceOp.createStateTendArr ();
    tracerTmp1       = spaceOp.createTracerArr    ();
    tracerTmp2       = spaceOp.createTracerArr    ();
    tracerTendArr    = spaceOp.createTracerTendArr();
    tracerTendArrTmp = spaceOp.createTracerTendArr();
  }


  void timeStep( StateArr &stateArr , TracerArr &tracerArr , real dt ) {
    // TODO: pass an MPI communicator to computeTendencies

    /////////////////////////////////////////////////////////////////
    // STAGE 1
    /////////////////////////////////////////////////////////////////
    memset( stateTendArr  , 0._fp );
    memset( tracerTendArr , 0._fp );
    // Accumulate tenedencies
    for (int spl = 0 ; spl < spaceOp.numSplit() ; spl++) {
      real dtloc = dt;
      spaceOp.computeTendencies( stateArr , stateTendArrTmp , tracerArr , tracerTendArrTmp , dtloc , spl );
      {
        auto &stateTendArr     = this->stateTendArr    ;
        auto &stateTendArrTmp  = this->stateTendArrTmp ;
        auto &tracerTendArr    = this->tracerTendArr   ;
        auto &tracerTendArrTmp = this->tracerTendArrTmp;
        parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
          for (int l=0; l < numState; l++) {
            stateTendArr (l,k,j,i) += stateTendArrTmp (l,k,j,i);
            tracerTendArr(l,k,j,i) += tracerTendArrTmp(l,k,j,i);
          }
        });
      }
    }
    // Apply accumulated tendencies
    {
      auto &stateTendArr    = this->stateTendArr   ;
      auto &stateTmp1      = this->stateTmp1     ;
      auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
        real &state0 = spaceOp.getState    (stateArr    ,loc  ,spl);
        real &state1 = spaceOp.getState    (stateTmp1   ,loc  ,spl);
        real &tend   = spaceOp.getStateTend(stateTendArr,loc,0,spl);
        state1 = state0 + dtloc * tend;
      };
      spaceOp.applyStateTendencies( applySingleTendency , spl );
    }
    {
      auto &tracerTendArr  = this->tracerTendArr ;
      auto &tracerTmp1     = this->tracerTmp1    ;
      auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
        real &tracer0 = spaceOp.getTracer    (tracerArr    ,loc  ,spl);
        real &tracer1 = spaceOp.getTracer    (tracerTmp1   ,loc  ,spl);
        real &tend    = spaceOp.getTracerTend(tracerTendArr,loc,0,spl);
        tracer1 = tracer0 + dtloc * tend;
      };
      spaceOp.applyTracerTendencies( applySingleTendency , spl );
    }

    /////////////////////////////////////////////////////////////////
    // STAGE 2
    /////////////////////////////////////////////////////////////////
    memset( stateTendArr  , 0._fp );
    memset( tracerTendArr , 0._fp );
    // Accumulate tenedencies
    for (int spl = 0 ; spl < spaceOp.numSplit() ; spl++) {
      real dtloc = dt;
      spaceOp.computeTendencies( stateTmp1 , stateTendArrTmp , tracerTmp1 , tracerTendArrTmp , dtloc , spl );
      {
        auto &stateTendArr     = this->stateTendArr    ;
        auto &stateTendArrTmp  = this->stateTendArrTmp ;
        auto &tracerTendArr    = this->tracerTendArr   ;
        auto &tracerTendArrTmp = this->tracerTendArrTmp;
        parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
          for (int l=0; l < numState; l++) {
            stateTendArr (l,k,j,i) += stateTendArrTmp (l,k,j,i);
            tracerTendArr(l,k,j,i) += tracerTendArrTmp(l,k,j,i);
          }
        });
      }
    }
    // Apply accumulated tendencies
    {
      auto &stateTendArr = this->stateTendArr;
      auto &stateTmp1    = this->stateTmp1   ;
      auto &stateTmp2    = this->stateTmp2   ;
      auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
        real &state0 = spaceOp.getState    (stateArr    ,loc  ,spl);
        real &state1 = spaceOp.getState    (stateTmp1   ,loc  ,spl);
        real &state2 = spaceOp.getState    (stateTmp2   ,loc  ,spl);
        real &tend   = spaceOp.getStateTend(stateTendArr,loc,0,spl);
        state2 = 3.*state0/4. + state1/4. + dtloc*tend/4.;
      };
      spaceOp.applyStateTendencies( applySingleTendency , spl );
    }
    {
      auto &tracerTendArr = this->tracerTendArr;
      auto &tracerTmp1    = this->tracerTmp1   ;
      auto &tracerTmp2    = this->tracerTmp2   ;
      auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
        real &tracer0 = spaceOp.getTracer    (tracerArr    ,loc  ,spl);
        real &tracer1 = spaceOp.getTracer    (tracerTmp1   ,loc  ,spl);
        real &tracer2 = spaceOp.getTracer    (tracerTmp2   ,loc  ,spl);
        real &tend    = spaceOp.getTracerTend(tracerTendArr,loc,0,spl);
        tracer2 = 3.*tracer0/4. + tracer1/4. + dtloc*tend/4.;
      };
      spaceOp.applyTracerTendencies( applySingleTendency , spl );
    }

    /////////////////////////////////////////////////////////////////
    // STAGE 3
    /////////////////////////////////////////////////////////////////
    memset( stateTendArr  , 0._fp );
    memset( tracerTendArr , 0._fp );
    // Accumulate tenedencies
    for (int spl = 0 ; spl < spaceOp.numSplit() ; spl++) {
      real dtloc = dt;
      spaceOp.computeTendencies( stateTmp2 , stateTendArrTmp , tracerTmp2 , tracerTendArrTmp , dtloc , spl );
      {
        auto &stateTendArr     = this->stateTendArr    ;
        auto &stateTendArrTmp  = this->stateTendArrTmp ;
        auto &tracerTendArr    = this->tracerTendArr   ;
        auto &tracerTendArrTmp = this->tracerTendArrTmp;
        parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
          for (int l=0; l < numState; l++) {
            stateTendArr (l,k,j,i) += stateTendArrTmp (l,k,j,i);
            tracerTendArr(l,k,j,i) += tracerTendArrTmp(l,k,j,i);
          }
        });
      }
    }
    // Apply acculumated tendencies
    {
      auto &stateTendArr = this->stateTendArr;
      auto &stateTmp1    = this->stateTmp1   ;
      auto &stateTmp2    = this->stateTmp2   ;
      auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
        real &state0 = spaceOp.getState    (stateArr    ,loc  ,spl);
        real &state1 = spaceOp.getState    (stateTmp1   ,loc  ,spl);
        real &state2 = spaceOp.getState    (stateTmp2   ,loc  ,spl);
        real &tend   = spaceOp.getStateTend(stateTendArr,loc,0,spl);
        state0 = state0/3. + 2.*state2/3. + 2.*dtloc*tend/3.;
      };

      spaceOp.applyStateTendencies( applySingleTendency , spl );
    }
    {
      auto &tracerTendArr = this->tracerTendArr;
      auto &tracerTmp1    = this->tracerTmp1   ;
      auto &tracerTmp2    = this->tracerTmp2   ;
      auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
        real &tracer0 = spaceOp.getTracer    (tracerArr    ,loc  ,spl);
        real &tracer1 = spaceOp.getTracer    (tracerTmp1   ,loc  ,spl);
        real &tracer2 = spaceOp.getTracer    (tracerTmp2   ,loc  ,spl);
        real &tend    = spaceOp.getTracerTend(tracerTendArr,loc,0,spl);
        tracer0 = tracer0/3. + 2.*tracer2/3. + 2.*dtloc*tend/3.;
      };

      spaceOp.applyTracerTendencies( applySingleTendency , spl );
    }
  }


  void finalize(StateArr &state, TracerArr &tracers) {
    spaceOp.finalize( state , tracers );
  }


  const char * getTemporalName() { return "ADER-DT"; }

};

