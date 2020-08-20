
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

  #ifndef TIME_RK_D2_S3
    static_assert(false,"ERROR: You included the wrong Temporal_*_defines.h file")
  #endif
  static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

  int  static constexpr nstages = 3;
  real2d A0;
  real2d A1;
  real1d b0;
  real1d b1;
  TendArr tend_s1;
  TendArr tend_s2;
  TendArr tend_s3;
  StateArr stateTmp;

  
  void init(std::string inFile) {
    Spatial::init(inFile);

    stateTmp = createStateArr();
    tend_s1 = createTendArr();
    tend_s2 = createTendArr();
    tend_s3 = createTendArr();

    A0 = real2d("A0",nstages,nstages);
    A1 = real2d("A1",nstages,nstages);
    b0 = real1d("b0",nstages);
    b1 = real1d("b1",nstages);

    A0(1,0) = 0.484545051359423;
    A0(2,0) = 0.484545051359423;
    A0(2,1) = 0.484545038970590;
    b0(0)   = 0.458432412971593;
    b0(1)   = 0.395204641671001;
    b0(2)   = 0.146362945357405;

    A1(1,0) = 0.117391953398452;
    A1(2,0) = 0.117391950396978;
    A1(2,1) = 0.117391950396979;
    b1(0)   = 0.095747226697978;
    b1(1)   = 0.035459719531028;
    b1(2)   = 0.035459720437660;
  }


  void timeStep( StateArr &stateArr , real dt ) {
    //////////////////////////////////////////////////////////////////////
    // Stage 1
    //////////////////////////////////////////////////////////////////////
    for (int spl = 0 ; spl < numSplit() ; spl++) {
      computeTendencies( stateArr , tend_s1 , dt , spl );

      auto &stateTmp = this->stateTmp;
      auto &tend_s1  = this->tend_s1 ;
      auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
        real &s    = get(stateTmp,loc  ,spl);

        real &s0   = get(stateArr,loc  ,spl);

        real &t1_0 = get(tend_s1 ,loc,0,spl);
        real &t1_1 = get(tend_s1 ,loc,1,spl);

        s = s0 + A0(1,0)*dt*t1_0 + A1(1,0)*dt*dt*t1_1;
      };

      applyTendencies( applySingleTendency , spl );
    }

    //////////////////////////////////////////////////////////////////////
    // Stage 2
    //////////////////////////////////////////////////////////////////////
    for (int spl = 0 ; spl < numSplit() ; spl++) {
      computeTendencies( stateTmp , tend_s2 , dt , spl );

      auto &stateTmp = this->stateTmp;
      auto &tend_s1  = this->tend_s1 ;
      auto &tend_s2  = this->tend_s2 ;
      auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
        real &s    = get(stateTmp,loc  ,spl);

        real &s0   = get(stateArr,loc  ,spl);

        real &t1_0 = get(tend_s1 ,loc,0,spl);
        real &t1_1 = get(tend_s1 ,loc,1,spl);

        real &t2_0 = get(tend_s2 ,loc,0,spl);
        real &t2_1 = get(tend_s2 ,loc,1,spl);

        s = s0 + A0(2,0)*dt*t1_0 + A1(2,0)*dt*dt*t1_1 + 
                 A0(2,1)*dt*t2_0 + A1(2,1)*dt*dt*t2_1;
      };

      applyTendencies( applySingleTendency , spl );
    }

    //////////////////////////////////////////////////////////////////////
    // Stage 2
    //////////////////////////////////////////////////////////////////////
    for (int spl = 0 ; spl < numSplit() ; spl++) {
      computeTendencies( stateTmp , tend_s3 , dt , spl );

      auto &tend_s1  = this->tend_s1;
      auto &tend_s2  = this->tend_s2;
      auto &tend_s3  = this->tend_s3;
      auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
        real &s    = get(stateArr,loc  ,spl);

        real &t1_0 = get(tend_s1 ,loc,0,spl);
        real &t1_1 = get(tend_s1 ,loc,1,spl);

        real &t2_0 = get(tend_s2 ,loc,0,spl);
        real &t2_1 = get(tend_s2 ,loc,1,spl);

        real &t3_0 = get(tend_s3 ,loc,0,spl);
        real &t3_1 = get(tend_s3 ,loc,1,spl);

        s += b0(0)*dt   *t1_0 + b0(1)*dt   *t2_0 + b0(2)*dt   *t3_0 +
             b1(0)*dt*dt*t1_1 + b1(1)*dt*dt*t2_1 + b1(2)*dt*dt*t3_1;
      };

      applyTendencies( applySingleTendency , spl );
    }
  }


  const char * getTemporalName() { return "ADER-DT"; }

};

