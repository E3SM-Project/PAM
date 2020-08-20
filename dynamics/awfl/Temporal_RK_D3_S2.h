
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

  #ifndef TIME_RK_D3_S2
    static_assert(false,"ERROR: You included the wrong Temporal_*_defines.h file");
  #endif
  static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

  int  static constexpr nstages = 2;
  real2d A0;
  real2d A1;
  real2d A2;
  real1d b0;
  real1d b1;
  real1d b2;
  TendArr tend_s1;
  TendArr tend_s2;
  StateArr stateTmp;

  
  void init(std::string inFile) {
    Spatial::init(inFile);

    stateTmp = createStateArr();
    tend_s1 = createTendArr();
    tend_s2 = createTendArr();

    A0 = real2d("A0",nstages,nstages);
    A1 = real2d("A1",nstages,nstages);
    A2 = real2d("A2",nstages,nstages);
    b0 = real1d("b0",nstages);
    b1 = real1d("b1",nstages);
    b2 = real1d("b2",nstages);

    A0(1,0) = 0.527276308170006;
    b0(0)=0.527276308170006;
    b0(1)=0.472723691829994;

    A1(1,0) = 0.139010152578696;
    b1(0)= 0.138304642985474;   
    b1(1)=0.112439354001911;

    A2(1,0) = 0.024432253349948;
    b2(0)=0.021904464176284;
    b2(1)=0.019762202490383;

    // A0(1,0) = 0.5;
    // b0(0)=0.5;
    // b0(1)=0.5;
       
    // A1(1,0) = 1./8.;
    // b1(0)=1./8.;   
    // b1(1)=1./8.;
       
    // A2(1,0) = 1./48.;
    // b2(0)=1./48.;
    // b2(1)=1./48.;
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
        real &t1_2 = get(tend_s1 ,loc,2,spl);

        s = s0 + A0(1,0)*dt*t1_0 + A1(1,0)*dt*dt*t1_1 + A2(1,0)*dt*dt*dt*t1_2;
      };

      applyTendencies( applySingleTendency , spl );
    }

    //////////////////////////////////////////////////////////////////////
    // Stage 2
    //////////////////////////////////////////////////////////////////////
    for (int spl = 0 ; spl < numSplit() ; spl++) {
      computeTendencies( stateTmp , tend_s2 , dt , spl );

      auto &tend_s1  = this->tend_s1 ;
      auto &tend_s2  = this->tend_s2 ;
      auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
        real &s    = get(stateArr,loc  ,spl);

        real &t1_0 = get(tend_s1 ,loc,0,spl);
        real &t1_1 = get(tend_s1 ,loc,1,spl);
        real &t1_2 = get(tend_s1 ,loc,2,spl);

        real &t2_0 = get(tend_s2 ,loc,0,spl);
        real &t2_1 = get(tend_s2 ,loc,1,spl);
        real &t2_2 = get(tend_s2 ,loc,2,spl);

        s += b0(0)*dt      *t1_0 + b0(1)*dt      *t2_0 +
             b1(0)*dt*dt   *t1_1 + b1(1)*dt*dt   *t2_1 +
             b2(0)*dt*dt*dt*t1_2 + b2(1)*dt*dt*dt*t2_2;
      };

      applyTendencies( applySingleTendency , spl );
    }
  }


  const char * getTemporalName() { return "ADER-DT"; }

};

