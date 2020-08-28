
#pragma once

#include "const.h"
#include "phys_params.h"
#include "TransformMatrices.h"
#include "WenoLimiter.h"
#include "Profiles.h"


/*********************************************
 ***** Required inside the Spatial class *****
 *********************************************
class Location;
    - Stores the indices of a single location on the grid

class StateArr; // OR
typedef [class] StateArr;
    - Declare a type for the model state

class TendArr; // OR
typedef [class] TendArr;
    - Declare a type for the model tendencies
    - It must have room for the nTimeDerivs dimension for the time integrator

StateArr createStateArr()
    - Create and return a new StateArr object

TendArr createTendArr(int nTimeDerivs)
    - Create and return a new TendArr object

YAKL_INLINE real &get(StateArr const &state , Location const &loc , int splitIndex)
    - Return the state value at the given location

YAKL_INLINE real &get(TendArr const &tend , Location const &loc , int timeDeriv , int splitIndex)
    - Return the tendency value at the given location

int numSplit()
    - Return the number of split components for this operator
    - The temporal operator will iterate through the splittings

real computeTimeStep(real cfl)
    - Return the time step in seconds from the cfl value

void init(int nTimeDerivs, bool timeAvg, std::string inFile)
    - Initialize internal data structures

void initState( StateArr &state )
    - Initialize the state

void computeTendencies( StateArr const &state , TendArr &tend , real dt , int splitIndex)
    - Compute tendency and time derivatives of the tendency if they are requested

template <class F> void applyTendencies( F const &applySingleTendency , int splitIndex )
    - Loop through the domain, and apply tendencies to the state

const char * getSpatialName()
    - Return the name and info for this spatial operator

void output(StateArr const &state, real etime)
    - Output to file
*/


template <int nTimeDerivs, bool timeAvg, int nAder>
class Spatial_euler3d_cons_expl_cart_fv_Agrid {
public:
  
  static_assert(nTimeDerivs == 1 , "ERROR: This Spatial class isn't setup to use nTimeDerivs > 1");

  int static constexpr hs = (ord-1)/2;
  int static constexpr numState = 5;

  // Stores a single index location
  struct Location {
    int l;
    int k;
    int j;
    int i;
  };

  typedef real4d StateArr;  // Spatial index
  typedef real4d TracerArr; // Spatial index

  typedef real4d StateTendArr;   // (time derivative & spatial index)
  typedef real4d TracerTendArr;  // (time derivative & spatial index)

  real5d stateLimits;
  real5d tracerLimits;
  real5d stateFluxLimits;
  real4d tracerFlux;
  // Hydrostatically balanced values for density, potential temperature, and pressure
  real1d hyDensCells;
  real1d hyPressureCells;
  real1d hyThetaCells;
  real1d hyDensThetaCells;
  real2d hyDensGLL;
  real2d hyPressureGLL;
  real2d hyThetaGLL;
  real2d hyDensThetaGLL;
  // Matrices to transform DOFs from one form to another
  SArray<real,2,ord,ngll> coefs_to_gll;
  SArray<real,2,ord,ngll> coefs_to_deriv_gll;
  SArray<real,2,ord,ngll> sten_to_gll;
  SArray<real,2,ord,ngll> sten_to_deriv_gll;
  SArray<real,3,ord,ord,ord> wenoRecon;
  SArray<real,1,hs+2> idl;
  real sigma;
  // For ADER spatial derivative computation
  SArray<real,2,ngll,ngll> derivMatrix;
  // For quadrature
  SArray<real,1,ord> gllWts_ord;
  SArray<real,1,ord> gllPts_ord;
  SArray<real,1,ngll> gllWts_ngll;
  SArray<real,1,ngll> gllPts_ngll;

  int static constexpr idR = 0;  // density perturbation
  int static constexpr idU = 1;  // u
  int static constexpr idV = 2;  // v
  int static constexpr idW = 3;  // w
  int static constexpr idT = 4;  // potential temperature perturbation

  int static constexpr BC_PERIODIC = 0;
  int static constexpr BC_WALL     = 1;

  int static constexpr DATA_SPEC_THERMAL       = 1;
  int static constexpr DATA_SPEC_THERMAL_MOIST = 2;

  bool sim2d;

  real dx;
  real dy;
  real dz;

  real dtInit;

  bool dimSwitch;

  int numTracers;
  int tracerCount;
  std::vector<std::string> tracerName;
  std::vector<std::string> tracerDesc;
  bool1d                   tracerPos;
  bool1d                   tracerAddsMass;

  // Values read from input file
  int         nx;
  int         ny;
  int         nz;
  real        xlen;
  real        ylen;
  real        zlen;
  int         bc_x;
  int         bc_y;
  int         bc_z;
  bool        weno_scalars;
  bool        weno_winds;
  std::string outFile;
  int         dataSpec;


  Spatial_euler3d_cons_expl_cart_fv_Agrid() {
    numTracers = 0;
    tracerCount = 0;
  }


  static_assert(ord%2 == 1,"ERROR: ord must be an odd integer");



  // Initialize a tracer
  int addTracer(std::string name , std::string desc , bool posDef , bool addsMass) {
    auto &tracerPos = this->tracerPos;

    int tr = tracerCount;
    tracerName.push_back(name);
    tracerDesc.push_back(desc);
    parallel_for( 1 , YAKL_LAMBDA (int i) {
      tracerPos     (tr) = posDef;
      tracerAddsMass(tr) = addsMass;
    });

    tracerCount++;
  }



  template <class F, class PHYS>
  void initTracerByLocation(std::string name , F const &initMass , TracerArr &tracers, PHYS const & physics) const {
    int tracerID = -1;
    for (int tr=0; tr < numTracers; tr++) {
      if (name == tracerName[tr]) tracerID = tr;
    }
    if (tracerID == -1) endrun("ERROR: Trying to init a tracer that has not been added");

    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      int tr = tracerID;
      tracers(tr,hs+k,hs+j,hs+i) = 0;
      for (int kk=0; kk<ord; kk++) {
        for (int jj=0; jj<ord; jj++) {
          for (int ii=0; ii<ord; ii++) {
            real zloc = (k+0.5_fp)*dz + gllPts_ord(kk)*dz;
            real yloc;
            if (sim2d) {
              yloc = ylen/2;
            } else {
              yloc = (j+0.5_fp)*dy + gllPts_ord(jj)*dy;
            }
            real xloc = (i+0.5_fp)*dx + gllPts_ord(ii)*dx;

            // Compute constant theta hydrostatic background state
            real th = 300;
            real rh = profiles::initConstTheta_density(th,zloc);

            typename PHYS::DynState dynState;
            dynState.rho       = rh;
            dynState.w         = 0;
            dynState.rho_theta = rh*th;

            real wt = gllWts_ord(kk) * gllWts_ord(jj) * gllWts_ord(ii);
            tracers(tr,hs+k,hs+j,hs+i) += initMass(xloc,yloc,zloc,xlen,ylen,zlen,dynState)*wt;
          }
        }
      }
    });
  }



  template <class PHYS>
  void adjustStateForMoisture(StateArr &state , TracerArr const &tracers , PHYS const &physics) const {
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // Add tracer density to dry density
      real rho_dry = state(idR,hs+k,hs+j,hs+i) + hyDensCells(hs+k);
      for (int tr=0; tr < numState; tr++) {
        if (tracerAddsMass(tr)) {
          state(idR,hs+k,hs+j,hs+i) += tracers(tr,hs+k,hs+j,hs+i);
        }
      }
      real rho_moist = state(idR,hs+k,hs+j,hs+i) + hyDensCells(hs+k);
      // Adjust momenta for moist density
      state(idU,hs+k,hs+j,hs+i) = state(idU,hs+k,hs+j,hs+i) / rho_dry * rho_moist;
      state(idV,hs+k,hs+j,hs+i) = state(idV,hs+k,hs+j,hs+i) / rho_dry * rho_moist;
      state(idW,hs+k,hs+j,hs+i) = state(idW,hs+k,hs+j,hs+i) / rho_dry * rho_moist;
      // Compute the dry temperature (same as the moist temperature)
      real rho_theta_dry = state(idT,hs+k,hs+j,hs+i) + hyDensThetaCells(hs+k);
      real temp = physics.tempFromRhoTheta(rho_dry,rho_theta_dry);
      // Compute moist theta from moist density, tracer densities, and temperature
      typename PHYS::MicroTracers tracersLoc;
      for (int tr=0; tr < numState; tr++) {
        tracersLoc(tr) = tracers(tr,hs+k,hs+j,hs+i);
      }
      real theta_moist = physics.thetaFromTemp(rho_moist, tracersLoc, temp);
      state(idT,hs+k,hs+j,hs+i) = rho_moist*theta_moist;
    });
  }



  StateArr createStateArr() const {
    return StateArr("stateArr",numState,nz+2*hs,ny+2*hs,nx+2*hs);
  }



  TracerArr createTracerArr() const {
    return StateArr("stateArr",numTracers,nz+2*hs,ny+2*hs,nx+2*hs);
  }



  StateTendArr createStateTendArr() const {
    return StateTendArr("stateTendArr",numState,nz,ny,nx);
  }



  TracerTendArr createTracerTendArr() const {
    return TracerTendArr("tracerTendArr",numTracers,nz,ny,nx);
  }



  YAKL_INLINE real &getState(StateArr const &state , Location const &loc , int splitIndex) {
    return state(loc.l , hs+loc.k , hs+loc.j , hs+loc.i);
  }



  YAKL_INLINE real &getTracer(TracerArr const &tracers , Location const &loc , int splitIndex) {
    return tracers(loc.l , hs+loc.k , hs+loc.j , hs+loc.i);
  }



  YAKL_INLINE real &getStateTend(StateTendArr const &tend , Location const &loc , int timeDeriv , int splitIndex) {
    return tend(loc.l , loc.k , loc.j , loc.i);
  }



  YAKL_INLINE real &getTracerTend(TracerTendArr const &tend , Location const &loc , int timeDeriv , int splitIndex) {
    return tend(loc.l , loc.k , loc.j , loc.i);
  }



  int numSplit() const {
    return 4;
  }



  real computeTimeStep(real cfl, StateArr const &state) {
    auto &hyDensCells          = this->hyDensCells         ;
    auto &hyDensThetaCells     = this->hyDensThetaCells    ;
    auto &dx                   = this->dx                  ;
    auto &dy                   = this->dy                  ;
    auto &dz                   = this->dz                  ;

    if (dtInit <= 0) {
      real maxwave = 0;
      real3d dt3d("dt3d",nz,ny,nx);
      parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        real r = state(idR,hs+k,hs+j,hs+i) + hyDensCells(hs+k);
        real u = state(idU,hs+k,hs+j,hs+i) / r;
        real v = state(idV,hs+k,hs+j,hs+i) / r;
        real w = state(idW,hs+k,hs+j,hs+i) / r;
        real t = ( state(idT,hs+k,hs+j,hs+i) + hyDensThetaCells(hs+k) ) / r;
        real p = C0*pow(r*t,GAMMA);
        real cs = sqrt(GAMMA*p/r);
        real udt = cfl * dx / max( abs(u-cs) , abs(u+cs) );
        real vdt = cfl * dy / max( abs(v-cs) , abs(v+cs) );
        real wdt = cfl * dz / max( abs(w-cs) , abs(w+cs) );
        dt3d(k,j,i) = min( min(udt,vdt) , wdt );
      });
      yakl::ParallelMin<real,memDevice> pmin(nz*nx*ny);
      dtInit = pmin(dt3d.data());
    }
    
    return dtInit;
  }



  // Initialize crap needed by recon()
  void init(std::string inFile, int numTracers) {
    this->numTracers = numTracers;
    tracerPos      = bool1d("tracerPos"     ,numTracers);
    tracerAddsMass = bool1d("tracerAddsMass",numTracers);

    dtInit = 0;
    dimSwitch = true;

    // Read the input file
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config                 ) { endrun("ERROR: Invalid YAML input file"); }
    if ( !config["nx"]           ) { endrun("ERROR: No nx in input file"); }
    if ( !config["ny"]           ) { endrun("ERROR: No ny in input file"); }
    if ( !config["nz"]           ) { endrun("ERROR: No nz in input file"); }
    if ( !config["xlen"]         ) { endrun("ERROR: No xlen in input file"); }
    if ( !config["ylen"]         ) { endrun("ERROR: No ylen in input file"); }
    if ( !config["zlen"]         ) { endrun("ERROR: No zlen in input file"); }
    if ( !config["bc_x"]         ) { endrun("ERROR: No bc_x in input file"); }
    if ( !config["bc_y"]         ) { endrun("ERROR: No bc_y in input file"); }
    if ( !config["bc_z"]         ) { endrun("ERROR: No bc_z in input file"); }
    if ( !config["weno_scalars"] ) { endrun("ERROR: No weno_scalars in input file"); }
    if ( !config["weno_winds"]   ) { endrun("ERROR: No weno_winds in input file"); }
    if ( !config["initData"]     ) { endrun("ERROR: No initData in input file"); }
    if ( !config["outFile"]      ) { endrun("ERROR: No outFile in input file"); }

    nx = config["nx"].as<int>();
    ny = config["ny"].as<int>();
    nz = config["nz"].as<int>();

    sim2d = ny == 1;

    xlen = config["xlen"].as<real>();
    ylen = config["ylen"].as<real>();
    zlen = config["zlen"].as<real>();

    weno_scalars = config["weno_scalars"].as<bool>();
    weno_winds   = config["weno_winds"].as<bool>();

    std::string dataStr = config["initData"].as<std::string>();
    if        (dataStr == "thermal") {
      dataSpec = DATA_SPEC_THERMAL;
    } else {
      endrun("ERROR: Invalid dataSpec");
    }

    std::string bc_x_str = config["bc_x"].as<std::string>();
    if        (bc_x_str == "periodic" ) {
      bc_x = BC_PERIODIC;
    } else if (bc_x_str == "wall"     ) {
      bc_x = BC_WALL;
    } else {
      endrun("Invalid bc_x");
    }

    std::string bc_y_str = config["bc_y"].as<std::string>();
    if        (bc_y_str == "periodic" ) {
      bc_y = BC_PERIODIC;
    } else if (bc_y_str == "wall"     ) {
      bc_y = BC_WALL;
    } else {
      endrun("Invalid bc_y");
    }

    std::string bc_z_str = config["bc_z"].as<std::string>();
    if        (bc_z_str == "periodic" ) {
      bc_z = BC_PERIODIC;
    } else if (bc_z_str == "wall"     ) {
      bc_z = BC_WALL;
    } else {
      endrun("Invalid bc_z");
    }

    outFile = config["outFile"].as<std::string>();

    dx = xlen/nx;
    dy = ylen/ny;
    dz = zlen/nz;

    TransformMatrices::weno_sten_to_coefs(this->wenoRecon);

    // Store to_gll and wenoRecon
    {
      SArray<real,2,ord,ord>  g2c;
      SArray<real,2,ord,ord>  s2c;
      SArray<real,2,ord,ngll> c2g_lower;
      SArray<real,2,ord,ord>  c2g;
      SArray<real,2,ord,ord>  c2d;

      TransformMatrices::gll_to_coefs      (g2c      );
      TransformMatrices::sten_to_coefs     (s2c      );
      TransformMatrices::coefs_to_gll_lower(c2g_lower);
      TransformMatrices::coefs_to_gll      (c2g      );
      TransformMatrices::coefs_to_deriv    (c2d      );

      this->coefs_to_gll       = c2g_lower;
      this->coefs_to_deriv_gll = c2g_lower * c2d;
      this->sten_to_gll        = c2g_lower       * s2c;
      this->sten_to_deriv_gll  = c2g_lower * c2d * s2c;

    }
    // Store ader derivMatrix
    {
      SArray<real,2,ngll,ngll> g2c;
      SArray<real,2,ngll,ngll> c2d;
      SArray<real,2,ngll,ngll> c2g;

      TransformMatrices::gll_to_coefs  (g2c);
      TransformMatrices::coefs_to_deriv(c2d);
      TransformMatrices::coefs_to_gll  (c2g);

      this->derivMatrix = c2g * c2d * g2c;
    }
    TransformMatrices::get_gll_points (this->gllPts_ord);
    TransformMatrices::get_gll_weights(this->gllWts_ord);
    TransformMatrices::get_gll_points (this->gllPts_ngll);
    TransformMatrices::get_gll_weights(this->gllWts_ngll);

    weno::wenoSetIdealSigma(this->idl,this->sigma);

    stateLimits     = real5d("stateLimits"    ,numState  ,2,nz+1,ny+1,nx+1);
    tracerLimits    = real5d("tracerLimits"   ,numTracers,2,nz+1,ny+1,nx+1);
    stateFluxLimits = real5d("stateFluxLimits",numState  ,2,nz+1,ny+1,nx+1);
    tracerFlux      = real4d("tracerFlux"     ,numTracers  ,nz+1,ny+1,nx+1);

    hyDensCells          = real1d("hyDensCells       ",nz+2*hs);
    hyPressureCells      = real1d("hyPressureCells   ",nz+2*hs);
    hyThetaCells         = real1d("hyThetaCells      ",nz+2*hs);
    hyDensThetaCells     = real1d("hyDensThetaCells  ",nz+2*hs);
    hyDensGLL            = real2d("hyDensGLL         ",nz,ngll);
    hyPressureGLL        = real2d("hyPressureGLL     ",nz,ngll);
    hyThetaGLL           = real2d("hyThetaGLL        ",nz,ngll);
    hyDensThetaGLL       = real2d("hyDensThetaGLL    ",nz,ngll);
  }



  // Initialize the state
  void initState( StateArr &state ) {
    auto nx                 = this->nx               ;
    auto ny                 = this->ny               ;
    auto nz                 = this->nz               ;
    auto dx                 = this->dx               ;
    auto dy                 = this->dy               ;
    auto dz                 = this->dz               ;
    auto gllPts_ord         = this->gllPts_ord       ;
    auto gllWts_ord         = this->gllWts_ord       ;
    auto gllPts_ngll        = this->gllPts_ngll      ;
    auto gllWts_ngll        = this->gllWts_ngll      ;
    auto &hyDensCells       = this->hyDensCells      ;
    auto &hyThetaCells      = this->hyThetaCells     ;
    auto &hyPressureCells   = this->hyPressureCells  ;
    auto &hyDensThetaCells  = this->hyDensThetaCells ;
    auto &hyDensGLL         = this->hyDensGLL        ;
    auto &hyThetaGLL        = this->hyThetaGLL       ;
    auto &hyPressureGLL     = this->hyPressureGLL    ;
    auto &hyDensThetaGLL    = this->hyDensThetaGLL   ;
    auto &dataSpec          = this->dataSpec         ;
    auto &sim2d             = this->sim2d            ;
    auto &xlen              = this->xlen             ;
    auto &ylen              = this->ylen             ;
    // Setup hydrostatic background state
    parallel_for( Bounds<1>(nz+2*hs) , YAKL_LAMBDA (int k) {
      // Compute cell averages
      hyDensCells     (k) = 0;
      hyPressureCells (k) = 0;
      hyThetaCells    (k) = 0;
      hyDensThetaCells(k) = 0;
      for (int kk=0; kk<ord; kk++) {
        real zloc = (k-hs+0.5_fp)*dz + gllPts_ord(kk)*dz;
        if        (dataSpec == DATA_SPEC_THERMAL || dataSpec == DATA_SPEC_THERMAL_MOIST) {
          // Compute constant theta hydrostatic background state
          real th  = 300;
          real rh = profiles::initConstTheta_density (th,zloc);
          real ph = profiles::initConstTheta_pressure(th,zloc);
          real wt = gllWts_ord(kk);
          hyDensCells     (k) += rh    * wt;
          hyThetaCells    (k) += th    * wt;
          hyDensThetaCells(k) += rh*th * wt;
          hyPressureCells (k) += ph    * wt;
        }
      }
    });
    parallel_for( Bounds<1>(nz) , YAKL_LAMBDA (int k) {
      // Compute ngll GLL points
      for (int kk=0; kk<ngll; kk++) {
        real zloc = (k+0.5_fp)*dz + gllPts_ngll(kk)*dz;
        if        (dataSpec == DATA_SPEC_THERMAL || dataSpec == DATA_SPEC_THERMAL_MOIST) {
          // Compute constant theta hydrostatic background state
          real th = 300;
          real rh = profiles::initConstTheta_density (th,zloc);
          real ph = profiles::initConstTheta_pressure(th,zloc);
          hyDensGLL     (k,kk) = rh;
          hyThetaGLL    (k,kk) = th;
          hyDensThetaGLL(k,kk) = rh*th;
          hyPressureGLL (k,kk) = ph;
        }
      }
    });

    // Compute the state
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      for (int l=0; l < numState; l++) {
        state(l,hs+k,hs+j,hs+i) = 0;
      }
      for (int kk=0; kk<ord; kk++) {
        for (int jj=0; jj<ord; jj++) {
          for (int ii=0; ii<ord; ii++) {
            real zloc = (k+0.5_fp)*dz + gllPts_ord(kk)*dz;
            real yloc;
            if (sim2d) {
              yloc = ylen/2;
            } else {
              yloc = (j+0.5_fp)*dy + gllPts_ord(jj)*dy;
            }
            real xloc = (i+0.5_fp)*dx + gllPts_ord(ii)*dx;
            real wt = gllWts_ord(kk) * gllWts_ord(jj) * gllWts_ord(ii);
            if        (dataSpec == DATA_SPEC_THERMAL || dataSpec == DATA_SPEC_THERMAL_MOIST) {
              // Compute constant theta hydrostatic background state
              real th = 300;
              real rh = profiles::initConstTheta_density(th,zloc);
              real tp = profiles::ellipsoid_linear(xloc, yloc, zloc, xlen/2, ylen/2, 2000, 2000, 2000, 2000, 2 );
              real t = th + tp;
              state(idT,hs+k,hs+j,hs+i) += (rh*t - rh*th) * wt;
            }
          }
        }
      }
    });
  }



  // Compute state and tendency time derivatives from the state
  void computeTendencies( StateArr  &state   , StateTendArr  &stateTend  ,
                          TracerArr &tracers , TracerTendArr &tracerTend ,
                          real &dt , int splitIndex ) {
    if (dimSwitch) {
      if        (splitIndex == 0) {
        computeTendenciesX( state , stateTend , tracers , tracerTend , dt );
      } else if (splitIndex == 1) {
        if (sim2d) {
          memset(stateTend  , 0._fp);
          memset(tracerTend , 0._fp);
        }
        else { computeTendenciesY( state , stateTend , tracers , tracerTend , dt ); }
      } else if (splitIndex == 2) {
        dt /= 2;
        computeTendenciesZ( state , stateTend , tracers , tracerTend , dt );
      } else if (splitIndex == 3) {
        dt /= 2;
        computeTendenciesZ( state , stateTend , tracers , tracerTend , dt );
      }
    } else {
      if        (splitIndex == 0) {
        dt /= 2;
        computeTendenciesZ( state , stateTend , tracers , tracerTend , dt );
      } else if (splitIndex == 1) {
        dt /= 2;
        computeTendenciesZ( state , stateTend , tracers , tracerTend , dt );
      } else if (splitIndex == 2) {
        if (sim2d) {
          memset(stateTend  , 0._fp);
          memset(tracerTend , 0._fp);
        }
        else { computeTendenciesY( state , stateTend , tracers , tracerTend , dt ); }
      } else if (splitIndex == 3) {
        computeTendenciesX( state , stateTend , tracers , tracerTend , dt );
      }
    }
    if (splitIndex == numSplit()-1) dimSwitch = ! dimSwitch;
  } // computeTendencies



  void computeTendenciesX( StateArr  &state   , StateTendArr  &stateTend  ,
                           TracerArr &tracers , TracerTendArr &tracerTend ,
                           real &dt ) {
    auto &nx                      = this->nx                     ;
    auto &weno_scalars            = this->weno_scalars           ;
    auto &weno_winds              = this->weno_winds             ;
    auto &c2g                     = this->coefs_to_gll           ;
    auto &s2g                     = this->sten_to_gll            ;
    auto &wenoRecon               = this->wenoRecon              ;
    auto &idl                     = this->idl                    ;
    auto &sigma                   = this->sigma                  ;
    auto &hyDensCells             = this->hyDensCells            ;
    auto &hyDensThetaCells        = this->hyDensThetaCells       ;
    auto &sim2d                   = this->sim2d                  ;
    auto &derivMatrix             = this->derivMatrix            ;
    auto &dx                      = this->dx                     ;
    auto &stateLimits             = this->stateLimits            ;
    auto &tracerLimits            = this->tracerLimits           ;
    auto &stateFluxLimits         = this->stateFluxLimits        ;
    auto &tracerFlux              = this->tracerFlux             ;
    auto &tracerPos               = this->tracerPos              ;
    auto &numTracers              = this->numTracers             ;
    auto &bc_x                    = this->bc_x                   ;

    // Pre-process the tracers by dividing by density inside the domain
    // After this, we can reconstruct tracers only (not rho * tracer)
    parallel_for( Bounds<4>(numTracers,nz,ny,nx) , YAKL_LAMBDA (int tr, int k, int j, int i) {
      tracers(tr,hs+k,hs+j,hs+i) /= (state(idR,hs+k,hs+j,hs+i) + hyDensCells(hs+k));
    });

    // Populate the halos
    if        (bc_x == BC_PERIODIC) {
      parallel_for( Bounds<3>(nz,ny,hs) , YAKL_LAMBDA(int k, int j, int ii) {
        for (int l=0; l < numState; l++) {
          state  (l,hs+k,hs+j,      ii) = state  (l,hs+k,hs+j,nx+ii);
          state  (l,hs+k,hs+j,hs+nx+ii) = state  (l,hs+k,hs+j,hs+ii);
        }
        for (int l=0; l < numTracers; l++) {
          tracers(l,hs+k,hs+j,      ii) = tracers(l,hs+k,hs+j,nx+ii);
          tracers(l,hs+k,hs+j,hs+nx+ii) = tracers(l,hs+k,hs+j,hs+ii);
        }
      });
    } else if (bc_x == BC_WALL) {
      parallel_for( Bounds<3>(nz,ny,hs) , YAKL_LAMBDA(int k, int j, int ii) {
        for (int l=0; l < numState; l++) {
          if (l == idU) {
            state(l,hs+k,hs+j,      ii) = 0;
            state(l,hs+k,hs+j,hs+nx+ii) = 0;
          } else {
            state  (l,hs+k,hs+j,      ii) = state  (l,hs+k,hs+j,hs     );
            state  (l,hs+k,hs+j,hs+nx+ii) = state  (l,hs+k,hs+j,hs+nx-1);
          }
        }
        for (int l=0; l < numTracers; l++) {
          tracers(l,hs+k,hs+j,      ii) = tracers(l,hs+k,hs+j,hs     );
          tracers(l,hs+k,hs+j,hs+nx+ii) = tracers(l,hs+k,hs+j,hs+nx-1);
        }
      });
    }

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // We need density and momentum to evolve the tracers with ADER
      SArray<real,2,nAder,ngll> r_DTs , ru_DTs;

      { // BEGIN: Reconstruct, time-average, and store the state and fluxes
        ////////////////////////////////////////////////////////////////
        // Reconstruct rho, u, v, w, theta
        ////////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> rv_DTs , rw_DTs , rt_DTs;
        { // BEGIN: Reconstruct the state
          SArray<real,1,ord> stencil;

          // Density
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idR,hs+k,hs+j,i+ii); }
          reconstruct_gll_values( stencil , r_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int ii=0; ii < ngll; ii++) { r_DTs(0,ii) += hyDensCells(hs+k); } // Add hydrostasis back on

          // u values and derivatives
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idU,hs+k,hs+j,i+ii); }
          reconstruct_gll_values( stencil , ru_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // v
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idV,hs+k,hs+j,i+ii); }
          reconstruct_gll_values( stencil , rv_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // w
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idW,hs+k,hs+j,i+ii); }
          reconstruct_gll_values( stencil , rw_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // theta
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idT,hs+k,hs+j,i+ii); }
          reconstruct_gll_values( stencil , rt_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) += hyDensThetaCells(hs+k); } // Add hydrostasis back on
        } // END: Reconstruct the state

        ///////////////////////////////////////////////////////////////
        // Compute other values needed for centered tendencies and DTs
        ///////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> ruu_DTs , ruv_DTs , ruw_DTs , rut_DTs , rt_gamma_DTs;
        for (int ii=0; ii < ngll; ii++) {
          real r = r_DTs (0,ii);
          real u = ru_DTs(0,ii) / r;
          real v = rv_DTs(0,ii) / r;
          real w = rw_DTs(0,ii) / r;
          real t = rt_DTs(0,ii) / r;
          ruu_DTs     (0,ii) = r*u*u;
          ruv_DTs     (0,ii) = r*u*v;
          ruw_DTs     (0,ii) = r*u*w;
          rut_DTs     (0,ii) = r*u*t;
          rt_gamma_DTs(0,ii) = pow(r*t,GAMMA);
        }

        //////////////////////////////////////////
        // Compute time derivatives if necessary
        //////////////////////////////////////////
        if (nAder > 1) {
          diffTransformEulerConsX( r_DTs , ru_DTs , rv_DTs , rw_DTs , rt_DTs , ruu_DTs , ruv_DTs , ruw_DTs ,
                                   rut_DTs , rt_gamma_DTs , derivMatrix , dx );
        }

        //////////////////////////////////////////
        // Time average if necessary
        //////////////////////////////////////////
        // density and momentum can't be overwritten because they will be used for tracers
        SArray<real,1,ngll> r_tavg, ru_tavg;
        if (timeAvg) {
          compute_timeAvg( r_DTs  , r_tavg  , dt );
          compute_timeAvg( ru_DTs , ru_tavg , dt );
          compute_timeAvg( rv_DTs           , dt );
          compute_timeAvg( rw_DTs           , dt );
          compute_timeAvg( rt_DTs           , dt );
          compute_timeAvg( ruu_DTs          , dt );
          compute_timeAvg( ruv_DTs          , dt );
          compute_timeAvg( ruw_DTs          , dt );
          compute_timeAvg( rut_DTs          , dt );
          compute_timeAvg( rt_gamma_DTs     , dt );
        } else {
          for (int ii=0; ii < ngll; ii++) {
            r_tavg (ii) = r_DTs (0,ii);
            ru_tavg(ii) = ru_DTs(0,ii);
          }
        }

        //////////////////////////////////////////
        // Store cell edge estimates of the state
        //////////////////////////////////////////
        // Left interface
        stateLimits(idR,1,k,j,i  ) = r_tavg  (0     );
        stateLimits(idU,1,k,j,i  ) = ru_tavg (0     );
        stateLimits(idV,1,k,j,i  ) = rv_DTs(0,0     );
        stateLimits(idW,1,k,j,i  ) = rw_DTs(0,0     );
        stateLimits(idT,1,k,j,i  ) = rt_DTs(0,0     );
        // Right interface
        stateLimits(idR,0,k,j,i+1) = r_tavg  (ngll-1);
        stateLimits(idU,0,k,j,i+1) = ru_tavg (ngll-1);
        stateLimits(idV,0,k,j,i+1) = rv_DTs(0,ngll-1);
        stateLimits(idW,0,k,j,i+1) = rw_DTs(0,ngll-1);
        stateLimits(idT,0,k,j,i+1) = rt_DTs(0,ngll-1);

        //////////////////////////////////////////
        // Store cell edge estimates of the flux
        //////////////////////////////////////////
        // Left interface
        stateFluxLimits(idR,1,k,j,i  ) = ru_tavg  (0     );
        stateFluxLimits(idU,1,k,j,i  ) = ruu_DTs(0,0     ) + C0*rt_gamma_DTs(0,0     );
        stateFluxLimits(idV,1,k,j,i  ) = ruv_DTs(0,0     );
        stateFluxLimits(idW,1,k,j,i  ) = ruw_DTs(0,0     );
        stateFluxLimits(idT,1,k,j,i  ) = rut_DTs(0,0     );
        // Right interface
        stateFluxLimits(idR,0,k,j,i+1) = ru_tavg  (ngll-1);
        stateFluxLimits(idU,0,k,j,i+1) = ruu_DTs(0,ngll-1) + C0*rt_gamma_DTs(0,ngll-1);
        stateFluxLimits(idV,0,k,j,i+1) = ruv_DTs(0,ngll-1);
        stateFluxLimits(idW,0,k,j,i+1) = ruw_DTs(0,ngll-1);
        stateFluxLimits(idT,0,k,j,i+1) = rut_DTs(0,ngll-1);
      } // END: Reconstruct, time-average, and store the state and fluxes

      // r_DTs and ru_DTs still exist and are computed
      { // BEGIN: Reconstruct, time-average, and store tracer fluxes
        // Only process one tracer at a time to save on local memory / register requirements
        for (int tr=0; tr < numTracers; tr++) {
          SArray<real,2,nAder,ngll> rt_DTs; // Density * tracer
          { // BEGIN: Reconstruct the tracer
            SArray<real,1,ord> stencil;
            for (int ii=0; ii < ord; ii++) { stencil(ii) = tracers(tr,hs+k,hs+j,i+ii); }
            reconstruct_gll_values( stencil , rt_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) *= r_DTs(0,ii); }
            if (tracerPos(tr)) {
              for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) = max( 0._fp , rt_DTs(0,ii) ); }
            }
          } // END: Reconstruct the tracer

          // Compute the tracer flux
          SArray<real,2,nAder,ngll> rut_DTs; // Density * uwind * tracer
          for (int ii=0; ii < ngll; ii++) {
            rut_DTs(0,ii) = rt_DTs(0,ii) * ru_DTs(0,ii) / r_DTs(0,ii);
          }

          //////////////////////////////////////////
          // Compute time derivatives if necessary
          //////////////////////////////////////////
          if (nAder > 1) {
            diffTransformTracer( r_DTs , ru_DTs , rt_DTs , rut_DTs , derivMatrix , dx );
          }

          //////////////////////////////////////////
          // Time average if necessary
          //////////////////////////////////////////
          if (timeAvg) {
            compute_timeAvg( rt_DTs  , dt );
            compute_timeAvg( rut_DTs , dt );
          }
          if (tracerPos(tr)) {
            for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) = max( 0._fp , rt_DTs(0,ii) ); }
          }

          ////////////////////////////////////////////////////////////
          // Store cell edge estimates of the tracer
          ////////////////////////////////////////////////////////////
          tracerLimits(tr,1,k,j,i  ) = rt_DTs (0,0     ); // Left interface
          tracerLimits(tr,0,k,j,i+1) = rt_DTs (0,ngll-1); // Right interface
        }
      } // END: Reconstruct, time-average, and store tracer fluxes

    });

    ////////////////////////////////////////////////
    // BCs for the state edge estimates
    ////////////////////////////////////////////////
    parallel_for( Bounds<2>(nz,ny) , YAKL_LAMBDA (int k, int j) {
      for (int l=0; l < numState; l++) {
        if        (bc_x == BC_PERIODIC) {
          stateLimits     (l,0,k,j,0 ) = stateLimits     (l,0,k,j,nx);
          stateFluxLimits (l,0,k,j,0 ) = stateFluxLimits (l,0,k,j,nx);
          stateLimits     (l,1,k,j,nx) = stateLimits     (l,1,k,j,0 );
          stateFluxLimits (l,1,k,j,nx) = stateFluxLimits (l,1,k,j,0 );
        } else if (bc_x == BC_WALL    ) {
          stateLimits     (l,0,k,j,0 ) = stateLimits     (l,1,k,j,0 );
          stateFluxLimits (l,0,k,j,0 ) = stateFluxLimits (l,1,k,j,0 );
          stateLimits     (l,1,k,j,nx) = stateLimits     (l,0,k,j,nx);
          stateFluxLimits (l,1,k,j,nx) = stateFluxLimits (l,0,k,j,nx);
        }
      }
      for (int l=0; l < numTracers; l++) {
        if        (bc_x == BC_PERIODIC) {
          tracerLimits    (l,0,k,j,0 ) = tracerLimits    (l,0,k,j,nx);
          tracerLimits    (l,1,k,j,nx) = tracerLimits    (l,1,k,j,0 );
        } else if (bc_x == BC_WALL    ) {
          tracerLimits    (l,0,k,j,0 ) = tracerLimits    (l,1,k,j,0 );
          tracerLimits    (l,1,k,j,nx) = tracerLimits    (l,0,k,j,nx);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<3>(nz,ny,nx+1) , YAKL_LAMBDA (int k, int j, int i) {
      // Get left and right state
      real r_L = stateLimits(idR,0,k,j,i)    ;   real r_R = stateLimits(idR,1,k,j,i)    ;
      real u_L = stateLimits(idU,0,k,j,i)/r_L;   real u_R = stateLimits(idU,1,k,j,i)/r_R;
      real v_L = stateLimits(idV,0,k,j,i)/r_L;   real v_R = stateLimits(idV,1,k,j,i)/r_R;
      real w_L = stateLimits(idW,0,k,j,i)/r_L;   real w_R = stateLimits(idW,1,k,j,i)/r_R;
      real t_L = stateLimits(idT,0,k,j,i)/r_L;   real t_R = stateLimits(idT,1,k,j,i)/r_R;
      // Compute average state
      real r = 0.5_fp * (r_L + r_R);
      real u = 0.5_fp * (u_L + u_R);
      real v = 0.5_fp * (v_L + v_R);
      real w = 0.5_fp * (w_L + w_R);
      real t = 0.5_fp * (t_L + t_R);
      real p = C0 * pow(r*t,GAMMA);
      real cs2 = GAMMA*p/r;
      real cs  = sqrt(cs2);

      // COMPUTE UPWIND STATE FLUXES
      // Get left and right fluxes
      real f1_L = stateFluxLimits(idR,0,k,j,i);   real f1_R = stateFluxLimits(idR,1,k,j,i);
      real f2_L = stateFluxLimits(idU,0,k,j,i);   real f2_R = stateFluxLimits(idU,1,k,j,i);
      real f3_L = stateFluxLimits(idV,0,k,j,i);   real f3_R = stateFluxLimits(idV,1,k,j,i);
      real f4_L = stateFluxLimits(idW,0,k,j,i);   real f4_R = stateFluxLimits(idW,1,k,j,i);
      real f5_L = stateFluxLimits(idT,0,k,j,i);   real f5_R = stateFluxLimits(idT,1,k,j,i);
      // Compute upwind characteristics
      // Waves 1-3, velocity: u
      real w1, w2, w3;
      if (u > 0) {
        w1 = f1_L - f5_L/t;
        w2 = f3_L - v*f5_L/t;
        w3 = f4_L - w*f5_L/t;
      } else {
        w1 = f1_R - f5_R/t;
        w2 = f3_R - v*f5_R/t;
        w3 = f4_R - w*f5_R/t;
      }
      // Wave 5, velocity: u-cs
      real w5 =  u*f1_R/(2*cs) - f2_R/(2*cs) + f5_R/(2*t);
      // Wave 6, velocity: u+cs
      real w6 = -u*f1_L/(2*cs) + f2_L/(2*cs) + f5_L/(2*t);
      // Use right eigenmatrix to compute upwind flux
      stateFluxLimits(idR,0,k,j,i) = w1 + w5 + w6;
      stateFluxLimits(idU,0,k,j,i) = u*w1 + (u-cs)*w5 + (u+cs)*w6;
      stateFluxLimits(idV,0,k,j,i) = w2 + v*w5 + v*w6;
      stateFluxLimits(idW,0,k,j,i) = w3 + w*w5 + w*w6;

      real massFlux = stateFluxLimits(idR,0,k,j,i);

      if (u > 0) {
        stateFluxLimits(idT,0,k,j,i) = massFlux * stateLimits(idT,0,k,j,i) / r_L;
      } else {
        stateFluxLimits(idT,0,k,j,i) = massFlux * stateLimits(idT,1,k,j,i) / r_R;
      }

      // COMPUTE UPWIND TRACER FLUXES
      // Handle it one tracer at a time
      for (int tr=0; tr < numTracers; tr++) {
        if (u > 0) {
          tracerFlux(tr,k,j,i) = massFlux * tracerLimits(tr,0,k,j,i) / r_L;
        } else {
          tracerFlux(tr,k,j,i) = massFlux * tracerLimits(tr,1,k,j,i) / r_R;
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Limit the tracer fluxes for positivity
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<4>(numTracers,nz,ny,nx+1) , YAKL_LAMBDA (int tr, int k, int j, int i) {
      real constexpr eps = 1.e-10;
      real u = 0.5_fp * ( stateLimits(idU,0,k,j,i) + stateLimits(idU,1,k,j,i) );
      // Solid wall BCs mean u == 0 at boundaries, so we assume periodic if u != 0
      if (tracerPos(tr)) {
        // Compute and apply the flux reduction factor of the upwind cell
        if      (u > 0) {
          // upwind is to the left of this interface
          int ind_i = i-1;
          if (ind_i == -1) ind_i = nx-1;
          real f1 = min( tracerFlux(tr,k,j,ind_i  ) , 0._fp );
          real f2 = max( tracerFlux(tr,k,j,ind_i+1) , 0._fp );
          real fluxOut = dt*(f2-f1)/dx;
          real dens = state(idR,hs+k,hs+j,hs+ind_i) + hyDensCells(hs+k);
          tracerFlux(tr,k,j,i) *= min( 1._fp , tracers(tr,hs+k,hs+j,hs+ind_i) * dens / (fluxOut + eps) );
        } else if (u < 0) {
          // upwind is to the right of this interface
          int ind_i = i;
          if (ind_i == nx) ind_i = 0;
          real f1 = min( tracerFlux(tr,k,j,ind_i  ) , 0._fp );
          real f2 = max( tracerFlux(tr,k,j,ind_i+1) , 0._fp );
          real fluxOut = dt*(f2-f1)/dx;
          real dens = state(idR,hs+k,hs+j,hs+ind_i) + hyDensCells(hs+k);
          tracerFlux(tr,k,j,i) *= min( 1._fp , tracers(tr,hs+k,hs+j,hs+ind_i) * dens / (fluxOut + eps) );
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA(int k, int j, int i) {
      for (int l = 0; l < numState; l++) {
        if (sim2d && l == idV) {
          stateTend(l,k,j,i) = 0;
        } else {
          stateTend(l,k,j,i) = - ( stateFluxLimits(l,0,k,j,i+1) - stateFluxLimits(l,0,k,j,i) ) / dx;
        }
      }
      for (int l = 0; l < numTracers; l++) {
        // Compute tracer tendency
        tracerTend(l,k,j,i) = - ( tracerFlux(l,k,j,i+1) - tracerFlux(l,k,j,i  ) ) / dx;
        // Multiply density back onto tracers
        tracers(l,hs+k,hs+j,hs+i) *= (state(idR,hs+k,hs+j,hs+i) + hyDensCells(hs+k));
      }
    });
  }



  void computeTendenciesY( StateArr  &state   , StateTendArr  &stateTend  ,
                           TracerArr &tracers , TracerTendArr &tracerTend ,
                           real &dt ) {
    auto &ny                      = this->ny                     ;
    auto &weno_scalars            = this->weno_scalars           ;
    auto &weno_winds              = this->weno_winds             ;
    auto &c2g                     = this->coefs_to_gll           ;
    auto &s2g                     = this->sten_to_gll            ;
    auto &wenoRecon               = this->wenoRecon              ;
    auto &idl                     = this->idl                    ;
    auto &sigma                   = this->sigma                  ;
    auto &hyDensCells             = this->hyDensCells            ;
    auto &hyDensThetaCells        = this->hyDensThetaCells       ;
    auto &sim2d                   = this->sim2d                  ;
    auto &derivMatrix             = this->derivMatrix            ;
    auto &dy                      = this->dy                     ;
    auto &stateLimits             = this->stateLimits            ;
    auto &stateFluxLimits         = this->stateFluxLimits        ;
    auto &tracerLimits            = this->tracerLimits           ;
    auto &tracerFlux              = this->tracerFlux             ;
    auto &tracerPos               = this->tracerPos              ;
    auto &numTracers              = this->numTracers             ;
    auto &bc_y                    = this->bc_y                   ;

    // Pre-process the tracers by dividing by density inside the domain
    // After this, we can reconstruct tracers only (not rho * tracer)
    parallel_for( Bounds<4>(numTracers,nz,ny,nx) , YAKL_LAMBDA (int tr, int k, int j, int i) {
      tracers(tr,hs+k,hs+j,hs+i) /= (state(idR,hs+k,hs+j,hs+i) + hyDensCells(hs+k));
    });

    // Populate the halos
    if        (bc_y == BC_PERIODIC) {
      parallel_for( Bounds<3>(nz,nx,hs) , YAKL_LAMBDA(int k, int i, int jj) {
        for (int l=0; l < numState; l++) {
          state(l,hs+k,      jj,hs+i) = state(l,hs+k,ny+jj,hs+i);
          state(l,hs+k,hs+ny+jj,hs+i) = state(l,hs+k,hs+jj,hs+i);
        }
        for (int l=0; l < numTracers; l++) {
          tracers(l,hs+k,      jj,hs+i) = tracers(l,hs+k,ny+jj,hs+i);
          tracers(l,hs+k,hs+ny+jj,hs+i) = tracers(l,hs+k,hs+jj,hs+i);
        }
      });
    } else if (bc_y == BC_WALL) {
      parallel_for( Bounds<3>(nz,nx,hs) , YAKL_LAMBDA(int k, int i, int jj) {
        for (int l=0; l < numState; l++) {
          if (l == idV) {
            state(l,hs+k,      jj,hs+i) = 0;
            state(l,hs+k,hs+ny+jj,hs+i) = 0;
          } else {
            state(l,hs+k,      jj,hs+i) = state(l,hs+k,hs     ,hs+i);
            state(l,hs+k,hs+ny+jj,hs+i) = state(l,hs+k,hs+ny-1,hs+i);
          }
        }
        for (int l=0; l < numTracers; l++) {
          tracers(l,hs+k,      jj,hs+i) = tracers(l,hs+k,hs     ,hs+i);
          tracers(l,hs+k,hs+ny+jj,hs+i) = tracers(l,hs+k,hs+ny-1,hs+i);
        }
      });
    }

    // Loop through all cells, reconstruct in y-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // These are needed by the tracers
      SArray<real,2,nAder,ngll> r_DTs , rv_DTs;

      { // BEGIN: Reconstruct, time-average, and store state and sate fluxes
        ////////////////////////////////////////////////////////////////
        // Reconstruct rho, u, v, w, theta
        ////////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> ru_DTs , rw_DTs , rt_DTs;
        {
          SArray<real,1,ord> stencil;

          // Density
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idR,hs+k,j+jj,hs+i); }
          reconstruct_gll_values( stencil , r_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int jj=0; jj < ngll; jj++) { r_DTs(0,jj) += hyDensCells(hs+k); } // Add hydrostasis back on

          // u values and derivatives
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idU,hs+k,j+jj,hs+i); }
          reconstruct_gll_values( stencil , ru_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // v
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idV,hs+k,j+jj,hs+i); }
          reconstruct_gll_values( stencil , rv_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // w
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idW,hs+k,j+jj,hs+i); }
          reconstruct_gll_values( stencil , rw_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // theta
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idT,hs+k,j+jj,hs+i); }
          reconstruct_gll_values( stencil , rt_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) += hyDensThetaCells(hs+k); } // Add hydrostasis back on
        }

        ///////////////////////////////////////////////////////////////
        // Compute other values needed for centered tendencies and DTs
        ///////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> rvu_DTs , rvv_DTs , rvw_DTs , rvt_DTs , rt_gamma_DTs;
        for (int jj=0; jj < ngll; jj++) {
          real r = r_DTs (0,jj);
          real u = ru_DTs(0,jj) / r;
          real v = rv_DTs(0,jj) / r;
          real w = rw_DTs(0,jj) / r;
          real t = rt_DTs(0,jj) / r;
          rvu_DTs     (0,jj) = r*v*u;
          rvv_DTs     (0,jj) = r*v*v;
          rvw_DTs     (0,jj) = r*v*w;
          rvt_DTs     (0,jj) = r*v*t;
          rt_gamma_DTs(0,jj) = pow(r*t,GAMMA);
        }

        //////////////////////////////////////////
        // Compute time derivatives if necessary
        //////////////////////////////////////////
        if (nAder > 1) {
          diffTransformEulerConsY( r_DTs , ru_DTs , rv_DTs , rw_DTs , rt_DTs , rvu_DTs , rvv_DTs , rvw_DTs ,
                                   rvt_DTs , rt_gamma_DTs , derivMatrix , dy );
        }

        //////////////////////////////////////////
        // Time average if necessary
        //////////////////////////////////////////
        // Don't overwrite r and rv because we need them for tracers
        SArray<real,1,ngll> r_tavg, rv_tavg;
        if (timeAvg) {
          compute_timeAvg( r_DTs  , r_tavg  , dt );
          compute_timeAvg( ru_DTs           , dt );
          compute_timeAvg( rv_DTs , rv_tavg , dt );
          compute_timeAvg( rw_DTs           , dt );
          compute_timeAvg( rt_DTs           , dt );
          compute_timeAvg( rvu_DTs          , dt );
          compute_timeAvg( rvv_DTs          , dt );
          compute_timeAvg( rvw_DTs          , dt );
          compute_timeAvg( rvt_DTs          , dt );
          compute_timeAvg( rt_gamma_DTs     , dt );
        } else {
          for (int jj=0; jj < ngll; jj++) {
            r_tavg (jj) = r_DTs (0,jj);
            rv_tavg(jj) = rv_DTs(0,jj);
          }
        }

        //////////////////////////////////////////
        // Store cell edge estimates of the state
        //////////////////////////////////////////
        // Left interface
        stateLimits(idR,1,k,j  ,i) = r_tavg  (0     );
        stateLimits(idU,1,k,j  ,i) = ru_DTs(0,0     );
        stateLimits(idV,1,k,j  ,i) = rv_tavg (0     );
        stateLimits(idW,1,k,j  ,i) = rw_DTs(0,0     );
        stateLimits(idT,1,k,j  ,i) = rt_DTs(0,0     );
        // Right interface       
        stateLimits(idR,0,k,j+1,i) = r_tavg  (ngll-1);
        stateLimits(idU,0,k,j+1,i) = ru_DTs(0,ngll-1);
        stateLimits(idV,0,k,j+1,i) = rv_tavg (ngll-1);
        stateLimits(idW,0,k,j+1,i) = rw_DTs(0,ngll-1);
        stateLimits(idT,0,k,j+1,i) = rt_DTs(0,ngll-1);

        //////////////////////////////////////////
        // Store cell edge estimates of the flux
        //////////////////////////////////////////
        // Left interface
        stateFluxLimits(idR,1,k,j  ,i) = rv_tavg  (0     );
        stateFluxLimits(idU,1,k,j  ,i) = rvu_DTs(0,0     );
        stateFluxLimits(idV,1,k,j  ,i) = rvv_DTs(0,0     ) + C0*rt_gamma_DTs(0,0     );
        stateFluxLimits(idW,1,k,j  ,i) = rvw_DTs(0,0     );
        stateFluxLimits(idT,1,k,j  ,i) = rvt_DTs(0,0     );
        // Right interface      
        stateFluxLimits(idR,0,k,j+1,i) = rv_tavg  (ngll-1);
        stateFluxLimits(idU,0,k,j+1,i) = rvu_DTs(0,ngll-1);
        stateFluxLimits(idV,0,k,j+1,i) = rvv_DTs(0,ngll-1) + C0*rt_gamma_DTs(0,ngll-1);
        stateFluxLimits(idW,0,k,j+1,i) = rvw_DTs(0,ngll-1);
        stateFluxLimits(idT,0,k,j+1,i) = rvt_DTs(0,ngll-1);
      } // END: Reconstruct, time-average, and store state and sate fluxes

      // r_DTs and rv_DTs still exist and are computed
      { // BEGIN: Reconstruct, time-average, and store tracer fluxes
        // Only process one tracer at a time to save on local memory / register requirements
        for (int tr=0; tr < numTracers; tr++) {
          SArray<real,2,nAder,ngll> rt_DTs; // Density * tracer
          { // BEGIN: Reconstruct the tracer
            SArray<real,1,ord> stencil;
            for (int jj=0; jj < ord; jj++) { stencil(jj) = tracers(tr,hs+k,j+jj,hs+i); }
            reconstruct_gll_values( stencil , rt_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) *= r_DTs(0,jj); }
            if (tracerPos(tr)) {
              for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) = max( 0._fp , rt_DTs(0,jj) ); }
            }
          } // END: Reconstruct the tracer

          // Compute the tracer flux
          SArray<real,2,nAder,ngll> rvt_DTs; // Density * vwind * tracer
          for (int jj=0; jj < ngll; jj++) {
            rvt_DTs(0,jj) = rt_DTs(0,jj) * rv_DTs(0,jj) / r_DTs(0,jj);
          }

          //////////////////////////////////////////
          // Compute time derivatives if necessary
          //////////////////////////////////////////
          if (nAder > 1) {
            diffTransformTracer( r_DTs , rv_DTs , rt_DTs , rvt_DTs , derivMatrix , dy );
          }

          //////////////////////////////////////////
          // Time average if necessary
          //////////////////////////////////////////
          if (timeAvg) {
            compute_timeAvg( rt_DTs  , dt );
            compute_timeAvg( rvt_DTs , dt );
          }
          if (tracerPos(tr)) {
            for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) = max( 0._fp , rt_DTs(0,jj) ); }
          }

          ////////////////////////////////////////////////////////////
          // Store cell edge estimates of the tracer
          ////////////////////////////////////////////////////////////
          tracerLimits(tr,1,k,j  ,i) = rt_DTs (0,0     ); // Left interface
          tracerLimits(tr,0,k,j+1,i) = rt_DTs (0,ngll-1); // Right interface
        }
      } // END: Reconstruct, time-average, and store tracer fluxes

    });

    ////////////////////////////////////////////////
    // BCs for the state edge estimates
    ////////////////////////////////////////////////
    parallel_for( Bounds<2>(nz,nx) , YAKL_LAMBDA (int k, int i) {
      for (int l=0; l < numState; l++) {
        if        (bc_y == BC_PERIODIC) {
          stateLimits     (l,0,k,0 ,i) = stateLimits     (l,0,k,ny,i);
          stateFluxLimits (l,0,k,0 ,i) = stateFluxLimits (l,0,k,ny,i);
          stateLimits     (l,1,k,ny,i) = stateLimits     (l,1,k,0 ,i);
          stateFluxLimits (l,1,k,ny,i) = stateFluxLimits (l,1,k,0 ,i);
        } else if (bc_y == BC_WALL    ) {
          stateLimits     (l,0,k,0 ,i) = stateLimits     (l,1,k,0 ,i);
          stateFluxLimits (l,0,k,0 ,i) = stateFluxLimits (l,1,k,0 ,i);
          stateLimits     (l,1,k,ny,i) = stateLimits     (l,0,k,ny,i);
          stateFluxLimits (l,1,k,ny,i) = stateFluxLimits (l,0,k,ny,i);
        }
      }
      for (int l=0; l < numTracers; l++) {
        if        (bc_y == BC_PERIODIC) {
          tracerLimits    (l,0,k,0 ,i) = tracerLimits    (l,0,k,ny,i);
          tracerLimits    (l,1,k,ny,i) = tracerLimits    (l,1,k,0 ,i);
        } else if (bc_y == BC_WALL    ) {
          tracerLimits    (l,0,k,0 ,i) = tracerLimits    (l,1,k,0 ,i);
          tracerLimits    (l,1,k,ny,i) = tracerLimits    (l,0,k,ny,i);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<3>(nz,ny+1,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // Get left and right state
      real r_L = stateLimits(idR,0,k,j,i)    ;   real r_R = stateLimits(idR,1,k,j,i)    ;
      real u_L = stateLimits(idU,0,k,j,i)/r_L;   real u_R = stateLimits(idU,1,k,j,i)/r_R;
      real v_L = stateLimits(idV,0,k,j,i)/r_L;   real v_R = stateLimits(idV,1,k,j,i)/r_R;
      real w_L = stateLimits(idW,0,k,j,i)/r_L;   real w_R = stateLimits(idW,1,k,j,i)/r_R;
      real t_L = stateLimits(idT,0,k,j,i)/r_L;   real t_R = stateLimits(idT,1,k,j,i)/r_R;
      // Compute average state
      real r = 0.5_fp * (r_L + r_R);
      real u = 0.5_fp * (u_L + u_R);
      real v = 0.5_fp * (v_L + v_R);
      real w = 0.5_fp * (w_L + w_R);
      real t = 0.5_fp * (t_L + t_R);
      real p = C0 * pow(r*t,GAMMA);
      real cs2 = GAMMA*p/r;
      real cs  = sqrt(cs2);

      // COMPUTE UPWIND STATE FLUXES
      // Get left and right fluxes
      real f1_L = stateFluxLimits(idR,0,k,j,i);   real f1_R = stateFluxLimits(idR,1,k,j,i);
      real f2_L = stateFluxLimits(idU,0,k,j,i);   real f2_R = stateFluxLimits(idU,1,k,j,i);
      real f3_L = stateFluxLimits(idV,0,k,j,i);   real f3_R = stateFluxLimits(idV,1,k,j,i);
      real f4_L = stateFluxLimits(idW,0,k,j,i);   real f4_R = stateFluxLimits(idW,1,k,j,i);
      real f5_L = stateFluxLimits(idT,0,k,j,i);   real f5_R = stateFluxLimits(idT,1,k,j,i);
      // Compute upwind characteristics
      // Waves 1-3, velocity: v
      real w1, w2, w3;
      if (v > 0) {
        w1 = f1_L - f5_L/t;
        w2 = f2_L - u*f5_L/t;
        w3 = f4_L - w*f5_L/t;
      } else {
        w1 = f1_R - f5_R/t;
        w2 = f2_R - u*f5_R/t;
        w3 = f4_R - w*f5_R/t;
      }
      // Wave 5, velocity: v-cs
      real w5 =  v*f1_R/(2*cs) - f3_R/(2*cs) + f5_R/(2*t);
      // Wave 6, velocity: v+cs
      real w6 = -v*f1_L/(2*cs) + f3_L/(2*cs) + f5_L/(2*t);
      // Use right eigenmatrix to compute upwind flux
      stateFluxLimits(idR,0,k,j,i) = w1 + w5 + w6;
      stateFluxLimits(idU,0,k,j,i) = w2 + u*w5 + u*w6;
      stateFluxLimits(idV,0,k,j,i) = v*w1 + (v-cs)*w5 + (v+cs)*w6;
      stateFluxLimits(idW,0,k,j,i) = w3 + w*w5 + w*w6;

      real massFlux = stateFluxLimits(idR,0,k,j,i);

      if (v > 0) {
        stateFluxLimits(idT,0,k,j,i) = massFlux * stateLimits(idT,0,k,j,i) / r_L;
      } else {
        stateFluxLimits(idT,0,k,j,i) = massFlux * stateLimits(idT,1,k,j,i) / r_R;
      }

      // COMPUTE UPWIND TRACER FLUXES
      // Handle it one tracer at a time
      for (int tr=0; tr < numTracers; tr++) {
        if (v > 0) {
          tracerFlux(tr,k,j,i) = massFlux * tracerLimits(tr,0,k,j,i) / r_L;
        } else {
          tracerFlux(tr,k,j,i) = massFlux * tracerLimits(tr,1,k,j,i) / r_R;
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Limit the tracer fluxes for positivity
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<4>(numTracers,nz,ny+1,nx) , YAKL_LAMBDA (int tr, int k, int j, int i) {
      real constexpr eps = 1.e-10;
      real v = 0.5_fp * ( stateLimits(idV,0,k,j,i) + stateLimits(idV,1,k,j,i) );
      // Solid wall BCs mean u == 0 at boundaries, so we assume periodic if u != 0
      if (tracerPos(tr)) {
        // Compute and apply the flux reduction factor of the upwind cell
        if      (v > 0) {
          // upwind is to the left of this interface
          int ind_j = j-1;
          if (ind_j == -1) ind_j = ny-1;
          real f1 = min( tracerFlux(tr,k,ind_j  ,i) , 0._fp );
          real f2 = max( tracerFlux(tr,k,ind_j+1,i) , 0._fp );
          real fluxOut = dt*(f2-f1)/dy;
          real dens = state(idR,hs+k,hs+ind_j,hs+i) + hyDensCells(hs+k);
          tracerFlux(tr,k,j,i) *= min( 1._fp , tracers(tr,hs+k,hs+ind_j,hs+i) * dens / (fluxOut + eps) );
        } else if (v < 0) {
          // upwind is to the right of this interface
          int ind_j = j;
          if (ind_j == ny) ind_j = 0;
          real f1 = min( tracerFlux(tr,k,j,ind_j  ) , 0._fp );
          real f2 = max( tracerFlux(tr,k,j,ind_j+1) , 0._fp );
          real fluxOut = dt*(f2-f1)/dy;
          real dens = state(idR,hs+k,hs+ind_j,hs+i) + hyDensCells(hs+k);
          tracerFlux(tr,k,j,i) *= min( 1._fp , tracers(tr,hs+k,hs+ind_j,hs+i) * dens / (fluxOut + eps) );
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA(int k, int j, int i) {
      for (int l=0; l < numState; l++) {
        stateTend(l,k,j,i) = - ( stateFluxLimits(l,0,k,j+1,i) - stateFluxLimits(l,0,k,j,i) ) / dy;
      }
      for (int l=0; l < numTracers; l++) {
        // Compute the tracer tendency
        tracerTend(l,k,j,i) = - ( tracerFlux(l,k,j+1,i) - tracerFlux(l,k,j,i) ) / dy;
        // Multiply density back onto the tracers
        tracers(l,hs+k,hs+j,hs+i) *= (state(idR,hs+k,hs+j,hs+i) + hyDensCells(hs+k));
      }
    });
  }



  void computeTendenciesZ( StateArr  &state   , StateTendArr  &stateTend  ,
                           TracerArr &tracers , TracerTendArr &tracerTend ,
                           real &dt ) {
    auto &nz                      = this->nz                     ;
    auto &weno_scalars            = this->weno_scalars           ;
    auto &weno_winds              = this->weno_winds             ;
    auto &c2g                     = this->coefs_to_gll           ;
    auto &s2g                     = this->sten_to_gll            ;
    auto &wenoRecon               = this->wenoRecon              ;
    auto &idl                     = this->idl                    ;
    auto &sigma                   = this->sigma                  ;
    auto &hyDensCells             = this->hyDensCells            ;
    auto &hyDensGLL               = this->hyDensGLL              ;
    auto &hyDensThetaGLL          = this->hyDensThetaGLL         ;
    auto &hyPressureGLL           = this->hyPressureGLL          ;
    auto &sim2d                   = this->sim2d                  ;
    auto &derivMatrix             = this->derivMatrix            ;
    auto &dz                      = this->dz                     ;
    auto &stateLimits             = this->stateLimits            ;
    auto &stateFluxLimits         = this->stateFluxLimits        ;
    auto &tracerLimits            = this->tracerLimits           ;
    auto &tracerFlux              = this->tracerFlux             ;
    auto &tracerPos               = this->tracerPos              ;
    auto &numTracers              = this->numTracers             ;
    auto &bc_z                    = this->bc_z                   ;
    auto &gllWts_ngll             = this->gllWts_ngll            ;

    // Pre-process the tracers by dividing by density inside the domain
    // After this, we can reconstruct tracers only (not rho * tracer)
    parallel_for( Bounds<4>(numTracers,nz,ny,nx) , YAKL_LAMBDA (int tr, int k, int j, int i) {
      tracers(tr,hs+k,hs+j,hs+i) /= (state(idR,hs+k,hs+j,hs+i) + hyDensCells(hs+k));
    });

    // Populate the halos
    if        (bc_z == BC_PERIODIC) {
      parallel_for( Bounds<3>(ny,nx,hs) , YAKL_LAMBDA(int j, int i, int kk) {
        for (int l=0; l < numState; l++) {
          state(l,      kk,hs+j,hs+i) = state(l,nz+kk,hs+j,hs+i);
          state(l,hs+nz+kk,hs+j,hs+i) = state(l,hs+kk,hs+j,hs+i);
        }
        for (int l=0; l < numTracers; l++) {
          tracers(l,      kk,hs+j,hs+i) = tracers(l,nz+kk,hs+j,hs+i);
          tracers(l,hs+nz+kk,hs+j,hs+i) = tracers(l,hs+kk,hs+j,hs+i);
        }
      });
    } else if (bc_z == BC_WALL) {
      parallel_for( Bounds<3>(ny,nx,hs) , YAKL_LAMBDA(int j, int i, int kk) {
        for (int l=0; l < numState; l++) {
          if (l == idW) {
            state(l,      kk,hs+j,hs+i) = 0;
            state(l,hs+nz+kk,hs+j,hs+i) = 0;
          } else {
            state(l,      kk,hs+j,hs+i) = state(l,hs     ,hs+j,hs+i);
            state(l,hs+nz+kk,hs+j,hs+i) = state(l,hs+nz-1,hs+j,hs+i);
          }
        }
        for (int l=0; l < numTracers; l++) {
          tracers(l,      kk,hs+j,hs+i) = tracers(l,hs     ,hs+j,hs+i);
          tracers(l,hs+nz+kk,hs+j,hs+i) = tracers(l,hs+nz-1,hs+j,hs+i);
        }
      });
    }

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // We need these to persist to evolve tracers with ADER
      SArray<real,2,nAder,ngll> r_DTs , rw_DTs;

      { // BEGIN: reconstruct, time-avg, and store state & state fluxes
        ////////////////////////////////////////////////////////////////
        // Reconstruct rho, u, v, w, theta
        ////////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> ru_DTs , rv_DTs , rt_DTs;
        {
          SArray<real,1,ord> stencil;

          // Density
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idR,k+kk,hs+j,hs+i); }
          reconstruct_gll_values( stencil , r_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int kk=0; kk < ngll; kk++) { r_DTs(0,kk) += hyDensGLL(k,kk); } // Add hydrostasis back on

          // u values and derivatives
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idU,k+kk,hs+j,hs+i); }
          reconstruct_gll_values( stencil , ru_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // v
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idV,k+kk,hs+j,hs+i); }
          reconstruct_gll_values( stencil , rv_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // w
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idW,k+kk,hs+j,hs+i); }
          reconstruct_gll_values( stencil , rw_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );
          if (bc_z == BC_WALL) {
            if (k == nz-1) rw_DTs(0,ngll-1) = 0;
            if (k == 0   ) rw_DTs(0,0     ) = 0;
          }

          // theta
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idT,k+kk,hs+j,hs+i); }
          reconstruct_gll_values( stencil , rt_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int kk=0; kk < ngll; kk++) { rt_DTs(0,kk) += hyDensThetaGLL(k,kk); } // Add hydrostasis back on
        }

        ///////////////////////////////////////////////////////////////
        // Compute other values needed for centered tendencies and DTs
        ///////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> rwu_DTs , rwv_DTs , rww_DTs , rwt_DTs , rt_gamma_DTs;
        for (int kk=0; kk < ngll; kk++) {
          real r = r_DTs (0,kk);
          real u = ru_DTs(0,kk) / r;
          real v = rv_DTs(0,kk) / r;
          real w = rw_DTs(0,kk) / r;
          real t = rt_DTs(0,kk) / r;
          rwu_DTs    (0,kk) = r*w*u;
          rwv_DTs    (0,kk) = r*w*v;
          rww_DTs    (0,kk) = r*w*w;
          rwt_DTs    (0,kk) = r*w*t;
          rt_gamma_DTs(0,kk) = pow(r*t,GAMMA);
        }

        //////////////////////////////////////////
        // Compute time derivatives if necessary
        //////////////////////////////////////////
        if (nAder > 1) {
          diffTransformEulerConsZ( r_DTs , ru_DTs , rv_DTs , rw_DTs , rt_DTs , rwu_DTs , rwv_DTs , rww_DTs ,
                                   rwt_DTs , rt_gamma_DTs , derivMatrix , hyPressureGLL , k , dz , bc_z , nz );
        }

        //////////////////////////////////////////
        // Time average if necessary
        //////////////////////////////////////////
        // We can't alter density and momentum because they're needed for tracers later
        SArray<real,1,ngll> r_tavg, rw_tavg;
        if (timeAvg) {
          compute_timeAvg( r_DTs  , r_tavg  , dt );
          compute_timeAvg( ru_DTs           , dt );
          compute_timeAvg( rv_DTs           , dt );
          compute_timeAvg( rw_DTs , rw_tavg , dt );
          compute_timeAvg( rt_DTs           , dt );
          compute_timeAvg( rwu_DTs          , dt );
          compute_timeAvg( rwv_DTs          , dt );
          compute_timeAvg( rww_DTs          , dt );
          compute_timeAvg( rwt_DTs          , dt );
          compute_timeAvg( rt_gamma_DTs     , dt );
        } else {
          for (int ii=0; ii < ngll; ii++) {
            r_tavg (ii) = r_DTs (0,ii);
            rw_tavg(ii) = rw_DTs(0,ii);
          }
        }

        //////////////////////////////////////////
        // Store cell edge estimates of the state
        //////////////////////////////////////////
        // Left interface
        stateLimits(idR,1,k  ,j,i) = r_tavg  (0     );
        stateLimits(idU,1,k  ,j,i) = ru_DTs(0,0     );
        stateLimits(idV,1,k  ,j,i) = rv_DTs(0,0     );
        stateLimits(idW,1,k  ,j,i) = rw_tavg (0     );
        stateLimits(idT,1,k  ,j,i) = rt_DTs(0,0     );
        // Right interface       
        stateLimits(idR,0,k+1,j,i) = r_tavg  (ngll-1);
        stateLimits(idU,0,k+1,j,i) = ru_DTs(0,ngll-1);
        stateLimits(idV,0,k+1,j,i) = rv_DTs(0,ngll-1);
        stateLimits(idW,0,k+1,j,i) = rw_tavg (ngll-1);
        stateLimits(idT,0,k+1,j,i) = rt_DTs(0,ngll-1);

        //////////////////////////////////////////
        // Store cell edge estimates of the flux
        //////////////////////////////////////////
        // Left interface
        stateFluxLimits(idR,1,k  ,j,i) = rw_tavg  (0     );
        stateFluxLimits(idU,1,k  ,j,i) = rwu_DTs(0,0     );
        stateFluxLimits(idV,1,k  ,j,i) = rwv_DTs(0,0     );
        stateFluxLimits(idW,1,k  ,j,i) = rww_DTs(0,0     ) + C0*rt_gamma_DTs(0,0     ) - hyPressureGLL(k,0     );
        stateFluxLimits(idT,1,k  ,j,i) = rwt_DTs(0,0     );
        // Right interface      
        stateFluxLimits(idR,0,k+1,j,i) = rw_tavg  (ngll-1);
        stateFluxLimits(idU,0,k+1,j,i) = rwu_DTs(0,ngll-1);
        stateFluxLimits(idV,0,k+1,j,i) = rwv_DTs(0,ngll-1);
        stateFluxLimits(idW,0,k+1,j,i) = rww_DTs(0,ngll-1) + C0*rt_gamma_DTs(0,ngll-1) - hyPressureGLL(k,ngll-1);
        stateFluxLimits(idT,0,k+1,j,i) = rwt_DTs(0,ngll-1);

        ////////////////////////////////////////////
        // Assign gravity source term
        ////////////////////////////////////////////
        real ravg = 0;
        for (int kk=0; kk < ngll; kk++) {
          ravg += (r_tavg(kk) - hyDensGLL(k,kk)) * gllWts_ngll(kk);
        }
        stateTend(idR,k,j,i) = 0;
        stateTend(idU,k,j,i) = 0;
        stateTend(idV,k,j,i) = 0;
        stateTend(idW,k,j,i) = -GRAV*ravg;
        stateTend(idT,k,j,i) = 0;
      } // END: reconstruct, time-avg, and store state & state fluxes

      // r_DTs and rw_DTs still exist and are computed
      { // BEGIN: Reconstruct, time-average, and store tracer fluxes
        // Only process one tracer at a time to save on local memory / register requirements
        for (int tr=0; tr < numTracers; tr++) {
          SArray<real,2,nAder,ngll> rt_DTs; // Density * tracer
          { // BEGIN: Reconstruct the tracer
            SArray<real,1,ord> stencil;
            for (int kk=0; kk < ord; kk++) { stencil(kk) = tracers(tr,k+kk,hs+j,hs+i); }
            reconstruct_gll_values( stencil , rt_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            for (int kk=0; kk < ngll; kk++) { rt_DTs(0,kk) *= r_DTs(0,kk); }
            if (tracerPos(tr)) {
              for (int kk=0; kk < ngll; kk++) { rt_DTs(0,kk) = max( 0._fp , rt_DTs(0,kk) ); }
            }
          } // END: Reconstruct the tracer

          // Compute the tracer flux
          SArray<real,2,nAder,ngll> rwt_DTs; // Density * wwind * tracer
          for (int kk=0; kk < ngll; kk++) {
            rwt_DTs(0,kk) = rt_DTs(0,kk) * rw_DTs(0,kk) / r_DTs(0,kk);
          }

          //////////////////////////////////////////
          // Compute time derivatives if necessary
          //////////////////////////////////////////
          if (nAder > 1) {
            diffTransformTracer( r_DTs , rw_DTs , rt_DTs , rwt_DTs , derivMatrix , dz );
          }

          //////////////////////////////////////////
          // Time average if necessary
          //////////////////////////////////////////
          if (timeAvg) {
            compute_timeAvg( rt_DTs  , dt );
            compute_timeAvg( rwt_DTs , dt );
          }
          if (tracerPos(tr)) {
            for (int kk=0; kk < ngll; kk++) { rt_DTs(0,kk) = max( 0._fp , rt_DTs(0,kk) ); }
          }
          if (bc_z == BC_WALL) {
            if (k == nz-1) rwt_DTs(0,ngll-1) = 0;
            if (k == 0   ) rwt_DTs(0,0     ) = 0;
          }

          ////////////////////////////////////////////////////////////
          // Store cell edge estimates of the tracer
          ////////////////////////////////////////////////////////////
          tracerLimits(tr,1,k  ,j,i) = rt_DTs (0,0     ); // Left interface
          tracerLimits(tr,0,k+1,j,i) = rt_DTs (0,ngll-1); // Right interface
        }
      } // END: Reconstruct, time-average, and store tracer fluxes

    });

    ////////////////////////////////////////////////
    // BCs for the state edge estimates
    ////////////////////////////////////////////////
    parallel_for( Bounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) {
      for (int l = 0; l < numState; l++) {
        if        (bc_z == BC_PERIODIC) {
          stateLimits     (l,0,0 ,j,i) = stateLimits     (l,0,nz,j,i);
          stateFluxLimits (l,0,0 ,j,i) = stateFluxLimits (l,0,nz,j,i);
          stateLimits     (l,1,nz,j,i) = stateLimits     (l,1,0 ,j,i);
          stateFluxLimits (l,1,nz,j,i) = stateFluxLimits (l,1,0 ,j,i);
        } else if (bc_z == BC_WALL    ) {
          stateLimits     (l,0,0 ,j,i) = stateLimits     (l,1,0 ,j,i);
          stateFluxLimits (l,0,0 ,j,i) = stateFluxLimits (l,1,0 ,j,i);
          stateLimits     (l,1,nz,j,i) = stateLimits     (l,0,nz,j,i);
          stateFluxLimits (l,1,nz,j,i) = stateFluxLimits (l,0,nz,j,i);
        }
      }
      for (int l = 0; l < numTracers; l++) {
        if        (bc_z == BC_PERIODIC) {
          tracerLimits(l,0,0 ,j,i) = tracerLimits(l,0,nz,j,i);
          tracerLimits(l,1,nz,j,i) = tracerLimits(l,1,0 ,j,i);
        } else if (bc_z == BC_WALL    ) {
          tracerLimits(l,0,0 ,j,i) = tracerLimits(l,1,0 ,j,i);
          tracerLimits(l,1,nz,j,i) = tracerLimits(l,0,nz,j,i);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<3>(nz+1,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // Get left and right state
      real r_L = stateLimits(idR,0,k,j,i)    ;   real r_R = stateLimits(idR,1,k,j,i)    ;
      real u_L = stateLimits(idU,0,k,j,i)/r_L;   real u_R = stateLimits(idU,1,k,j,i)/r_R;
      real v_L = stateLimits(idV,0,k,j,i)/r_L;   real v_R = stateLimits(idV,1,k,j,i)/r_R;
      real w_L = stateLimits(idW,0,k,j,i)/r_L;   real w_R = stateLimits(idW,1,k,j,i)/r_R;
      real t_L = stateLimits(idT,0,k,j,i)/r_L;   real t_R = stateLimits(idT,1,k,j,i)/r_R;
      // Compute average state
      real r = 0.5_fp * (r_L + r_R);
      real u = 0.5_fp * (u_L + u_R);
      real v = 0.5_fp * (v_L + v_R);
      real w = 0.5_fp * (w_L + w_R);
      real t = 0.5_fp * (t_L + t_R);
      real p = C0 * pow(r*t,GAMMA);
      real cs2 = GAMMA*p/r;
      real cs  = sqrt(cs2);
      // Get left and right fluxes
      real f1_L = stateFluxLimits(idR,0,k,j,i);   real f1_R = stateFluxLimits(idR,1,k,j,i);
      real f2_L = stateFluxLimits(idU,0,k,j,i);   real f2_R = stateFluxLimits(idU,1,k,j,i);
      real f3_L = stateFluxLimits(idV,0,k,j,i);   real f3_R = stateFluxLimits(idV,1,k,j,i);
      real f4_L = stateFluxLimits(idW,0,k,j,i);   real f4_R = stateFluxLimits(idW,1,k,j,i);
      real f5_L = stateFluxLimits(idT,0,k,j,i);   real f5_R = stateFluxLimits(idT,1,k,j,i);
      // Compute upwind characteristics
      // Waves 1-3, velocity: w
      real w1, w2, w3;
      if (w > 0) {
        w1 = f1_L - f5_L/t;
        w2 = f2_L - u*f5_L/t;
        w3 = f3_L - v*f5_L/t;
      } else {
        w1 = f1_R - f5_R/t;
        w2 = f2_R - u*f5_R/t;
        w3 = f3_R - v*f5_R/t;
      }
      // Wave 5, velocity: w-cs
      real w5 =  w*f1_R/(2*cs) - f4_R/(2*cs) + f5_R/(2*t);
      // Wave 6, velocity: w+cs
      real w6 = -w*f1_L/(2*cs) + f4_L/(2*cs) + f5_L/(2*t);
      // Use right eigenmatrix to compute upwind flux
      stateFluxLimits(idR,0,k,j,i) = w1 + w5 + w6;
      stateFluxLimits(idU,0,k,j,i) = w2 + u*w5 + u*w6;
      stateFluxLimits(idV,0,k,j,i) = w3 + v*w5 + v*w6;
      stateFluxLimits(idW,0,k,j,i) = w*w1 + (w-cs)*w5 + (w+cs)*w6;

      real massFlux = stateFluxLimits(idR,0,k,j,i);

      if (w > 0) {
        stateFluxLimits(idT,0,k,j,i) = massFlux * stateLimits(idT,0,k,j,i) / r_L;
      } else {
        stateFluxLimits(idT,0,k,j,i) = massFlux * stateLimits(idT,1,k,j,i) / r_R;
      }

      // COMPUTE UPWIND TRACER FLUXES
      // Handle it one tracer at a time
      for (int tr=0; tr < numTracers; tr++) {
        if (w > 0) {
          tracerFlux(tr,k,j,i) = massFlux * tracerLimits(tr,0,k,j,i) / r_L;
        } else {
          tracerFlux(tr,k,j,i) = massFlux * tracerLimits(tr,1,k,j,i) / r_R;
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Limit the tracer fluxes for positivity
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<4>(numTracers,nz+1,ny,nx) , YAKL_LAMBDA (int tr, int k, int j, int i) {
      real constexpr eps = 1.e-10;
      real w = 0.5_fp * ( stateLimits(idW,0,k,j,i) + stateLimits(idW,1,k,j,i) );
      // Solid wall BCs mean w == 0 at boundaries
      if (tracerPos(tr)) {
        // Compute and apply the flux reduction factor of the upwind cell
        if      (w > 0) {
          int ind_k = k-1;
          // upwind is to the left of this interface
          real f1 = min( tracerFlux(tr,ind_k  ,j,i) , 0._fp );
          real f2 = max( tracerFlux(tr,ind_k+1,j,i) , 0._fp );
          real fluxOut = dt*(f2-f1)/dz;
          real dens = state(idR,hs+ind_k,hs+j,hs+i) + hyDensCells(hs+ind_k);
          tracerFlux(tr,k,j,i) *= min( 1._fp , tracers(tr,hs+ind_k,hs+j,hs+i) * dens / (fluxOut + eps) );
        } else if (w < 0) {
          int ind_k = k;
          // upwind is to the right of this interface
          real f1 = min( tracerFlux(tr,ind_k  ,j,i) , 0._fp );
          real f2 = max( tracerFlux(tr,ind_k+1,j,i) , 0._fp );
          real fluxOut = dt*(f2-f1)/dz;
          real dens = state(idR,hs+ind_k,hs+j,hs+i) + hyDensCells(hs+ind_k);
          tracerFlux(tr,k,j,i) *= min( 1._fp , tracers(tr,hs+ind_k,hs+j,hs+i) * dens / (fluxOut + eps) );
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA(int k, int j, int i) {
      for (int l=0; l < numState; l++) {
        if (sim2d && l == idV) {
          stateTend(l,k,j,i) = 0;
        } else {
          stateTend(l,k,j,i) += - ( stateFluxLimits(l,0,k+1,j,i) - stateFluxLimits(l,0,k,j,i) ) / dz;
        }
      }
      for (int l=0; l < numTracers; l++) {
        // Compute tracer tendency
        tracerTend(l,k,j,i) = - ( tracerFlux(l,k+1,j,i) - tracerFlux(l,k,j,i) ) / dz;
        // Multiply density back onto the tracers
        tracers(l,hs+k,hs+j,hs+i) *= (state(idR,hs+k,hs+j,hs+i) + hyDensCells(hs+k));
      }
    });
  }



  template <class F> void applyStateTendencies( F const &applySingleTendency , int splitIndex ) {
    parallel_for( Bounds<4>(numState,nz,ny,nx) , YAKL_LAMBDA (int l, int k, int j, int i) {
      applySingleTendency({l,k,j,i});
    });
  }



  template <class F> void applyTracerTendencies( F const &applySingleTendency , int splitIndex ) {
    parallel_for( Bounds<4>(numTracers,nz,ny,nx) , YAKL_LAMBDA (int l, int k, int j, int i) {
      applySingleTendency({l,k,j,i});
    });
  }



  const char * getName() { return "1-D Uniform Transport with upwind FV on A-grid"; }



  void output(StateArr const &state, TracerArr const &tracers, real etime) const {
    auto &dx                    = this->dx                   ;
    auto &dy                    = this->dy                   ;
    auto &dz                    = this->dz                   ;
    auto &hyDensCells           = this->hyDensCells          ;
    auto &hyThetaCells          = this->hyThetaCells         ;
    auto &hyDensThetaCells      = this->hyDensThetaCells     ;
    auto &hyPressureCells       = this->hyPressureCells      ;

    yakl::SimpleNetCDF nc;
    int ulIndex = 0; // Unlimited dimension index to place this data at
    // Create or open the file
    if (etime == 0.) {
      nc.create(outFile);
      // x-coordinate
      real1d xloc("xloc",nx);
      parallel_for( nx , YAKL_LAMBDA (int i) { xloc(i) = (i+0.5)*dx; });
      nc.write(xloc.createHostCopy(),"x",{"x"});
      // y-coordinate
      real1d yloc("yloc",ny);
      parallel_for( ny , YAKL_LAMBDA (int i) { yloc(i) = (i+0.5)*dy; });
      nc.write(yloc.createHostCopy(),"y",{"y"});
      // z-coordinate
      real1d zloc("zloc",nz);
      parallel_for( nz , YAKL_LAMBDA (int i) { zloc(i) = (i+0.5)*dz; });
      nc.write(zloc.createHostCopy(),"z",{"z"});
      // hydrostatic density, theta, and pressure
      parallel_for( nz , YAKL_LAMBDA (int k) { zloc(k) = hyDensCells(hs+k); });
      nc.write(zloc.createHostCopy(),"hyDens"    ,{"z"});
      parallel_for( nz , YAKL_LAMBDA (int k) { zloc(k) = hyPressureCells(hs+k); });
      nc.write(zloc.createHostCopy(),"hyPressure",{"z"});
      parallel_for( nz , YAKL_LAMBDA (int k) { zloc(k) = hyThetaCells(hs+k); });
      nc.write(zloc.createHostCopy(),"hyTheta"   ,{"z"});
      // Create time variable
      nc.write1(0._fp,"t",0,"t");
    } else {
      nc.open(outFile,yakl::NETCDF_MODE_WRITE);
      ulIndex = nc.getDimSize("t");
      // Write the elapsed time
      nc.write1(etime,"t",ulIndex,"t");
    }
    real3d data("data",nz,ny,nx);
    // rho'
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idR,hs+k,hs+j,hs+i); });
    nc.write1(data.createHostCopy(),"dens_pert",{"z","y","x"},ulIndex,"t");
    // u
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idU,hs+k,hs+j,hs+i); });
    nc.write1(data.createHostCopy(),"u",{"z","y","x"},ulIndex,"t");
    // v
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idV,hs+k,hs+j,hs+i); });
    nc.write1(data.createHostCopy(),"v",{"z","y","x"},ulIndex,"t");
    // w
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idW,hs+k,hs+j,hs+i); });
    nc.write1(data.createHostCopy(),"w",{"z","y","x"},ulIndex,"t");
    // theta'
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      real r =   state(idR,hs+k,hs+j,hs+i) + hyDensCells     (hs+k);
      real t = ( state(idT,hs+k,hs+j,hs+i) + hyDensThetaCells(hs+k) ) / r;
      data(k,j,i) = t - hyThetaCells(hs+k);
    });
    nc.write1(data.createHostCopy(),"pot_temp_pert",{"z","y","x"},ulIndex,"t");
    // pressure'
    parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      real r =   state(idR,hs+k,hs+j,hs+i) + hyDensCells     (hs+k);
      real t = ( state(idT,hs+k,hs+j,hs+i) + hyDensThetaCells(hs+k) ) / r;
      real p = C0*pow(r*t,GAMMA);
      data(k,j,i) = p - hyPressureCells(hs+k);
    });
    nc.write1(data.createHostCopy(),"pressure_pert",{"z","y","x"},ulIndex,"t");

    for (int tr=0; tr < numTracers; tr++) {
      parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        real r = state(idR,hs+k,hs+j,hs+i) + hyDensCells(hs+k);
        data(k,j,i) = tracers(tr,hs+k,hs+j,hs+i)/r;
      });
      nc.write1(data.createHostCopy(),std::string("tracer_")+tracerName[tr],{"z","y","x"},ulIndex,"t");
    }

    // Close the file
    nc.close();
  }



  void finalize(StateArr const &state , TracerArr const &tracers) {}



  // ord stencil values to ngll GLL values; store in DTs
  YAKL_INLINE void reconstruct_gll_values( SArray<real,1,ord> const stencil ,
                                           SArray<real,2,nAder,ngll> &DTs ,
                                           SArray<real,2,ord,ngll> const &coefs_to_gll ,
                                           SArray<real,2,ord,ngll> const &sten_to_gll  ,
                                           SArray<real,3,ord,ord,ord> const &wenoRecon ,
                                           SArray<real,1,hs+2> const &idl              ,
                                           real sigma, bool doweno ) {
    if (doweno) {

      // Reconstruct values
      SArray<real,1,ord> wenoCoefs;
      weno::compute_weno_coefs( wenoRecon , stencil , wenoCoefs , idl , sigma );
      // Transform ord weno coefficients into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real tmp = 0;
        for (int s=0; s < ord; s++) {
          tmp += coefs_to_gll(s,ii) * wenoCoefs(s);
        }
        DTs(0,ii) = tmp;
      }

    } else {

      // Transform ord stencil cell averages into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real tmp = 0;
        for (int s=0; s < ord; s++) {
          tmp += sten_to_gll(s,ii) * stencil(s);
        }
        DTs(0,ii) = tmp;
      }

    } // if doweno
  }



  YAKL_INLINE void diffTransformEulerConsX( SArray<real,2,nAder,ngll> &r  ,
                                            SArray<real,2,nAder,ngll> &ru ,
                                            SArray<real,2,nAder,ngll> &rv ,
                                            SArray<real,2,nAder,ngll> &rw ,
                                            SArray<real,2,nAder,ngll> &rt ,
                                            SArray<real,2,nAder,ngll> &ruu ,
                                            SArray<real,2,nAder,ngll> &ruv ,
                                            SArray<real,2,nAder,ngll> &ruw ,
                                            SArray<real,2,nAder,ngll> &rut ,
                                            SArray<real,2,nAder,ngll> &rt_gamma ,
                                            SArray<real,2,ngll,ngll> const &deriv ,
                                            real dx ) {
    // zero out the non-linear DTs
    for (int kt=1; kt < nAder; kt++) {
      for (int ii=0; ii < ngll; ii++) {
        ruu     (kt,ii) = 0;
        ruv     (kt,ii) = 0;
        ruw     (kt,ii) = 0;
        rut     (kt,ii) = 0;
        rt_gamma(kt,ii) = 0;
      }
    }

    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the state at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real df1_dx = 0;
        real df2_dx = 0;
        real df3_dx = 0;
        real df4_dx = 0;
        real df5_dx = 0;
        for (int s=0; s<ngll; s++) {
          df1_dx += deriv(s,ii) * ( ru (kt,s) );
          if (kt == 0) { df2_dx += deriv(s,ii) * ( ruu(kt,s) + C0*rt_gamma(kt,s)   ); }
          else         { df2_dx += deriv(s,ii) * ( ruu(kt,s) + C0*rt_gamma(kt,s)/2 ); }
          df3_dx += deriv(s,ii) * ( ruv(kt,s) );
          df4_dx += deriv(s,ii) * ( ruw(kt,s) );
          df5_dx += deriv(s,ii) * ( rut(kt,s) );
        }
        r (kt+1,ii) = -df1_dx/dx/(kt+1._fp);
        ru(kt+1,ii) = -df2_dx/dx/(kt+1._fp);
        rv(kt+1,ii) = -df3_dx/dx/(kt+1._fp);
        rw(kt+1,ii) = -df4_dx/dx/(kt+1._fp);
        rt(kt+1,ii) = -df5_dx/dx/(kt+1._fp);
      }

      // Compute ru* at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real tot_ruu = 0;
        real tot_ruv = 0;
        real tot_ruw = 0;
        real tot_rut = 0;
        for (int ir=0; ir<=kt+1; ir++) {
          tot_ruu += ru(ir,ii) * ru(kt+1-ir,ii) - r(ir,ii) * ruu(kt+1-ir,ii);
          tot_ruv += ru(ir,ii) * rv(kt+1-ir,ii) - r(ir,ii) * ruv(kt+1-ir,ii);
          tot_ruw += ru(ir,ii) * rw(kt+1-ir,ii) - r(ir,ii) * ruw(kt+1-ir,ii);
          tot_rut += ru(ir,ii) * rt(kt+1-ir,ii) - r(ir,ii) * rut(kt+1-ir,ii);
        }
        ruu(kt+1,ii) = tot_ruu / r(0,ii);
        ruv(kt+1,ii) = tot_ruv / r(0,ii);
        ruw(kt+1,ii) = tot_ruw / r(0,ii);
        rut(kt+1,ii) = tot_rut / r(0,ii);

        // Compute rt_gamma at the next time level
        real tot_rt_gamma = 0;
        for (int ir=0; ir<=kt; ir++) {
          tot_rt_gamma += (kt+1._fp -ir) * ( GAMMA*rt_gamma(ir,ii)*rt(kt+1-ir,ii) - rt(ir,ii)*rt_gamma(kt+1-ir,ii) );
        }
        rt_gamma(kt+1,ii) = ( GAMMA*rt_gamma(0,ii)*rt(kt+1,ii) + tot_rt_gamma / (kt+1._fp) ) / rt(0,ii);
      }
    }

    // Fix the rt_gamma
    for (int kt=1; kt < nAder; kt++) {
      for (int ii=0; ii < ngll; ii++) {
        rt_gamma(kt,ii) /= 2;
      }
    }
  }



  YAKL_INLINE void diffTransformEulerConsY( SArray<real,2,nAder,ngll> &r  ,
                                            SArray<real,2,nAder,ngll> &ru ,
                                            SArray<real,2,nAder,ngll> &rv ,
                                            SArray<real,2,nAder,ngll> &rw ,
                                            SArray<real,2,nAder,ngll> &rt ,
                                            SArray<real,2,nAder,ngll> &rvu ,
                                            SArray<real,2,nAder,ngll> &rvv ,
                                            SArray<real,2,nAder,ngll> &rvw ,
                                            SArray<real,2,nAder,ngll> &rvt ,
                                            SArray<real,2,nAder,ngll> &rt_gamma ,
                                            SArray<real,2,ngll,ngll> const &deriv , 
                                            real dy ) {
    // zero out the non-linear DTs
    for (int kt=1; kt < nAder; kt++) {
      for (int ii=0; ii < ngll; ii++) {
        rvu     (kt,ii) = 0;
        rvv     (kt,ii) = 0;
        rvw     (kt,ii) = 0;
        rvt     (kt,ii) = 0;
        rt_gamma(kt,ii) = 0;
      }
    }

    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the state at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real drv_dy    = 0;
        real drvu_dy   = 0;
        real drvv_p_dy = 0;
        real drvw_dy   = 0;
        real drvt_dy   = 0;
        for (int s=0; s<ngll; s++) {
          drv_dy    += deriv(s,ii) * rv(kt,s);
          drvu_dy   += deriv(s,ii) * rvu(kt,s);
          if (kt == 0) { drvv_p_dy += deriv(s,ii) * ( rvv(kt,s) + C0*rt_gamma(kt,s)   ); }
          else         { drvv_p_dy += deriv(s,ii) * ( rvv(kt,s) + C0*rt_gamma(kt,s)/2 ); }
          drvw_dy   += deriv(s,ii) * rvw(kt,s);
          drvt_dy   += deriv(s,ii) * rvt(kt,s);
        }
        r (kt+1,ii) = -drv_dy   /dy/(kt+1);
        ru(kt+1,ii) = -drvu_dy  /dy/(kt+1);
        rv(kt+1,ii) = -drvv_p_dy/dy/(kt+1);
        rw(kt+1,ii) = -drvw_dy  /dy/(kt+1);
        rt(kt+1,ii) = -drvt_dy  /dy/(kt+1);
      }

      // Compute ru* at the next time level
      for (int ii=0; ii<ngll; ii++) {
        // Compute the non-linear differential transforms
        real tot_rvu = 0;
        real tot_rvv = 0;
        real tot_rvw = 0;
        real tot_rvt = 0;
        for (int l=0; l<=kt+1; l++) {
          tot_rvu += rv(l,ii) * ru(kt+1-l,ii) - r(l,ii) * rvu(kt+1-l,ii);
          tot_rvv += rv(l,ii) * rv(kt+1-l,ii) - r(l,ii) * rvv(kt+1-l,ii);
          tot_rvw += rv(l,ii) * rw(kt+1-l,ii) - r(l,ii) * rvw(kt+1-l,ii);
          tot_rvt += rv(l,ii) * rt(kt+1-l,ii) - r(l,ii) * rvt(kt+1-l,ii);
        }
        rvu(kt+1,ii) = tot_rvu / r(0,ii);
        rvv(kt+1,ii) = tot_rvv / r(0,ii);
        rvw(kt+1,ii) = tot_rvw / r(0,ii);
        rvt(kt+1,ii) = tot_rvt / r(0,ii);

        // Compute rt_gamma at the next time level
        real tot_rt_gamma = 0;
        for (int l=0; l<=kt; l++) {
          tot_rt_gamma += (kt+1._fp -l) * ( GAMMA*rt_gamma(l,ii)*rt(kt+1-l,ii) - rt(l,ii)*rt_gamma(kt+1-l,ii) );
        }
        rt_gamma(kt+1,ii) = ( GAMMA*rt_gamma(0,ii)*rt(kt+1,ii) + tot_rt_gamma / (kt+1._fp) ) / rt(0,ii);
      }
    }

    // Fix the rt_gamma
    for (int kt=1; kt < nAder; kt++) {
      for (int ii=0; ii < ngll; ii++) {
        rt_gamma(kt,ii) /= 2;
      }
    }
  }



  YAKL_INLINE void diffTransformEulerConsZ( SArray<real,2,nAder,ngll> &r  ,
                                            SArray<real,2,nAder,ngll> &ru ,
                                            SArray<real,2,nAder,ngll> &rv ,
                                            SArray<real,2,nAder,ngll> &rw ,
                                            SArray<real,2,nAder,ngll> &rt ,
                                            SArray<real,2,nAder,ngll> &rwu ,
                                            SArray<real,2,nAder,ngll> &rwv ,
                                            SArray<real,2,nAder,ngll> &rww ,
                                            SArray<real,2,nAder,ngll> &rwt ,
                                            SArray<real,2,nAder,ngll> &rt_gamma ,
                                            SArray<real,2,ngll,ngll> const &deriv , 
                                            real2d const &hyPressureGLL , 
                                            int k , real dz , int bc_z , int nz ) {
    // zero out the non-linear DTs
    for (int kt=1; kt < nAder; kt++) {
      for (int ii=0; ii < ngll; ii++) {
        rwu     (kt,ii) = 0;
        rwv     (kt,ii) = 0;
        rww     (kt,ii) = 0;
        rwt     (kt,ii) = 0;
        rt_gamma(kt,ii) = 0;
      }
    }

    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the state at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real drw_dz    = 0;
        real drwu_dz   = 0;
        real drwv_dz   = 0;
        real drww_p_dz = 0;
        real drwt_dz   = 0;
        for (int s=0; s<ngll; s++) {
          drw_dz    += deriv(s,ii) * rw(kt,s);
          drwu_dz   += deriv(s,ii) * rwu(kt,s);
          drwv_dz   += deriv(s,ii) * rwv(kt,s);
          if (kt == 0) { drww_p_dz += deriv(s,ii) * ( rww(kt,s) + C0*rt_gamma(kt,s) - hyPressureGLL(k,s) ); }
          else         { drww_p_dz += deriv(s,ii) * ( rww(kt,s) + C0*rt_gamma(kt,s)/2                    ); }
          drwt_dz   += deriv(s,ii) * rwt(kt,s);
        }
        r (kt+1,ii) = -drw_dz   /dz/(kt+1);
        ru(kt+1,ii) = -drwu_dz  /dz/(kt+1);
        rv(kt+1,ii) = -drwv_dz  /dz/(kt+1);
        rw(kt+1,ii) = -drww_p_dz/dz/(kt+1);
        rt(kt+1,ii) = -drwt_dz  /dz/(kt+1);
        if (bc_z == BC_WALL) {
          if (k == nz-1) rw(kt+1,ngll-1) = 0;
          if (k == 0   ) rw(kt+1,0     ) = 0;
        }
      }

      // Compute ru* at the next time level
      for (int ii=0; ii<ngll; ii++) {
        // Compute the non-linear differential transforms
        real tot_rwu = 0;
        real tot_rwv = 0;
        real tot_rww = 0;
        real tot_rwt = 0;
        for (int l=0; l<=kt+1; l++) {
          tot_rwu += rw(l,ii) * ru(kt+1-l,ii) - r(l,ii) * rwu(kt+1-l,ii);
          tot_rwv += rw(l,ii) * rv(kt+1-l,ii) - r(l,ii) * rwv(kt+1-l,ii);
          tot_rww += rw(l,ii) * rw(kt+1-l,ii) - r(l,ii) * rww(kt+1-l,ii);
          tot_rwt += rw(l,ii) * rt(kt+1-l,ii) - r(l,ii) * rwt(kt+1-l,ii);
        }
        rwu(kt+1,ii) = tot_rwu / r(0,ii);
        rwv(kt+1,ii) = tot_rwv / r(0,ii);
        rww(kt+1,ii) = tot_rww / r(0,ii);
        rwt(kt+1,ii) = tot_rwt / r(0,ii);

        // Compute rt_gamma at the next time level
        real tot_rt_gamma = 0;
        for (int l=0; l<=kt; l++) {
          tot_rt_gamma += (kt+1-l) * ( GAMMA*rt_gamma(l,ii)*rt(kt+1-l,ii) - rt(l,ii)*rt_gamma(kt+1-l,ii) );
        }
        rt_gamma(kt+1,ii) = ( GAMMA*rt_gamma(0,ii)*rt(kt+1,ii) + tot_rt_gamma / (kt+1) ) / rt(0,ii);
      }
    }

    // Fix the rt_gamma
    for (int kt=1; kt < nAder; kt++) {
      for (int ii=0; ii < ngll; ii++) {
        rt_gamma(kt,ii) /= 2;
      }
    }
  }



  YAKL_INLINE void diffTransformTracer( SArray<real,2,nAder,ngll> const &r  ,
                                        SArray<real,2,nAder,ngll> const &ru ,
                                        SArray<real,2,nAder,ngll> &rt ,
                                        SArray<real,2,nAder,ngll> &rut ,
                                        SArray<real,2,ngll,ngll> const &deriv , 
                                        real dx ) {
    // zero out the non-linear DT
    for (int kt=1; kt < nAder; kt++) {
      for (int ii=0; ii < ngll; ii++) {
        rut(kt,ii) = 0;
      }
    }
    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the rho*tracer at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real df_dx = 0;
        for (int s=0; s<ngll; s++) {
          df_dx += deriv(s,ii) * rut(kt,s);
        }
        rt(kt+1,ii) = -df_dx/dx/(kt+1._fp);
      }
      // Compute rut at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real tot_rut = 0;
        for (int ir=0; ir<=kt+1; ir++) {
          tot_rut += ru(ir,ii) * rt(kt+1-ir,ii) - r(ir,ii) * rut(kt+1-ir,ii);
        }
        rut(kt+1,ii) = tot_rut / r(0,ii);
      }
    }
  }



  YAKL_INLINE void compute_timeAvg( SArray<real,3,numState,nAder,ngll> &dts , real dt ) {
    real dtmult = dt;
    for (int kt=1; kt<nAder; kt++) {
      for (int l=0; l<numState; l++) {
        for (int ii=0; ii<ngll; ii++) {
          dts(l,0,ii) += dts(l,kt,ii) * dtmult / (kt+1._fp);
        }
      }
      dtmult *= dt;
    }
  }



  YAKL_INLINE void compute_timeAvg( SArray<real,2,nAder,ngll> &dts , real dt ) {
    real dtmult = dt;
    for (int kt=1; kt<nAder; kt++) {
      for (int ii=0; ii<ngll; ii++) {
        dts(0,ii) += dts(kt,ii) * dtmult / (kt+1._fp);
      }
      dtmult *= dt;
    }
  }



  YAKL_INLINE void compute_timeAvg( SArray<real,2,nAder,ngll> const &dts , SArray<real,1,ngll> &tavg , real dt ) {
    for (int ii=0; ii<ngll; ii++) {
      tavg(ii) = dts(0,ii);
    }
    real dtmult = dt;
    for (int kt=1; kt<nAder; kt++) {
      for (int ii=0; ii<ngll; ii++) {
        tavg(ii) += dts(kt,ii) * dtmult / (kt+1._fp);
      }
      dtmult *= dt;
    }
  }

};


