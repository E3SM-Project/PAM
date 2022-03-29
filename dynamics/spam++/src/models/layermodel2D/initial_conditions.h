#ifndef _ICS_H_
#define _ICS_H_

#include "common.h"
#include "geometry.h"
#include "topology.h"
#include "variable_sets.h"
#include "params.h"
#include <math.h>


struct dbv_constants {
real const g = 9.80616;
real const Lx = 5000. * 1000.;
real const Ly = 5000. * 1000.;
real const coriolis = 0.00006147;
real const H0 = 750.0;
real const ox = 0.1;
real const oy = 0.1;
real const sigmax = 3./40.*Lx;
real const sigmay = 3./40.*Ly;
real const dh = 75.0;
real const xc1 = (0.5-ox) * Lx;
real const yc1 = (0.5-oy) * Ly;
real const xc2 = (0.5+ox) * Lx;
real const yc2 = (0.5+oy) * Ly;
real const xc = 0.5 * Lx;
real const yc = 0.5 * Ly;
real const c = 0.05;
real const a = 1.0/3.0;
real const D = 0.5 * Lx;
};
dbv_constants dbl_vortex_constants;


struct smallbubble_constants {
real const g = 9.80616;
real const Lx = 1000.;
real const Ly = 1500.;
real const xc = 0.5 * Lx;
real const yc = 0.5 * Ly;
real const theta0 = 300.0;
real const zc = 350.;
real const dss = 0.5;
real const rc = 250.;
real const rh0 = 0.8;
};
smallbubble_constants rb_constants;

// These settings seem to work well
// nx = 100 200
// ny = 150 300
// dt = 0.025 0.0125
// Nsteps = 36000 72000
// Nout = 4000 8000





struct largebubble_constants  {
  real const g = 9.80616;
  real const Lx = 20000.;
  real const Ly = 20000.;
  real const xc = 0.5 * Lx;
  real const yc = 0.5 * Ly;
  real const theta0 = 300.0;
  real const zc = 3000.;
  real const xrad = 2000.;
  real const zrad = 2000.;
  real const amp = 2.0;
  real const Cpv = 1859.;
  real const Cpd = 1003.;
};
largebubble_constants lrb_constants;

// These settings seem to work well
// nx = 200 400
// ny = 200 400
// dt = 0.2 0.1
// Nsteps = 10000 20000
// Nout = 500 1000


struct dw_constants {
real const A = 0.98;
real const v1 = 0.1;
real const v2 = 0.2;
real const p0 = 20.0;
};
dw_constants density_wave_constants;








// SOME SUBROUTINE OR SET OF SUBROUTINES HERE THAT INITIALIZES VARIABLES, HAMILTONIANS AND SETS IC SPECIFIC PARAMETERS

void set_ic_specific_params(std::string inFile, ModelParameters &params)
{
  
  std::string strDataInit = "";
  std::string strTracerInit1 = "";
  std::string strTracerInit2 = "";
  std::string strTracerInit3 = "";
  std::string strTracerFCTInit1 = "";
  std::string strTracerFCTInit2 = "";
  std::string strTracerFCTInit3 = "";
  
  // Read in equals-separated key = value file line by line
  std::ifstream fInStream(inFile);
  std::string line;
  while (std::getline(fInStream, line)) {
    // Remove spaces and tabs from the line
    line.erase (std::remove(line.begin(), line.end(), ' '), line.end());
    line.erase (std::remove(line.begin(), line.end(), '\t'), line.end());

    // If the line isn't empty and doesn't begin with a comment specifier, split it based on the colon
    if (!line.empty() && line.find("//",0) != 0) {
      // Find the colon
      uint splitloc = line.find('=',0);
      // Store the key and value strings
      std::string key   = line.substr(0,splitloc);
      std::string value = line.substr(splitloc+1,line.length()-splitloc);

      // Transform the value into a string stream for convenience
      std::stringstream ssVal(value);

      // Match the key, and store the value
           if ( !strcmp( "dataInit"   , key.c_str() ) ) { ssVal >> strDataInit       ;}
           else if ( !strcmp( "tracerInit1"   , key.c_str() ) ) { ssVal >> strTracerInit1       ;}
           else if ( !strcmp( "tracerInit2"   , key.c_str() ) ) { ssVal >> strTracerInit2       ;}
           else if ( !strcmp( "tracerInit3"   , key.c_str() ) ) { ssVal >> strTracerInit3       ;}
           else if ( !strcmp( "tracerFCTInit1"   , key.c_str() ) ) { ssVal >> strTracerFCTInit1       ;}
           else if ( !strcmp( "tracerFCTInit2"   , key.c_str() ) ) { ssVal >> strTracerFCTInit2       ;}
           else if ( !strcmp( "tracerFCTInit3"   , key.c_str() ) ) { ssVal >> strTracerFCTInit3       ;}
      //else {
      //  std::cout << "Error: key " << key << " not understood in file " << inFile << "\n";
      //}
    }
  }

  if (!strcmp("",strDataInit.c_str())) { std::cout << "Error: key " << "dataInit" << " not set.\n"; exit(-1); }

    size_t splitloc = strDataInit.find("//",0);
    std::string sub_str;
    if (splitloc != std::string::npos){
      sub_str = strDataInit.substr(0,splitloc);
    } else {
      sub_str = strDataInit;
    }
    if      ( !strcmp(sub_str.c_str(),"doublevortex" ) ) { params.data_init_cond = DATA_INIT::DOUBLEVORTEX  ; }
    else if ( !strcmp(sub_str.c_str(),"rb" ) ) { params.data_init_cond = DATA_INIT::RB  ; }
    else if ( !strcmp(sub_str.c_str(),"mrb" ) ) { params.data_init_cond = DATA_INIT::MRB  ; }
    else if ( !strcmp(sub_str.c_str(),"lrb"   ) ) { params.data_init_cond = DATA_INIT::LRB    ; }
    else if ( !strcmp(sub_str.c_str(),"mlrb"   ) ) { params.data_init_cond = DATA_INIT::MLRB    ; }
    else if ( !strcmp(sub_str.c_str(),"densitywave"   ) ) { params.data_init_cond = DATA_INIT::DENSITYWAVE    ; }
    else  {
      std::cout << "Error: unrecognized dataInit " << strDataInit << "\n";
      exit(-1);
    }

    splitloc = strTracerInit1.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strTracerInit1.substr(0,splitloc);
    } else {
      sub_str = strTracerInit1;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.tracer_init_cond[0] = TRACER_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"doublesquare" ) ) { params.tracer_init_cond[0] = TRACER_INIT::DOUBLESQUARE  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.tracer_init_cond[0] = TRACER_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized tracerInit " << strTracerInit1 << "\n";
      exit(-1);
    }

    splitloc = strTracerInit2.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strTracerInit2.substr(0,splitloc);
    } else {
      sub_str = strTracerInit2;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.tracer_init_cond[1] = TRACER_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"doublesquare" ) ) { params.tracer_init_cond[1] = TRACER_INIT::DOUBLESQUARE  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.tracer_init_cond[1] = TRACER_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized tracerInit " << strTracerInit2 << "\n";
      exit(-1);
    }
    
    splitloc = strTracerInit3.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strTracerInit3.substr(0,splitloc);
    } else {
      sub_str = strTracerInit3;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.tracer_init_cond[2] = TRACER_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"doublesquare" ) ) { params.tracer_init_cond[2] = TRACER_INIT::DOUBLESQUARE  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.tracer_init_cond[2] = TRACER_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized tracerInit " << strTracerInit3 << "\n";
      exit(-1);
    }
        
    splitloc = strTracerFCTInit1.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strTracerFCTInit1.substr(0,splitloc);
    } else {
      sub_str = strTracerFCTInit1;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.tracer_init_cond[0+ntracers_nofct] = TRACER_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"doublesquare" ) ) { params.tracer_init_cond[0+ntracers_nofct] = TRACER_INIT::DOUBLESQUARE  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.tracer_init_cond[0+ntracers_nofct] = TRACER_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized tracerInit " << strTracerFCTInit1 << "\n";
      exit(-1);
    }

    splitloc = strTracerFCTInit2.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strTracerFCTInit2.substr(0,splitloc);
    } else {
      sub_str = strTracerFCTInit2;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.tracer_init_cond[1+ntracers_nofct] = TRACER_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"doublesquare" ) ) { params.tracer_init_cond[1+ntracers_nofct] = TRACER_INIT::DOUBLESQUARE  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.tracer_init_cond[1+ntracers_nofct] = TRACER_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized tracerInit " << strTracerFCTInit2 << "\n";
      exit(-1);
    }
    
    splitloc = strTracerFCTInit3.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strTracerFCTInit3.substr(0,splitloc);
    } else {
      sub_str = strTracerFCTInit3;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.tracer_init_cond[2+ntracers_nofct] = TRACER_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"doublesquare" ) ) { params.tracer_init_cond[2+ntracers_nofct] = TRACER_INIT::DOUBLESQUARE  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.tracer_init_cond[2+ntracers_nofct] = TRACER_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized tracerInit " << strTracerFCTInit3 << "\n";
      exit(-1);
    }
    
//THIS SHOULD REALLY BE GENERALIZABLE IE ASK THE RELEVANT CONSTANTS STRUCT FOR LX/XC...
    
  if (params.data_init_cond == DATA_INIT::DOUBLEVORTEX)
  {
  params.xlen = dbl_vortex_constants.Lx;
  params.xc = dbl_vortex_constants.Lx/2.;
  params.ylen = dbl_vortex_constants.Ly;
  params.yc = dbl_vortex_constants.Ly/2.;
  params.g = dbl_vortex_constants.g;
  }

  if (params.data_init_cond == DATA_INIT::RB || params.data_init_cond == DATA_INIT::MRB)
  {
  params.xlen = rb_constants.Lx;
  params.xc = rb_constants.Lx/2.;
  params.ylen = rb_constants.Ly;
  params.yc = rb_constants.Ly/2.;
  params.g = rb_constants.g;
  }
  
  if (params.data_init_cond == DATA_INIT::LRB || params.data_init_cond == DATA_INIT::MLRB)
  {
  params.xlen = lrb_constants.Lx;
  params.xc = lrb_constants.Lx/2.;
  params.ylen = lrb_constants.Ly;
  params.yc = lrb_constants.Ly/2.;
  params.g = lrb_constants.g;
  }

  if (params.data_init_cond == DATA_INIT::DENSITYWAVE)
  {
  params.xlen = 2.0;
  params.xc = 0.0;
  params.ylen = 2.0;
  params.yc = 0.0;
  params.g = 0.0;
  }
  
}



// *******   Initial Conditions   ***********//

real YAKL_INLINE double_vortex_coriolis(real x, real y)
{
  return dbl_vortex_constants.coriolis;
}

real YAKL_INLINE double_vortex_h(real x, real y)
{
  real xprime1 = dbl_vortex_constants.Lx / (M_PI * dbl_vortex_constants.sigmax) * sin(M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc1));
  real yprime1 = dbl_vortex_constants.Ly / (M_PI * dbl_vortex_constants.sigmay) * sin(M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc1));
  real xprime2 = dbl_vortex_constants.Lx / (M_PI * dbl_vortex_constants.sigmax) * sin(M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc2));
  real yprime2 = dbl_vortex_constants.Ly / (M_PI * dbl_vortex_constants.sigmay) * sin(M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc2));
  real xprimeprime1 = dbl_vortex_constants.Lx / (2.0 * M_PI * dbl_vortex_constants.sigmax) * sin(2 * M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc1));
  real yprimeprime1 = dbl_vortex_constants.Ly / (2.0 * M_PI * dbl_vortex_constants.sigmay) * sin(2 * M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc1));
  real xprimeprime2 = dbl_vortex_constants.Lx / (2.0 * M_PI * dbl_vortex_constants.sigmax) * sin(2 * M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc2));
  real yprimeprime2 = dbl_vortex_constants.Ly / (2.0 * M_PI * dbl_vortex_constants.sigmay) * sin(2 * M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc2));

  return dbl_vortex_constants.H0 - dbl_vortex_constants.dh * (exp(-0.5 * (xprime1 * xprime1 + yprime1 * yprime1)) + exp(-0.5 * (xprime2 * xprime2 + yprime2 * yprime2)) - 4. * M_PI * dbl_vortex_constants.sigmax * dbl_vortex_constants.sigmay / dbl_vortex_constants.Lx / dbl_vortex_constants.Ly);
}

vec<2> YAKL_INLINE double_vortex_v(real x, real y) {
vec<2> vvec;

real xprime1 = dbl_vortex_constants.Lx / (M_PI * dbl_vortex_constants.sigmax) * sin(M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc1));
real yprime1 = dbl_vortex_constants.Ly / (M_PI * dbl_vortex_constants.sigmay) * sin(M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc1));
real xprime2 = dbl_vortex_constants.Lx / (M_PI * dbl_vortex_constants.sigmax) * sin(M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc2));
real yprime2 = dbl_vortex_constants.Ly / (M_PI * dbl_vortex_constants.sigmay) * sin(M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc2));
real xprimeprime1 = dbl_vortex_constants.Lx / (2.0 * M_PI * dbl_vortex_constants.sigmax) * sin(2 * M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc1));
real yprimeprime1 = dbl_vortex_constants.Ly / (2.0 * M_PI * dbl_vortex_constants.sigmay) * sin(2 * M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc1));
real xprimeprime2 = dbl_vortex_constants.Lx / (2.0 * M_PI * dbl_vortex_constants.sigmax) * sin(2 * M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc2));
real yprimeprime2 = dbl_vortex_constants.Ly / (2.0 * M_PI * dbl_vortex_constants.sigmay) * sin(2 * M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc2));

vvec.u = - dbl_vortex_constants.g * dbl_vortex_constants.dh / dbl_vortex_constants.coriolis / dbl_vortex_constants.sigmay * (yprimeprime1 * exp(-0.5*(xprime1 * xprime1 + yprime1 * yprime1)) + yprimeprime2 * exp(-0.5*(xprime2 * xprime2 + yprime2 * yprime2)));
vvec.v = dbl_vortex_constants.g * dbl_vortex_constants.dh / dbl_vortex_constants.coriolis / dbl_vortex_constants.sigmax * (xprimeprime1 * exp(-0.5*(xprime1 * xprime1 + yprime1 * yprime1)) + xprimeprime2 * exp(-0.5*(xprime2 * xprime2 + yprime2 * yprime2)));
return vvec;
}

real YAKL_INLINE double_vortex_S(real x, real y)
{
  //real sval = g * (1. + c * sin(2. * M_PI / Lx * (x - xc)) * sin(2. * M_PI / Ly * (y - yc)) * exp(-((x-xc)*(x-xc) + (y-yc)*(y-yc))/(a*a*D*D)));
  real sval = dbl_vortex_constants.g * (1. + dbl_vortex_constants.c * exp(-((x-dbl_vortex_constants.xc)*(x-dbl_vortex_constants.xc) + (y-dbl_vortex_constants.yc)*(y-dbl_vortex_constants.yc))/(dbl_vortex_constants.a*dbl_vortex_constants.a*dbl_vortex_constants.D*dbl_vortex_constants.D)));
  //real sval = g * (1. + c * sin(2. * M_PI / Lx * (x- xc)));
  //real sval = g;
  //real sval = g * (1. + c * ((x > 0.35 * Lx && x < 0.65 * Lx && y > 0.35 * Ly && y < 0.65 * Ly ) ? 1. : 0.));
  return sval * double_vortex_h(x,y);
}

// CAN WE GENERALIZE THESE? HOW? NEED TO CALL OUT TO THE CORRECT HEIGHT FUNCTION...
// ALSO NEED TO SET VARIOUS LX/LY SIZES
// MAYBE THIS LATTER IS DOABLE VIA PARAMS?
// SIMILAR WITH GAUSSIANS
// MAYBE WHAT WE DO IS HAVE A TRACER STRUCT, AND THEN SET THE VALUES FOR THE TRACER STRUCTURE ACCORDINGLY?
// IN FACT, I THINK ALL TEST CASES STRUCTURES SHOULD INHERIT FROM A MASTER STRUCTURE THAT HAS LX/LY/ETC ALREADY SET UP?
// SOMETHING LIKE THIS...

real YAKL_INLINE tracer_square_cent(real x, real y)         {return (x > 0.35*dbl_vortex_constants.Lx && x < 0.65*dbl_vortex_constants.Lx && y > 0.35*dbl_vortex_constants.Ly && y < 0.65*dbl_vortex_constants.Ly                        ) ? 0.005 : 0.;}
real YAKL_INLINE tracer_square_ur(real x, real y)         {return (x > 0.6*dbl_vortex_constants.Lx && x < 0.9*dbl_vortex_constants.Lx && y > 0.6*dbl_vortex_constants.Ly && y < 0.9*dbl_vortex_constants.Ly                        ) ? 0.005 : 0.;}
real YAKL_INLINE tracer_square_ll(real x, real y)         {return (x > 0.1*dbl_vortex_constants.Lx && x < 0.4*dbl_vortex_constants.Lx && y > 0.1*dbl_vortex_constants.Ly && y < 0.4*dbl_vortex_constants.Ly                        ) ? 0.005 : 0.;}
real YAKL_INLINE tracer_square_urpll(real x, real y)         {return tracer_square_ur(x,y) + tracer_square_ll(x,y);}

real YAKL_INLINE double_vortex_tracer_square_cent (real x, real y) {return double_vortex_h(x,y) * tracer_square_cent(x, y);}
real YAKL_INLINE double_vortex_tracer_square_urpll (real x, real y) {return double_vortex_h(x,y) * tracer_square_urpll(x, y);}

real YAKL_INLINE double_vortex_tracer_gaussian(real x, real y)     {return  double_vortex_h(x,y) * 0.005 * exp(-((x-dbl_vortex_constants.xc)*(x-dbl_vortex_constants.xc) + (y-dbl_vortex_constants.yc)*(y-dbl_vortex_constants.yc))/(dbl_vortex_constants.a*dbl_vortex_constants.a*dbl_vortex_constants.D*dbl_vortex_constants.D));}
//{ return 0.005 *double_vortex_h(x,y) * exp(-100. * pow((x-xc)/Lx,2.)) * exp(-100. * pow((y-yc)/Ly,2.)); }

// ADD MORE ICs HERE!!!

//wavespeed = sqrt(g * H0)
//dt = Constant(get_dt(wavespeed, cval, order, variant, Lx, Ly, nx, ny))


// Universal

real YAKL_INLINE isentropic_T(real x, real z, real theta0, real g)
{
  return theta0 - z * g / thermo.cst.Cpd;
}

real YAKL_INLINE isentropic_p(real x, real z, real theta0, real g)
{
  return thermo.cst.pr * pow(isentropic_T(x, z, theta0, g) / theta0, 1./thermo.cst.kappa_d);
}

real YAKL_INLINE isentropic_rho(real x, real z, real theta0, real g) {
  real p = isentropic_p(x, z, theta0, g);
  real T = isentropic_T(x, z, theta0, g);
  real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
  return 1./alpha;
}


real YAKL_INLINE linear_ellipsoid(real x, real z, real x0, real z0, real xrad, real zrad, real amp)
{
  real xn = (x-x0)/xrad;
  real zn = (z-z0)/zrad;
  real dist = sqrt( xn*xn + zn*zn );
  return amp * std::max( 1._fp - dist , 0._fp );  
}

real YAKL_INLINE flat_geop(real x, real z, real g)
{
  return g * z;
}

// Returns saturation vapor pressure
 real YAKL_INLINE saturation_vapor_pressure(real temp) {
   real tc = temp - 273.15;
   return 610.94 * exp( 17.625*tc / (243.04+tc) );
 }
 

// (Moist) Rising Bubble (small-scale)

real YAKL_INLINE rb_entropicvar(real x, real z) {
  real p = isentropic_p(x, z, rb_constants.theta0, rb_constants.g);
  real T = isentropic_T(x, z, rb_constants.theta0, rb_constants.g);
  real r = sqrt((x-rb_constants.xc)*(x-rb_constants.xc) + (z-rb_constants.zc)*(z-rb_constants.zc));
  real dtheta = (r<rb_constants.rc) ? rb_constants.dss/2. * (1. + cos(M_PI * r/rb_constants.rc)) : 0.;
  real dT = dtheta * pow(p/thermo.cst.pr, thermo.cst.kappa_d);
  return thermo.compute_entropic_var(p, T+dT, 0, 0, 0, 0);
}

real YAKL_INLINE rb_rho(real x, real z)
{
  return isentropic_rho(x, z, rb_constants.theta0, rb_constants.g);
}

real YAKL_INLINE rb_entropicdensity(real x, real z) {
  return rb_entropicvar(x,z) * rb_rho(x,z);
}

real YAKL_INLINE rb_rho_acousticbalance(real x, real z) {
  real rho_b = rb_rho(x,z);
  real theta = rb_entropicvar(x,z);
  return rho_b * rb_constants.theta0 / theta;

}

real YAKL_INLINE rb_entropicdensity_acousticbalance(real x, real z) {
  return rb_entropicvar(x,z) * rb_rho_acousticbalance(x,z);
}

real YAKL_INLINE rb_geop(real x, real z)
{
  return flat_geop(x,z,rb_constants.g);
}

// The moist bubble is just the dry bubble with an added moisture perturbation
// There is no attempt at balancing anything, just rho_d and theta_h are in balance!
// This is not quite hydrostatic balance even...


// We assume a formula here for SVP that might not be consistent with the thermodynamics
real YAKL_INLINE mrb_rho_v(real x, real z) {
  real r = sqrt((x-rb_constants.xc)*(x-rb_constants.xc) + (z-rb_constants.zc)*(z-rb_constants.zc));
  real rh = (r<rb_constants.rc) ? rb_constants.rh0 * (1. + cos(M_PI * r/rb_constants.rc)) : 0.;
  real Th = isentropic_T(x, z, rb_constants.theta0, rb_constants.g);
  real svp = saturation_vapor_pressure(Th);
  real pv = svp * rh;
  return pv / (thermo.cst.Rv * Th);
}

real YAKL_INLINE mrb_rho_d(real x, real z) {
  real p = isentropic_p(x, z, rb_constants.theta0, rb_constants.g);
  real T = isentropic_T(x, z, rb_constants.theta0, rb_constants.g);
  real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
  return 1./alpha;
}

real YAKL_INLINE mrb_rho(real x, real z) {
  real rhod = mrb_rho_d(x,z);
  real rhov = mrb_rho_v(x,z);
  return rhod + rhov;
}

real YAKL_INLINE mrb_entropicdensity(real x, real z) {
  real p = isentropic_p(x, z, rb_constants.theta0, rb_constants.g);
  real T = isentropic_T(x, z, rb_constants.theta0, rb_constants.g);
  real r = sqrt((x-rb_constants.xc)*(x-rb_constants.xc) + (z-rb_constants.zc)*(z-rb_constants.zc));
  real dtheta = (r<rb_constants.rc) ? rb_constants.dss/2. * (1. + cos(M_PI * r/rb_constants.rc)) : 0.;
  real dT = dtheta * pow(p/thermo.cst.pr, thermo.cst.kappa_d);
  real theta = thermo.compute_entropic_var(p, T+dT, 1, 0, 0, 0);
  return theta * mrb_rho(x,z);
}


// (Moist) Rising Bubble (large-scale)

real YAKL_INLINE lrb_entropicvar(real x, real z) {
  
  real p = isentropic_p(x, z, lrb_constants.theta0, lrb_constants.g);
  real T0 = isentropic_T(x, z, lrb_constants.theta0, lrb_constants.g);
  real dtheta = linear_ellipsoid(x, z, lrb_constants.xc, lrb_constants.zc, lrb_constants.xrad, lrb_constants.zrad, lrb_constants.amp);
  real dT = dtheta * pow(p/thermo.cst.pr, thermo.cst.kappa_d);
  return thermo.compute_entropic_var(p, T0+dT, 0, 0, 0, 0);

}

real YAKL_INLINE lrb_rho(real x, real z)
{
  return isentropic_rho(x, z, lrb_constants.theta0, lrb_constants.g);
}

real YAKL_INLINE lrb_entropicdensity(real x, real z) {
  return lrb_entropicvar(x,z) * lrb_rho(x,z);
}

real YAKL_INLINE lrb_rho_acousticbalance(real x, real z) {
  real rho_b = lrb_rho(x,z);
  real theta = lrb_entropicvar(x,z);
  return rho_b * lrb_constants.theta0 / theta;

}

real YAKL_INLINE lrb_entropicdensity_acousticbalance(real x, real z) {
  return lrb_entropicvar(x,z) * lrb_rho_acousticbalance(x,z);
}

real YAKL_INLINE lrb_geop(real x, real z)
{
  return flat_geop(x,z,lrb_constants.g);
}


real YAKL_INLINE mlrb_rho(real x, real z) {}
real YAKL_INLINE mlrb_entropicdensity(real x, real z) {}
real YAKL_INLINE mlrb_rho_v(real x, real z) {}


// DENSITY WAVE
real YAKL_INLINE density_wave_rho(real x, real z)
{
  return 1. + density_wave_constants.A * sin(2. * M_PI * (x+z));
}

real YAKL_INLINE density_wave_entropicvar(real x, real z)
{
  real p = density_wave_constants.p0;
  real T = density_wave_constants.p0 / (density_wave_rho(x,z) * thermo.cst.Rd);
  return thermo.compute_entropic_var(p, T, 0, 0, 0, 0);
}

real YAKL_INLINE density_wave_entropicdensity(real x, real z) {
  return density_wave_entropicvar(x,z) * density_wave_rho(x,z);
}

vec<2> YAKL_INLINE density_wave_v(real x, real z) {
vec<2> vvec;

vvec.u = density_wave_constants.v1;
vvec.v = density_wave_constants.v2;
return vvec;
}










template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<nprog> &progvars, VariableSet<nconst> &constvars, 
Geometry<nquadx, nquady, nquadz> &primal_geom, Geometry<nquadx, nquady, nquadz> &dual_geom)
{

  if (params.data_init_cond == DATA_INIT::DOUBLEVORTEX)
  {
    std::cout << "IC: double vortex " << "\n";

      dual_geom.set_2form_values(double_vortex_h, progvars.fields_arr[DENSVAR], 0);
#ifdef _TSWE
      dual_geom.set_2form_values(double_vortex_S, progvars.fields_arr[DENSVAR], 1);
#endif
  primal_geom.set_1form_values(double_vortex_v, progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
  primal_geom.set_2form_values(double_vortex_coriolis, constvars.fields_arr[CORIOLISVAR], 0);

// HOW DO GENERALIZE THESE?
// WANT TO SCALE TRACER FIELDS BY ACTUAL HEIGHT FIELDS...
// SHOULD BE USABLE FOR ANY IC!
for (int i=0; i<ntracers; i++)
{
if (params.tracer_init_cond[i] == TRACER_INIT::GAUSSIAN) {dual_geom.set_2form_values(double_vortex_tracer_gaussian, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
if (params.tracer_init_cond[i] == TRACER_INIT::SQUARE) {dual_geom.set_2form_values(double_vortex_tracer_square_cent, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
if (params.tracer_init_cond[i] == TRACER_INIT::DOUBLESQUARE) {dual_geom.set_2form_values(double_vortex_tracer_square_urpll, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
}

  }

  if (params.data_init_cond == DATA_INIT::RB)
  {
    std::cout << "IC: small rising bubble " << "\n";
    if (params.acoustic_balance)
    {
      std::cout << "acoustically balanced " << "\n";
    dual_geom.set_2form_values(rb_rho_acousticbalance, progvars.fields_arr[DENSVAR], 0);
    dual_geom.set_2form_values(rb_entropicdensity_acousticbalance, progvars.fields_arr[DENSVAR], 1);
    }
    else
    {
      dual_geom.set_2form_values(rb_rho, progvars.fields_arr[DENSVAR], 0);
      dual_geom.set_2form_values(rb_entropicdensity, progvars.fields_arr[DENSVAR], 1);      
    }
    dual_geom.set_2form_values(rb_geop, constvars.fields_arr[HSVAR], 0);
  }

  
  // ADD ACOUSTIC BALANCING HERE
  // FIX THIS UP
  if (params.data_init_cond == DATA_INIT::MRB)
  {
    std::cout << "IC: small moist rising bubble " << "\n";
    
    #ifdef _USE_RHOD_VARIANT
    dual_geom.set_2form_values(mrb_rho_d, progvars.fields_arr[DENSVAR], 0);
    #else
    dual_geom.set_2form_values(mrb_rho, progvars.fields_arr[DENSVAR], 0);    
    #endif
    dual_geom.set_2form_values(mrb_entropicdensity, progvars.fields_arr[DENSVAR], 1);
    dual_geom.set_2form_values(mrb_rho_v, progvars.fields_arr[DENSVAR], 2);
    dual_geom.set_2form_values(rb_geop, constvars.fields_arr[HSVAR], 0);
  }
  
// MERGE WITH MLRB below?
// ie if CE dont set vapor, but if MCE do?
// acoustic balancing/rho change though between dry and moist so maybe not...
  if (params.data_init_cond == DATA_INIT::LRB)
  {
    std::cout << "IC: large rising bubble " << "\n";
    thermo.cst.Cpv = lrb_constants.Cpv;
    thermo.cst.Cpd = lrb_constants.Cpd;
    thermo.cst.Cvd = thermo.cst.Cpd - thermo.cst.Rd;
    thermo.cst.Cvv = thermo.cst.Cpv - thermo.cst.Rv;
    
    if (params.acoustic_balance)
    {
      std::cout << "acoustically balanced " << "\n";
    dual_geom.set_2form_values(lrb_rho_acousticbalance, progvars.fields_arr[DENSVAR], 0);
    dual_geom.set_2form_values(lrb_entropicdensity_acousticbalance, progvars.fields_arr[DENSVAR], 1);
    }
    else
    {
      dual_geom.set_2form_values(lrb_rho, progvars.fields_arr[DENSVAR], 0);
      dual_geom.set_2form_values(lrb_entropicdensity, progvars.fields_arr[DENSVAR], 1);      
    }
    dual_geom.set_2form_values(lrb_geop, constvars.fields_arr[HSVAR], 0);
  }

  
  // ADD ACOUSTIC BALANCING HERE
  // FIX THIS UP
  if (params.data_init_cond == DATA_INIT::MLRB)
  {
    std::cout << "IC: large moist rising bubble " << "\n";
    thermo.cst.Cpv = lrb_constants.Cpv;
    thermo.cst.Cpd = lrb_constants.Cpd;
    thermo.cst.Cvd = thermo.cst.Cpd - thermo.cst.Rd;
    thermo.cst.Cvv = thermo.cst.Cpv - thermo.cst.Rv;
    
    #ifdef _USE_RHOD_VARIANT
    dual_geom.set_2form_values(mlrb_rhod, progvars.fields_arr[DENSVAR], 0);
    #else
    dual_geom.set_2form_values(mlrb_rho, progvars.fields_arr[DENSVAR], 0);    
    #endif
    dual_geom.set_2form_values(mlrb_entropicdensity, progvars.fields_arr[DENSVAR], 1);
    dual_geom.set_2form_values(mlrb_rho_v, progvars.fields_arr[DENSVAR], 2);
    dual_geom.set_2form_values(lrb_geop, constvars.fields_arr[HSVAR], 0);
  }
 

//NEEDS ADJUSTING FOR CEFCT CASE ie all FCT!
//THIS SHOULD PROBABLY BE A COMPILE FLAG?
//IE AVOID CODE DUPLICATION?
  if (params.data_init_cond == DATA_INIT::DENSITYWAVE)
  {
    std::cout << "IC: densitywave " << "\n";
      dual_geom.set_2form_values(density_wave_rho, progvars.fields_arr[DENSVAR], 0);
      dual_geom.set_2form_values(density_wave_entropicdensity, progvars.fields_arr[DENSVAR], 1);      
      primal_geom.set_1form_values(density_wave_v, progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
  }
}
#endif
