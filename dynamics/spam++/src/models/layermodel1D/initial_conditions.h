#ifndef _ICS_H_
#define _ICS_H_

#include "common.h"
#include "geometry.h"
#include "topology.h"
#include "variable_sets.h"
#include "params.h"
#include <math.h>


struct gaussian_constants {
real const g = 9.80616;
real const Lx = 5000. * 1000.;
real const H0 = 750.0;
real const dh = 75.0;
real const dv = 50.0;
real const xc = 0.5 * Lx;
real const sigmax = 3./40.*Lx;
real const c = 0.05;
real const a = 1.0/3.0;
real const D = 0.5 * Lx;
};
gaussian_constants gauss_constants;



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
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.data_init_cond = DATA_INIT::GAUSSIAN  ; }
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
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.tracerFCT_init_cond[0] = TRACER_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"doublesquare" ) ) { params.tracerFCT_init_cond[0] = TRACER_INIT::DOUBLESQUARE  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.tracerFCT_init_cond[0] = TRACER_INIT::SQUARE    ; }
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
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.tracerFCT_init_cond[1] = TRACER_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"doublesquare" ) ) { params.tracerFCT_init_cond[1] = TRACER_INIT::DOUBLESQUARE  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.tracerFCT_init_cond[1] = TRACER_INIT::SQUARE    ; }
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
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.tracerFCT_init_cond[2] = TRACER_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"doublesquare" ) ) { params.tracerFCT_init_cond[2] = TRACER_INIT::DOUBLESQUARE  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.tracerFCT_init_cond[2] = TRACER_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized tracerInit " << strTracerFCTInit3 << "\n";
      exit(-1);
    }
    
    
  if (params.data_init_cond == DATA_INIT::GAUSSIAN)
  {
  params.xlen = gauss_constants.Lx;
  params.xc = gauss_constants.Lx/2.;
  params.g = gauss_constants.g;
  }


  
}



// *******   Initial Conditions   ***********//

real YAKL_INLINE gaussian_h(real x)
{
  real xprime1 = gauss_constants.Lx / (M_PI * gauss_constants.sigmax) * sin(M_PI / gauss_constants.Lx * (x - gauss_constants.xc));
  return gauss_constants.H0 + gauss_constants.dh * exp(-0.5 * (xprime1 * xprime1));
}

real YAKL_INLINE gaussian_v(real x)
{
  //real xprime1 = gauss_constants.Lx / (M_PI * gauss_constants.sigmax) * sin(M_PI / gauss_constants.Lx * (x - gauss_constants.xc));
  //return gauss_constants.dv * exp(-0.5 * (xprime1 * xprime1));
  return gauss_constants.dv;
}


real YAKL_INLINE gaussian_S(real x)
{
  real sval = gauss_constants.g * (1. + gauss_constants.c * exp(-((x-gauss_constants.xc)*(x-gauss_constants.xc))/(gauss_constants.a*gauss_constants.a*gauss_constants.D*gauss_constants.D)));
  return sval * gaussian_h(x);
}

//FIX THESE!

// CAN WE GENERALIZE THESE? HOW? NEED TO CALL OUT TO THE CORRECT HEIGHT FUNCTION...
// ALSO NEED TO SET VARIOUS LX/LY SIZES
// MAYBE THIS LATTER IS DOABLE VIA PARAMS?
// SIMILAR WITH GAUSSIANS
// MAYBE WHAT WE DO IS HAVE A TRACER STRUCT, AND THEN SET THE VALUES FOR THE TRACER STRUCTURE ACCORDINGLY?
// IN FACT, I THINK ALL TEST CASES STRUCTURES SHOULD INHERIT FROM A MASTER STRUCTURE THAT HAS LX/LY/ETC ALREADY SET UP?
// SOMETHING LIKE THIS...

real YAKL_INLINE _tracer_square_cent(real x)         {return (x > 0.35*gauss_constants.Lx && x < 0.65*gauss_constants.Lx) ? 0.005 : 0.;}
real YAKL_INLINE _tracer_square_r(real x)         {return (x > 0.6*gauss_constants.Lx && x < 0.9*gauss_constants.Lx) ? 0.005 : 0.;}
real YAKL_INLINE _tracer_square_l(real x)         {return (x > 0.1*gauss_constants.Lx && x < 0.4*gauss_constants.Lx) ? 0.005 : 0.;}
real YAKL_INLINE _tracer_square_lr(real x)         {return _tracer_square_l(x) + _tracer_square_r(x);}

real YAKL_INLINE tracer_square_cent (real x) {return gaussian_h(x) * _tracer_square_cent(x);}
real YAKL_INLINE tracer_square_lr (real x) {return gaussian_h(x) * _tracer_square_lr(x);}

real YAKL_INLINE tracer_gaussian(real x)     {return  gaussian_h(x) * 0.005 * exp(-((x-gauss_constants.xc)*(x-gauss_constants.xc))/(gauss_constants.a*gauss_constants.a*gauss_constants.D*gauss_constants.D));}
//{ return 0.005 *double_vortex_h(x,y) * exp(-100. * pow((x-xc)/Lx,2.)) * exp(-100. * pow((y-yc)/Ly,2.)); }

// ADD MORE ICs HERE!!!

//wavespeed = sqrt(g * H0)
//dt = Constant(get_dt(wavespeed, cval, order, variant, Lx, Ly, nx, ny))


template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<nprog> &progvars, VariableSet<nconst> &constvars, 
Geometry<1, nquadx, nquady, nquadz> &primal_geom, Geometry<1, nquadx, nquady, nquadz> &dual_geom)
{

  if (params.data_init_cond == DATA_INIT::GAUSSIAN)
  {
    std::cout << "IC: gaussian " << "\n";
      dual_geom.set_1form_values(gaussian_h, progvars.fields_arr[DENSVAR], 0);
#ifdef _TSWE
      dual_geom.set_1form_values(gaussian_S, progvars.fields_arr[DENSVAR], 1);
#endif

// HOW DO GENERALIZE THESE?
// WANT TO SCALE TRACER FIELDS BY ACTUAL HEIGHT FIELDS...
// SHOULD BE USABLE FOR ANY IC!
for (int i=0; i<ntracers; i++)
{
if (params.tracer_init_cond[i] == TRACER_INIT::GAUSSIAN) {dual_geom.set_1form_values(tracer_gaussian, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
if (params.tracer_init_cond[i] == TRACER_INIT::SQUARE) {dual_geom.set_1form_values(tracer_square_cent, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
if (params.tracer_init_cond[i] == TRACER_INIT::DOUBLESQUARE) {dual_geom.set_1form_values(tracer_square_lr, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
}
for (int i=0; i<ntracers_fct; i++)
{
if (params.tracerFCT_init_cond[i] == TRACER_INIT::GAUSSIAN) {dual_geom.set_1form_values(tracer_gaussian, progvars.fields_arr[DENSFCTVAR], i+ndensityfct-ntracers_fct);}
if (params.tracerFCT_init_cond[i] == TRACER_INIT::SQUARE) {dual_geom.set_1form_values(tracer_square_cent, progvars.fields_arr[DENSFCTVAR], i+ndensityfct-ntracers_fct);}
if (params.tracerFCT_init_cond[i] == TRACER_INIT::DOUBLESQUARE) {dual_geom.set_1form_values(tracer_square_lr, progvars.fields_arr[DENSFCTVAR], i+ndensityfct-ntracers_fct);}
}

primal_geom.set_1form_values(gaussian_v, progvars.fields_arr[VVAR], 0);

  }
  
}
#endif