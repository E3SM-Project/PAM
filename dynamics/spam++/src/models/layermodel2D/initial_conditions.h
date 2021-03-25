#ifndef _ICS_H_
#define _ICS_H_

#include "common.h"
#include "geometry.h"
#include "topology.h"
#include "variable_sets.h"
#include "params.h"

struct dbv_constants {
real const g = 9.80616;
real const Lx = 5000. * 1000.;
real const Ly = 5000. * 1000.;
real const coriolis = 0.00006147;
real const H0 = 750.0;
real const ox = 0.1;
real const oy = 0.1;
real const sigmax = 3./40.*Lx;
real const sigmay = 3./40.*Lx;
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


struct drybubble_constants {
real const g = 9.80616;
real const Lx = 1000.;
real const Ly = 1500.;
real const xc = 0.5 * Lx;
real const yc = 0.5 * Ly;
real const theta0 = 300.0;
real const zc = 350.;
real const dss = 0.5;
real const rc = 250.;
};
drybubble_constants rb_constants;

struct moistbubble_constants : drybubble_constants {
};
moistbubble_constants mrb_constants;













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
    else if ( !strcmp(sub_str.c_str(),"mrb"   ) ) { params.data_init_cond = DATA_INIT::MRB    ; }
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
    
    
  if (params.data_init_cond == DATA_INIT::DOUBLEVORTEX)
  {
  params.xlen = dbl_vortex_constants.Lx;
  params.xc = dbl_vortex_constants.Lx/2.;
  params.ylen = dbl_vortex_constants.Ly;
  params.yc = dbl_vortex_constants.Ly/2.;
  params.g = dbl_vortex_constants.g;
  }

  if (params.data_init_cond == DATA_INIT::RB)
  {
  params.xlen = rb_constants.Lx;
  params.xc = rb_constants.Lx/2.;
  params.ylen = rb_constants.Ly;
  params.yc = rb_constants.Ly/2.;
  params.g = rb_constants.g;
  }
  
  if (params.data_init_cond == DATA_INIT::MRB)
  {
  params.xlen = mrb_constants.Lx;
  params.xc = mrb_constants.Lx/2.;
  params.ylen = mrb_constants.Ly;
  params.yc = mrb_constants.Ly/2.;
  params.g = mrb_constants.g;
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

// CAN WE GENERALIZE?
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



real YAKL_INLINE rb_T(real x, real z)
{
  return rb_constants.theta0 - z * rb_constants.g / thermo.cst.Cpd;
}

real YAKL_INLINE rb_p(real x, real z)
{
  return thermo.cst.pr * pow(rb_T(x,z) / rb_constants.theta0, 1./thermo.cst.kappa_d);
}

real YAKL_INLINE rb_rho(real x, real z) {
  real p = rb_p(x,z);
  real T = rb_T(x,z);
  real alpha = thermo.compute_alpha(p,T,0,0,0,0);
  return 1./alpha;
}

real YAKL_INLINE rb_entropicdensity(real x, real z) {
  real p = rb_p(x,z);
  real T0 = rb_T(x,z);
  real r = sqrt((x-rb_constants.xc)*(x-rb_constants.xc) + (z-rb_constants.zc)*(z-rb_constants.zc));
  real dtheta = (r<rb_constants.rc) ? rb_constants.dss/2. * (1. + cos(M_PI * r/rb_constants.rc)) : 0.;
  real dT = dtheta * pow(p/thermo.cst.pr,thermo.cst.kappa_d);
// add DT PERTURBATION INTO THIS!
  real entropic_var = thermo.compute_entropic_var(p,T0+dT,0,0,0,0);
  return entropic_var * rb_rho(x,z);
}

real YAKL_INLINE rb_geop(real x, real z)
{
  return rb_rho(x,z) * rb_constants.g * z;
}



real YAKL_INLINE mrb_rho(real x, real y) {}
real YAKL_INLINE mrb_entropicdensity(real x, real y) {}
real YAKL_INLINE mrb_rho_v(real x, real y) {}


template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<nprog> &progvars, VariableSet<nconst> &constvars, 
Geometry<2, nquadx, nquady, nquadz> &primal_geom, Geometry<2, nquadx, nquady, nquadz> &dual_geom)
{

  if (params.data_init_cond == DATA_INIT::DOUBLEVORTEX)
  {
      dual_geom.set_2form_values(double_vortex_h, progvars.fields_arr[DENSVAR], 0);
#ifdef _TSWE
      dual_geom.set_2form_values(double_vortex_S, progvars.fields_arr[DENSVAR], 1);
#endif

// HOW DO GENERALIZE THESE?
// WANT TO SCALE TRACER FIELDS BY ACTUAL HEIGHT FIELDS...
for (int i=0; i<ntracers; i++)
{
if (params.tracer_init_cond[i] == TRACER_INIT::GAUSSIAN) {dual_geom.set_2form_values(double_vortex_tracer_gaussian, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
if (params.tracer_init_cond[i] == TRACER_INIT::SQUARE) {dual_geom.set_2form_values(double_vortex_tracer_square_cent, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
if (params.tracer_init_cond[i] == TRACER_INIT::DOUBLESQUARE) {dual_geom.set_2form_values(double_vortex_tracer_square_urpll, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
}
for (int i=0; i<ntracers_fct; i++)
{
if (params.tracerFCT_init_cond[i] == TRACER_INIT::GAUSSIAN) {dual_geom.set_2form_values(double_vortex_tracer_gaussian, progvars.fields_arr[DENSFCTVAR], i);}
if (params.tracerFCT_init_cond[i] == TRACER_INIT::SQUARE) {dual_geom.set_2form_values(double_vortex_tracer_square_cent, progvars.fields_arr[DENSFCTVAR], i);}
if (params.tracerFCT_init_cond[i] == TRACER_INIT::DOUBLESQUARE) {dual_geom.set_2form_values(double_vortex_tracer_square_urpll, progvars.fields_arr[DENSFCTVAR], i);}
}
      primal_geom.set_1form_values(double_vortex_v, progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
      primal_geom.set_2form_values(double_vortex_coriolis, constvars.fields_arr[CORIOLISVAR], 0);
  }

  if (params.data_init_cond == DATA_INIT::RB)
  {
    dual_geom.set_2form_values(rb_rho, progvars.fields_arr[DENSVAR], 0);
    dual_geom.set_2form_values(rb_entropicdensity, progvars.fields_arr[DENSVAR], 1);
    dual_geom.set_2form_values(rb_geop, constvars.fields_arr[HSVAR], 0);
  }
  
  if (params.data_init_cond == DATA_INIT::MRB)
  {
    #ifdef _MCERHO
    dual_geom.set_2form_values(mrb_rho, progvars.fields_arr[DENSVAR], 0);
    #endif
    #ifdef _MCERHOD
    dual_geom.set_2form_values(mrb_rhod, progvars.fields_arr[DENSVAR], 0);
    #endif
    dual_geom.set_2form_values(mrb_entropicdensity, progvars.fields_arr[DENSVAR], 1);
    dual_geom.set_2form_values(mrb_rho_v, progvars.fields_arr[DENSFCTVAR], 0);
  }
  
}
#endif