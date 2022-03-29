#ifndef _ICS_H_
#define _ICS_H_

#include "common.h"
#include "geometry.h"
#include "topology.h"
#include "variable_sets.h"
#include "params.h"
#include <math.h>

// SOME SUBROUTINE OR SET OF SUBROUTINES HERE THAT INITIALIZES VARIABLES, HAMILTONIANS AND SETS IC SPECIFIC PARAMETERS

void set_ic_specific_params(std::string inFile, ModelParameters &params)
{
  
  std::string strDataInit1 = "";
  std::string strDataInit2 = "";
  std::string strDataInit3 = "";
  std::string strDataFCTInit1 = "";
  std::string strDataFCTInit2 = "";
  std::string strDataFCTInit3 = "";
  std::string strDataQInit1 = "";
  std::string strDataQInit2 = "";
  std::string strDataQInit3 = "";
  std::string strWindInit = "";
  
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
      // Match the key, and store the value
           if ( !strcmp( "dataInit1"   , key.c_str() ) ) { ssVal >> strDataInit1       ;}
      else if ( !strcmp( "dataInit2"   , key.c_str() ) ) { ssVal >> strDataInit2       ;}
      else if ( !strcmp( "dataInit3"   , key.c_str() ) ) { ssVal >> strDataInit3       ;}
      else if ( !strcmp( "dataFCTInit1"   , key.c_str() ) ) { ssVal >> strDataFCTInit1       ;}
      else if ( !strcmp( "dataFCTInit2"   , key.c_str() ) ) { ssVal >> strDataFCTInit2       ;}
      else if ( !strcmp( "dataFCTInit3"   , key.c_str() ) ) { ssVal >> strDataFCTInit3       ;}
      else if ( !strcmp( "dataQInit1"   , key.c_str() ) ) { ssVal >> strDataQInit1       ;}
      else if ( !strcmp( "dataQInit2"   , key.c_str() ) ) { ssVal >> strDataQInit2       ;}
      else if ( !strcmp( "dataQInit3"   , key.c_str() ) ) { ssVal >> strDataQInit3       ;}
      else if ( !strcmp( "windInit"    , key.c_str() ) ) { ssVal >> strWindInit        ;}
      //else {
      //  std::cout << "Error: key " << key << " not understood in file " << inFile << "\n";
      //}
    }
  }

  if (!strcmp("",strDataInit1.c_str())) { std::cout << "Error: key " << "dataInit1" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataInit2.c_str())) { std::cout << "Error: key " << "dataInit2" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataInit3.c_str())) { std::cout << "Error: key " << "dataInit3" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataFCTInit1.c_str())) { std::cout << "Error: key " << "dataFCTInit1" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataFCTInit2.c_str())) { std::cout << "Error: key " << "dataFCTInit2" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataFCTInit3.c_str())) { std::cout << "Error: key " << "dataFCTInit3" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataQInit1.c_str())) { std::cout << "Error: key " << "dataQInit1" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataQInit2.c_str())) { std::cout << "Error: key " << "dataQInit2" << " not set.\n"; exit(-1); }
  if (!strcmp("",strDataQInit3.c_str())) { std::cout << "Error: key " << "dataQInit3" << " not set.\n"; exit(-1); }
  if (!strcmp("",strWindInit.c_str()))  { std::cout << "Error: key " << "windInit"  << " not set.\n"; exit(-1); }
  size_t splitloc = strDataInit1.find("//",0);
  std::string sub_str;
  if (splitloc != std::string::npos){
    sub_str = strDataInit1.substr(0,splitloc);
  } else {
    sub_str = strDataInit1;
  }
  if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.data_init_cond[0] = DATA_INIT::GAUSSIAN  ; }
  else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.data_init_cond[0] = DATA_INIT::VORTICES  ; }
  else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.data_init_cond[0] = DATA_INIT::SQUARE    ; }
  else if ( !strcmp(sub_str.c_str(),"doublesquare"   ) ) { params.data_init_cond[0] = DATA_INIT::DOUBLESQUARE    ; }
  else  {
    std::cout << "Error: unrecognized dataInit " << strDataInit1 << "\n";
    exit(-1);
  }

  splitloc = strDataInit2.find("//",0);
  if (splitloc != std::string::npos){
    sub_str = strDataInit2.substr(0,splitloc);
  } else {
    sub_str = strDataInit2;
  }
  if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.data_init_cond[1] = DATA_INIT::GAUSSIAN  ; }
  else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.data_init_cond[1] = DATA_INIT::VORTICES  ; }
  else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.data_init_cond[1] = DATA_INIT::SQUARE    ; }
  else if ( !strcmp(sub_str.c_str(),"doublesquare"   ) ) { params.data_init_cond[1] = DATA_INIT::DOUBLESQUARE    ; }
  else  {
    std::cout << "Error: unrecognized dataInit " << strDataInit2 << "\n";
    exit(-1);
  }

  splitloc = strDataInit3.find("//",0);
  if (splitloc != std::string::npos){
    sub_str = strDataInit3.substr(0,splitloc);
  } else {
    sub_str = strDataInit3;
  }
  if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.data_init_cond[2] = DATA_INIT::GAUSSIAN  ; }
  else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.data_init_cond[2] = DATA_INIT::VORTICES  ; }
  else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.data_init_cond[2] = DATA_INIT::SQUARE    ; }
  else if ( !strcmp(sub_str.c_str(),"doublesquare"   ) ) { params.data_init_cond[2] = DATA_INIT::DOUBLESQUARE    ; }
  else  {
    std::cout << "Error: unrecognized dataInit " << strDataInit3 << "\n";
    exit(-1);
  }

  splitloc = strDataFCTInit1.find("//",0);
  if (splitloc != std::string::npos){
    sub_str = strDataFCTInit1.substr(0,splitloc);
  } else {
    sub_str = strDataFCTInit1;
  }
  if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.data_init_cond[0+ntnofctdofs] = DATA_INIT::GAUSSIAN  ; }
  else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.data_init_cond[0+ntnofctdofs] = DATA_INIT::VORTICES  ; }
  else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.data_init_cond[0+ntnofctdofs] = DATA_INIT::SQUARE    ; }
  else if ( !strcmp(sub_str.c_str(),"doublesquare"   ) ) { params.data_init_cond[0+ntnofctdofs] = DATA_INIT::DOUBLESQUARE    ; }
  else  {
    std::cout << "Error: unrecognized dataInit " << strDataFCTInit1 << "\n";
    exit(-1);
  }

  splitloc = strDataFCTInit2.find("//",0);
  if (splitloc != std::string::npos){
    sub_str = strDataFCTInit2.substr(0,splitloc);
  } else {
    sub_str = strDataFCTInit2;
  }
  if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.data_init_cond[1+ntnofctdofs] = DATA_INIT::GAUSSIAN  ; }
  else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.data_init_cond[1+ntnofctdofs] = DATA_INIT::VORTICES  ; }
  else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.data_init_cond[1+ntnofctdofs] = DATA_INIT::SQUARE    ; }
  else if ( !strcmp(sub_str.c_str(),"doublesquare"   ) ) { params.data_init_cond[1+ntnofctdofs] = DATA_INIT::DOUBLESQUARE    ; }
  else  {
    std::cout << "Error: unrecognized dataInit " << strDataFCTInit2 << "\n";
    exit(-1);
  }

  splitloc = strDataFCTInit3.find("//",0);
  if (splitloc != std::string::npos){
    sub_str = strDataFCTInit3.substr(0,splitloc);
  } else {
    sub_str = strDataFCTInit3;
  }
  if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.data_init_cond[2+ntnofctdofs] = DATA_INIT::GAUSSIAN  ; }
  else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.data_init_cond[2+ntnofctdofs] = DATA_INIT::VORTICES  ; }
  else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.data_init_cond[2+ntnofctdofs] = DATA_INIT::SQUARE    ; }
  else if ( !strcmp(sub_str.c_str(),"doublesquare"   ) ) { params.data_init_cond[2+ntnofctdofs] = DATA_INIT::DOUBLESQUARE    ; }
  else  {
    std::cout << "Error: unrecognized dataInit " << strDataFCTInit3 << "\n";
    exit(-1);
  }

  splitloc = strDataQInit1.find("//",0);
  if (splitloc != std::string::npos){
    sub_str = strDataQInit1.substr(0,splitloc);
  } else {
    sub_str = strDataQInit1;
  }
  if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.dataQ_init_cond[0] = DATA_INIT::GAUSSIAN  ; }
  else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.dataQ_init_cond[0] = DATA_INIT::VORTICES  ; }
  else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.dataQ_init_cond[0] = DATA_INIT::SQUARE    ; }
  else if ( !strcmp(sub_str.c_str(),"doublesquare"   ) ) { params.dataQ_init_cond[0] = DATA_INIT::DOUBLESQUARE    ; }
  else  {
    std::cout << "Error: unrecognized dataInit " << strDataQInit1 << "\n";
    exit(-1);
  }

  splitloc = strDataQInit2.find("//",0);
  if (splitloc != std::string::npos){
    sub_str = strDataQInit2.substr(0,splitloc);
  } else {
    sub_str = strDataQInit2;
  }
  if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.dataQ_init_cond[1] = DATA_INIT::GAUSSIAN  ; }
  else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.dataQ_init_cond[1] = DATA_INIT::VORTICES  ; }
  else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.dataQ_init_cond[1] = DATA_INIT::SQUARE    ; }
  else if ( !strcmp(sub_str.c_str(),"doublesquare"   ) ) { params.dataQ_init_cond[1] = DATA_INIT::DOUBLESQUARE    ; }
  else  {
    std::cout << "Error: unrecognized dataInit " << strDataQInit2 << "\n";
    exit(-1);
  }

  splitloc = strDataQInit3.find("//",0);
  if (splitloc != std::string::npos){
    sub_str = strDataQInit3.substr(0,splitloc);
  } else {
    sub_str = strDataQInit3;
  }
  if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.dataQ_init_cond[2] = DATA_INIT::GAUSSIAN  ; }
  else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.dataQ_init_cond[2] = DATA_INIT::VORTICES  ; }
  else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.dataQ_init_cond[2] = DATA_INIT::SQUARE    ; }
  else if ( !strcmp(sub_str.c_str(),"doublesquare"   ) ) { params.dataQ_init_cond[2] = DATA_INIT::DOUBLESQUARE    ; }
  else  {
    std::cout << "Error: unrecognized dataInit " << strDataQInit3 << "\n";
    exit(-1);
  }

  splitloc = strWindInit.find("//",0);
  if (splitloc != std::string::npos){
    sub_str = strWindInit.substr(0,splitloc);
  } else {
    sub_str = strWindInit;
  }
  if      ( !strcmp(sub_str.c_str(),"uniform_x"     ) ) { params.wind_init_cond = WIND_INIT::UNIFORM_X     ; }
  else if ( !strcmp(sub_str.c_str(),"uniform_y" ) ) { params.wind_init_cond = WIND_INIT::UNIFORM_Y ; }
  else if ( !strcmp(sub_str.c_str(),"uniform_xy" ) ) { params.wind_init_cond = WIND_INIT::UNIFORM_XY ; }
  else if ( !strcmp(sub_str.c_str(),"deformational" ) ) { params.wind_init_cond = WIND_INIT::DEFORMATIONAL ; }
  else if ( !strcmp(sub_str.c_str(),"doublevortex" ) ) { params.wind_init_cond = WIND_INIT::DOUBLEVORTEX ; }
  else  {
    std::cout << "Error: unrecognized windInit " << strWindInit << "\n";
    exit(-1);
  }

params.etime = 0.0;


params.xlen = 1.0;
params.xc = 0.5;
params.zlen = 1.0;
params.zc = 0.5;


}



// *******   Initial Conditions   ***********//

// Initial conditions related variables and functions
#define gaussian_2d(x,y)   (1. * exp(-100. * pow(x-0.5,2.)) * exp(-100. * pow(y-0.5,2.)))
real YAKL_INLINE gaussian(real x, real y)         { return gaussian_2d(x,y); }

#define vortex1_2d(x,y)   (1. *  exp(-100. * pow(x-0.75,2.)) * exp(-100. * pow(y-0.75,2.)))
#define vortex2_2d(x,y)   (0.5 * exp(-50.  * pow(x-0.25,2.)) * exp(-75.  * pow(y-0.25,2.)))
real YAKL_INLINE vortices(real x, real y)         { return vortex1_2d(x,y)   + vortex2_2d(x,y); }

real YAKL_INLINE square(real x, real y)         {return (x > 0.35 && x < 0.65 && y > 0.35 && y < 0.65                        ) ? 1. : 0.;}

real YAKL_INLINE square_ur(real x, real y)         {return (x > 0.6 && x < 0.9 && y > 0.6 && y < 0.9                        ) ? 1. : 0.;}
real YAKL_INLINE square_ll(real x, real y)         {return (x > 0.1 && x < 0.4 && y > 0.1 && y < 0.4                        ) ? 1. : 0.;}
real YAKL_INLINE doublesquare(real x, real y)         {return square_ur(x,y) + square_ll(x,y);}




#define C_UNIFORM_WIND 1.

vecext<2> YAKL_INLINE uniform_x_wind(real x, real y) {
  vecext<2> vvec;
  vvec.u = C_UNIFORM_WIND;
  return vvec;
}

vecext<2> YAKL_INLINE uniform_y_wind(real x, real y) {
  vecext<2> vvec;
  vvec.w = -C_UNIFORM_WIND;
  return vvec;
}

vecext<2> YAKL_INLINE uniform_xy_wind(real x, real y) {
  vecext<2> vvec;
  vvec.u = -C_UNIFORM_WIND/sqrt(2.);
  vvec.w = C_UNIFORM_WIND/sqrt(2.);
  return vvec;
}

// FIX THIS
vecext<2> YAKL_INLINE deformational_wind(real x, real y) {
  vecext<2> vvec;
  vvec.u = C_UNIFORM_WIND;
  return vvec;
}

vecext<2> YAKL_INLINE doublevortex_wind(real x, real y) {
vecext<2> vvec;

real xprime1 = 1.0 / (M_PI * 3./40.) * sin(M_PI / 1.0 * (x - 0.4));
real yprime1 = 1.0 / (M_PI * 3./40.) * sin(M_PI / 1.0 * (y - 0.4));
real xprime2 = 1.0 / (M_PI * 3./40.) * sin(M_PI / 1.0 * (x - 0.6));
real yprime2 = 1.0 / (M_PI * 3./40.) * sin(M_PI / 1.0 * (y - 0.6));
real xprimeprime1 = 1.0 / (2.0 * M_PI * 3./40.) * sin(2 * M_PI / 1.0 * (x - 0.4));
real yprimeprime1 = 1.0 / (2.0 * M_PI * 3./40.) * sin(2 * M_PI / 1.0 * (y - 0.4));
real xprimeprime2 = 1.0 / (2.0 * M_PI * 3./40.) * sin(2 * M_PI / 1.0 * (x - 0.6));
real yprimeprime2 = 1.0 / (2.0 * M_PI * 3./40.) * sin(2 * M_PI / 1.0 * (y - 0.6));

vvec.u = -1.0 * (yprimeprime1 * exp(-0.5*(xprime1 * xprime1 + yprime1 * yprime1)) + yprimeprime2 * exp(-0.5*(xprime2 * xprime2 + yprime2 * yprime2)));
vvec.w = 1.0 * (xprimeprime1 * exp(-0.5*(xprime1 * xprime1 + yprime1 * yprime1)) + xprimeprime2 * exp(-0.5*(xprime2 * xprime2 + yprime2 * yprime2)));
return vvec;
}

template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<nprog> &progvars, VariableSet<nconst> &constvars, 
Geometry<nquadx, nquady, nquadz> &pgeom, Geometry<nquadx, nquady, nquadz> &dgeom)
{
    for (int i=0; i<ntdofs; i++)
    {
    if (params.data_init_cond[i] == DATA_INIT::GAUSSIAN) {dgeom.set_11form_values(gaussian, progvars.fields_arr[DENSVAR], i);}
    if (params.data_init_cond[i] == DATA_INIT::VORTICES) {dgeom.set_11form_values(vortices, progvars.fields_arr[DENSVAR], i);}
    if (params.data_init_cond[i] == DATA_INIT::SQUARE)   {dgeom.set_11form_values(square,   progvars.fields_arr[DENSVAR], i);}
    if (params.data_init_cond[i] == DATA_INIT::DOUBLESQUARE)   {dgeom.set_11form_values(doublesquare,   progvars.fields_arr[DENSVAR], i);}
    }
    for (int i=0; i<nQdofs; i++)
    {
    if (params.dataQ_init_cond[i] == DATA_INIT::GAUSSIAN) {pgeom.set_11form_values(gaussian, progvars.fields_arr[QXZVAR], i);}
    if (params.dataQ_init_cond[i] == DATA_INIT::VORTICES) {pgeom.set_11form_values(vortices, progvars.fields_arr[QXZVAR], i);}
    if (params.dataQ_init_cond[i] == DATA_INIT::SQUARE)   {pgeom.set_11form_values(square,   progvars.fields_arr[QXZVAR], i);}
    if (params.dataQ_init_cond[i] == DATA_INIT::DOUBLESQUARE)   {pgeom.set_11form_values(doublesquare,   progvars.fields_arr[QXZVAR], i);}
    }
    if (params.wind_init_cond == WIND_INIT::UNIFORM_X    ) {
      pgeom.set_01form_values(uniform_x_wind,     constvars.fields_arr[WVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
      pgeom.set_10form_values(uniform_x_wind,     constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
    }
    if (params.wind_init_cond == WIND_INIT::UNIFORM_Y    ) {
      pgeom.set_01form_values(uniform_y_wind,     constvars.fields_arr[WVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
      pgeom.set_10form_values(uniform_y_wind,     constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
    }
    if (params.wind_init_cond == WIND_INIT::UNIFORM_XY    ) {
      pgeom.set_01form_values(uniform_xy_wind,     constvars.fields_arr[WVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
      pgeom.set_10form_values(uniform_xy_wind,     constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
    }
    if (params.wind_init_cond == WIND_INIT::DEFORMATIONAL) {
      pgeom.set_01form_values(deformational_wind, constvars.fields_arr[WVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
      pgeom.set_10form_values(deformational_wind, constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
    }
    if (params.wind_init_cond == WIND_INIT::DOUBLEVORTEX) {
      pgeom.set_01form_values(doublevortex_wind, constvars.fields_arr[WVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
      pgeom.set_10form_values(doublevortex_wind, constvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
    }

}

#endif