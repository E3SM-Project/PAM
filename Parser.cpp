
#include "Parser.h"


void readParamsFile(std::string fNameIn, Domain &dom, Parallel &par, FileIO &io) {

  std::string strTimeMethod;
  std::string strEqnSet;
  std::string strDataInit;

  // Initialize all REQUIRED read-in values to -999
  dom.nx_glob   = 0;
  dom.ny_glob   = 0;
  dom.nz_glob   = 0;
  dom.simLength = -999;
  par.nproc_x   = -999;
  par.nproc_y   = -999;
  outFreq       = -999;

  // Create some default parameters
  dom.xlen      = 20000;      // 20km default
  dom.ylen      = 20000;      // 20km default
  dom.zlen      = 10000;      // 10km default
  dom.doWeno    = 0;          // WENO off by default
  dom.cfl       = 0.7;        // Allow for any reasonable range of winds
  timeMethod    = TIME_ADER;
  strTimeMethod = "ADER";
  dom.eqnSet    = EQN_THETA_CONS;
  strEqnSet     = "theta_cons";
  dom.dataInit  = DATA_INIT_THERMAL;
  strDataInit   = "thermal";
  strakaVisc    = 0;          // Turn off Straka viscosity by default

  // Read in colon-separated key: value file line by line
  std::ifstream fInStream(fNameIn);
  std::string line;
  while (std::getline(fInStream, line)) {
    // Remove spaces and tabs from the line
    line.erase (std::remove(line.begin(), line.end(), ' '), line.end());
    line.erase (std::remove(line.begin(), line.end(), '\t'), line.end());

    // If the line isn't empty and doesn't begin with a comment specifier, split it based on the colon
    if (!line.empty() && line.find("//",0) != 0) {
      // Find the colon
      uint splitloc = line.find(':',0);
      // Store the key and value strings
      std::string key   = line.substr(0,splitloc);
      std::string value = line.substr(splitloc+1,line.length()-splitloc);

      // Transform the value into a string stream for convenience
      std::stringstream ssVal(value);

      // Match the key, and store the value
      if      ( !strcmp( "nx"        , key.c_str() ) ) { ssVal >> dom.nx_glob  ; }
      else if ( !strcmp( "ny"        , key.c_str() ) ) { ssVal >> dom.ny_glob  ; }
      else if ( !strcmp( "nz"        , key.c_str() ) ) { ssVal >> dom.nz_glob  ; }
      else if ( !strcmp( "xlen"      , key.c_str() ) ) { ssVal >> dom.xlen     ; }
      else if ( !strcmp( "ylen"      , key.c_str() ) ) { ssVal >> dom.ylen     ; }
      else if ( !strcmp( "zlen"      , key.c_str() ) ) { ssVal >> dom.zlen     ; }
      else if ( !strcmp( "cfl"       , key.c_str() ) ) { ssVal >> dom.cfl      ; }
      else if ( !strcmp( "simLength" , key.c_str() ) ) { ssVal >> dom.simLength; }
      else if ( !strcmp( "parNx"     , key.c_str() ) ) { ssVal >> par.nproc_x  ; }
      else if ( !strcmp( "parNy"     , key.c_str() ) ) { ssVal >> par.nproc_y  ; }
      else if ( !strcmp( "outFreq"   , key.c_str() ) ) { ssVal >> outFreq      ; }
      else if ( !strcmp( "doWeno"    , key.c_str() ) ) { ssVal >> dom.doWeno   ; }
      else if ( !strcmp( "strakaVisc", key.c_str() ) ) { ssVal >> strakaVisc   ; }
      else if ( !strcmp( "timeMethod", key.c_str() ) ) { ssVal >> strTimeMethod; handleTimeMethod(strTimeMethod,fNameIn); }
      else if ( !strcmp( "eqnSet"    , key.c_str() ) ) { ssVal >> strEqnSet    ; handleEqnSet    (strEqnSet    ,fNameIn, dom); }
      else if ( !strcmp( "dataInit"  , key.c_str() ) ) { ssVal >> strDataInit  ; handleDataInit  (strDataInit  ,fNameIn, dom); }
      else {
        std::cout << "Error: key " << key << " not understood in file " << fNameIn << "\n";
      }
    }
  }

  // Test to make sure all required values were initialized
  if (dom.nx_glob   == 0   ) { std::cout << "Error: key " << "nx"        << " not set."; exit(-1); }
  if (dom.ny_glob   == 0   ) { std::cout << "Error: key " << "ny"        << " not set."; exit(-1); }
  if (dom.nz_glob   == 0   ) { std::cout << "Error: key " << "nz"        << " not set."; exit(-1); }
  if (dom.simLength == -999) { std::cout << "Error: key " << "simLength" << " not set."; exit(-1); }
  if (par.nproc_x   == -999) { std::cout << "Error: key " << "parNx"     << " not set."; exit(-1); }
  if (par.nproc_y   == -999) { std::cout << "Error: key " << "parNy"     << " not set."; exit(-1); }
  if (outFreq       == -999) { std::cout << "Error: key " << "outFreq"   << " not set."; exit(-1); }

  // Print out the values
  if (par.masterproc) {
    std::cout << "nx:         " << dom.nx_glob   << "\n";
    std::cout << "ny:         " << dom.ny_glob   << "\n";
    std::cout << "nz:         " << dom.nz_glob   << "\n";
    std::cout << "xlen:       " << dom.xlen      << "\n";
    std::cout << "ylen:       " << dom.ylen      << "\n";
    std::cout << "zlen:       " << dom.zlen      << "\n";
    std::cout << "cfl:        " << dom.cfl       << "\n";
    std::cout << "simLength:  " << dom.simLength << "\n";
    std::cout << "parNx:      " << par.nproc_x   << "\n";
    std::cout << "parNy:      " << par.nproc_y   << "\n";
    std::cout << "outFreq:    " << outFreq       << "\n";
    std::cout << "doWeno:     " << dom.doWeno    << "\n";
    std::cout << "timeMethod: " << strTimeMethod << "\n";
    std::cout << "eqnSet:     " << strEqnSet     << "\n";
    std::cout << "dataInit:   " << strDataInit   << "\n";
    std::cout << "strakaVisc: " << strakaVisc    << "\n";
  }

}

void handleTimeMethod(std::string &str, std::string &fNameIn) {
  size_t splitloc = str.find("//",0);
  std::string strloc;
  if (splitloc != std::string::npos){
    strloc = str.substr(0,splitloc);
  } else {
    strloc = str;
  }
  if      ( !strcmp(strloc.c_str(),"SSPRK3") ) { timeMethod = TIME_SSPRK3; }
  else if ( !strcmp(strloc.c_str(),"ADER"  ) ) { timeMethod = TIME_ADER  ; }
  else  {
    std::cout << "Error: unrecognized timeMethod " << str << " in file " << fNameIn << "\n";
    exit(-1);
  }
}

void handleEqnSet(std::string &str, std::string &fNameIn, Domain &dom) {
  size_t splitloc = str.find("//",0);
  std::string strloc;
  if (splitloc != std::string::npos){
    strloc = str.substr(0,splitloc);
  } else {
    strloc = str;
  }
  if      ( !strcmp(strloc.c_str(),"theta_cons") ) { dom.eqnSet = EQN_THETA_CONS; }
  else if ( !strcmp(strloc.c_str(),"theta_prim") ) { dom.eqnSet = EQN_THETA_PRIM; }
  else  {
    std::cout << "Error: unrecognized eqnSet " << str << " in file " << fNameIn << "\n";
    exit(-1);
  }
}

void handleDataInit(std::string &str, std::string &fNameIn, Domain &dom) {
  size_t splitloc = str.find("//",0);
  std::string strloc;
  if (splitloc != std::string::npos){
    strloc = str.substr(0,splitloc);
  } else {
    strloc = str;
  }
  if      ( !strcmp(strloc.c_str(),"thermal"  ) ) { dom.dataInit = DATA_INIT_THERMAL  ; }
  else if ( !strcmp(strloc.c_str(),"collision") ) { dom.dataInit = DATA_INIT_COLLISION; }
  else if ( !strcmp(strloc.c_str(),"straka"   ) ) { dom.dataInit = DATA_INIT_STRAKA   ; }
  else  {
    std::cout << "Error: unrecognized dataInit " << str << " in file " << fNameIn << "\n";
    exit(-1);
  }
}


