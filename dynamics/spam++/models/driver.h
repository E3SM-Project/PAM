#ifndef _DRIVER_H_
#define _DRIVER_H_

#include "common.h"
#include "tendencies.h"
#include "variable_set.h"
#include "util.h"
#include "topology.h"
#include <string>

class ModelTendencies : Tendencies {};

class PrognosticVars : VariableSet {};
class ConstantVars : VariableSet {};
class DiagnosticVars : VariableSet {};

void set_initial_conditions(VariableSet &constant_vars, VariableSet &prognostic_vars);
void readParamsFile(std::string inFile) {};

// parameters
int nx, ny, nz;
int Nsteps, nout;
real dt;
real xc, yc, zc;
int xlen, ylen, zlen;

#endif
