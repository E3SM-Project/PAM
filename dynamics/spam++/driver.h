#ifndef _DRIVER_H_
#define _DRIVER_H_

#include "common.h"
#include "tendencies.h"
#include "finitevolume.h"
#include "divergence.h"
#include "RKSimple.h"
#include "variable_sets.h"
#include "util.h"
#include "topology.h"
#include "geometry.h"
#include <cmath>
#include <iostream>

#include "driver.h"

// HOW DO WE MAKE THIS UNIVERSAL? VIA PRE-PROCESSOR FLAGS, OF COURSE...
#include "advection.h"


// Time scheme types
int constexpr TIME_TYPE_KGRK = 1;
int constexpr TIME_TYPE_ADER = 2;

// Reconstruction types
int constexpr RECONSTRUCTION_TYPE_FV   = 1;
int constexpr RECONSTRUCTION_TYPE_WENO = 2;


void readParamsFile(std::string inFile) {};

void set_initial_conditions(VariableSet &prognostic_vars, VariableSet &constant_vars, VariableSet &diagnostic_vars);
void initialize_variables(VariableSet &prognostic_vars, VariableSet &constant_vars, VariableSet &diagnostic_vars);

// parameters
int nx, ny, nz;
int Nsteps, nout;
real dt, time;

real xc, yc, zc; // UNIFORM RECT SPECIFIC...
int xlen, ylen, zlen; // UNIFORM RECT SPECIFIC...

#endif
