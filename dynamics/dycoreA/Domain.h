
#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "const.h"

class Domain {

public:

  ulong nx_glob;
  ulong ny_glob;
  ulong nz_glob;

  int nx;
  int ny;
  int nz;

  int run2d;

  int doWeno;

  real xlen;
  real ylen;
  real zlen;

  real dx;
  real dy;
  real dz;

  real cfl;
  real simLength;

  real etime;

  int dataInit;

  real dt;

  // Hydrostatic background state cell averages (no halos)
  realArr hyDensCells;
  realArr hyDensThetaCells;
  realArr hyThetaCells;
  realArr hyPressureCells;
  realArr hyEnergyCells;

  // Hydrostatic background state at tord GLL points within cells (no halos)
  realArr hyDensGLL;
  realArr hyDensThetaGLL;
  realArr hyThetaGLL;
  realArr hyPressureGLL;
  realArr hyEnergyGLL;
};

#endif
