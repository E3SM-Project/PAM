
#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "grid.h"

class Domain {

public:

  real   dx;          // x-grid spacing
  real   dy;          // y-grid spacing
  real   dtn;         // Kurant-Sub-cycled time step
  real   dtfactor;    // Kurant sub-cycling factor
  real   at, bt, ct;  // Adams-Bashforth constants
  int    na, nb, nc;  // Indices for the A-B steps
  int    nstep;       // Current number of performed time steps
  int    ncycle;      // Number of kurant subcycles over the dynamical timestep
  int    icycle;      // Current subcycle stage
  int    nstop;       // Time step number to stop the integration
  real   dt;          // Model time step

  real1d dz;          // Z-grid spacing factor
  real2d z;           // Height of the pressure levels above surface,m
  real2d pres;        // Pressure,mb at scalar levels
  real2d zi;          // Height of the interface levels
  real2d presi;       // Pressure,mb at interface levels
  real2d adz;         // Ratio of the thickness of scalar levels to dz
  real2d adzw;        // Ratio of the thinckness of w levels to dz
  real1d dt3;         // Last three time stepps dims(3)

  real1d fcorz;       // Vertical Coriolis parameter
  real1d fcor;        // Coriolis parameter
  real1d longitude0;  // Latitude of the domain's center
  real1d latitude0;   // Longitude of the domain's center
  real1d z0;          // Roughness length
  real1d uhl;         // Current large-scale velocity in x near sfc
  real1d vhl;         // Current large-scale velocity in y near sfc
  real1d taux0;       // Surface stress in x, m2/s2
  real1d tauy0;       // Surface stress in y, m2/s2

  Domain(int ncrms) {
    z          = real2d("z"         ,ncrms,nz );
    pres       = real2d("pres"      ,ncrms,nzm);
    zi         = real2d("zi "       ,ncrms,nz );
    presi      = real2d("presi"     ,ncrms,nz );
    adz        = real2d("adz"       ,ncrms,nzm);
    adzw       = real2d("adzw"      ,ncrms,nz );
    dz         = real1d("dz"        ,ncrms    );
    dt3        = real1d("dt3"       ,3        );
    fcor       = real1d("fcor"      ,ncrms    );
    fcorz      = real1d("fcorz"     ,ncrms    );
    longitude0 = real1d("longitude0",ncrms    );
    latitude0  = real1d("latitude0" ,ncrms    );
    z0         = real1d("z0"        ,ncrms    );
    uhl        = real1d("uhl"       ,ncrms    );
    vhl        = real1d("vhl"       ,ncrms    );
    taux0      = real1d("taux0"     ,ncrms    );
    tauy0      = real1d("tauy0"     ,ncrms    );
  }

};

#endif

