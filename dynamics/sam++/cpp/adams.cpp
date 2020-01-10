
#include "adams.h"

// Uses Adams-Bashforth time stepping to advance u, v, and w
void adams(Domain const &dom, real1d const &dz, real2d const &rho, real2d const &rhow,
           real5d &dudt, real5d &dvdt, real5d dwdt, real4d &u, real4d &v, real4d &w, real4d &misc) {
  // Adams-Bashforth scheme
  real dtdx = dtn/dx;
  real dtdy = dtn/dy;

  for (int k=0; k<dom.nzm; k++) {
    for (int j=0; j<dom.ny; j++) {
      for (int i=0; i<dom.nx; i++) {
        for (int icrm=0; icrm<dom.ncrms; icrm++) {
          real dtdz = dtn / dom.dz(icrm);
          real rhox = rho (icrm,k) * dtdx;
          real rhoy = rho (icrm,k) * dtdy;
          real rhoz = rhow(icrm,k) * dtdz;
          real utend = ( dom.at*dudt(icrm,i,j,k,na) + dom.bt*dudt(icrm,i,j,k,nb) + dom.ct*dudt(icrm,i,j,k,nc) );
          real vtend = ( dom.at*dvdt(icrm,i,j,k,na) + dom.bt*dvdt(icrm,i,j,k,nb) + dom.ct*dvdt(icrm,i,j,k,nc) );
          real wtend = ( dom.at*dwdt(icrm,i,j,k,na) + dom.bt*dwdt(icrm,i,j,k,nb) + dom.ct*dwdt(icrm,i,j,k,nc) );
          dudt(icrm,i,j,k,nc) = u(icrm,i,j,k) + dom.dt3(na) * utend;
          dvdt(icrm,i,j,k,nc) = v(icrm,i,j,k) + dom.dt3(na) * vtend;
          dwdt(icrm,i,j,k,nc) = w(icrm,i,j,k) + dom.dt3(na) * wtend;
          u   (icrm,i,j,k) = 0.5 * ( u(icrm,i,j,k) + dudt(icrm,i,j,k,nc) ) * rhox;
          v   (icrm,i,j,k) = 0.5 * ( v(icrm,i,j,k) + dvdt(icrm,i,j,k,nc) ) * rhoy;
          w   (icrm,i,j,k) = 0.5 * ( w(icrm,i,j,k) + dwdt(icrm,i,j,k,nc) ) * rhoz;
          misc(icrm,i,j,k) = w(icrm,i,j,k) / rhoz;
        }
      }
    }
  }

}

