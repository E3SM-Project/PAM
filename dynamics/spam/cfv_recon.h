#ifndef _CFV_RECON_H_
#define _CFV_RECON_H_

#include "common.h"


template<uint ndofs, uint nd> void YAKL_INLINE cfv(SArray<real,3,ndofs,nd,2> &edgerecon, SArray<real,3,ndofs,nd,1> const &dens) {
    for (int l=0; l<ndofs; l++) {
      for (int d=0; d<nd; d++) {
          edgerecon(l,d,0) = dens(l,d,0);
          edgerecon(l,d,1) = dens(l,d,0);
        }
      }
}


template<uint ndofs, uint nd> void YAKL_INLINE cfv(SArray<real,3,ndofs,nd,2> &edgerecon, SArray<real,3,ndofs,nd,3> const &dens) {
    real er;
    for (int l=0; l<ndofs; l++) {
      for (int d=0; d<nd; d++) {
        er = (8.0/6.0)*dens(l,d,1) -(1.0/6.0)*(dens(l,d,0) + dens(l,d,2));
          edgerecon(l,d,0) = er;
          edgerecon(l,d,1) = er;
        }
      }
}


template<uint ndofs, uint nd> void YAKL_INLINE cfv(SArray<real,3,ndofs,nd,2> &edgerecon, SArray<real,3,ndofs,nd,5> const &dens) {
    real er;
    for (int l=0; l<ndofs; l++) {
      for (int d=0; d<nd; d++) {
        er = (46.0/30.0)*dens(l,d,2) -(9./30.0)*(dens(l,d,1) + dens(l,d,3)) +(1./30.0)*(dens(l,d,0) + dens(l,d,4));
          edgerecon(l,d,0) = er;
          edgerecon(l,d,1) = er;
        }
      }
}

template<uint ndofs, uint nd> void YAKL_INLINE cfv(SArray<real,3,ndofs,nd,2> &edgerecon, SArray<real,3,ndofs,nd,7> const &dens) {
    real er;
    for (int l=0; l<ndofs; l++) {
      for (int d=0; d<nd; d++) {
        er = (704.0/420.0)*dens(l,d,3) -(171./420.0)*(dens(l,d,2) + dens(l,d,4)) +(32./420.0)*(dens(l,d,1) + dens(l,d,5)) -(3./420.0)*(dens(l,d,0) + dens(l,d,6));
          edgerecon(l,d,0) = er;
          edgerecon(l,d,1) = er;
        }
      }
}


template<uint ndofs, uint nd> void YAKL_INLINE cfv(SArray<real,3,ndofs,nd,2> &edgerecon, SArray<real,3,ndofs,nd,9> const &dens) {
    real er;
    for (int l=0; l<ndofs; l++) {
      for (int d=0; d<nd; d++) {
        er = (2252.0/1260.0)*dens(l,d,4) -(625./1260.0)*(dens(l,d,3) + dens(l,d,5)) +(152./1260.0)*(dens(l,d,2) + dens(l,d,6)) -(25./1260.0)*(dens(l,d,1) + dens(l,d,7)) +(2./1260.0)*(dens(l,d,0) + dens(l,d,8));
          edgerecon(l,d,0) = er;
          edgerecon(l,d,1) = er;
        }
      }
}


#endif
