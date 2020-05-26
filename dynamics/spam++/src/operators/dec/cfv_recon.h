#ifndef _CFV_RECON_H_
#define _CFV_RECON_H_

#include "common.h"


template<uint ndofs> void YAKL_INLINE cfv(SArray<real,ndofs,ndims,2> &edgerecon, SArray<real,ndofs,ndims,1> const &dens) {
    for (int l=0; l<ndofs; l++) {
      for (int k=0; k<ndims; k++) {
          edgerecon(l,k,0) = dens(l,k,0);
          edgerecon(l,k,1) = dens(l,k,0);
        }
      }
}


template<uint ndofs> void YAKL_INLINE cfv(SArray<real,ndofs,ndims,2> &edgerecon, SArray<real,ndofs,ndims,3> const &dens) {
    real er;
    for (int l=0; l<ndofs; l++) {
      for (int k=0; k<ndims; k++) {
        er = (8.0/6.0)*dens(l,k,1) -(1.0/6.0)*(dens(l,k,0) + dens(l,k,2));
          edgerecon(l,k,0) = er;
          edgerecon(l,k,1) = er;
        }
      }
}


template<uint ndofs> void YAKL_INLINE cfv(SArray<real,ndofs,ndims,2> &edgerecon, SArray<real,ndofs,ndims,5> const &dens) {
    real er;
    for (int l=0; l<ndofs; l++) {
      for (int k=0; k<ndims; k++) {
        er = (46.0/30.0)*dens(l,k,2) -(9./30.0)*(dens(l,k,1) + dens(l,k,3)) +(1./30.0)*(dens(l,k,0) + dens(l,k,4));
          edgerecon(l,k,0) = er;
          edgerecon(l,k,1) = er;
        }
      }
}

template<uint ndofs> void YAKL_INLINE cfv(SArray<real,ndofs,ndims,2> &edgerecon, SArray<real,ndofs,ndims,7> const &dens) {
    real er;
    for (int l=0; l<ndofs; l++) {
      for (int k=0; k<ndims; k++) {
        er = (704.0/420.0)*dens(l,k,3) -(171./420.0)*(dens(l,k,2) + dens(l,k,4)) +(32./420.0)*(dens(l,k,1) + dens(l,k,5)) -(3./420.0)*(dens(l,k,0) + dens(l,k,6));
          edgerecon(l,k,0) = er;
          edgerecon(l,k,1) = er;
        }
      }
}


template<uint ndofs> void YAKL_INLINE cfv(SArray<real,ndofs,ndims,2> &edgerecon, SArray<real,ndofs,ndims,9> const &dens) {
    real er;
    for (int l=0; l<ndofs; l++) {
      for (int k=0; k<ndims; k++) {
        er = (2252.0/1260.0)*dens(l,k,4) -(625./1260.0)*(dens(l,k,3) + dens(l,k,5)) +(152./1260.0)*(dens(l,k,2) + dens(l,k,6)) -(25./1260.0)*(dens(l,k,1) + dens(l,k,7)) +(2./1260.0)*(dens(l,k,0) + dens(l,k,8));
          edgerecon(l,k,0) = er;
          edgerecon(l,k,1) = er;
        }
      }
}


#endif
