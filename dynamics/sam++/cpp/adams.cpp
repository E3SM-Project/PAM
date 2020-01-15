
#include "adams.h"

extern "C" void adams() {
  // Adams-Bashforth scheme
  real dtdx = dtn/dx;
  real dtdy = dtn/dy;

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  Kokkos::parallel_for( "adams" , nzm*ny*nx*ncrms , KOKKOS_LAMBDA (int iGlob) {
    int k, j, i, icrm;
    yakl::unpackIndices( iGlob , nzm , ny , nx , ncrms , k , j , i , icrm );
    real dtdz = dtn/dz(icrm);
    real rhox = rho (k,icrm)*dtdx;
    real rhoy = rho (k,icrm)*dtdy;
    real rhoz = rhow(k,icrm)*dtdz;
    dudt(nc-1,k,j,i,icrm) = u(k,j+offy_u,i+offx_u,icrm) + dt3(na-1) * ( at*dudt(na-1,k,j,i,icrm) + bt*dudt(nb-1,k,j,i,icrm) + ct*dudt(nc-1,k,j,i,icrm) );
    dvdt(nc-1,k,j,i,icrm) = v(k,j+offy_v,i+offx_v,icrm) + dt3(na-1) * ( at*dvdt(na-1,k,j,i,icrm) + bt*dvdt(nb-1,k,j,i,icrm) + ct*dvdt(nc-1,k,j,i,icrm) );
    dwdt(nc-1,k,j,i,icrm) = w(k,j+offy_w,i+offx_w,icrm) + dt3(na-1) * ( at*dwdt(na-1,k,j,i,icrm) + bt*dwdt(nb-1,k,j,i,icrm) + ct*dwdt(nc-1,k,j,i,icrm) );
    u   (k,j+offy_u,i+offx_u,icrm) = 0.5 * ( u(k,j+offy_u,i+offx_u,icrm) + dudt(nc-1,k,j,i,icrm) ) * rhox;
    v   (k,j+offy_v,i+offx_v,icrm) = 0.5 * ( v(k,j+offy_v,i+offx_v,icrm) + dvdt(nc-1,k,j,i,icrm) ) * rhoy;
    w   (k,j+offy_w,i+offx_w,icrm) = 0.5 * ( w(k,j+offy_w,i+offx_w,icrm) + dwdt(nc-1,k,j,i,icrm) ) * rhoz;
    misc(k,j       ,i       ,icrm) = 0.5 * ( w(k,j+offy_w,i+offx_w,icrm) + dwdt(nc-1,k,j,i,icrm) );
  });
}

