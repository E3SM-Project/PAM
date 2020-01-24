
#include "sgs.h"


extern "C" void kurant_sgs(real &cfl) {
  real2d tkhmax("tkhmax",nzm,ncrms);
  auto sgs_field_diag = :: sgs_field_diag;
  auto dz             = :: dz;
  auto dy             = :: dy;
  auto dx             = :: dx;
  auto dt             = :: dt;
  auto adzw           = :: adzw;
  auto grdf_x         = :: grdf_x;
  auto grdf_y         = :: grdf_y;
  auto grdf_z         = :: grdf_z;
  auto ncrms          = :: ncrms;

  // for (int k=0; k<nzm; k++) {
  //   for (int icrm=0; icrm<ncrms; icrm++) {
  yakl::parallel_for( nzm*ncrms , YAKL_LAMBDA (int iGlob) {
    int k, icrm;
    yakl::unpackIndices( iGlob , nzm,ncrms , k,icrm );

    tkhmax(k,icrm) = 0.;
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  yakl::parallel_for( nzm*ny*nx*ncrms , YAKL_LAMBDA (int iGlob) {
    int k, j, i, icrm;
    yakl::unpackIndices( iGlob , nzm,ny,nx,ncrms , k,j,i,icrm );

    yakl::atomicMax( tkhmax(k,icrm) , sgs_field_diag(1,k,offy_d+j,offx_d+i,icrm) );
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int icrm=0; icrm < ncrms; icrm++) {
  yakl::parallel_for( nzm*ncrms , YAKL_LAMBDA (int iGlob) {
    int k, icrm;
    yakl::unpackIndices( iGlob , nzm,ncrms , k,icrm );

    real dztmp = dz(icrm)*adzw(k,icrm);
    real xdir = 0.5*tkhmax(k,icrm)*grdf_x(k,icrm)*dt/(dx*dx);
    real ydir = 0.5*tkhmax(k,icrm)*grdf_y(k,icrm)*dt/(dy*dy)*YES3D;
    real zdir = 0.5*tkhmax(k,icrm)*grdf_z(k,icrm)*dt/(dztmp*dztmp);
    tkhmax(k,icrm) = max( max( xdir , ydir ) , zdir );
  });

  // Perform a max reduction over tkhmax
  yakl::ParallelMax<real,yakl::memDevice> pmax( nzm*ncrms );
  real cfl_loc = pmax( tkhmax.data() );
  cfl = max(cfl , cfl_loc);

}


