
#include "Exchange.h"


void Exchange::allocate(Domain &dom) {
  haloSendBufS = realArr("haloSendBufS",maxPack*dom.nz*hs*dom.nx);
  haloSendBufN = realArr("haloSendBufN",maxPack*dom.nz*hs*dom.nx);
  haloSendBufW = realArr("haloSendBufW",maxPack*dom.nz*dom.ny*hs);
  haloSendBufE = realArr("haloSendBufE",maxPack*dom.nz*dom.ny*hs);
  haloRecvBufS = realArr("haloRecvBufS",maxPack*dom.nz*hs*dom.nx);
  haloRecvBufN = realArr("haloRecvBufN",maxPack*dom.nz*hs*dom.nx);
  haloRecvBufW = realArr("haloRecvBufW",maxPack*dom.nz*dom.ny*hs);
  haloRecvBufE = realArr("haloRecvBufE",maxPack*dom.nz*dom.ny*hs);

  edgeSendBufS = realArr("edgeSendBufS",maxPack*dom.nz*dom.nx);
  edgeSendBufN = realArr("edgeSendBufN",maxPack*dom.nz*dom.nx);
  edgeSendBufW = realArr("edgeSendBufW",maxPack*dom.nz*dom.ny);
  edgeSendBufE = realArr("edgeSendBufE",maxPack*dom.nz*dom.ny);
  edgeRecvBufS = realArr("edgeRecvBufS",maxPack*dom.nz*dom.nx);
  edgeRecvBufN = realArr("edgeRecvBufN",maxPack*dom.nz*dom.nx);
  edgeRecvBufW = realArr("edgeRecvBufW",maxPack*dom.nz*dom.ny);
  edgeRecvBufE = realArr("edgeRecvBufE",maxPack*dom.nz*dom.ny);

  haloSendBufS_host = haloSendBufS.createHostCopy(); 
  haloSendBufN_host = haloSendBufN.createHostCopy(); 
  haloSendBufW_host = haloSendBufW.createHostCopy(); 
  haloSendBufE_host = haloSendBufE.createHostCopy(); 
  haloRecvBufS_host = haloRecvBufS.createHostCopy(); 
  haloRecvBufN_host = haloRecvBufN.createHostCopy(); 
  haloRecvBufW_host = haloRecvBufW.createHostCopy(); 
  haloRecvBufE_host = haloRecvBufE.createHostCopy(); 

  edgeSendBufS_host = edgeSendBufS.createHostCopy(); 
  edgeSendBufN_host = edgeSendBufN.createHostCopy(); 
  edgeSendBufW_host = edgeSendBufW.createHostCopy(); 
  edgeSendBufE_host = edgeSendBufE.createHostCopy(); 
  edgeRecvBufS_host = edgeRecvBufS.createHostCopy(); 
  edgeRecvBufN_host = edgeRecvBufN.createHostCopy(); 
  edgeRecvBufW_host = edgeRecvBufW.createHostCopy(); 
  edgeRecvBufE_host = edgeRecvBufE.createHostCopy(); 
}


void Exchange::haloInit() {
  nPack   = 0;
  nUnpack = 0;
}


void Exchange::haloPackN_x(Domain const &dom, realArr const &a, int const n) {
  auto &haloSendBufW = this->haloSendBufW;
  auto &haloSendBufE = this->haloSendBufE;
  auto &nPack        = this->nPack       ;
  yakl::parallel_for( n*dom.nz*dom.ny*hs , YAKL_LAMBDA (int iGlob) {
    int v, k, j, ii;
    yakl::unpackIndices(iGlob,n,dom.nz,dom.ny,hs,v,k,j,ii);
    int nGlob = dom.nz*dom.ny*hs;
    haloSendBufW(nPack*nGlob+iGlob) = a(v,hs+k,hs+j,hs    +ii);
    haloSendBufE(nPack*nGlob+iGlob) = a(v,hs+k,hs+j,dom.nx+ii);
  });
  nPack = nPack + n;
}


void Exchange::haloPackN_y(Domain const &dom, realArr const &a, int const n) {
  auto &haloSendBufS = this->haloSendBufS;
  auto &haloSendBufN = this->haloSendBufN;
  auto &nPack        = this->nPack       ;
  yakl::parallel_for( n*dom.nz*hs*dom.nx , YAKL_LAMBDA (int iGlob) {
    int v, k, ii, i;
    yakl::unpackIndices(iGlob,n,dom.nz,hs,dom.nx,v,k,ii,i);
    int nGlob = dom.nz*hs*dom.nx;
    haloSendBufS(nPack*nGlob+iGlob) = a(v,hs+k,hs    +ii,hs+i);
    haloSendBufN(nPack*nGlob+iGlob) = a(v,hs+k,dom.ny+ii,hs+i);
  });
  nPack = nPack + n;
}


void Exchange::haloUnpackN_x(Domain const &dom, realArr &a, int const n) {
  auto &haloRecvBufW = this->haloRecvBufW;
  auto &haloRecvBufE = this->haloRecvBufE;
  auto &nUnpack      = this->nUnpack     ;
  yakl::parallel_for( n*dom.nz*dom.ny*hs , YAKL_LAMBDA (int iGlob) {
    int v, k, j, ii;
    yakl::unpackIndices(iGlob,n,dom.nz,dom.ny,hs,v,k,j,ii);
    int nGlob = dom.nz*dom.ny*hs;
    a(v,hs+k,hs+j,          ii) = haloRecvBufW(nUnpack*nGlob+iGlob);
    a(v,hs+k,hs+j,dom.nx+hs+ii) = haloRecvBufE(nUnpack*nGlob+iGlob);
  });
  nUnpack = nUnpack + n;
}


void Exchange::haloUnpackN_y(Domain const &dom, realArr &a, int const n) {
  auto &haloRecvBufS = this->haloRecvBufS;
  auto &haloRecvBufN = this->haloRecvBufN;
  auto &nUnpack      = this->nUnpack     ;
  yakl::parallel_for( n*dom.nz*hs*dom.nx , YAKL_LAMBDA (int iGlob) {
    int v, k, ii, i;
    yakl::unpackIndices(iGlob,n,dom.nz,hs,dom.nx,v,k,ii,i);
    int nGlob = dom.nz*hs*dom.nx;
    a(v,hs+k,          ii,hs+i) = haloRecvBufS(nUnpack*nGlob+iGlob);
    a(v,hs+k,dom.ny+hs+ii,hs+i) = haloRecvBufN(nUnpack*nGlob+iGlob);
  });
  nUnpack = nUnpack + n;
}


void Exchange::haloExchange_x(Domain const &dom, Parallel const &par) {
  int ierr;

  if (par.nproc_x > 0) {
    yakl::fence();

    //Pre-post the receives
    ierr = MPI_Irecv( haloRecvBufW_host.data() , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( haloRecvBufE_host.data() , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &rReq[1] );

    haloSendBufW.deep_copy_to(haloSendBufW_host);
    haloSendBufE.deep_copy_to(haloSendBufE_host);
    yakl::fence();

    //Send the data
    ierr = MPI_Isend( haloSendBufW_host.data() , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( haloSendBufE_host.data() , nPack*dom.nz*dom.ny*hs , MPI_FLOAT , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);

    haloRecvBufW_host.deep_copy_to(haloRecvBufW);
    haloRecvBufE_host.deep_copy_to(haloRecvBufE);

  } else {
    haloExchange_x_loc(dom, haloSendBufW, haloSendBufE, haloRecvBufW, haloRecvBufE);
  }
}
void Exchange::haloExchange_x_loc(Domain const &dom, realArr &haloSendBufW, realArr &haloSendBufE, realArr &haloRecvBufW, realArr &haloRecvBufE) {
  yakl::parallel_for( nPack*dom.nz*dom.ny*hs , YAKL_LAMBDA (int iGlob) {
    haloRecvBufW(iGlob) = haloSendBufE(iGlob);
    haloRecvBufE(iGlob) = haloSendBufW(iGlob);
  });
}


void Exchange::haloExchange_y(Domain const &dom, Parallel const &par) {
  int ierr;

  if (par.nproc_y > 0) {
    yakl::fence();

    //Pre-post the receives
    ierr = MPI_Irecv( haloRecvBufS_host.data() , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( haloRecvBufN_host.data() , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &rReq[1] );

    haloSendBufS.deep_copy_to(haloSendBufS_host);
    haloSendBufN.deep_copy_to(haloSendBufN_host);
    yakl::fence();

    //Send the data
    ierr = MPI_Isend( haloSendBufS_host.data() , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( haloSendBufN_host.data() , nPack*dom.nz*hs*dom.nx , MPI_FLOAT , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);

    haloRecvBufS_host.deep_copy_to(haloRecvBufS);
    haloRecvBufN_host.deep_copy_to(haloRecvBufN);

  } else {
    haloExchange_y_loc(dom, haloSendBufS, haloSendBufN, haloRecvBufS, haloRecvBufN);
  }
}
void Exchange::haloExchange_y_loc(Domain const &dom, realArr &haloSendBufS, realArr &haloSendBufN, realArr &haloRecvBufS, realArr &haloRecvBufN) {
  yakl::parallel_for( nPack*dom.nz*hs*dom.nx , YAKL_LAMBDA (int iGlob) {
    haloRecvBufS(iGlob) = haloSendBufN(iGlob);
    haloRecvBufN(iGlob) = haloSendBufS(iGlob);
  });
}


void Exchange::edgePackN_x(Domain const &dom, realArr const &a, int const n) {
  auto &edgeSendBufW = this->edgeSendBufW;
  auto &edgeSendBufE = this->edgeSendBufE;
  auto &nPack        = this->nPack       ;
  // for (int v=0; v<n; v++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  yakl::parallel_for( n*dom.nz*dom.ny , YAKL_LAMBDA (int iGlob) {
    int v, k, j;
    yakl::unpackIndices(iGlob,n,dom.nz,dom.ny,v,k,j);
    int nGlob = dom.nz*dom.ny;
    edgeSendBufW(nPack*nGlob+iGlob) = a(v,1,k,j,0     );
    edgeSendBufE(nPack*nGlob+iGlob) = a(v,0,k,j,dom.nx);
  });
  nPack = nPack + n;
}


void Exchange::edgePackN_y(Domain const &dom, realArr const &a, int const n) {
  auto &edgeSendBufS = this->edgeSendBufS;
  auto &edgeSendBufN = this->edgeSendBufN;
  auto &nPack        = this->nPack       ;
  // for (int v=0; v<n; v++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( n*dom.nz*dom.nx , YAKL_LAMBDA (int iGlob) {
    int v, k, i;
    yakl::unpackIndices(iGlob,n,dom.nz,dom.nx,v,k,i);
    int nGlob = dom.nz*dom.nx;
    edgeSendBufS(nPack*nGlob+iGlob) = a(v,1,k,0     ,i);
    edgeSendBufN(nPack*nGlob+iGlob) = a(v,0,k,dom.ny,i);
  });
  nPack = nPack + n;
}


void Exchange::edgeUnpackN_x(Domain const &dom, realArr &a, int const n) {
  auto &edgeRecvBufW = this->edgeRecvBufW;
  auto &edgeRecvBufE = this->edgeRecvBufE;
  auto &nUnpack      = this->nUnpack     ;
  // for (int v=0; v<n; v++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int j=0; j<dom.ny; j++) {
  yakl::parallel_for( n*dom.nz*dom.ny , YAKL_LAMBDA (int iGlob) {
    int v, k, j;
    yakl::unpackIndices(iGlob,n,dom.nz,dom.ny,v,k,j);
    int nGlob = dom.nz*dom.ny;
    a(v,0,k,j,0     ) = edgeRecvBufW(nUnpack*nGlob+iGlob);
    a(v,1,k,j,dom.nx) = edgeRecvBufE(nUnpack*nGlob+iGlob);
  });
  nUnpack = nUnpack + n;
}


void Exchange::edgeUnpackN_y(Domain const &dom, realArr &a, int const n) {
  auto &edgeRecvBufS = this->edgeRecvBufS;
  auto &edgeRecvBufN = this->edgeRecvBufN;
  auto &nUnpack      = this->nUnpack     ;
  // for (int v=0; v<n; v++) {
  //   for (int k=0; k<dom.nz; k++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( n*dom.nz*dom.nx , YAKL_LAMBDA (int iGlob) {
    int v, k, i;
    yakl::unpackIndices(iGlob,n,dom.nz,dom.nx,v,k,i);
    int nGlob = dom.nz*dom.nx;
    a(v,0,k,0     ,i) = edgeRecvBufS(nUnpack*nGlob+iGlob);
    a(v,1,k,dom.ny,i) = edgeRecvBufN(nUnpack*nGlob+iGlob);
  });
  nUnpack = nUnpack + n;
}


void Exchange::edgeExchange_x(Domain const &dom, Parallel const &par) {
  int ierr;

  if (par.nproc_x > 0) {
    yakl::fence();

    //Pre-post the receives
    ierr = MPI_Irecv( edgeRecvBufW_host.data() , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( edgeRecvBufE_host.data() , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &rReq[1] );

    edgeSendBufW.deep_copy_to(edgeSendBufW_host);
    edgeSendBufE.deep_copy_to(edgeSendBufE_host);
    yakl::fence();

    //Send the data
    ierr = MPI_Isend( edgeSendBufW_host.data() , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( edgeSendBufE_host.data() , nPack*dom.nz*dom.ny , MPI_FLOAT , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);

    edgeRecvBufW_host.deep_copy_to(edgeRecvBufW);
    edgeRecvBufE_host.deep_copy_to(edgeRecvBufE);

  } else {
    edgeExchange_x_loc(dom, edgeSendBufW, edgeSendBufE, edgeRecvBufW, edgeRecvBufE);
  }
}
void Exchange::edgeExchange_x_loc(Domain const &dom, realArr &edgeSendBufW, realArr &edgeSendBufE, realArr &edgeRecvBufW, realArr &edgeRecvBufE) {
  yakl::parallel_for( nPack*dom.nz*dom.ny , YAKL_LAMBDA (int iGlob) {
    edgeRecvBufW(iGlob) = edgeSendBufE(iGlob);
    edgeRecvBufE(iGlob) = edgeSendBufW(iGlob);
  });
}


void Exchange::edgeExchange_y(Domain const &dom, Parallel const &par) {
  int ierr;

  if (par.nproc_y > 0) {
    yakl::fence();

    //Pre-post the receives
    ierr = MPI_Irecv( edgeRecvBufS_host.data() , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( edgeRecvBufN_host.data() , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &rReq[1] );

    edgeSendBufS.deep_copy_to(edgeSendBufS_host);
    edgeSendBufN.deep_copy_to(edgeSendBufN_host);
    yakl::fence();

    //Send the data
    ierr = MPI_Isend( edgeSendBufS_host.data() , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( edgeSendBufN_host.data() , nPack*dom.nz*dom.nx , MPI_FLOAT , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);

    edgeRecvBufS_host.deep_copy_to(edgeRecvBufS);
    edgeRecvBufN_host.deep_copy_to(edgeRecvBufN);

  } else {
    edgeExchange_y_loc(dom, edgeSendBufS, edgeSendBufN, edgeRecvBufS, edgeRecvBufN);
  }
}
void Exchange::edgeExchange_y_loc(Domain const &dom, realArr &edgeSendBufS, realArr &edgeSendBufN, realArr &edgeRecvBufS, realArr &edgeRecvBufN) {
  yakl::parallel_for( nPack*dom.nz*dom.nx, YAKL_LAMBDA (int iGlob) {
    edgeRecvBufS(iGlob) = edgeSendBufN(iGlob);
    edgeRecvBufN(iGlob) = edgeSendBufS(iGlob);
  });
}

