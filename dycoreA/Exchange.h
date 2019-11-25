
#ifndef _EXCHANGE_H_
#define _EXCHANGE_H_

#include "const.h"
#include "mpi.h"
#include "Domain.h"
#include "Parallel.h"


class Exchange {

protected:

  int const maxPack = numState*2;

  MPI_Request sReq [2];
  MPI_Request rReq [2];

  MPI_Status  sStat[2];
  MPI_Status  rStat[2];

  int nPack;
  int nUnpack;

  realArr haloSendBufS;
  realArr haloSendBufN;
  realArr haloSendBufW;
  realArr haloSendBufE;
  realArr haloRecvBufS;
  realArr haloRecvBufN;
  realArr haloRecvBufW;
  realArr haloRecvBufE;

  realArr edgeRecvBufE;
  realArr edgeRecvBufW;
  realArr edgeSendBufE;
  realArr edgeSendBufW;
  realArr edgeRecvBufN;
  realArr edgeRecvBufS;
  realArr edgeSendBufN;
  realArr edgeSendBufS;

  realArrHost haloSendBufS_host;
  realArrHost haloSendBufN_host;
  realArrHost haloSendBufW_host;
  realArrHost haloSendBufE_host;
  realArrHost haloRecvBufS_host;
  realArrHost haloRecvBufN_host;
  realArrHost haloRecvBufW_host;
  realArrHost haloRecvBufE_host;

  realArrHost edgeRecvBufE_host;
  realArrHost edgeRecvBufW_host;
  realArrHost edgeSendBufE_host;
  realArrHost edgeSendBufW_host;
  realArrHost edgeRecvBufN_host;
  realArrHost edgeRecvBufS_host;
  realArrHost edgeSendBufN_host;
  realArrHost edgeSendBufS_host;

public:


  void allocate(Domain &dom);


  void haloInit();


  void haloPackN_x(Domain const &dom, realArr const &a, int const n);


  void haloPackN_y(Domain const &dom, realArr const &a, int const n);


  void haloUnpackN_x(Domain const &dom, realArr &a, int const n);


  void haloUnpackN_y(Domain const &dom, realArr &a, int const n);


  void haloExchange_x(Domain const &dom, Parallel const &par);
  void haloExchange_x_loc(Domain const &dom, realArr &haloSendBufW, realArr &haloSendBufE, realArr &haloRecvBufW, realArr &haloRecvBufE);


  void haloExchange_y(Domain const &dom, Parallel const &par);
  void haloExchange_y_loc(Domain const &dom, realArr &haloSendBufS, realArr &haloSendBufN, realArr &haloRecvBufS, realArr &haloRecvBufN);


  void edgePackN_x(Domain const &dom, realArr const &a, int const n);


  void edgePackN_y(Domain const &dom, realArr const &a, int const n);


  void edgeUnpackN_x(Domain const &dom, realArr &a, int const n);


  void edgeUnpackN_y(Domain const &dom, realArr &a, int const n);


  void edgeExchange_x(Domain const &dom, Parallel const &par);
  void edgeExchange_x_loc(Domain const &dom, realArr &edgeSendBufW, realArr &edgeSendBufE, realArr &edgeRecvBufW, realArr &edgeRecvBufE);


  void edgeExchange_y(Domain const &dom, Parallel const &par);
  void edgeExchange_y_loc(Domain const &dom, realArr &edgeSendBufS, realArr &edgeSendBufN, realArr &edgeRecvBufS, realArr &edgeRecvBufN);

};

#endif
