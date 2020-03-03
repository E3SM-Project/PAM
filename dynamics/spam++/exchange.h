#ifndef _EXCHANGE_H_
#define _EXCHANGE_H_


#include "common.h"
#include "field.h"
#include "topology.h"


class Exchange {
public:
  void clone(Exchange &exch);
  void initialize(Topology &topo, int ndof0, int ndof1, int ndof2 = 0, int ndof3 = 0);
  void pack(Field &field);
  void exchange();
  void unpack(Field &field);

};

class PeriodicExchange : Exchange {

public:
  int ndofs0, ndofs1, ndofs2, ndofs3;
  int bufsize_x, bufsize_y, bufsize_z;

  Topology topology;

  realArrHost haloSendBuf_Xm_host;
  realArrHost haloRecvBuf_Xm_host;
  realArr haloSendBuf_Xm;
  realArr haloRecvBuf_Xm;
  realArrHost haloSendBuf_Xp_host;
  realArrHost haloRecvBuf_Xp_host;
  realArr haloSendBuf_Xp;
  realArr haloRecvBuf_Xp;

  realArrHost haloSendBuf_Ym_host;
  realArrHost haloRecvBuf_Ym_host;
  realArr haloSendBuf_Ym;
  realArr haloRecvBuf_Ym;
  realArrHost haloSendBuf_Yp_host;
  realArrHost haloRecvBuf_Yp_host;
  realArr haloSendBuf_Yp;
  realArr haloRecvBuf_Yp;

  realArrHost haloSendBuf_Zm_host;
  realArrHost haloRecvBuf_Zm_host;
  realArr haloSendBuf_Zm;
  realArr haloRecvBuf_Zm;
  realArrHost haloSendBuf_Zp_host;
  realArrHost haloRecvBuf_Zp_host;
  realArr haloSendBuf_Zp;
  realArr haloRecvBuf_Zp;

};

#endif
