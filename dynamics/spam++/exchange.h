#ifndef _EXCHANGE_H_
#define _EXCHANGE_H_


#include "common.h"
#include "field.h"
#include "topology.h"


template<int ndims, int ndof0, int ndof1, int ndof2, int ndof3> class Exchange {
public:

  int bufsize_x, bufsize_y, bufsize_z;

  Topology topology;
  int total_dofs;

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

  void clone(Exchange &exch);
  void pack(Field &field);
  void exchange();
  void unpack(Field &field);
  void initialize(Topology &topo);

};

class PeriodicExchange : Exchange {

public:


};

#endif
