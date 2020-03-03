

#include "exchange.h"

void Exchange::clone(Exchange &exch) {
  exch.initialize(topology, ndofs0, ndofs1, ndofs2, ndofs3);
}


// ALL BROKEN FOR PARALLEL
  void PeriodicExchange::initialize(Topology &topo, int ndof0, int ndof1, int ndof2 = 0, int ndof3 = 0) {
    topology = topo;
    ndofs0 = ndof0;
    ndofs1 = ndof1;
    ndofs2 = ndof2;
    ndofs3 = ndof3;

    int total_dofs;
    if (ndims == 1)
    {
      total_dofs = ndofs0 + ndofs1;
      bufsize_x = topology.halosize_x*total_dofs;

    }

    // This duplicates corner elements- probably okay
    if (ndims == 2)
    {
      total_dofs = ndofs0 + 2*ndofs1 + ndofs2;
      bufsize_x = topology.halosize_x*total_dofs*topology.n_cell_y;
      bufsize_y = topology.halosize_y*total_dofs*topology.n_cell_x;
    }

    // This duplicates corner elements- probably okay
    if (ndims == 3)
    {
      total_dofs = ndofs0 + 3*ndofs1 + 3*ndofs2 + ndofs3;
      bufsize_x = topology.halosize_x*total_dofs*topology.n_cell_y*topology.n_cell_z;
      bufsize_y = topology.halosize_y*total_dofs*topology.n_cell_x*topology.n_cell_z;
      bufsize_z = topology.halosize_z*total_dofs*topology.n_cell_x*topology.n_cell_y;
    }

    haloSendBuf_Xm = realArr("haloSendBuf_Xm", bufsize_x);
    haloRecvBuf_Xm = realArr("haloRecvBuf_Xm", bufsize_x);
    haloSendBuf_Xm_host = haloSendBuf_Xm.createHostCopy();
    haloRecvBuf_Xm_host = haloRecvBuf_Xm.createHostCopy();
    haloSendBuf_Xp = realArr("haloSendBuf_Xp", bufsize_x);
    haloRecvBuf_Xp = realArr("haloRecvBuf_Xp", bufsize_x);
    haloSendBuf_Xp_host = haloSendBuf_Xp.createHostCopy();
    haloRecvBuf_Xp_host = haloRecvBuf_Xp.createHostCopy();

    if (ndims >= 2) {
      haloSendBuf_Ym = realArr("haloSendBuf_Ym", bufsize_y);
      haloRecvBuf_Ym = realArr("haloRecvBuf_Ym", bufsize_y);
      haloSendBuf_Ym_host = haloSendBuf_Ym.createHostCopy();
      haloRecvBuf_Ym_host = haloRecvBuf_Ym.createHostCopy();
      haloSendBuf_Yp = realArr("haloSendBuf_Yp", bufsize_y);
      haloRecvBuf_Yp = realArr("haloRecvBuf_Yp", bufsize_y);
      haloSendBuf_Yp_host = haloSendBuf_Yp.createHostCopy();
      haloRecvBuf_Yp_host = haloRecvBuf_Yp.createHostCopy();
    }

    if (ndims == 3) {
      haloSendBuf_Zm = realArr("haloSendBuf_Zm", bufsize_z);
      haloRecvBuf_Zm = realArr("haloRecvBuf_Zm", bufsize_z);
      haloSendBuf_Zm_host = haloSendBuf_Zm.createHostCopy();
      haloRecvBuf_Zm_host = haloRecvBuf_Zm.createHostCopy();
      haloSendBuf_Zp = realArr("haloSendBuf_Zp", bufsize_z);
      haloRecvBuf_Zp = realArr("haloRecvBuf_Zp", bufsize_z);
      haloSendBuf_Zp_host = haloSendBuf_Zp.createHostCopy();
      haloRecvBuf_Zp_host = haloRecvBuf_Zp.createHostCopy();
    }

  }

  void PeriodicExchange::pack(Field &field) {

    int is = topology.is;
    int js = topology.js;
    int ks = topology.ks;

    // pack left (x-) and  right (x+)
    yakl::parallel_for("PackLeftRight", bufsize_x, YAKL_LAMBDA (int iGlob) {
      int ndof, ii, k, j;
      yakl::unpackIndices(iGlob, total_dofs, topology.halosize_x, topology.n_cells_z, topology.n_cells_y, ndof, ii, k, j);
      haloSendBuf_Xm(iGlob) = field.data(ndof, k+ks, j+js, ii+is+topology.n_cells_x);
      haloSendBuf_Xp(iGlob) = field.data(ndof, k+ks, j+js, ii+is);
    });

    //pack down (y-) and up (y+)
    if (ndims >=2) {
      yakl::parallel_for("PackDownUp", bufsize_y, YAKL_LAMBDA (int iGlob) {
        int ndof, i, k, jj;
        yakl::unpackIndices(iGlob, total_dofs, topology.halosize_y, topology.n_cells_z, topology.n_cells_x, ndof, jj, k, i);
        haloSendBuf_Ym(iGlob) = field.data(ndof, k+ks, jj+js+topology.n_cells_y, i+is);
        haloSendBuf_Yp(iGlob) = field.data(ndof, k+ks, jj+js, i+is);
      });
    }

    //pack bottom (z-) and top (z+)
    if (ndims ==3) {
      yakl::parallel_for("PackBottomTop", bufsize_z, YAKL_LAMBDA (int iGlob) {
        int ndof, i, kk, j;
        yakl::unpackIndices(iGlob, total_dofs, topology.halosize_z, topology.n_cells_y, topology.n_cells_x, ndof, kk, j, i);
        haloSendBuf_Zm(iGlob) = field.data(ndof, kk+ks+topology.n_cells_z, j+js, i+is);
        haloSendBuf_Zp(iGlob) = field.data(ndof, kk+ks, j+js, i+is);
      });
    }


  }

  void PeriodicExchange::unpack(Field &field) {

    int is = topology.is;
    int js = topology.js;
    int ks = topology.ks;

    // pack left (x-) and  right (x+)
    yakl::parallel_for("UnpackLeftRight", bufsize_x, YAKL_LAMBDA (int iGlob) {
      int ndof, ii, k, j;
      yakl::unpackIndices(iGlob, total_dofs, topology.halosize_x, topology.n_cells_z, topology.n_cells_y, ndof, ii, k, j);
      field.data(ndof, k+ks, j+js, ii+is+topology.n_cells_x) = haloRecvBuf_Xm(iGlob);
      field.data(ndof, k+ks, j+js, ii+is) = haloRecvBuf_Xp(iGlob);
    });

    //pack down (y-) and up (y+)
    if (ndims >=2) {
      yakl::parallel_for("UnpackDownUp", bufsize_y, YAKL_LAMBDA (int iGlob) {
        int ndof, i, k, jj;
        yakl::unpackIndices(iGlob, total_dofs, topology.halosize_y, topology.n_cells_z, topology.n_cells_x, ndof, jj, k, i);
        field.data(ndof, k+ks, jj+js+topology.n_cells_y, i+is) = haloRecvBuf_Ym(iGlob);
        field.data(ndof, k+ks, jj+js, i+is) = haloRecvBuf_Yp(iGlob);
      });
    }

    //pack bottom (z-) and top (z+)
    if (ndims ==3) {
      yakl::parallel_for("UnpackBottomTop", bufsize_z, YAKL_LAMBDA (int iGlob) {
        int ndof, i, kk, j;
        yakl::unpackIndices(iGlob, total_dofs, topology.halosize_z, topology.n_cells_y, topology.n_cells_x , ndof, kk, j, i);
        field.data(ndof, kk+ks+topology.n_cells_z, j+js, i+is) = haloRecvBuf_Zm(iGlob);
        field.data(ndof, kk+ks, j+js, i+is) = haloRecvBuf_Zp(iGlob);
      });
    }
  }



  void PeriodicExchange::exchange() {
    yakl::parallel_for( bufsize , YAKL_LAMBDA (int iGlob) {
      haloRecvBuf(iGlob) = haloSendBuf(iGlob);
    });

  }
