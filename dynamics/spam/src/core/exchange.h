#pragma once

#include "common.h"
#include "topology.h"

// Xm SendBuf holds the cells closest to x- ie the left side -> Go into Xp
// RecvBuf Xp SendBuf holds the cells closest to x+ ie the right side -> Go into
// Xm RecvBuf

// Xm RecvBuf holds the halo cells closest to x- ie the left side
// Xp SendBuf holds the halo cells closest to x+ ie the right side

class Exchange {
public:
  Topology topology;
  bool single_process;

  int bufsize_x, bufsize_y, bufsize_xy;
  int nens;
  int total_dofs; // total number of dofs in a "horiz" slice ie at each layer or
                  // interface
  // basedof = 0,1,2 in 2D or 0,1 in 1D ie vertices,edges,cells or
  // vertices,cells extdof = 0,1 (0 = interfaces, 1=layers)
  int basedof, extdof, ndofs;
  int _nz; //, _nloop, _nloop_halo;
  int mirror_size;

  realHost1d haloSendBuf_Xm_host;
  realHost1d haloRecvBuf_Xm_host;
  real1d haloSendBuf_Xm;
  real1d haloRecvBuf_Xm;
  realHost1d haloSendBuf_Xp_host;
  realHost1d haloRecvBuf_Xp_host;
  real1d haloSendBuf_Xp;
  real1d haloRecvBuf_Xp;

  realHost1d haloSendBuf_Ym_host;
  realHost1d haloRecvBuf_Ym_host;
  real1d haloSendBuf_Ym;
  real1d haloRecvBuf_Ym;
  realHost1d haloSendBuf_Yp_host;
  realHost1d haloRecvBuf_Yp_host;
  real1d haloSendBuf_Yp;
  real1d haloRecvBuf_Yp;

  realHost1d haloSendBuf_XYll_host;
  realHost1d haloRecvBuf_XYll_host;
  real1d haloSendBuf_XYll;
  real1d haloRecvBuf_XYll;
  realHost1d haloSendBuf_XYul_host;
  realHost1d haloRecvBuf_XYul_host;
  real1d haloSendBuf_XYul;
  real1d haloRecvBuf_XYul;
  realHost1d haloSendBuf_XYlr_host;
  realHost1d haloRecvBuf_XYlr_host;
  real1d haloSendBuf_XYlr;
  real1d haloRecvBuf_XYlr;
  realHost1d haloSendBuf_XYur_host;
  realHost1d haloRecvBuf_XYur_host;
  real1d haloSendBuf_XYur;
  real1d haloRecvBuf_XYur;

  // 4 corners in 2D
  MPI_Request sReq[4];
  MPI_Request rReq[4];

  MPI_Status sStat[4];
  MPI_Status rStat[4];

  bool is_initialized;

  Exchange();
  void printinfo();
  void initialize(const Exchange &exch);
  void initialize(const Topology &topo, int nd0, int nd1, int nd2);
  void pack(const real5d &data);
  void unpack(real5d &data);
  void exchange();
  void exchange_data(real5d &data);
  void exchange_x();
  void exchange_y();
  void exchange_corners();
  void exchange_mirror(real5d &data);
  void exchange_direct_x(real5d &data);
  void exchange_direct_y(real5d &data);
  void exchange_direct(real5d &data);
};

Exchange::Exchange() { this->is_initialized = false; }

void Exchange::printinfo() {
  std::cout << "exchange info\n" << std::flush;
  std::cout << "bufsize_x " << this->bufsize_x << " bufsize_y "
            << this->bufsize_y << "\n"
            << std::flush;
  std::cout << "total_dofs " << this->total_dofs << " basedof " << this->basedof
            << " extdof " << this->extdof << " ndofs " << this->ndofs << "\n"
            << std::flush;
}

void Exchange::initialize(const Exchange &exch) {
  initialize(exch.topology, exch.basedof, exch.extdof, exch.ndofs);
}

void Exchange::initialize(const Topology &topo, int bdof, int edof, int nd) {
  this->topology = topo;
  this->single_process = (topo.nprocx == 1 && topo.nprocy == 1);

  this->basedof = bdof;
  this->extdof = edof;
  this->ndofs = nd;

  this->total_dofs = this->ndofs;
  if (ndims == 2 && this->basedof == 1) {
    this->total_dofs = 2 * this->ndofs;
  } // 2 edges per cell in 2D

  if (this->extdof == 0) {
    this->_nz = this->topology.ni;
    // this->_nloop = this->topology.n_cells_interfaces;
    // this->_nloop_halo = this->topology.n_cells_interfaces_with_halo;
  }
  if (this->extdof == 1) {
    this->_nz = this->topology.nl;
    // this->_nloop = this->topology.n_cells_layers;
    // this->_nloop_halo = this->topology.n_cells_layers_with_halo;
  }
  this->nens = this->topology.nens;

  this->mirror_size = this->topology.n_cells_x * this->topology.n_cells_y *
                      this->topology.mirror_halo * this->nens;

  if (!single_process) {
    this->bufsize_x = this->topology.halosize_x * this->total_dofs *
                      this->topology.n_cells_y * this->_nz * this->nens;
    this->bufsize_y = this->topology.halosize_y * this->total_dofs *
                      this->topology.n_cells_x * this->_nz * this->nens;
    this->bufsize_xy = this->topology.halosize_y * this->total_dofs *
                       this->topology.halosize_x * this->_nz * this->nens;

    this->haloSendBuf_Xm = real1d("haloSendBuf_Xm", this->bufsize_x);
    this->haloRecvBuf_Xm = real1d("haloRecvBuf_Xm", this->bufsize_x);
    this->haloSendBuf_Xm_host = this->haloSendBuf_Xm.createHostCopy();
    this->haloRecvBuf_Xm_host = this->haloRecvBuf_Xm.createHostCopy();
    this->haloSendBuf_Xp = real1d("haloSendBuf_Xp", this->bufsize_x);
    this->haloRecvBuf_Xp = real1d("haloRecvBuf_Xp", this->bufsize_x);
    this->haloSendBuf_Xp_host = this->haloSendBuf_Xp.createHostCopy();
    this->haloRecvBuf_Xp_host = this->haloRecvBuf_Xp.createHostCopy();

    if (ndims == 2) {
      this->haloSendBuf_Ym = real1d("haloSendBuf_Ym", this->bufsize_y);
      this->haloRecvBuf_Ym = real1d("haloRecvBuf_Ym", this->bufsize_y);
      this->haloSendBuf_Ym_host = this->haloSendBuf_Ym.createHostCopy();
      this->haloRecvBuf_Ym_host = this->haloRecvBuf_Ym.createHostCopy();
      this->haloSendBuf_Yp = real1d("haloSendBuf_Yp", this->bufsize_y);
      this->haloRecvBuf_Yp = real1d("haloRecvBuf_Yp", this->bufsize_y);
      this->haloSendBuf_Yp_host = this->haloSendBuf_Yp.createHostCopy();
      this->haloRecvBuf_Yp_host = this->haloRecvBuf_Yp.createHostCopy();
    }

    if (ndims == 2) {

      this->haloSendBuf_XYll = real1d("haloSendBuf_XYll", this->bufsize_xy);
      this->haloRecvBuf_XYll = real1d("haloRecvBuf_XYll", this->bufsize_xy);
      this->haloSendBuf_XYll_host = this->haloSendBuf_XYll.createHostCopy();
      this->haloRecvBuf_XYll_host = this->haloRecvBuf_XYll.createHostCopy();

      this->haloSendBuf_XYul = real1d("haloSendBuf_XYul", this->bufsize_xy);
      this->haloRecvBuf_XYul = real1d("haloRecvBuf_XYul", this->bufsize_xy);
      this->haloSendBuf_XYul_host = this->haloSendBuf_XYul.createHostCopy();
      this->haloRecvBuf_XYul_host = this->haloRecvBuf_XYul.createHostCopy();

      this->haloSendBuf_XYlr = real1d("haloSendBuf_XYlr", this->bufsize_xy);
      this->haloRecvBuf_XYlr = real1d("haloRecvBuf_XYlr", this->bufsize_xy);
      this->haloSendBuf_XYlr_host = this->haloSendBuf_XYlr.createHostCopy();
      this->haloRecvBuf_XYlr_host = this->haloRecvBuf_XYlr.createHostCopy();

      this->haloSendBuf_XYur = real1d("haloSendBuf_XYur", this->bufsize_xy);
      this->haloRecvBuf_XYur = real1d("haloRecvBuf_XYur", this->bufsize_xy);
      this->haloSendBuf_XYur_host = this->haloSendBuf_XYur.createHostCopy();
      this->haloRecvBuf_XYur_host = this->haloRecvBuf_XYur.createHostCopy();
    }
  }

  this->is_initialized = true;
}

void Exchange::pack(const real5d &data) {

  int is = this->topology.is;
  int js = this->topology.js;
  int ks = this->topology.ks;

  YAKL_SCOPE(_nz, this->_nz);
  YAKL_SCOPE(n_cells_x, this->topology.n_cells_x);
  YAKL_SCOPE(n_cells_y, this->topology.n_cells_y);
  YAKL_SCOPE(halosize_x, this->topology.halosize_x);
  YAKL_SCOPE(halosize_y, this->topology.halosize_y);

  YAKL_SCOPE(haloSendBuf_Xp, this->haloSendBuf_Xp);
  YAKL_SCOPE(haloSendBuf_Xm, this->haloSendBuf_Xm);

  // pack left (x-) and  right (x+)
  parallel_for(
      "pack x",
      SimpleBounds<5>(this->total_dofs, this->topology.halosize_x, this->_nz,
                      this->topology.n_cells_y, this->nens),
      YAKL_LAMBDA(int ndof, int ii, int k, int j, int n) {
        int iGlob =
            n + nens * (j + n_cells_y * (k + _nz * (ii + halosize_x * ndof)));
        haloSendBuf_Xp(iGlob) =
            data(ndof, k + ks, j + js, ii + is + n_cells_x - halosize_x, n);
        haloSendBuf_Xm(iGlob) = data(ndof, k + ks, j + js, ii + is, n);
      });

  // pack down (y-) and up (y+)
  if (ndims == 2) {
    YAKL_SCOPE(haloSendBuf_Yp, this->haloSendBuf_Yp);
    YAKL_SCOPE(haloSendBuf_Ym, this->haloSendBuf_Ym);

    parallel_for(
        "pack y",
        SimpleBounds<5>(this->total_dofs, this->topology.halosize_y, this->_nz,
                        this->topology.n_cells_x, this->nens),
        YAKL_LAMBDA(int ndof, int jj, int k, int i, int n) {
          int iGlob =
              n + nens * (i + n_cells_x * (k + _nz * (jj + halosize_y * ndof)));
          haloSendBuf_Yp(iGlob) =
              data(ndof, k + ks, jj + js + n_cells_y - halosize_y, i + is, n);
          haloSendBuf_Ym(iGlob) = data(ndof, k + ks, jj + js, i + is, n);
        });
  }

  if (ndims == 2) {
    // pack corners

    YAKL_SCOPE(haloSendBuf_XYll, this->haloSendBuf_XYll);
    YAKL_SCOPE(haloSendBuf_XYur, this->haloSendBuf_XYur);
    YAKL_SCOPE(haloSendBuf_XYul, this->haloSendBuf_XYul);
    YAKL_SCOPE(haloSendBuf_XYlr, this->haloSendBuf_XYlr);

    parallel_for(
        "pack corners",
        SimpleBounds<5>(this->total_dofs, this->topology.halosize_y,
                        this->topology.halosize_x, this->_nz, this->nens),
        YAKL_LAMBDA(int ndof, int jj, int ii, int k, int n) {
          int iGlob =
              n +
              nens * (k + _nz * (ii + halosize_x * (jj + halosize_y * ndof)));
          haloSendBuf_XYll(iGlob) = data(ndof, k + ks, jj + js, ii + is, n);
          haloSendBuf_XYur(iGlob) =
              data(ndof, k + ks, jj + js + n_cells_y - halosize_y,
                   ii + is + n_cells_x - halosize_x, n);
          haloSendBuf_XYul(iGlob) =
              data(ndof, k + ks, jj + js + n_cells_y - halosize_y, ii + is, n);
          haloSendBuf_XYlr(iGlob) =
              data(ndof, k + ks, jj + js, ii + is + n_cells_x - halosize_x, n);
        });
  }
}

void Exchange::unpack(real5d &data) {

  int is = this->topology.is;
  int js = this->topology.js;
  int ks = this->topology.ks;

  YAKL_SCOPE(_nz, this->_nz);
  YAKL_SCOPE(n_cells_x, this->topology.n_cells_x);
  YAKL_SCOPE(n_cells_y, this->topology.n_cells_y);
  YAKL_SCOPE(halosize_x, this->topology.halosize_x);
  YAKL_SCOPE(halosize_y, this->topology.halosize_y);

  YAKL_SCOPE(haloRecvBuf_Xp, this->haloRecvBuf_Xp);
  YAKL_SCOPE(haloRecvBuf_Xm, this->haloRecvBuf_Xm);

  // unpack left (x-) and  right (x+)
  parallel_for(
      "unpack x",
      SimpleBounds<5>(this->total_dofs, this->topology.halosize_x, this->_nz,
                      this->topology.n_cells_y, this->nens),
      YAKL_LAMBDA(int ndof, int ii, int k, int j, int n) {
        int iGlob =
            n + nens * (j + n_cells_y * (k + _nz * (ii + halosize_x * ndof)));
        data(ndof, k + ks, j + js, ii + is - halosize_x, n) =
            haloRecvBuf_Xm(iGlob);
        data(ndof, k + ks, j + js, ii + is + n_cells_x, n) =
            haloRecvBuf_Xp(iGlob);
      });

  // unpack down (y-) and up (y+)
  if (ndims == 2) {
    YAKL_SCOPE(haloRecvBuf_Yp, this->haloRecvBuf_Yp);
    YAKL_SCOPE(haloRecvBuf_Ym, this->haloRecvBuf_Ym);

    parallel_for(
        "unpack y",
        SimpleBounds<5>(this->total_dofs, this->topology.halosize_y, this->_nz,
                        this->topology.n_cells_x, this->nens),
        YAKL_LAMBDA(int ndof, int jj, int k, int i, int n) {
          int iGlob =
              n + nens * (i + n_cells_x * (k + _nz * (jj + halosize_y * ndof)));
          data(ndof, k + ks, jj + js - halosize_y, i + is, n) =
              haloRecvBuf_Ym(iGlob);
          data(ndof, k + ks, jj + js + n_cells_y, i + is, n) =
              haloRecvBuf_Yp(iGlob);
        });
  }

  if (ndims == 2) {
    YAKL_SCOPE(haloRecvBuf_XYll, this->haloRecvBuf_XYll);
    YAKL_SCOPE(haloRecvBuf_XYur, this->haloRecvBuf_XYur);
    YAKL_SCOPE(haloRecvBuf_XYul, this->haloRecvBuf_XYul);
    YAKL_SCOPE(haloRecvBuf_XYlr, this->haloRecvBuf_XYlr);

    // unpack corners
    parallel_for(
        "unpack corners",
        SimpleBounds<5>(this->total_dofs, this->topology.halosize_y,
                        this->topology.halosize_x, this->_nz, this->nens),
        YAKL_LAMBDA(int ndof, int jj, int ii, int k, int n) {
          int iGlob =
              n +
              nens * (k + _nz * (ii + halosize_x * (jj + halosize_y * ndof)));
          data(ndof, k + ks, jj + js - halosize_y, ii + is - halosize_x, n) =
              haloRecvBuf_XYll(iGlob);
          data(ndof, k + ks, jj + js + n_cells_y, ii + is + n_cells_x, n) =
              haloRecvBuf_XYur(iGlob);
          data(ndof, k + ks, jj + js - halosize_y, ii + is + n_cells_x, n) =
              haloRecvBuf_XYlr(iGlob);
          data(ndof, k + ks, jj + js + n_cells_y, ii + is - halosize_x, n) =
              haloRecvBuf_XYul(iGlob);
        });
  }
}

void Exchange::exchange_x() {
  int ierr;

  if (this->topology.nprocx > 1) {
    yakl::fence();

    // Pre-post the receives
    ierr =
        MPI_Irecv(haloRecvBuf_Xm_host.data(), this->bufsize_x, REAL_MPI,
                  this->topology.x_neigh(0), 0, MPI_COMM_WORLD, &this->rReq[0]);
    ierr =
        MPI_Irecv(haloRecvBuf_Xp_host.data(), this->bufsize_x, REAL_MPI,
                  this->topology.x_neigh(1), 1, MPI_COMM_WORLD, &this->rReq[1]);

    // Copy send buffers to host
    haloSendBuf_Xm.deep_copy_to(haloSendBuf_Xm_host);
    haloSendBuf_Xp.deep_copy_to(haloSendBuf_Xp_host);
    yakl::fence();

    // Send the data
    ierr =
        MPI_Isend(haloSendBuf_Xm_host.data(), this->bufsize_x, REAL_MPI,
                  this->topology.x_neigh(0), 1, MPI_COMM_WORLD, &this->sReq[0]);
    ierr =
        MPI_Isend(haloSendBuf_Xp_host.data(), this->bufsize_x, REAL_MPI,
                  this->topology.x_neigh(1), 0, MPI_COMM_WORLD, &this->sReq[1]);

    // Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, this->sReq, this->sStat);
    ierr = MPI_Waitall(2, this->rReq, this->rStat);

    // Copy recv buffers to device
    haloRecvBuf_Xm_host.deep_copy_to(haloRecvBuf_Xm);
    haloRecvBuf_Xp_host.deep_copy_to(haloRecvBuf_Xp);
  }

  else {

    YAKL_SCOPE(haloSendBuf_Xp, this->haloSendBuf_Xp);
    YAKL_SCOPE(haloSendBuf_Xm, this->haloSendBuf_Xm);
    YAKL_SCOPE(haloRecvBuf_Xp, this->haloRecvBuf_Xp);
    YAKL_SCOPE(haloRecvBuf_Xm, this->haloRecvBuf_Xm);

    parallel_for(
        "exchange buffers x", SimpleBounds<1>(this->bufsize_x),
        YAKL_LAMBDA(int iGlob) {
          haloRecvBuf_Xp(iGlob) = haloSendBuf_Xm(iGlob);
          haloRecvBuf_Xm(iGlob) = haloSendBuf_Xp(iGlob);
        });
  }
}

void Exchange::exchange_direct_x(real5d &data) {
  int is = this->topology.is;
  int js = this->topology.js;
  int ks = this->topology.ks;

  YAKL_SCOPE(_nz, this->_nz);
  YAKL_SCOPE(n_cells_x, this->topology.n_cells_x);
  YAKL_SCOPE(halosize_x, this->topology.halosize_x);

  parallel_for(
      "exchange direct x",
      SimpleBounds<5>(this->total_dofs, this->topology.halosize_x, this->_nz,
                      this->topology.n_cells_y, this->nens),
      YAKL_LAMBDA(int ndof, int ii, int k, int j, int n) {
        data(ndof, k + ks, j + js, ii + is - halosize_x, n) =
            data(ndof, k + ks, j + js, ii + is + n_cells_x - halosize_x, n);
        data(ndof, k + ks, j + js, ii + is + n_cells_x, n) =
            data(ndof, k + ks, j + js, ii + is, n);
      });
}
void Exchange::exchange_direct_y(real5d &data) {
  int is = this->topology.is;
  int js = this->topology.js;
  int ks = this->topology.ks;

  YAKL_SCOPE(n_cells_y, this->topology.n_cells_y);
  YAKL_SCOPE(halosize_y, this->topology.halosize_y);

  parallel_for(
      "exchange direct y",
      SimpleBounds<5>(this->total_dofs, this->topology.halosize_y, this->_nz,
                      this->topology.n_cells_x + 2 * this->topology.halosize_x,
                      this->nens),
      YAKL_LAMBDA(int ndof, int jj, int k, int i, int n) {
        data(ndof, k + ks, jj + js - halosize_y, i, n) =
            data(ndof, k + ks, jj + js + n_cells_y - halosize_y, i, n);
        data(ndof, k + ks, jj + js + n_cells_y, i, n) =
            data(ndof, k + ks, jj + js, i, n);
      });
}

void Exchange::exchange_direct(real5d &data) {
  exchange_direct_x(data);
  if (ndims == 2) {
    exchange_direct_y(data);
  }
}

void Exchange::exchange_y() {
  int ierr;

  if (this->topology.nprocy > 1) {

    yakl::fence();

    // Pre-post the receives
    ierr =
        MPI_Irecv(haloRecvBuf_Ym_host.data(), this->bufsize_y, REAL_MPI,
                  this->topology.y_neigh(0), 0, MPI_COMM_WORLD, &this->rReq[0]);
    ierr =
        MPI_Irecv(haloRecvBuf_Yp_host.data(), this->bufsize_y, REAL_MPI,
                  this->topology.y_neigh(1), 1, MPI_COMM_WORLD, &this->rReq[1]);

    // Copy send buffers to host
    haloSendBuf_Ym.deep_copy_to(haloSendBuf_Ym_host);
    haloSendBuf_Yp.deep_copy_to(haloSendBuf_Yp_host);
    yakl::fence();

    // Send the data
    ierr =
        MPI_Isend(haloSendBuf_Ym_host.data(), this->bufsize_y, REAL_MPI,
                  this->topology.y_neigh(0), 1, MPI_COMM_WORLD, &this->sReq[0]);
    ierr =
        MPI_Isend(haloSendBuf_Yp_host.data(), this->bufsize_y, REAL_MPI,
                  this->topology.y_neigh(1), 0, MPI_COMM_WORLD, &this->sReq[1]);

    // Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, this->sReq, this->sStat);
    ierr = MPI_Waitall(2, this->rReq, this->rStat);

    // Copy recv buffers to device
    haloRecvBuf_Ym_host.deep_copy_to(haloRecvBuf_Ym);
    haloRecvBuf_Yp_host.deep_copy_to(haloRecvBuf_Yp);
  }

  else {
    YAKL_SCOPE(haloSendBuf_Yp, this->haloSendBuf_Yp);
    YAKL_SCOPE(haloSendBuf_Ym, this->haloSendBuf_Ym);
    YAKL_SCOPE(haloRecvBuf_Yp, this->haloRecvBuf_Yp);
    YAKL_SCOPE(haloRecvBuf_Ym, this->haloRecvBuf_Ym);
    parallel_for(
        "exchange buffers y", SimpleBounds<1>(this->bufsize_y),
        YAKL_LAMBDA(int iGlob) {
          haloRecvBuf_Yp(iGlob) = haloSendBuf_Ym(iGlob);
          haloRecvBuf_Ym(iGlob) = haloSendBuf_Yp(iGlob);
        });
  }
}

void Exchange::exchange_corners() {
  int ierr;

  if (this->topology.nprocx > 1 || this->topology.nprocy > 1) {
    yakl::fence();

    // Pre-post the receives
    ierr =
        MPI_Irecv(haloRecvBuf_XYll_host.data(), this->bufsize_xy, REAL_MPI,
                  this->topology.ur_neigh, 0, MPI_COMM_WORLD, &this->rReq[0]);
    ierr =
        MPI_Irecv(haloRecvBuf_XYur_host.data(), this->bufsize_xy, REAL_MPI,
                  this->topology.ll_neigh, 1, MPI_COMM_WORLD, &this->rReq[1]);
    ierr =
        MPI_Irecv(haloRecvBuf_XYul_host.data(), this->bufsize_xy, REAL_MPI,
                  this->topology.lr_neigh, 2, MPI_COMM_WORLD, &this->rReq[2]);
    ierr =
        MPI_Irecv(haloRecvBuf_XYlr_host.data(), this->bufsize_xy, REAL_MPI,
                  this->topology.ul_neigh, 3, MPI_COMM_WORLD, &this->rReq[3]);

    // Copy send buffers to host
    haloSendBuf_XYll.deep_copy_to(haloSendBuf_XYll_host);
    haloSendBuf_XYur.deep_copy_to(haloSendBuf_XYur_host);
    haloSendBuf_XYul.deep_copy_to(haloSendBuf_XYul_host);
    haloSendBuf_XYlr.deep_copy_to(haloSendBuf_XYlr_host);
    yakl::fence();

    // Send the data
    ierr =
        MPI_Isend(haloSendBuf_XYll_host.data(), this->bufsize_xy, REAL_MPI,
                  this->topology.ur_neigh, 1, MPI_COMM_WORLD, &this->sReq[0]);
    ierr =
        MPI_Isend(haloSendBuf_XYur_host.data(), this->bufsize_xy, REAL_MPI,
                  this->topology.ll_neigh, 0, MPI_COMM_WORLD, &this->sReq[1]);
    ierr =
        MPI_Isend(haloSendBuf_XYul_host.data(), this->bufsize_xy, REAL_MPI,
                  this->topology.lr_neigh, 3, MPI_COMM_WORLD, &this->sReq[2]);
    ierr =
        MPI_Isend(haloSendBuf_XYlr_host.data(), this->bufsize_xy, REAL_MPI,
                  this->topology.ul_neigh, 2, MPI_COMM_WORLD, &this->sReq[3]);

    // Wait for the sends and receives to finish
    ierr = MPI_Waitall(4, this->sReq, this->sStat);
    ierr = MPI_Waitall(4, this->rReq, this->rStat);

    // Copy recv buffers to device
    haloRecvBuf_XYll_host.deep_copy_to(haloRecvBuf_XYll);
    haloRecvBuf_XYur_host.deep_copy_to(haloRecvBuf_XYur);
    haloRecvBuf_XYul_host.deep_copy_to(haloRecvBuf_XYul);
    haloRecvBuf_XYlr_host.deep_copy_to(haloRecvBuf_XYlr);
  }

  else {
    YAKL_SCOPE(haloSendBuf_XYll, this->haloSendBuf_XYll);
    YAKL_SCOPE(haloSendBuf_XYur, this->haloSendBuf_XYur);
    YAKL_SCOPE(haloSendBuf_XYul, this->haloSendBuf_XYul);
    YAKL_SCOPE(haloSendBuf_XYlr, this->haloSendBuf_XYlr);
    YAKL_SCOPE(haloRecvBuf_XYll, this->haloRecvBuf_XYll);
    YAKL_SCOPE(haloRecvBuf_XYur, this->haloRecvBuf_XYur);
    YAKL_SCOPE(haloRecvBuf_XYul, this->haloRecvBuf_XYul);
    YAKL_SCOPE(haloRecvBuf_XYlr, this->haloRecvBuf_XYlr);

    parallel_for(
        "exchage buffers corners", SimpleBounds<1>(this->bufsize_xy),
        YAKL_LAMBDA(int iGlob) {
          haloRecvBuf_XYll(iGlob) = haloSendBuf_XYur(iGlob);
          haloRecvBuf_XYur(iGlob) = haloSendBuf_XYll(iGlob);
          haloRecvBuf_XYul(iGlob) = haloSendBuf_XYlr(iGlob);
          haloRecvBuf_XYlr(iGlob) = haloSendBuf_XYul(iGlob);
        });
  }
}

void Exchange::exchange_mirror(real5d &data) {
  int is = this->topology.is;
  int js = this->topology.js;
  int ks = this->topology.ks;

  YAKL_SCOPE(_nz, this->_nz);
  // vertical layers
  if (this->extdof == 1) {
    parallel_for(
        "exchange mirror",
        SimpleBounds<5>(this->total_dofs, this->topology.mirror_halo,
                        this->topology.n_cells_x, this->topology.n_cells_y,
                        this->nens),
        YAKL_LAMBDA(int ndof, int kk, int i, int j, int n) {
          // Mirror Top
          data(ndof, _nz + ks + kk, j + js, i + is, n) =
              data(ndof, _nz + ks - kk - 1, j + js, i + is, n);
          // Mirror Bottom
          data(ndof, ks - kk - 1, j + js, i + is, n) =
              data(ndof, ks + kk, j + js, i + is, n);
        });
  }

  // vertical interfaces
  if (this->extdof == 0) {
    parallel_for(
        "exchange mirror",
        SimpleBounds<5>(this->total_dofs, this->topology.mirror_halo,
                        this->topology.n_cells_x, this->topology.n_cells_y,
                        this->nens),
        YAKL_LAMBDA(int ndof, int kk, int i, int j, int n) {
          // Mirror Top
          data(ndof, _nz + ks + kk, j + js, i + is, n) =
              data(ndof, _nz + ks - kk - 2, j + js, i + is, n);
          // Mirror Bottom
          data(ndof, ks - kk - 1, j + js, i + is, n) =
              data(ndof, ks + kk + 1, j + js, i + is, n);
        });
  }
}

// EVENTUALLY WE SHOULD BE MORE CLEVER HERE IE GROUP ALL THE SEND/RECVS IN X/Y
// TOGETHER, ETC.
void Exchange::exchange() {
  exchange_x();
  if (ndims == 2) {
    exchange_y();
  }
  if (ndims == 2) {
    exchange_corners();
  }
}

void Exchange::exchange_data(real5d &data) {
  if (this->single_process) {
    this->exchange_direct(data);
  } else {
    yakl::timer_start("pack");
    this->pack(data);
    yakl::timer_stop("pack");
    yakl::timer_start("exchange");
    this->exchange();
    yakl::timer_stop("exchange");
    yakl::timer_start("unpack");
    this->unpack(data);
    yakl::timer_stop("unpack");
  }

#ifdef _EXTRUDED
  yakl::timer_start("exchange_mirror");
  exchange_mirror(data);
  yakl::timer_stop("exchange_mirror");
#endif
}
