#pragma once

#include "common.h"
#include "params.h"

int ij_to_l(int i, int j, int nx) { return j * nx + i; }

int wrap(int i, int nx) {
  if (i < 0) {
    return i + nx;
  } else if (i > nx - 1) {
    return i - nx;
  } else {
    return i;
  }
}

class Parallel {
public:
  int nx = -1;
  int ny = -1;
  int nz = -1;
  int nens = -1;

  int nprocx = -1;
  int nprocy = -1;

  int nranks = -1;
  int myrank = -1;
  int actualrank = -1;
  bool masterproc;
  int px = -1;
  int py = -1;

  int i_beg = -1;
  int i_end = -1;
  int j_beg = -1;
  int j_end = -1;

  int halox = -1;
  int haloy = -1;

  SArray<int, 1, 2> x_neigh;
  SArray<int, 1, 2> y_neigh;

  int ll_neigh = -1;
  int ur_neigh = -1;
  int ul_neigh = -1;
  int lr_neigh = -1;

  int nx_glob = -1;
  int ny_glob = -1;

  BND_TYPE xbnd;
  BND_TYPE ybnd;
};
