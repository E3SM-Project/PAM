#pragma once

#include "common.h"

// A LOT OF THIS CAN/SHOULD BE MOVED TO EITHER TOPOLOGY AND/OR GEOMETRY

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
