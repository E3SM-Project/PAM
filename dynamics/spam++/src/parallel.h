#ifndef _PARALLEL_H_
#define _PARALLEL_H_

#include "mpi.h"

class Parallel
{
public:
  int nx = -1;
  int ny = -1;
  int nz = -1;

  int nprocx = -1;
  int nprocy = -1;

  int nranks = -1;
  int myrank = -1;
  int masterproc = -1;
  int px = -1;
  int py = -1;

  int i_beg = -1;
  int i_end = -1;
  int j_beg = -1;
  int j_end = -1;

  int halox = -1;
  int haloy = -1;

  SArray<int,2> x_neigh;
  SArray<int,2> y_neigh;
  SArray<int,2> z_neigh;
  
  int ll_neigh, ur_neigh, ul_neigh, lr_neigh;

  int nx_glob = -1;
  int ny_glob = -1;

  BND_TYPE xbnd;
  BND_TYPE ybnd;

};


void parallel_logging()
{

}

#endif
