
#ifndef __GRID_H__
#define __GRID_H__

#include "const.h"

int constexpr nx = nx_gl/nsubdomains_x;
int constexpr ny = ny_gl/nsubdomains_y;
int constexpr nz = nz_gl+1;        // note that nz_gl = crm_nz
int constexpr nzm = nz-1;          // note that nzm   = crm_nz
int constexpr nsubdomains = nsubdomains_x * nsubdomains_y;
int constexpr RUN3D = ny_gl > 1;
int constexpr RUN2D = !RUN3D;

int constexpr nxp1 = nx + 1;
int constexpr nyp1 = ny + 1 * YES3D;
int constexpr nxp2 = nx + 2;
int constexpr nyp2 = ny + 2 * YES3D;
int constexpr nxp3 = nx + 3;
int constexpr nyp3 = ny + 3 * YES3D;
int constexpr nxp4 = nx + 4;
int constexpr nyp4 = ny + 4 * YES3D;

int constexpr dimx1_u = -1         ;      
int constexpr dimx2_u = nxp3       ;     
int constexpr dimy1_u = 1-(2)*YES3D;
int constexpr dimy2_u = nyp2       ;
int constexpr dimx1_v = -1         ;
int constexpr dimx2_v = nxp2       ;
int constexpr dimy1_v = 1-2*YES3D  ;    
int constexpr dimy2_v = nyp3       ;   
int constexpr dimx1_w = -1         ;
int constexpr dimx2_w = nxp2       ;
int constexpr dimy1_w = 1-(2)*YES3D;
int constexpr dimy2_w = nyp2       ;
int constexpr dimx1_s = -2         ;
int constexpr dimx2_s = nxp3       ;
int constexpr dimy1_s = 1-(3)*YES3D;
int constexpr dimy2_s = nyp3       ;

int constexpr ncols = nx*ny;
int constexpr nadams = 3;

#endif

