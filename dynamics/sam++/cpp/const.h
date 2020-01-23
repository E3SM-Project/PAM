
#ifndef __CONST_H__
#define __CONST_H__

#include "YAKL.h"
#include "Array.h"

typedef double real;

int  constexpr ncrms      = NCRMS;
int  constexpr crm_nx     = CRM_NX;
int  constexpr crm_ny     = CRM_NY;
int  constexpr crm_nz     = CRM_NZ;
int  constexpr crm_nx_rad = CRM_NX_RAD;
int  constexpr crm_ny_rad = CRM_NY_RAD;
real constexpr crm_dx     = CRM_DX;
real constexpr crm_dy     = crm_dx;
real constexpr crm_dt     = CRM_DT;

int  constexpr YES3D = YES3DVAL  ; // Domain dimensionality: 1 - 3D, 0 - 2D
int  constexpr nx_gl = crm_nx    ; // Number of grid points in X
int  constexpr ny_gl = crm_ny    ; // Number of grid points in Y
int  constexpr nz_gl = crm_nz    ; // Number of pressure (scalar) levels
int  constexpr nsubdomains_x  = 1; // No of subdomains in x
int  constexpr nsubdomains_y  = 1; // No of subdomains in y
int  constexpr navgmom_x = -1    ; 
int  constexpr navgmom_y = -1    ; 

int  constexpr nx = nx_gl/nsubdomains_x;
int  constexpr ny = ny_gl/nsubdomains_y;
int  constexpr nz = nz_gl+1            ;                     //  note that nz_gl = crm_nz
int  constexpr nzm = nz-1              ;                     //  note that nzm   = crm_nz
int  constexpr nsubdomains = nsubdomains_x * nsubdomains_y;
int  constexpr RUN3D = ny_gl > 1;
int  constexpr RUN2D = ! RUN3D;
int  constexpr nxp1 = nx + 1;
int  constexpr nyp1 = ny + 1 * YES3D;
int  constexpr nxp2 = nx + 2;
int  constexpr nyp2 = ny + 2 * YES3D;
int  constexpr nxp3 = nx + 3;
int  constexpr nyp3 = ny + 3 * YES3D;
int  constexpr nxp4 = nx + 4;
int  constexpr nyp4 = ny + 4 * YES3D;
int  constexpr dimx1_u = -1          ;
int  constexpr dimx2_u = nxp3        ;
int  constexpr dimy1_u = 1-(2)*YES3D ;
int  constexpr dimy2_u = nyp2        ;
int  constexpr dimx1_v = -1          ;
int  constexpr dimx2_v = nxp2        ;
int  constexpr dimy1_v = 1-2*YES3D   ;
int  constexpr dimy2_v = nyp3        ;
int  constexpr dimx1_w = -1          ;
int  constexpr dimx2_w = nxp2        ;
int  constexpr dimy1_w = 1-(2)*YES3D ;
int  constexpr dimy2_w = nyp2        ;
int  constexpr dimx1_s = -2          ;
int  constexpr dimx2_s = nxp3        ;
int  constexpr dimy1_s = 1-(3)*YES3D ;
int  constexpr dimy2_s = nyp3        ;
int  constexpr dimx1_d = 0           ;
int  constexpr dimx2_d = nxp1        ;
int  constexpr dimy1_d = 1-YES3D     ;
int  constexpr dimy2_d = nyp1        ;
int  constexpr ncols = nx*ny;
int  constexpr nadams = 3;

int  constexpr dimx_u     = dimx2_u - dimx1_u   + 1;
int  constexpr dimx_v     = dimx2_v - dimx1_v   + 1;
int  constexpr dimx_w     = dimx2_w - dimx1_w   + 1;
int  constexpr dimx_s     = dimx2_s - dimx1_s   + 1;
int  constexpr dimx_d     = dimx2_d - dimx1_d   + 1;

int  constexpr dimy_u     = dimy2_u - dimy1_u   + 1;
int  constexpr dimy_v     = dimy2_v - dimy1_v   + 1;
int  constexpr dimy_w     = dimy2_w - dimy1_w   + 1;
int  constexpr dimy_s     = dimy2_s - dimy1_s   + 1;
int  constexpr dimy_d     = dimy2_d - dimy1_d   + 1;
int  constexpr dimy_p     = ny      - (1-YES3D) + 1;
int  constexpr dimy_tk2   = nyp1    - (1-YES3D) + 1;
int  constexpr dimy_sstxy = ny      - (1-YES3D) + 1;
int  constexpr dimy_fcory = ny      - 0         + 1;

int  constexpr offx_u     = 1 - dimx1_u  ;
int  constexpr offx_v     = 1 - dimx1_v  ;
int  constexpr offx_w     = 1 - dimx1_w  ;
int  constexpr offx_s     = 1 - dimx1_s  ;
int  constexpr offx_d     = 1 - dimx1_d  ;

int  constexpr offy_u     = 1 - dimy1_u  ;
int  constexpr offy_v     = 1 - dimy1_v  ;
int  constexpr offy_w     = 1 - dimy1_w  ;
int  constexpr offy_s     = 1 - dimy1_s  ;
int  constexpr offy_d     = 1 - dimy1_d  ;
int  constexpr offy_p     = 1 - (1-YES3D);
int  constexpr offy_tk2   = 1 - (1-YES3D);
int  constexpr offy_sstxy = 1 - (1-YES3D);
int  constexpr offy_fcory = 1 - 0        ;


real constexpr cp    = 1004.          ;  // Specific heat of air, J/kg/K
real constexpr ggr   = 9.81           ;  // Gravity acceleration, m/s2
real constexpr lcond = 2.5104e+06     ;  // Latent heat of condensation, J/kg
real constexpr lfus  = 0.3336e+06     ;  // Latent heat of fusion, J/kg
real constexpr lsub  = 2.8440e+06     ;  // Latent heat of sublimation, J/kg
real constexpr rv    = 461.           ;  // Gas constant for water vapor, J/kg/K
real constexpr rgas  = 287.           ;  // Gas constant for dry air, J/kg/K
real constexpr diffelq = 2.21e-05     ;  // Diffusivity of water vapor, m2/s
real constexpr therco = 2.40e-02      ;  // Thermal conductivity of air, J/m/s/K
real constexpr muelq = 1.717e-05      ;  // Dynamic viscosity of air
real constexpr fac_cond = lcond/cp    ;
real constexpr fac_fus  = lfus/cp     ;
real constexpr fac_sub  = lsub/cp     ;
real constexpr pi = 3.141592653589793 ;  // sine, cosine, cosine, sine, 3.14159 !


int  constexpr nsgs_fields = 1;         // total number of prognostic sgs vars
int  constexpr nsgs_fields_diag = 2;    // total number of diagnostic sgs vars
bool constexpr do_sgsdiag_bound = true; // exchange boundaries for diagnostics fields

typedef yakl::Array<real,yakl::memDevice> umgReal1d;
typedef yakl::Array<real,yakl::memDevice> umgReal2d;
typedef yakl::Array<real,yakl::memDevice> umgReal3d;
typedef yakl::Array<real,yakl::memDevice> umgReal4d;
typedef yakl::Array<real,yakl::memDevice> umgReal5d;
typedef yakl::Array<real,yakl::memDevice> umgReal6d;
typedef yakl::Array<real,yakl::memDevice> umgReal7d;

typedef yakl::Array<real,yakl::memDevice> real1d;
typedef yakl::Array<real,yakl::memDevice> real2d;
typedef yakl::Array<real,yakl::memDevice> real3d;
typedef yakl::Array<real,yakl::memDevice> real4d;
typedef yakl::Array<real,yakl::memDevice> real5d;
typedef yakl::Array<real,yakl::memDevice> real6d;
typedef yakl::Array<real,yakl::memDevice> real7d;

typedef yakl::Array<real,yakl::memHost> realHost1d;
typedef yakl::Array<real,yakl::memHost> realHost2d;
typedef yakl::Array<real,yakl::memHost> realHost3d;
typedef yakl::Array<real,yakl::memHost> realHost4d;
typedef yakl::Array<real,yakl::memHost> realHost5d;
typedef yakl::Array<real,yakl::memHost> realHost6d;
typedef yakl::Array<real,yakl::memHost> realHost7d;


#endif

