
#include "supercell_init.h"


extern "C" void supercell_init_fortran( double *vert_interface_p , double *rho_d_col_p , double *uvel_col_p ,
                                        double *vvel_col_p , double *wvel_col_p , double *temp_col_p ,
                                        double *rho_v_col_p , double &Rd , double &Rv , double &grav, int &nz) {
  real1d vert_interface("vert_interface",vert_interface_p,nz+1);
  real1d rho_d_col     ("rho_d_col"     ,rho_d_col_p     ,nz  );
  real1d uvel_col      ("uvel_col"      ,uvel_col_p      ,nz  );
  real1d vvel_col      ("vvel_col"      ,vvel_col_p      ,nz  );
  real1d wvel_col      ("wvel_col"      ,wvel_col_p      ,nz  );
  real1d temp_col      ("temp_col"      ,temp_col_p      ,nz  );
  real1d rho_v_col     ("rho_v_col"     ,rho_v_col_p     ,nz  );

  supercell_init(vert_interface, rho_d_col, uvel_col, vvel_col, wvel_col, temp_col, rho_v_col, Rd, Rv, grav);
}


