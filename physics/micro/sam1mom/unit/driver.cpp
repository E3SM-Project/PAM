
#include "YAKL.h"
#include "YAKL_netcdf.h"

using yakl::fortran::parallel_for;
using yakl::fortran::Bounds;

typedef double real;

typedef yakl::Array<real,2,yakl::memDevice,yakl::styleFortran> real2d;

real2d transpose( real2d const &data );

extern "C"
void sam1mom_main_fortran(double &dt, int &ncol, int &nz, double *zint, double *rho, double *rhow, double *pres, 
                          double *tabs, double *qv, double *qn, double *qp);


int main() {
  yakl::init();
  {
    real2d data;
    real   dt;
    yakl::SimpleNetCDF nc;
    nc.open("sam1mom_data.nc");
    nc.read(data,"qv"         );    real2d qv          = transpose(data);
    nc.read(data,"qn"         );    real2d qn          = transpose(data);
    nc.read(data,"qp"         );    real2d qp          = transpose(data);
    nc.read(data,"zint"       );    real2d zint        = transpose(data);
    nc.read(data,"pressure"   );    real2d pressure    = transpose(data);
    nc.read(data,"temp"       );    real2d temp        = transpose(data);
    nc.read(data,"density"    );    real2d density     = transpose(data);
    nc.read(data,"density_int");    real2d density_int = transpose(data);
    nc.read(dt  ,"dt"         );

    int ncol = qv.dimension[0];
    int nz   = qv.dimension[1];

    auto qv_host          = qv         .createHostCopy();
    auto qn_host          = qn         .createHostCopy();
    auto qp_host          = qp         .createHostCopy();
    auto zint_host        = zint       .createHostCopy();
    auto pressure_host    = pressure   .createHostCopy();
    auto temp_host        = temp       .createHostCopy();
    auto density_host     = density    .createHostCopy();
    auto density_int_host = density_int.createHostCopy();

    sam1mom_main_fortran( dt , ncol , nz , zint_host.data() , density_host.data() , density_int_host.data() ,
                          pressure_host.data() , temp_host.data() , qv_host.data() , qn_host.data() , qp_host.data() );

    auto qv_fortran   = qv_host  .createDeviceCopy();
    auto qn_fortran   = qn_host  .createDeviceCopy();
    auto qp_fortran   = qp_host  .createDeviceCopy();
    auto temp_fortran = temp_host.createDeviceCopy();
  }
  yakl::finalize();
}



real2d transpose( real2d const &data ) {
  int nz   = data.dimension[0];
  int ncol = data.dimension[1];
  real2d ret("ret",ncol,nz);
  parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
    ret(i,k) = data(k,i);
  });
  return ret;
}

