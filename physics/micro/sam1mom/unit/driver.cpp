
#include "YAKL.h"
#include "YAKL_netcdf.h"

using yakl::fortran::parallel_for;
using yakl::fortran::Bounds;
using yakl::intrinsics::size;
using yakl::intrinsics::shape;

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
    nc.read(data,"qv"         );    real2d qv          = transpose(data);   data = real2d();
    nc.read(data,"qn"         );    real2d qn          = transpose(data);   data = real2d();
    nc.read(data,"qp"         );    real2d qp          = transpose(data);   data = real2d();
    nc.read(data,"zint"       );    real2d zint        = transpose(data);   data = real2d();
    nc.read(data,"pressure"   );    real2d pressure    = transpose(data);   data = real2d();
    nc.read(data,"temp"       );    real2d temp        = transpose(data);   data = real2d();
    nc.read(data,"density"    );    real2d density     = transpose(data);   data = real2d();
    nc.read(data,"density_int");    real2d density_int = transpose(data);   data = real2d();
    nc.read(dt  ,"dt"         );

    int ncol = size(qv,1);
    int nz   = size(qv,2);

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
  int dim1 = size(data,1);
  int dim2 = size(data,2);
  real2d ret("ret",dim2,dim1);
  parallel_for( Bounds<2>(dim1,dim2) , YAKL_LAMBDA (int i1, int i2) {
    ret(i2,i1) = data(i1,i2);
  });
  return ret;
}



