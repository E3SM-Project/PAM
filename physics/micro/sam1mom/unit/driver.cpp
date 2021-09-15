
#include "YAKL.h"
#include "YAKL_netcdf.h"
#include "Sam1mom.h"

using yakl::intrinsics::size;
using yakl::intrinsics::shape;

typedef sam1mom::real1d real1d_f;
typedef sam1mom::real2d real2d_f;
typedef sam1mom::real3d real3d_f;


real2d_f transpose( real2d_f const &data );

extern "C"
void sam1mom_main_fortran(double &dt, int &ncol, int &nz, double *zint, double *rho, double *rhow, double *pres, 
                          double *tabs, double *qv, double *qn, double *qp);



int main() {
  yakl::init();
  {
    real2d_f qv_fortran  ;
    real2d_f qn_fortran  ;
    real2d_f qp_fortran  ;
    real2d_f temp_fortran;

    real2d_f qv_cxx  ;
    real2d_f qn_cxx  ;
    real2d_f qp_cxx  ;
    real2d_f temp_cxx;

    {
      real2d_f data;
      real   dt;
      yakl::SimpleNetCDF nc;
      nc.open("sam1mom_data.nc");
      nc.read(data,"qv"         );    real2d_f qv          = transpose(data);   data = real2d_f();
      nc.read(data,"qn"         );    real2d_f qn          = transpose(data);   data = real2d_f();
      nc.read(data,"qp"         );    real2d_f qp          = transpose(data);   data = real2d_f();
      nc.read(data,"zint"       );    real2d_f zint        = transpose(data);   data = real2d_f();
      nc.read(data,"pressure"   );    real2d_f pressure    = transpose(data);   data = real2d_f();
      nc.read(data,"temp"       );    real2d_f temp        = transpose(data);   data = real2d_f();
      nc.read(data,"density"    );    real2d_f density     = transpose(data);   data = real2d_f();
      nc.read(data,"density_int");    real2d_f density_int = transpose(data);   data = real2d_f();
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

      qv_fortran   = qv_host  .createDeviceCopy();
      qn_fortran   = qn_host  .createDeviceCopy();
      qp_fortran   = qp_host  .createDeviceCopy();
      temp_fortran = temp_host.createDeviceCopy();
    }

    {
      real2d_f data;
      real   dt;
      yakl::SimpleNetCDF nc;
      nc.open("sam1mom_data.nc");
      nc.read(data,"qv"         );    real2d_f qv          = transpose(data);   data = real2d_f();
      nc.read(data,"qn"         );    real2d_f qn          = transpose(data);   data = real2d_f();
      nc.read(data,"qp"         );    real2d_f qp          = transpose(data);   data = real2d_f();
      nc.read(data,"zint"       );    real2d_f zint        = transpose(data);   data = real2d_f();
      nc.read(data,"pressure"   );    real2d_f pressure    = transpose(data);   data = real2d_f();
      nc.read(data,"temp"       );    real2d_f temp        = transpose(data);   data = real2d_f();
      nc.read(data,"density"    );    real2d_f density     = transpose(data);   data = real2d_f();
      nc.read(data,"density_int");    real2d_f density_int = transpose(data);   data = real2d_f();
      nc.read(dt  ,"dt"         );

      int ncol = size(qv,1);
      int nz   = size(qv,2);

      sam1mom::Sam1mom micro;
      micro.main( dt , zint , density , density_int , pressure , temp , qv , qn , qp );

      qv_cxx   = qv  ;
      qn_cxx   = qn  ;
      qp_cxx   = qp  ;
      temp_cxx = temp;
    }

    int ncol = size(qv_cxx,1);
    int nz   = size(qv_cxx,2);

    real2d_f qv_adiff  ("qv_adiff  ",ncol,nz);
    real2d_f qn_adiff  ("qn_adiff  ",ncol,nz);
    real2d_f qp_adiff  ("qp_adiff  ",ncol,nz);
    real2d_f temp_adiff("temp_adiff",ncol,nz);
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int icol) {
      qv_adiff  (icol,k) = abs( qv_cxx  (icol,k) - qv_fortran  (icol,k) );
      qn_adiff  (icol,k) = abs( qn_cxx  (icol,k) - qn_fortran  (icol,k) );
      qp_adiff  (icol,k) = abs( qp_cxx  (icol,k) - qp_fortran  (icol,k) );
      temp_adiff(icol,k) = abs( temp_cxx(icol,k) - temp_fortran(icol,k) );
    });
    std::cout << "Relative diff qv  : " << std::scientific << std::setprecision(16) << yakl::intrinsics::sum(qv_adiff  ) / yakl::intrinsics::sum( qv_fortran   ) << "\n";
    std::cout << "Relative diff qn  : " << std::scientific << std::setprecision(16) << yakl::intrinsics::sum(qn_adiff  ) / yakl::intrinsics::sum( qn_fortran   ) << "\n";
    std::cout << "Relative diff qp  : " << std::scientific << std::setprecision(16) << yakl::intrinsics::sum(qp_adiff  ) / yakl::intrinsics::sum( qp_fortran   ) << "\n";
    std::cout << "Relative diff temp: " << std::scientific << std::setprecision(16) << yakl::intrinsics::sum(temp_adiff) / yakl::intrinsics::sum( temp_fortran ) << "\n";

  }
  yakl::finalize();
}



real2d_f transpose( real2d_f const &data ) {
  int dim1 = size(data,1);
  int dim2 = size(data,2);
  real2d_f ret(data.label(),dim2,dim1);
  yakl::fortran::parallel_for( yakl::fortran::Bounds<2>(dim1,dim2) , YAKL_LAMBDA (int i1, int i2) {
    ret(i2,i1) = data(i1,i2);
  });
  return ret;
}



