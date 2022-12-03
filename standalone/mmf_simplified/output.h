
#pragma once

#include "pam_coupler.h"
#include "YAKL_netcdf.h"

inline void output( pam::PamCoupler const &coupler , std::string out_prefix , real etime ) {
  int nranks;
  int myrank;
  MPI_Comm_size( MPI_COMM_WORLD , &nranks );
  MPI_Comm_rank( MPI_COMM_WORLD , &myrank );

  auto dx = coupler.get_dx();
  auto dy = coupler.get_dy();
  auto nx = coupler.get_nx();
  auto ny = coupler.get_ny();
  auto nz = coupler.get_nz();

  auto &dm = coupler.get_data_manager_device_readonly();

  MPI_Barrier(MPI_COMM_WORLD);
  for (int rr=0; rr < nranks; rr++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (rr == myrank) {



      std::string fname = out_prefix + std::string("_") + std::to_string(myrank) + std::string(".nc");

      yakl::SimpleNetCDF nc;
      int ulIndex = 0; // Unlimited dimension index to place this data at
      // Create or open the file
      if (etime == 0.) {
        nc.create(fname);

        // x-coordinate
        real1d xloc("xloc",nx);
        parallel_for( "Spatial.h output 1" , nx , YAKL_LAMBDA (int i) { xloc(i) = (i+0.5)*dx; });
        nc.write(xloc.createHostCopy(),"x",{"x"});

        // y-coordinate
        real1d yloc("yloc",ny);
        parallel_for( "Spatial.h output 2" , ny , YAKL_LAMBDA (int i) { yloc(i) = (i+0.5)*dy; });
        nc.write(yloc.createHostCopy(),"y",{"y"});

        // z-coordinate
        auto zint = dm.get<real const,2>("vertical_interface_height");
        real1d zmid("zmid",nz);
        parallel_for( "Spatial.h output 3" , nz , YAKL_LAMBDA (int i) {
          zmid(i) = ( zint(i,0) + zint(i+1,0) ) / 2;
        });
        nc.write(zmid.createHostCopy(),"z",{"z"});

        // Create time variable
        nc.write1(0._fp,"t",0,"t");
      } else {
        nc.open(fname,yakl::NETCDF_MODE_WRITE);
        ulIndex = nc.getDimSize("t");

        // Write the elapsed time
        nc.write1(etime,"t",ulIndex,"t");
      }

      std::vector<std::string> tracer_names = coupler.get_tracer_names();
      int num_tracers = coupler.get_num_tracers();
      // Create MultiField of all state and tracer full variables, since we're doing the same operation on each
      pam::MultiField<real const,4> fields;
      fields.add_field( dm.get<real const,4>("density_dry") );
      fields.add_field( dm.get<real const,4>("uvel"       ) );
      fields.add_field( dm.get<real const,4>("vvel"       ) );
      fields.add_field( dm.get<real const,4>("wvel"       ) );
      fields.add_field( dm.get<real const,4>("temp"       ) );
      for (int tr=0; tr < num_tracers; tr++) {
        fields.add_field( dm.get<real const,4>(tracer_names[tr]) );
      }

      // First, write out standard coupler state
      real3d data("data",nz,ny,nx);
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = fields(0,k,j,i,0); });
      nc.write1(data,"density"    ,{"z","y","x"},ulIndex,"t");
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = fields(1,k,j,i,0); });
      nc.write1(data,"uvel"       ,{"z","y","x"},ulIndex,"t");
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = fields(2,k,j,i,0); });
      nc.write1(data,"vvel"       ,{"z","y","x"},ulIndex,"t");
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = fields(3,k,j,i,0); });
      nc.write1(data,"wvel"       ,{"z","y","x"},ulIndex,"t");
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = fields(4,k,j,i,0); });
      nc.write1(data,"temperature",{"z","y","x"},ulIndex,"t");
      for (int tr=0; tr < num_tracers; tr++) {
        parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
          // data(k,j,i) = fields(5+tr,k,j,i,0) / (state(idR,hs+k,hs+j,hs+i,0) + hyDensCells(k,0));
          data(k,j,i) = fields(5+tr,k,j,i,0);
        });
        nc.write1(data,tracer_names[tr],{"z","y","x"},ulIndex,"t");
      }

      // Close the file
      nc.close();



    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}


