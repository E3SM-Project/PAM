
#pragma once

#include "mpi.h"
#include "pam_coupler.h"
#include "YAKL_netcdf.h"

inline void output( pam::PamCoupler const &coupler , std::string out_prefix , real etime ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;

  int nranks;
  int myrank;
  MPI_Comm_size( MPI_COMM_WORLD , &nranks );
  MPI_Comm_rank( MPI_COMM_WORLD , &myrank );

  auto dx = coupler.get_dx();
  auto dy = coupler.get_dy();
  auto nx = coupler.get_nx();
  auto ny = coupler.get_ny();
  auto nz = coupler.get_nz();
  auto nens = coupler.get_nens();

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
        real1d xp1loc("xp1loc",nx+1);
        parallel_for( "Spatial.h output 1" , nx , YAKL_LAMBDA (int i) { xloc(i) = (i+0.5)*dx; });
        parallel_for( "Spatial.h output 2" , nx+1 , YAKL_LAMBDA (int i) { xp1loc(i) = (i)*dx; });
        nc.write(xloc.createHostCopy(),"x",{"x"});
        nc.write(xp1loc.createHostCopy(),"xp1",{"xp1"});

        // y-coordinate
        real1d yloc("yloc",ny);
        real1d yp1loc("yp1loc",ny+1);
        parallel_for( "Spatial.h output 3" , ny , YAKL_LAMBDA (int i) { yloc(i) = (i+0.5)*dy; });
        parallel_for( "Spatial.h output 4" , ny+1 , YAKL_LAMBDA (int i) { yp1loc(i) = (i)*dy; });
        nc.write(yloc.createHostCopy(),"y",{"y"});
        nc.write(yp1loc.createHostCopy(),"yp1",{"yp1"});

        // z-coordinate
        auto zint = dm.get<real const,2>("vertical_interface_height");
        auto dz = dm.get<real const,2>("vertical_cell_dz");
        auto zmid = dm.get<real const,2>("vertical_midpoint_height");
        nc.write(zmid.createHostCopy(),"z", {"z","nens"});
        nc.write(dz.createHostCopy(),"dz", {"z","nens"});
        nc.write(zint.createHostCopy(),"zp1", {"zp1","nens"});

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
      real4d data("data",nz,ny,nx,nens);
      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int n) { data(k,j,i,n) = fields(0,k,j,i,n); });
      nc.write1(data,"density"    ,{"z","y","x","nens"},ulIndex,"t");
      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int n) { data(k,j,i,n) = fields(1,k,j,i,n); });
      nc.write1(data,"uvel"       ,{"z","y","x","nens"},ulIndex,"t");
      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int n) { data(k,j,i,n) = fields(2,k,j,i,n); });
      nc.write1(data,"vvel"       ,{"z","y","x","nens"},ulIndex,"t");
      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int n) { data(k,j,i,n) = fields(3,k,j,i,n); });
      nc.write1(data,"wvel"       ,{"z","y","x","nens"},ulIndex,"t");
      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int n) { data(k,j,i,n) = fields(4,k,j,i,n); });
      nc.write1(data,"temperature",{"z","y","x","nens"},ulIndex,"t");
      for (int tr=0; tr < num_tracers; tr++) {
        parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int n) {
          // data(k,j,i) = fields(5+tr,k,j,i,0) / (state(idR,hs+k,hs+j,hs+i,0) + hyDensCells(k,0));
          data(k,j,i,n) = fields(5+tr,k,j,i,n);
        });
        nc.write1(data,tracer_names[tr],{"z","y","x","nens"},ulIndex,"t");
      }

    auto uvel_stag = dm.get<real const, 4>("uvel_stag");
    auto vvel_stag = dm.get<real const, 4>("vvel_stag");
    auto wvel_stag = dm.get<real const, 4>("wvel_stag");

//WRITE OUT OTHER VARIABLES!

      // Close the file
      nc.close();



    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
