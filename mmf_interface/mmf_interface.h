
#pragma once

#include "pam_coupler.h"
#include "DataManager.h"
#include "Options.h"
#include <string>
#include <vector>

// This is a vector because CPU threaded regions will require a different PamCoupler instance for each thread
extern std::vector<PamCoupler> pam_interface_couplers;

namespace mmf_interface {

  inline void check_thread_id(int thread_id) {
    if (thread_id >= pam_interface_couplers.size()) {
      endrun("ERROR: thread_id is larger than pam_interface_couplers.size()");
    }
    if (thread_id < 0) {
      endrun("ERROR: thread_id is  less than zero");
    }
  }


  // This is intended to be called at the beginning of the simulation, not at the beginning of every GCM time step
  inline void init(int nthreads=1) {
    if (nthreads < 0) endrun("ERROR: nthreads is less than zero");
    pam_interface_couplers = std::vector<PamCoupler>(nthreads);
  }


  // This is intended to be called at the end of the simulation, not at the end of every GCM time step
  inline void finalize() {
    pam_interface_couplers = std::vector<PamCoupler>();
  }


  // This is intended to be called at the beginning of the simulation, not at the beginning of every GCM time step
  inline void allocate_coupler_state(int nz, int ny, int nx, int nens, int thread_id=0) {
    check_thread_id(thread_id);
    pam_interface_couplers[thread_id].allocate_coupler_state( nz , ny , nx , nens );
  }


  inline PamCoupler & get_coupler(int thread_id=0) {
    check_thread_id(thread_id);
    return pam_interface_couplers[thread_id];
  }


  // This is intended to be called at the beginning of the simulation, not at the beginning of every GCM time step
  inline void set_phys_constants(real R_d, real R_v, real cp_d, real cp_v, real grav, real p0, int thread_id=0) {
    check_thread_id(thread_id);
    pam_interface_couplers[thread_id].set_phys_constants( R_d ,  R_v ,  cp_d ,  cp_v ,  grav ,  p0 );
  }


  inline void set_grid(real xlen, real ylen, realConst2d zint_in, int thread_id=0) {
    check_thread_id(thread_id);
    pam_interface_couplers[thread_id].set_grid(xlen, ylen, zint_in);
  }


  inline void set_grid(real xlen, real ylen, realConst1d zint_in, int thread_id=0) {
    check_thread_id(thread_id);
    pam_interface_couplers[thread_id].set_grid(xlen, ylen, zint_in);
  }


  inline void register_dimension( std::string name , int len , int thread_id=0 ) {
    check_thread_id(thread_id);
    pam_interface_couplers[thread_id].dm.add_dimension( name , len );
  }


  // Allocate room and name the field to be retrieved later
  template <class T>
  inline void register_and_allocate_array( std::string name , std::string desc , std::vector<int> dims , int thread_id=0 ) {
    check_thread_id(thread_id);
    auto &dm = pam_interface_couplers[thread_id].dm;
    dm.register_and_allocate<T>( name , desc , dims );
  }


  // Allocate room and name the field to be retrieved later
  template <class T>
  inline void unregister_and_deallocate( std::string name , int thread_id=0 ) {
    check_thread_id(thread_id);
    auto &dm = pam_interface_couplers[thread_id].dm;
    dm.unregister_and_deallocate<T>( name );
  }


  // Retrieve the field for writing or reading
  template <class T, int N>
  inline Array<T,N,memDevice,styleC> get_array( std::string name , int thread_id=0 ) {
    check_thread_id(thread_id);
    auto &dm = pam_interface_couplers[thread_id].dm;
    return dm.get<T,N>( name );
  }


  // Set a configuration option
  template <class T>
  inline void set_option( std::string name , T value , int thread_id=0 ) {
    check_thread_id(thread_id);
    pam_interface_couplers[thread_id].set_option<T>(name,value);
  }


  // Set a configuration option
  template <class T>
  inline T get_option( std::string name , int thread_id=0 ) {
    check_thread_id(thread_id);
    return pam_interface_couplers[thread_id].get_option<T>(name);
  }


  // Set a configuration option
  inline bool option_is_set( std::string name , int thread_id=0 ) {
    check_thread_id(thread_id);
    return pam_interface_couplers[thread_id].option_exists(name);
  }


  // Set a configuration option
  inline void remove_option( std::string name , int thread_id=0 ) {
    check_thread_id(thread_id);
    pam_interface_couplers[thread_id].delete_option(name);
  }

}


