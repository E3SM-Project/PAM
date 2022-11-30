
#pragma once

#include "pam_coupler.h"
#include "Options.h"
#include <string>
#include <vector>

// IMPORTANT: The pam_interface routines only deal with host-side data. These are the GCM-facing routines,
// and only some of them are callable from Fortran bindings. All routines with Fortran bindings have a comment
// that says "THIS HAS FORTRAN BINDINGS". To use this data on the GPU, it must be replicated to the GPU data manager
// from a C++-side interface routines

namespace pam_interface {

  // This is a vector because CPU threaded regions will require a different PamCoupler instance for each thread
  extern std::vector<pam::PamCoupler> couplers;


  // User provided function to perform initialization as called by the GCM
  // This is called once at the beginning of the simulation, not every time step
  // THIS HAS FORTRAN BINDINGS
  extern std::function<void()> gcm_initialize;


  // User provided function to perform a tendency calculation as called by the GCM
  // This is called every time step
  // THIS HAS FORTRAN BINDINGS
  extern std::function<void()> gcm_tendency;


  // User provided function to finalize as called by the GCM
  // This is called once at the beginning of the simulation, not every time step
  // THIS HAS FORTRAN BINDINGS
  extern std::function<void()> gcm_finalize;


  // This is intended to be called at the end of the simulation, not at the end of every GCM time step
  // THIS HAS FORTRAN BINDINGS
  inline void finalize() { couplers.clear(); }


  // Obtains the coupler for this thread ID: std::this_thread::get_id()
  inline pam::PamCoupler & get_coupler() {
    std::thread::id tid = std::this_thread::get_id();
    for (int i=0; i < couplers.size(); i++) { if (tid == couplers[i].get_thread_id()) return couplers[i]; }
    // If we got here, there isn't a coupler for this thread yet, so let's create one and return it
    couplers.emplace_back();
    return couplers.back();
  }


  // Register a dimension name in this thread's coupler
  // When you register and allocate an array, the dimensions are named for you. Therefore, this exists to predefine
  // a name for a dimension of a given size so that you have the dimensions names you want.
  // This is here because it's not easy to pass arrays of strings from Fortran to C++ and vice versa
  // THIS HAS FORTRAN BINDINGS
  inline void register_dimension(std::string name, int len) {
    get_coupler().get_data_manager_host_readwrite().add_dimension(name,len);
  }


  // Allocate an array, and store its metadata in the host-side data manager "dm_host" for this thread's coupler
  // THIS HAS FORTRAN BINDINGS
  template <class T>
  inline void register_and_allocate_array(std::string name, std::string desc, std::vector<int> dims) {
    get_coupler().get_data_manager_host_readwrite().register_and_allocate<T>( name , desc , dims );
  }


  // Allocate an array, and store its metadata in the host-side data manager "dm_host" for this thread's coupler
  // THIS HAS FORTRAN BINDINGS
  template <class T>
  inline void register_existing_array(std::string name, std::string desc, std::vector<int> dims, T * ptr) {
    get_coupler().get_data_manager_host_readwrite().register_existing<T>( name , desc , dims , ptr );
  }


  // Deallocate an array, and erase its metadata from the host-side data manager "dm_host" for this thread's coupler
  // THIS HAS FORTRAN BINDINGS
  inline void unregister_and_deallocate(std::string name) {
    get_coupler().get_data_manager_host_readwrite().unregister_and_deallocate(name);
  }


  // Deallocate an array, and erase its metadata from the host-side data manager "dm_host" for this thread's coupler
  // THIS HAS FORTRAN BINDINGS
  inline void unregister(std::string name) {
    get_coupler().get_data_manager_host_readwrite().unregister(name);
  }


  // Get an array from the host-side data manager "dm_host" for this thread's coupler
  // THIS HAS FORTRAN BINDINGS
  template <class T, int N>
  inline Array<T,N,memHost,styleC> get_array(std::string name) {
    return get_coupler().get_data_manager_host_readwrite().get<T,N>(name);
  }


  // Check if an array of this name has been registered and allocated in this thread's coupler
  // THIS HAS FORTRAN BINDINGS
  inline bool array_exists(std::string name) {
    return get_coupler().get_data_manager_host_readwrite().entry_exists(name);
  }


  // Set a {key,value} option pair in this thread's coupler
  // THIS HAS FORTRAN BINDINGS
  template <class T>
  inline void set_option(std::string name, T value) { get_coupler().set_option<T>(name,value); }


  // Get the option value for the given name from this thread's coupler
  // THIS HAS FORTRAN BINDINGS
  template <class T>
  inline T get_option(std::string name) { return get_coupler().get_option<T>(name); }


  // Check if the option of this name has been set in this thread's coupler
  // THIS HAS FORTRAN BINDINGS
  inline bool option_is_set(std::string name) { return get_coupler().option_exists(name); }


  // Remove this option from this thread's coupler
  // THIS HAS FORTRAN BINDINGS
  inline void remove_option(std::string name) { get_coupler().delete_option(name); }


  // Remove this option from this thread's coupler
  // THIS HAS FORTRAN BINDINGS
  inline void make_readonly(std::string name) { get_coupler().get_data_manager_host_readwrite().make_readonly(name); }
}


