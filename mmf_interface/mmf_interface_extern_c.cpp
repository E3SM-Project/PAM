
#include "mmf_interface.h"


extern "C" mmf_interface_init( int thread_id ) {
  mmf_interface::init(thread_id);
}


extern "C" mmf_interface_finalize() {
  mmf_interface::finalize();
}


extern "C" mmf_interface_allocate_coupler_state(int nz, int ny, int nx, int nens, int thread_id) {
  mmf_interface::allocate_coupler_state(nz,ny,nx,nens,thread_id);
}


extern "C" mmf_interface_set_phys_constants(double R_d, double R_v, double cp_d, double cp_v, double grav, double p0,
                                            int thread_id) {
  mmf_interface::mmf_interface_set_phys_constants(R_d, R_v, cp_d, cp_v, grav, p0, thread_id);
}


// Set grid with 2-D z-interface array dimensioned nz,nens with nens the fastest varying index
extern "C" mmf_interface_set_grid_zgrid2d(double xlen, double ylen, double *zint_p, int thread_id) {
  mmf_interface::check_thread_id(thread_id);
  int nz   = pam_interface_couplers[thread_id].get_nz  ();
  int nens = pam_interface_couplers[thread_id].get_nens();
  real2d zint("zint",zint_p,nz,nens);
  mmf_interface::set_grid(xlen,ylen,zint,thread_id);
}


// Set grid with 1-D z-interface array dimensioned nz
extern "C" mmf_interface_set_grid_zgrid1d(double xlen, double ylen, double *zint_p, int thread_id) {
  mmf_interface::check_thread_id(thread_id);
  int nens = pam_interface_couplers[thread_id].get_nens();
  real1d zint("zint",zint_p,nz);
  mmf_interface::set_grid(xlen,ylen,zint,thread_id);
}


// Set grid with 1-D z-interface array dimensioned nz
extern "C" mmf_interface_register_dimension(char const *name, int len, int thread_id) {
  mmf_interface::register_dimension(std::string(name) , len , thread_id);
}


// Allocate room and name the field to be retrieved later. Dimensions should be in C ordering
extern "C" void mmf_interface_register_and_allocate_array_double( char const * name , char const * desc , int ndims ,
                                                                  int *dims_p , int thread_id ) {
  std::vector<int> dims(ndims);   for (int i=0; i < ndims; i++) { dims[i] = dims_p[i]; }
  mmf_interface::register_and_allocate_array<double>( std::string(name) , std::string(desc) , dims , int thread_id );
}
extern "C" void mmf_interface_register_and_allocate_array_float ( char const * name , char const * desc , int ndims ,
                                                                  int *dims_p , int thread_id ) {
  std::vector<int> dims(ndims);   for (int i=0; i < ndims; i++) { dims[i] = dims_p[i]; }
  mmf_interface::register_and_allocate_array<float>( std::string(name) , std::string(desc) , dims , int thread_id );
}
extern "C" void mmf_interface_register_and_allocate_array_int   ( char const * name , char const * desc , int ndims ,
                                                                  int *dims_p , int thread_id ) {
  std::vector<int> dims(ndims);   for (int i=0; i < ndims; i++) { dims[i] = dims_p[i]; }
  mmf_interface::register_and_allocate_array<int>( std::string(name) , std::string(desc) , dims , int thread_id );
}
extern "C" void mmf_interface_register_and_allocate_array_bool  ( char const * name , char const * desc , int ndims ,
                                                                  int *dims_p , int thread_id ) {
  std::vector<int> dims(ndims);   for (int i=0; i < ndims; i++) { dims[i] = dims_p[i]; }
  mmf_interface::register_and_allocate_array<bool>( std::string(name) , std::string(desc) , dims , int thread_id );
}


// Allocate room and name the field to be retrieved later. Dimensions should be in C ordering
extern "C" void mmf_interface_unregister_and_deallocate_double( char const * name , int thread_id ) {
  mmf_interface::unregister_and_deallocate<double>( std::string(name) , thread_id );
}
extern "C" void mmf_interface_unregister_and_deallocate_float ( char const * name , int thread_id ) {
  mmf_interface::unregister_and_deallocate<float>( std::string(name) , thread_id );
}
extern "C" void mmf_interface_unregister_and_deallocate_int   ( char const * name , int thread_id ) {
  mmf_interface::unregister_and_deallocate<int>( std::string(name) , thread_id );
}
extern "C" void mmf_interface_unregister_and_deallocate_bool  ( char const * name , int thread_id ) {
  mmf_interface::unregister_and_deallocate<bool>( std::string(name) , thread_id );
}



