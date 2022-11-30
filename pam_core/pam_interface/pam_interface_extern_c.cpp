
#include "pam_interface.h"


extern "C" void pam_interface_finalize() {
  pam_interface::finalize();
}


extern "C" void pam_interface_mirror_array_readonly_bool(char const *name, char const *desc, int *dims, int ndims, bool *ptr) {
  auto reg = [] (std::string name , std::string desc , std::vector<int> dims , bool *ptr ) {
    pam_interface::register_existing_array<bool const>(name,desc,dims,const_cast<bool const *>(ptr));
  };
  if (ndims == 1) reg(name,desc,{dims[0]},ptr);
  if (ndims == 2) reg(name,desc,{dims[0],dims[1]},ptr);
  if (ndims == 3) reg(name,desc,{dims[0],dims[1],dims[2]},ptr);
  if (ndims == 4) reg(name,desc,{dims[0],dims[1],dims[2],dims[3]},ptr);
  if (ndims == 5) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4]},ptr);
  if (ndims == 6) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]},ptr);
  if (ndims == 7) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5],dims[6]},ptr);
}
extern "C" void pam_interface_mirror_array_readonly_int(char const *name, char const *desc, int *dims, int ndims, int *ptr) {
  auto reg = [] (std::string name , std::string desc , std::vector<int> dims , int *ptr ) {
    pam_interface::register_existing_array<int const>(name,desc,dims,const_cast<int const *>(ptr));
  };
  if (ndims == 1) reg(name,desc,{dims[0]},ptr);
  if (ndims == 2) reg(name,desc,{dims[0],dims[1]},ptr);
  if (ndims == 3) reg(name,desc,{dims[0],dims[1],dims[2]},ptr);
  if (ndims == 4) reg(name,desc,{dims[0],dims[1],dims[2],dims[3]},ptr);
  if (ndims == 5) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4]},ptr);
  if (ndims == 6) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]},ptr);
  if (ndims == 7) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5],dims[6]},ptr);
}
extern "C" void pam_interface_mirror_array_readonly_float(char const *name, char const *desc, int *dims, int ndims, float *ptr) {
  auto reg = [] (std::string name , std::string desc , std::vector<int> dims , float *ptr ) {
    pam_interface::register_existing_array<float const>(name,desc,dims,const_cast<float const *>(ptr));
  };
  if (ndims == 1) reg(name,desc,{dims[0]},ptr);
  if (ndims == 2) reg(name,desc,{dims[0],dims[1]},ptr);
  if (ndims == 3) reg(name,desc,{dims[0],dims[1],dims[2]},ptr);
  if (ndims == 4) reg(name,desc,{dims[0],dims[1],dims[2],dims[3]},ptr);
  if (ndims == 5) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4]},ptr);
  if (ndims == 6) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]},ptr);
  if (ndims == 7) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5],dims[6]},ptr);
}
extern "C" void pam_interface_mirror_array_readonly_double(char const *name, char const *desc, int *dims, int ndims, double *ptr) {
  auto reg = [] (std::string name , std::string desc , std::vector<int> dims , double *ptr ) {
    pam_interface::register_existing_array<double const>(name,desc,dims,const_cast<double const *>(ptr));
  };
  if (ndims == 1) reg(name,desc,{dims[0]},ptr);
  if (ndims == 2) reg(name,desc,{dims[0],dims[1]},ptr);
  if (ndims == 3) reg(name,desc,{dims[0],dims[1],dims[2]},ptr);
  if (ndims == 4) reg(name,desc,{dims[0],dims[1],dims[2],dims[3]},ptr);
  if (ndims == 5) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4]},ptr);
  if (ndims == 6) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]},ptr);
  if (ndims == 7) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5],dims[6]},ptr);
}


extern "C" void pam_interface_mirror_array_readwrite_bool(char const *name, char const *desc, int *dims, int ndims, bool *ptr) {
  auto reg = [] (std::string name , std::string desc , std::vector<int> dims , bool *ptr ) {
    pam_interface::register_existing_array<bool>(name,desc,dims,ptr);
  };
  if (ndims == 1) reg(name,desc,{dims[0]},ptr);
  if (ndims == 2) reg(name,desc,{dims[0],dims[1]},ptr);
  if (ndims == 3) reg(name,desc,{dims[0],dims[1],dims[2]},ptr);
  if (ndims == 4) reg(name,desc,{dims[0],dims[1],dims[2],dims[3]},ptr);
  if (ndims == 5) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4]},ptr);
  if (ndims == 6) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]},ptr);
  if (ndims == 7) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5],dims[6]},ptr);
}
extern "C" void pam_interface_mirror_array_readwrite_int(char const *name, char const *desc, int *dims, int ndims, int *ptr) {
  auto reg = [] (std::string name , std::string desc , std::vector<int> dims , int *ptr ) {
    pam_interface::register_existing_array<int>(name,desc,dims,ptr);
  };
  if (ndims == 1) reg(name,desc,{dims[0]},ptr);
  if (ndims == 2) reg(name,desc,{dims[0],dims[1]},ptr);
  if (ndims == 3) reg(name,desc,{dims[0],dims[1],dims[2]},ptr);
  if (ndims == 4) reg(name,desc,{dims[0],dims[1],dims[2],dims[3]},ptr);
  if (ndims == 5) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4]},ptr);
  if (ndims == 6) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]},ptr);
  if (ndims == 7) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5],dims[6]},ptr);
}
extern "C" void pam_interface_mirror_array_readwrite_float(char const *name, char const *desc, int *dims, int ndims, float *ptr) {
  auto reg = [] (std::string name , std::string desc , std::vector<int> dims , float *ptr ) {
    pam_interface::register_existing_array<float>(name,desc,dims,ptr);
  };
  if (ndims == 1) reg(name,desc,{dims[0]},ptr);
  if (ndims == 2) reg(name,desc,{dims[0],dims[1]},ptr);
  if (ndims == 3) reg(name,desc,{dims[0],dims[1],dims[2]},ptr);
  if (ndims == 4) reg(name,desc,{dims[0],dims[1],dims[2],dims[3]},ptr);
  if (ndims == 5) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4]},ptr);
  if (ndims == 6) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]},ptr);
  if (ndims == 7) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5],dims[6]},ptr);
}
extern "C" void pam_interface_mirror_array_readwrite_double(char const *name, char const *desc, int *dims, int ndims, double *ptr) {
  auto reg = [] (std::string name , std::string desc , std::vector<int> dims , double *ptr ) {
    pam_interface::register_existing_array<double>(name,desc,dims,ptr);
  };
  if (ndims == 1) reg(name,desc,{dims[0]},ptr);
  if (ndims == 2) reg(name,desc,{dims[0],dims[1]},ptr);
  if (ndims == 3) reg(name,desc,{dims[0],dims[1],dims[2]},ptr);
  if (ndims == 4) reg(name,desc,{dims[0],dims[1],dims[2],dims[3]},ptr);
  if (ndims == 5) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4]},ptr);
  if (ndims == 6) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]},ptr);
  if (ndims == 7) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5],dims[6]},ptr);
}



extern "C" void pam_interface_set_option_bool(char const *name , bool value ) {
  pam_interface::set_option( name , value );
}
extern "C" void pam_interface_set_option_int(char const *name , int value ) {
  pam_interface::set_option( name , value );
}
extern "C" void pam_interface_set_option_string(char const *name , char const *value ) {
  pam_interface::set_option( name , std::string(value) );
}
extern "C" void pam_interface_set_option_float(char const *name , float value ) {
  pam_interface::set_option( name , value );
}
extern "C" void pam_interface_set_option_double(char const *name , double value ) {
  pam_interface::set_option( name , value );
}


extern "C" void pam_interface_get_option_bool(char const *name , bool &value ) {
  value = pam_interface::get_option<bool>( name );
}
extern "C" void pam_interface_get_option_int(char const *name , int &value ) {
  value = pam_interface::get_option<int>( name );
}
extern "C" void pam_interface_get_option_stringlen(char const *name , int &len ) {
  std::string value_loc = pam_interface::get_option<std::string>( name );
  len = value_loc.size();
}
extern "C" void pam_interface_get_option_string(char const *name , char *value ) {
  std::string value_loc = pam_interface::get_option<std::string>( name );
  for (int i=0; i < value_loc.size(); i++) { value[i] = value_loc[i]; }
}
extern "C" void pam_interface_get_option_float(char const *name , float &value ) {
  value = pam_interface::get_option<float>( name );
}
extern "C" void pam_interface_get_option_double(char const *name , double &value ) {
  value = pam_interface::get_option<double>( name );
}


extern "C" void pam_interface_option_exists(char const *name , bool &exists ) {
  exists = pam_interface::option_is_set( name );
}


extern "C" void pam_interface_remove_option(char const *name ) {
  pam_interface::remove_option( name );
}


extern "C" void pam_interface_make_readonly(char const *name) {
  pam_interface::make_readonly( name );
}


extern "C" void pam_interface_create_array_bool(char const *name, char const *desc, int *dims, int ndims) {
  auto reg = [] (std::string name , std::string desc , std::vector<int> dims ) {
    pam_interface::register_and_allocate_array<bool>(name,desc,dims);
  };
  if (ndims == 1) reg(name,desc,{dims[0]});
  if (ndims == 2) reg(name,desc,{dims[0],dims[1]});
  if (ndims == 3) reg(name,desc,{dims[0],dims[1],dims[2]});
  if (ndims == 4) reg(name,desc,{dims[0],dims[1],dims[2],dims[3]});
  if (ndims == 5) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4]});
  if (ndims == 6) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]});
  if (ndims == 7) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5],dims[6]});
}
extern "C" void pam_interface_create_array_int(char const *name, char const *desc, int *dims, int ndims) {
  auto reg = [] (std::string name , std::string desc , std::vector<int> dims ) {
    pam_interface::register_and_allocate_array<int>(name,desc,dims);
  };
  if (ndims == 1) reg(name,desc,{dims[0]});
  if (ndims == 2) reg(name,desc,{dims[0],dims[1]});
  if (ndims == 3) reg(name,desc,{dims[0],dims[1],dims[2]});
  if (ndims == 4) reg(name,desc,{dims[0],dims[1],dims[2],dims[3]});
  if (ndims == 5) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4]});
  if (ndims == 6) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]});
  if (ndims == 7) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5],dims[6]});
}
extern "C" void pam_interface_create_array_float(char const *name, char const *desc, int *dims, int ndims) {
  auto reg = [] (std::string name , std::string desc , std::vector<int> dims ) {
    pam_interface::register_and_allocate_array<float>(name,desc,dims);
  };
  if (ndims == 1) reg(name,desc,{dims[0]});
  if (ndims == 2) reg(name,desc,{dims[0],dims[1]});
  if (ndims == 3) reg(name,desc,{dims[0],dims[1],dims[2]});
  if (ndims == 4) reg(name,desc,{dims[0],dims[1],dims[2],dims[3]});
  if (ndims == 5) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4]});
  if (ndims == 6) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]});
  if (ndims == 7) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5],dims[6]});
}
extern "C" void pam_interface_create_array_double(char const *name, char const *desc, int *dims, int ndims) {
  auto reg = [] (std::string name , std::string desc , std::vector<int> dims ) {
    pam_interface::register_and_allocate_array<double>(name,desc,dims);
  };
  if (ndims == 1) reg(name,desc,{dims[0]});
  if (ndims == 2) reg(name,desc,{dims[0],dims[1]});
  if (ndims == 3) reg(name,desc,{dims[0],dims[1],dims[2]});
  if (ndims == 4) reg(name,desc,{dims[0],dims[1],dims[2],dims[3]});
  if (ndims == 5) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4]});
  if (ndims == 6) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]});
  if (ndims == 7) reg(name,desc,{dims[0],dims[1],dims[2],dims[3],dims[4],dims[5],dims[6]});
}


extern "C" void pam_interface_destroy_array(char const *name) {
  pam_interface::unregister_and_deallocate(name);
}


// TODO: This doesn't work becuase you need a double pointer or something to actually change the pointer address
extern "C" void pam_interface_get_array_bool(char const *name, bool * &ptr, int *dims, int ndims) {
  using pam_interface::get_array;
  if (ndims==1) { auto var=get_array<bool,1>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==2) { auto var=get_array<bool,2>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==3) { auto var=get_array<bool,3>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==4) { auto var=get_array<bool,4>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==5) { auto var=get_array<bool,5>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==6) { auto var=get_array<bool,6>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==7) { auto var=get_array<bool,7>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
}
extern "C" void pam_interface_get_array_int(char const *name, int * &ptr, int *dims, int ndims) {
  using pam_interface::get_array;
  if (ndims==1) { auto var=get_array<int,1>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==2) { auto var=get_array<int,2>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==3) { auto var=get_array<int,3>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==4) { auto var=get_array<int,4>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==5) { auto var=get_array<int,5>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==6) { auto var=get_array<int,6>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==7) { auto var=get_array<int,7>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
}
extern "C" void pam_interface_get_array_float(char const *name, float * &ptr, int *dims, int ndims) {
  using pam_interface::get_array;
  if (ndims==1) { auto var=get_array<float,1>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==2) { auto var=get_array<float,2>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==3) { auto var=get_array<float,3>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==4) { auto var=get_array<float,4>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==5) { auto var=get_array<float,5>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==6) { auto var=get_array<float,6>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==7) { auto var=get_array<float,7>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
}
extern "C" void pam_interface_get_array_double(char const *name, double * &ptr, int *dims, int ndims) {
  using pam_interface::get_array;
  if (ndims==1) { auto var=get_array<double,1>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==2) { auto var=get_array<double,2>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==3) { auto var=get_array<double,3>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==4) { auto var=get_array<double,4>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==5) { auto var=get_array<double,5>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==6) { auto var=get_array<double,6>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
  if (ndims==7) { auto var=get_array<double,7>(name); for (int i=0; i<ndims; i++) { dims[i]=var.dimension[i]; }; ptr=var.data(); }
}


extern "C" void pam_interface_array_exists(char const *name, bool &exists) {
  exists = pam_interface::array_exists(name);
}


extern "C" void pam_interface_register_dimension(char const *name, int len) {
  pam_interface::register_dimension(name,len);
}




