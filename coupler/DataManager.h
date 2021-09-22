
#pragma once

#include "pam_const.h"
#include <typeinfo>

using yakl::Array;



// Aggregate multiple fields into a single field that makes it easier to operate
// on them together inside the same kernel. used mostly for tracers
template <int MAX_FIELDS, class T>
class MultipleFields {
public:
  yakl::SArray<T,1,MAX_FIELDS> fields;
  int num_fields;

  YAKL_INLINE MultipleFields() { num_fields = 0; }

  YAKL_INLINE MultipleFields(MultipleFields const &rhs) {
    this->num_fields = rhs.num_fields;
    for (int i=0; i < num_fields; i++) {
      this->fields(i) = rhs.fields(i);
    }
  }

  YAKL_INLINE MultipleFields & operator=(MultipleFields const &rhs) {
    this->num_fields = rhs.num_fields;
    for (int i=0; i < num_fields; i++) {
      this->fields(i) = rhs.fields(i);
    }
    return *this;
  }

  YAKL_INLINE MultipleFields(MultipleFields &&rhs) {
    this->num_fields = rhs.num_fields;
    for (int i=0; i < num_fields; i++) {
      this->fields(i) = rhs.fields(i);
    }
  }

  YAKL_INLINE MultipleFields& operator=(MultipleFields &&rhs) {
    this->num_fields = rhs.num_fields;
    for (int i=0; i < num_fields; i++) {
      this->fields(i) = rhs.fields(i);
    }
    return *this;
  }

  YAKL_INLINE void add_field( T &field ) {
    this->fields(num_fields) = field;
    num_fields++;
  }

  YAKL_INLINE real &operator() (int tr, int k, int j, int i, int iens) const {
    return this->fields(tr)(k,j,i,iens);
  }
};




class DataManager {
public:

  struct Entry {
    std::string              name;
    std::string              desc;
    size_t                   type_hash;
    void *                   ptr;
    std::vector<int>         dims;
    std::vector<std::string> dim_names;
    bool                     positive;
  };

  struct Dimension {
    std::string name;
    int         len;
  };

  std::vector<Entry>       entries;
  std::vector<Dimension>   dimensions;
  bool check_data;
  bool die_on_failed_check;


  DataManager() {
    entries             = std::vector<Entry>();
    dimensions          = std::vector<Dimension>();
    check_data          = false;
    die_on_failed_check = false;
  }


  ~DataManager() {
    finalize();
  }


  template <class T>
  void register_and_allocate( std::string name ,
                              std::string desc ,
                              std::vector<int> dims ,
                              std::vector<std::string> dim_names ,
                              bool positive = false ) {
    // Make sure we don't have a duplicate entry
    if ( find_entry(name) != -1) {
      endrun("ERROR: Duplicate entry name");
    }
    if (dims.size() != dim_names.size()) {
      endrun("ERROR: Must have the same number of dims and dim_names");
    }

    // Make sure the dimensions are the same size as existing ones of the same name
    for (int i=0; i < dim_names.size(); i++) {
      int dimid = find_dimension(dim_names[i]);
      if (dimid == -1) {
        Dimension loc;
        loc.name = dim_names[i];
        loc.len  = dims     [i];
        dimensions.push_back(loc);
      } else {
        if (dimensions[dimid].len != dims[i]) {
          endrun("ERROR: Dimension already exists but has a different length");
        }
      }
    }


    Entry loc;
    loc.name      = name;
    loc.desc      = desc;
    loc.type_hash = get_type_hash<T>();
    loc.ptr       = yakl::yaklAllocDevice( get_data_size(dims)*sizeof(T) , name.c_str() );
    loc.dims      = dims;
    loc.dim_names = dim_names;
    loc.positive  = positive;

    entries.push_back( loc );
  }


  template <class T, int N>
  Array<T,N,memDevice,styleC> get( std::string name ) const {
    // Make sure we have this name as an entry
    int id = find_entry_or_error( name );
    // Make sure it's the right type and dimensionality
    validate_type<T>(id);
    validate_dims<N>(id);
    Array<T,N,memDevice,styleC> ret( name.c_str() , (T *) entries[id].ptr , entries[id].dims );
    return ret;
  }


  template <class T>
  Array<T,2,memDevice,styleC> get_lev_col( std::string name ) const {
    // Make sure we have this name as an entry
    int id = find_entry_or_error( name );
    // Make sure it's the right type
    validate_type<T>(id);
    validate_dims_lev_col(id);
    int nlev = entries[id].dims[0];
    int ncol = 1;
    for (int i=1; i < entries[id].dims.size(); i++) {
      ncol *= entries[id].dims[i];
    }
    Array<T,2,memDevice,styleC> ret( name.c_str() , (T *) entries[id].ptr , nlev , ncol );
    return ret;
  }


  template <class T>
  Array<T,1,memDevice,styleC> get_collapsed( std::string name ) const {
    // Make sure we have this name as an entry
    int id = find_entry_or_error( name );
    // Make sure it's the right type
    validate_type<T>(id);
    int ncells = entries[id].dims[0];
    for (int i=1; i < entries[id].dims.size(); i++) {
      ncells *= entries[id].dims[i];
    }
    Array<T,1,memDevice,styleC> ret( name.c_str() , (T *) entries[id].ptr , ncells );
    return ret;
  }


  void validate_all() const {
    // Check for NaNs
    for (int id = 0; id < entries.size(); id++) {
      if      (entry_type_is_same<short int>             (id)) { validate_single_nan<short int>             (id); }
      else if (entry_type_is_same<int>                   (id)) { validate_single_nan<int>                   (id); }
      else if (entry_type_is_same<long int>              (id)) { validate_single_nan<long int>              (id); }
      else if (entry_type_is_same<long long int>         (id)) { validate_single_nan<long long int>         (id); }
      else if (entry_type_is_same<unsigned short int>    (id)) { validate_single_nan<unsigned short int>    (id); }
      else if (entry_type_is_same<unsigned int>          (id)) { validate_single_nan<unsigned int>          (id); }
      else if (entry_type_is_same<unsigned long int>     (id)) { validate_single_nan<unsigned long int>     (id); }
      else if (entry_type_is_same<unsigned long long int>(id)) { validate_single_nan<unsigned long long int>(id); }
      else if (entry_type_is_same<float>                 (id)) { validate_single_nan<float>                 (id); }
      else if (entry_type_is_same<double>                (id)) { validate_single_nan<double>                (id); }
      else if (entry_type_is_same<long double>           (id)) { validate_single_nan<long double>           (id); }

      // Check for inf
      if      (entry_type_is_same<float>                 (id)) { validate_single_inf<float>                 (id); }
      else if (entry_type_is_same<double>                (id)) { validate_single_inf<double>                (id); }
      else if (entry_type_is_same<long double>           (id)) { validate_single_inf<long double>           (id); }

      // Check for negative values in positive-definite variables
      if      (entry_type_is_same<short int>             (id)) { validate_single_neg<short int>             (id); }
      else if (entry_type_is_same<int>                   (id)) { validate_single_neg<int>                   (id); }
      else if (entry_type_is_same<long int>              (id)) { validate_single_neg<long int>              (id); }
      else if (entry_type_is_same<long long int>         (id)) { validate_single_neg<long long int>         (id); }
      else if (entry_type_is_same<float>                 (id)) { validate_single_neg<float>                 (id); }
      else if (entry_type_is_same<double>                (id)) { validate_single_neg<double>                (id); }
      else if (entry_type_is_same<long double>           (id)) { validate_single_neg<long double>           (id); }
    }
  }


  template <class T>
  void validate_single_nan(int id) const {
    T *arr_dev = (T *) entries[id].ptr;
    size_t nelems = get_num_elems(id);
    T *arr = (T *) yakl::yaklAllocHost( nelems*sizeof(T) , "Nan check");
    yakl::memcpy_device_to_host(arr,arr_dev,nelems*sizeof(T));
    for (int i=0; i < nelems; i++) {
      if ( std::isnan( arr[i] ) ) {
        std::cerr << "WARNING: NaN discovered in: " << entries[id].name << " at global index: " << i << "\n";
        if (die_on_failed_check) {
          endrun("");
        }
      }
    }
  }


  template <class T>
  void validate_single_inf(int id) const {
    T *arr_dev = (T *) entries[id].ptr;
    size_t nelems = get_num_elems(id);
    T *arr = (T *) yakl::yaklAllocHost( nelems*sizeof(T) , "Inf check");
    yakl::memcpy_device_to_host(arr,arr_dev,nelems*sizeof(T));
    for (int i=0; i < nelems; i++) {
      if ( std::isinf( arr[i] ) ) {
        std::cerr << "WARNING: inf discovered in: " << entries[id].name << " at global index: " << i << "\n";
        if (die_on_failed_check) {
          endrun("");
        }
      }
    }
  }


  template <class T>
  void validate_single_neg(int id) const {
    if (entries[id].positive) {
    T *arr_dev = (T *) entries[id].ptr;
    size_t nelems = get_num_elems(id);
    T *arr = (T *) yakl::yaklAllocHost( nelems*sizeof(T) , "Neg check");
    yakl::memcpy_device_to_host(arr,arr_dev,nelems*sizeof(T));
      for (int i=0; i < nelems; i++) {
        if ( arr[i] < 0. ) {
          std::cerr << "WARNING: negative value discovered in positive-definite entry: " << entries[id].name
                    << " at global index: " << i << "\n";
          if (die_on_failed_check) {
            endrun("");
          }
        }
      }
    }
  }


  size_t get_num_elems(int id) const {
    size_t nelems = entries[id].dims[0];
    for (int i=1; i < entries[id].dims.size(); i++) { nelems *= entries[id].dims[i]; }
    return nelems;
  }


  int find_entry( std::string name ) const {
    for (int i=0; i < entries.size(); i++) {
      if (entries[i].name == name) return i;
    }
    return -1;
  }


  int find_dimension( std::string name ) const {
    for (int i=0; i < dimensions.size(); i++) {
      if (dimensions[i].name == name) return i;
    }
    return -1;
  }


  int find_entry_or_error( std::string name ) const {
    int id = find_entry( name );
    if (id >= 0) return id;
    endrun("ERROR: Could not find entry in coupler data");
    return -1;
  }


  int get_data_size( std::vector<int> dims ) const {
    int size = 1;
    for (int i=0; i < dims.size(); i++) { size *= dims[i]; }
    return size;
  }


  int get_dimension_size( std::string name ) const {
    int id = find_dimension( name );
    if (id == -1) { endrun("ERROR: Could not find dimension."); }
    return dimensions[id].len;
  }


  template <class T> size_t get_type_hash() const {
    return typeid(T).hash_code();
  }


  template <class T> size_t entry_type_is_same(int id) const {
    return entries[id].type_hash == typeid(T).hash_code();
  }


  template <class T>
  void validate_type(int id) const {
    if ( entries[id].type_hash != get_type_hash<T>() ) {
      endrun("ERROR: Requested Array type does not match entry type");
    }
  }


  template <int N>
  void validate_dims(int id) const {
    if ( N != entries[id].dims.size() ) {
      endrun("ERROR: Requested dimensions is different from the entry dimensions");
    }
  }


  void validate_dims_lev_col(int id) const {
    if ( entries[id].dims.size() < 2 ) {
      endrun("ERROR: Requested data is only one-dimensional");
    }
  }


  void finalize() {
    for (int i=0; i < entries.size(); i++) {
      yakl::yaklFreeDevice( entries[i].ptr , entries[i].name.c_str() );
    }
    entries              = std::vector<Entry>();
    dimensions           = std::vector<Dimension>();
  }


};

