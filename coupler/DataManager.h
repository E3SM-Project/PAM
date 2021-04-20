
#pragma once

#include "pam_const.h"
#include <typeinfo>

using yakl::Array;



// Aggregate multiple fields into a single field that makes it easier to operate
// on them together inside the same kernel. used mosetly for tracers
template <int MAX_FIELDS, class T>
class MultipleFields {
public:
  yakl::SArray<T,1,MAX_FIELDS> fields;
  int num_fields;

  MultipleFields() { num_fields = 0; }

  MultipleFields(MultipleFields const &rhs) {
    this->num_fields = rhs.num_fields;
    for (int i=0; i < num_fields; i++) {
      this->fields(i) = rhs.fields(i);
    }
  }

  MultipleFields & operator=(MultipleFields const &rhs) {
    this->num_fields = rhs.num_fields;
    for (int i=0; i < num_fields; i++) {
      this->fields(i) = rhs.fields(i);
    }
    return *this;
  }

  MultipleFields(MultipleFields &&rhs) {
    this->num_fields = rhs.num_fields;
    for (int i=0; i < num_fields; i++) {
      this->fields(i) = rhs.fields(i);
    }
  }

  MultipleFields& operator=(MultipleFields &&rhs) {
    this->num_fields = rhs.num_fields;
    for (int i=0; i < num_fields; i++) {
      this->fields(i) = rhs.fields(i);
    }
    return *this;
  }

  void add_tracer( T &tracer ) {
    this->fields(num_fields) = tracer;
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
    std::string              type;
    void *                   ptr;
    std::vector<int>         dims;
    std::vector<std::string> dim_names;
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
                              std::vector<std::string> dim_names ) {
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
    loc.type      = get_type_string<T>();
    loc.ptr       = yakl::yaklAllocDevice( get_data_size(dims)*sizeof(T) , name.c_str() );
    loc.dims      = dims;
    loc.dim_names = dim_names;

    entries.push_back( loc );
  }


  template <class T, int N>
  Array<T,N,memDevice,styleC> get( std::string name ) {
    // Make sure we have this name as an entry
    int id = find_entry_or_error( name );
    // Make sure it's the right type
    if ( entries[id].type != get_type_string<T>() ) {
      endrun("ERROR: Requested Array type does not match entry type");
    }
    // Make sure the dimensionality is correct
    if ( N != entries[id].dims.size() ) {
      endrun("ERROR: Requested dimensions is different from the entry dimensions");
    }
    Array<T,N,memDevice,styleC> ret( name.c_str() , (T *) entries[id].ptr , entries[id].dims );
    if (check_data) { validate_array( ret , name ); }
    return ret;
  }


  template <class T>
  Array<T,2,memDevice,styleC> get_lev_col( std::string name ) {
    // Make sure we have this name as an entry
    int id = find_entry_or_error( name );
    // Make sure it's the right type
    if ( entries[id].type != get_type_string<T>() ) {
      endrun("ERROR: Requested Array type does not match entry type");
    }
    // Make sure the dimensionality is correct
    if ( entries[id].dims.size() < 2 ) {
      endrun("ERROR: Requested data is only one-dimensional");
    }
    int nlev = entries[id].dims[0];
    int ncol = 1;
    for (int i=1; i < entries[id].dims.size(); i++) {
      ncol *= entries[id].dims[i];
    }
    Array<T,2,memDevice,styleC> ret( name.c_str() , (T *) entries[id].ptr , nlev , ncol );
    if (check_data) { validate_array( ret , name ); }
    return ret;
  }


  template <class T>
  Array<T,1,memDevice,styleC> get_collapsed( std::string name ) {
    // Make sure we have this name as an entry
    int id = find_entry_or_error( name );
    // Make sure it's the right type
    if ( entries[id].type != get_type_string<T>() ) {
      endrun("ERROR: Requested Array type does not match entry type");
    }
    int ncells = entries[id].dims[0];
    for (int i=1; i < entries[id].dims.size(); i++) {
      ncells *= entries[id].dims[i];
    }
    Array<T,1,memDevice,styleC> ret( name.c_str() , (T *) entries[id].ptr , ncells );
    if (check_data) { validate_array( ret , name ); }
    return ret;
  }


  template <class T, int N>
  void validate_array( Array<T,N,memDevice,styleC> const &arr , std::string name ) const {
    auto arrHost = arr.createHostCopy();
    for (unsigned i=0; i < arr.get_elem_count(); i++) {
      T val = arr.myData[i];
      if ( std::isnan(val) ) {
        std::cerr << "ERROR: NaN discovered in: " << name << " at global index: " << i << "\n";
        if (die_on_failed_check) {
          endrun("");
        }
      }
      if ( std::isinf(val) ) {
        std::cerr << "ERROR: inf discovered in: " << name << " at global index: " << i << "\n";
        if (die_on_failed_check) {
          endrun("");
        }
      }
    }
  }


  int find_entry( std::string name ) const {
    for (int i=0; i < entries.size(); i++) {
      if (entries[i].name == name) return i;
    }
    return -1;
  }


  int find_dimension( std::string name ) const {
    for (int i=0; i < entries.size(); i++) {
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


  template <class T> std::string get_type_string() const {
    return std::string( typeid(T).name() );
  }


  void finalize() {
    for (int i=0; i < entries.size(); i++) {
      yakl::yaklFreeDevice( entries[i].ptr , entries[i].name.c_str() );
    }
    entries              = std::vector<Entry>();
    dimensions           = std::vector<Dimension>();
  }


};

