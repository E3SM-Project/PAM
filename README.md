# SCAMPAM: SCAlable, Modular, and Portable Atmospheric Model

PAM is a portable atmospheric model written in C++ with performance portability in mind using a kernel launching approach. It works on CPUs, Nvidia GPUs, and AMD GPUs. The goal for PAM is to allow easily exchangeable dynamical core and physics options through a simple and clear interface layer that's common to each. An emphasis is placed on algorithms that give good hardware uitilization on accelerated architectures and MPI patterns that give good scalability.

* [Directory Structure](#directory-structure)
* [Running standalone PAM](#running-standalone-pam)
* [Data types](#data-types)
* [`DataManager` class](#datamanager-class)
* [Driver, coupler state, and dimensions](#driver-coupler-state-and-dimensions)
* [`Dycore` class](#dycore-class)
* [`Microphysics` class](#microphysics-class)
* [Debugging utilities](#debugging-utilities)
* [`MultipleFields` class](#multiplefields-class)

## Directory Structure

* `coupler`: DataManager class and coupled state allocator function
* `dynamics`: Dynamical cores, currently "AWFL"
  * `awfl`: ADER WENO Finite-Volume (Collocated)
  * `spam`: Structure Preserving Atmospheric Model (Staggered)
* `externals`: git submodule dependencies
* `include`: PAM include files for the whole project
* `physics`: physics modules, currently micro only
  * `micro`: The various microphysics options for PAM
    * `none`: No microphysics
    * `kessler`: Kessler 3-species 1-moment microphysics
    * `p3`: P3 5-species, 2-moment microphysics
* `standalone`: standalone driver for PAM
* `utils`: Various PAM utilities, mostly in python

## Running standalone PAM

1. Clone the repo
2. `cd SCAMPAM && git submodule update --init --recursive`
3. `cd standalone/build`
4. Edit `cmakescript.sh` to choose your dycore and microphysics option
5. `source [machine].sh`
6. `./cmakescript`
7. Use `SCAMPAM/utils/generate_vertical_levels.py` to generate a vertical coordinates file (see the file's documentation for help)
8. Edit `../inputs/input_euler3d.yaml` to whatever parameters you want
9. `./driver ../inputs/input_euler3d.yaml`

## Data types

In `include/pam_const.h`, some convenient `typedef`s are created to hide the ugly templated `Array` syntax. For instance, a `real3d` is a 3-D array of type `real` that can only be accessed on the device in kernels. An `int2d` is a 2-D array of integers that can only be accessed on the device in kernels. A `realHost4d` is a 4-D array of real values that can only be accessed on the host outside kernels.

## `DataManager` class

The `DataManager` class exists as a convenient way store and manage persistent data that needs to exist for an entire CRM run (an entire "GCM physics time step"). To allocate data, one uses:
```C++
template <class T>
void register_and_allocate( std::string name ,
                            std::string desc ,
                            std::vector<int> dims ,
                            std::vector<std::string> dim_names ,
                            bool positive = false )
```

For convenience, you can pass an initializer list like `{nz,ny,nx,nens}` for the `dims` and `{"z","y","x","nens"}` for the `dim_names`. The DataManager class will automatically check to make sure you have not registered another variable by the same name or changed the length of a dimension name. It will also store the number of dimensions and the type `T` of the data. An example of registering a field is:
```C++
dm.register_and_allocate<real>( "water_vapor" , "Water Vapor Density" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
```

Only C-style array objects are used and supported in the `DataManager`.

All DataManager `get*` functions are YAKL Arrays that wrap internal DataManager pointers into contiguous memory. **Important**: DataManager calls only return YAKL `Array` objects that wrap pointers -- meaning there is virtually no cost to accessing already registered entries in the DataManager. It only copies metadata and then wraps an already existing pointer. 

To access this data later, you have three options available to you:
1. `template <class T, int N> Array<T,N,memDevice,styleC> get( std::string name )`
  * When you use `auto zint = dm.get<real,2>("vertical_midpoint_height")`, the DataManager will make sure that the entry with that name indeed has two dimensions and is of type `real`.
2. `template <class T> Array<T,2,memDevice,styleC> get_lev_col( std::string name )`
  * When you use `auto var = dm.get_lev_col<real>("density")`, the DataManager assumes the first dimension is levels, and it then collapses all other dimensions into a single column dimension; thus all returns from `get_lev_col` are 2-D YAKL `Array` objects.
3. `template <class T>  Array<T,1,memDevice,styleC> get_collapsed( std::string name )`
  * When you use `auto var = dm.get_collapsed<real>("density")`, the DataManager collapses all dimensions into a single dimensions and always returns a 1-D YAKL `Array` object.

There is also a `real get_dimension_size(std::string name)` function to return the size of the dimension of that name.

Finally, there is a `void validate_all()` function that will validate all entries in the DataManager to see if they are inf, NaN, or negative (only for entries registered as positive definite).

## Driver, coupler state, and dimensions

The driver will have entries into the main `DataManager` object called `vertical_midpoint_height` and `vertical_interface_height` initialized and the coupler dynamics state allocated before the dycore and micro `init(...)` functions are called.

The main dimensions in the model are `nz`, `ny`, `nx`, and `nens` for the z, y, and x dimensions and the ensemble dimension, respectively. The ensemble direction exists because typically many independent CRMs are simulated on a given MPI task.

The coupler dynamics state consists of the following:

```C++
dm.register_and_allocate<real>( "density_dry"  , "dry density"          , {nz,ny,nx,nens} , {"z","y","x","nens"} );
dm.register_and_allocate<real>( "uvel"         , "x-direction velocity" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
dm.register_and_allocate<real>( "vvel"         , "y-direction velocity" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
dm.register_and_allocate<real>( "wvel"         , "z-direction velocity" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
dm.register_and_allocate<real>( "temp"         , "temperature"          , {nz,ny,nx,nens} , {"z","y","x","nens"} );
dm.register_and_allocate<real>( "pressure_dry" , "dry pressure"         , {nz,ny,nx,nens} , {"z","y","x","nens"} );
```

To complete the coupler state, the microphysics and other parameterizations will register any tracers or other persistent data with the DataManager. 

The coupled state is assumed to be collocated on an A-grid on static height-based coordinates.

## `Dycore` class

The `Dycore` class needs to have the following defined to fit in the PAM interface. The typical calling sequence is:
* `init(...)`
* `add_tracer(...)` called multiple times by microphysics and potentially other classes
* `init_state_and_tracers(...)`
  * Inside the main time stepping loop, `compute_time_step()` and `timestep()` are called
* `finalize()`

```C++
class Dycore {
  // Allocate and initialize internal data that doesn't depend on the state
  // All other calls in this class assume this has previously been called
  void init(std::string inFile, int ny, int nx, int nens, real xlen, real ylen, int num_tracers, DataManager &dm);
  
  // Allocate and initialize internal data that depends on the state or micro class
  // The micro class is largely there to give the tracer index for water vapor
  template <class MICRO> void init_state_and_tracers( DataManager &dm , MICRO const &micro );
  
  // Compute the maximum stable time step given a CFL value
  // All necessary data is obtained from the DataManager object
  // Micro mainly exists to give the water vapor tracer index
  template <class MICRO> real compute_time_step(real cfl, DataManager &dm, MICRO const &micro);
  
  // Perform a dynamics time step:
  // (1) Convert coupler state to internal dynamics state
  // (2) Integrate space-time PDEs to compute and apply tendencies
  // (3) Convert internal dynamics state to coupler state
  template <class MICRO> void timeStep( DataManager &dm , MICRO const &micro , real dt );
  
  // Deallocate internal data structures and clean up
  void finalize(DataManager &dm);
  
  // Register a tracer in the dycore, and store whether it is positive definite
  // and whether it potentially contributes to the full mass (moisture loading)
  // The dycore doesn't actually have to moisture load; it can ignore the last bool
  // The return value is the index the dycore uses for the given tracer
  int add_tracer(DataManager &dm , std::string name , std::string desc , bool pos_def , bool mass_var);
};
```

## `Microphysics` class

The `Microphysics` class needs to have the following defined in order to fit in the PAM interface:

```C++
class Microphysics {
  // Allocate and initialize internal data structures
  // DC (Dycore) object exists mostly for the add_tracer(...) function
  // Register tracers and any other persistent variables in the DataManager
  // Need to store the water vapor tracer index (returned by add_tracer)
  // All other functions in this class assume this has been called first
  template <class DC> void init(std::string infile , int ny, int nx, int nens , DC &dycore , DataManager &dm);
  
  // Perform a microphysics time step
  // (1) Get coupler state, and transform it into whatever the micro needs
  // (2) Perform the microphysics step
  // (3) Transform results back to the coupler state
  void timeStep( DataManager &dm , real dt 
  
  // Return the tracer index for water vapor, originally returned by dycore.add_tracer()
  real get_water_vapor_index();
  
  // This struct needs to have at least the following defined
  struct Constants {
    real R_d    ;   
    real cp_d   ;   
    real cv_d   ;   
    real gamma_d;
    real kappa_d;
    real R_v    ;   
    real cp_v   ;   
    real cv_v   ;   
    real p0     ;   
  };  
  
  // Need to have an instance of the Constants struct named constants
  Constants constants;
};
```

It is assumed that at least one tracer (water vapor) will be defined by the microphysics.

## Debugging utilities

In `include/pam_const.h`, there are various validation routines to determine if an array contains inf, NaN, or negative values:
* `template <class T, int N, int MEM, int STYLE> void validate_array_nan( yakl::Array<T,N,MEM,STYLE> const &arr);`
* `template <class T, int N, int MEM, int STYLE> void validate_array_inf( yakl::Array<T,N,MEM,STYLE> const &arr);`
* `template <class T, int N, int MEM, int STYLE> void validate_array_inf_nan( yakl::Array<T,N,MEM,STYLE> const &arr);`
* `template <class T, int N, int MEM, int STYLE> void validate_array_positive( yakl::Array<T,N,MEM,STYLE> const &arr);`

These should be used in `Microphysics` and `Dycore` classes, particularly when the `PAM_DEBUG` macro variable is defined.

## `MultipleFields` class

There is a `MultipleFields` class that allows you to "glue together" multiple YAKL `Array` objects of the same dimensionality to make them appear to have an extra dimension. This is useful when performing the same function on multiple `Array` objects with potentially arbitrary number and name -- operations such as column-averaging, zeroing the memory, aggregating all tracers, etc. After adding fields, the object can then be accessed in kernels as if it were an array of one dimension higher. An example of using this is as follows:

```C++
MultipleFields<max_tracers,real4d> dm_tracers;
for (int tr = 0; tr < num_tracers; tr++) {
  auto trac = dm.get<real,4>( tracer_name[tr] );
  dm_tracers.add_field( trac );
}
parallel_for( Bounds<5>(num_tracers,nz,ny,nx,nens) , YAKL_LAMBDA (int l, int k, int j, int i, int iens) {
  dm_tracers(l,k,j,i,iens) = 0;
});
```

