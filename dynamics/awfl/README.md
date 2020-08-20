# PFAST (Performance-portable Framework for Arbitrary Space-Time PDE Modeling)

* **Performance Portable**: Uses the Yet Another Kernel Launcher (YAKL) framework for C++ performance portability
* **Flexible**: The C++ framework allows different dimensionalities, grids, splitting techniques, and PDEs. What all models must have in common is stepping a model state forward in time using spatial and temporal operators.
* **Clear**: Rather than significant levels of abstraction, which can make it unclear what is going on, abstraction is used only when it is clearly needed for generality. Most of the work in PFAST will go into the `Spatial` classes, which determing the dimensionality, grid, splitting, and PDEs being used. Because of this, PFAST encourages a common set of routines for code reusability rather than abstraction. Also, the input files are in easy-to-read, easy-to-write YAML files

## Driver

The driver is a simple C++ program that includes a spatial operator, then includes a temporal operator, and then steps forward in time using time steps given by the spatial operator.

## Model Flow

**Driver**: The driver makes calls to file I/O, time step computation, and the time stepping routine.

**Temporal**:
* The `Temporal` class inherits the `Spatial` class
* `Temporal` calls the `Spatial` class's `computeTendencies` function in order to get the tendencies for a single time step.
* PFAST supports arbitrary operator splitting (defined entirely by the `Spatial` class, multiple time derivatives, and time-averaged tendencies for certain ADER temporal operators.
* The `Temporal` operator creates a C++ lambda called `applySingleTendency`, which is passes to the `Spatial` class's `applyTendencies` function to tell the `Spatial` class what to do with the tendencies
* The `Temporal` class also cycles through all of the split stages each time step, and the `Spatial` class decides how many there are and what to do in each stage

**Spatial**:
* Most of the code lives in the `Spatial` class.
* The `Spatial` class defines its own `StateArr` and `TendArr` classes to hold the model's state and tendencies, respectively, allowing for any desired representation of the model state.
* The `Spatial` class initializes the model state, handles file I/O, defines how many stages of splitting there are, computes the tendencies for each splitting stage, and applies the tendencies for each splitting stage.
* The `Spatial` class implements a `Location` class to define a single location on the grid for a given spltting stage, as well as `get` functions to obtain the state and tendency at a given `Location` from `stateArr` and `tendArr`.
* The `Spatial` class iterates over the `applySingleTendency` from the `Temporal` object to apply tendencies to whatever portions of its grid are appropritate for the given splitting stage.
* The splitting capability adds significant flexibility to the models. Together with the arbitrary definitions of `stateArr` and `tendArr` as well as the `get` functions and `Location` class, it can handle PDE operator splitting, or even adaptive multi-grid implementations.

## Interfaces

### Temporal

```C++
class Temporal {
  /* REQUIRED:
  int  static constexpr nTimeDerivs = [# tendency time derivatives needed by the time stepping scheme];
  bool static constexpr timeAvg     = [whether the spatial operator should return a time-averaged tendency];

  static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

  void init(std::string inFile):
      - Process input file with filename "inFile"
      - Allocate and initialize internal stuff

  void timeStep( StateArr &state , real dt ): 
      - Perform a single time step

  const char * getTemporalName():
      - Return the name and info about this temporal operator
  */
};
```

### Spatial

```C++
class Spatial {
  /*********************************************
   ***** Required inside the Spatial class *****
   *********************************************
  class Location;
      - Stores the indices of a single location on the grid

  class StateArr; // OR
  typedef [class] StateArr;
      - Declare a type for the model state

  class TendArr; // OR
  typedef [class] TendArr;
      - Declare a type for the model tendencies
      - It must have room for the nTimeDerivs dimension for the time integrator

  StateArr createStateArr()
      - Create and return a new StateArr object

  TendArr createTendArr(int nTimeDerivs)
      - Create and return a new TendArr object

  YAKL_INLINE real &get(StateArr const &state , Location const &loc , int splitIndex)
      - Return the state value at the given location for this split stage

  YAKL_INLINE real &get(TendArr const &tend , Location const &loc , int timeDeriv , int splitIndex)
      - Return the tendency value at the given location for this split stage

  int numSplit()
      - Return the number of split components for this operator
      - The temporal operator will iterate through the splittings

  real computeTimeStep(real cfl)
      - Return the time step in seconds based on the cfl value

  void init(int nTimeDerivs, bool timeAvg, std::string inFile)
      - Initialize internal data structures
      - Read YAML input file for any needed parameters

  void initState( StateArr &state )
      - Initialize the state

  void computeTendencies( StateArr const &state , TendArr &tend , real dt , int splitIndex)
      - Compute tendency and time derivatives of the tendency if they are requested for this split stage

  template <class F> void applyTendencies( F const &applySingleTendency , int splitIndex )
      - Loop through the domain, and apply tendencies to the state for this split stage

  const char * getSpatialName()
      - Return the name and info for this spatial operator

  void output(StateArr const &state, real etime)
      - Output to file
  */
}
```

