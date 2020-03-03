


//THIS IS ACTUALLY AN IO CLASS FOR the advection equation
// Can it be generalized?
// Yes- outputInit does constant and prognostic variables
// output does prognostic variables

// CAREFUL HERE THOUGH- WE ALSO HAVE DIAGNOSTIC QUANTITIES that will be useful later on...
// although maybe we can punt on this for now?

  void FileIO::initalize() {}

  void FileIO::output(VariableSet &vars) {
    // open netcdf file
    // open variables for constants and vars
    // write variables for constants and vars
    // write elapsed time/time step
  }

  void FileIO::outputInit(VariableSet &const_vars, VariableSet &vars) {
    // create netcdf file
    // create dimensions
    // create variables for constants and vars
    // write variables for constants and vars
    // write elapsed time/time step
  }
