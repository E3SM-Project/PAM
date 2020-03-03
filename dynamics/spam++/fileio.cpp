


//THIS IS ACTUALLY AN IO CLASS FOR the advection equation
// Can it be generalized?
// Yes- outputInit does constant and prognostic variables
// output does prognostic variables

// CAREFUL HERE THOUGH- WE ALSO HAVE DIAGNOSTIC QUANTITIES that will be useful later on...
// although maybe we can punt on this for now?

  void FileIO::initalize() {
     file.open ("output.txt");
  }

  void FileIO::output(VariableSet &vars, int nstep, real time) {
    // open netcdf file
    // open variables for constants and vars
    // write variables for constants and vars
    // write elapsed time/time step

    file << "Variables at step" << nstep << " and t=" << time << "\n";
    for (int i=0; i<vars.field_arr.size(); i++)
    {
      file << vars.field_arr[i].name << "\n";
    }

  }

  void FileIO::outputInit(VariableSet &vars, VariableSet &const_vars, real time) {

    file << "Constants at step" << 0 << " and t=" << time << "\n";
    for (int i=0; i<const_vars.field_arr.size(); i++)
    {
      file << const_vars.field_arr[i].name << "\n";
    }

    file << "Variables at step" << 0 << " and t=" << time << "\n";
    for (int i=0; i<vars.field_arr.size(); i++)
    {
      file << vars.field_arr[i].name << "\n";
    }
    // create netcdf file
    // create dimensions
    // create variables for constants and vars
    // write variables for constants and vars
    // write elapsed time/time step
  }

  void FileIO::close() {
    file.close();
  }
