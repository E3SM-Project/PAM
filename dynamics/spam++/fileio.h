
#ifndef _FILEIO_H_
#define _FILEIO_H_


#include "common.h"
#include <iostream>
#include <fstream>

class FileIO {
  ofstream file;

public:

  void initalize();
  void output(VariableSet &vars, int nstep, int time);
  void outputInit(VariableSet &vars, VariableSet &const_vars, real time);

};
