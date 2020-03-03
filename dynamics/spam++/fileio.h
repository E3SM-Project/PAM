
#ifndef _FILEIO_H_
#define _FILEIO_H_


#include "common.h"
#include "STDLIB.h"

class FileIO {

public:

  void initalize();
  void output(VariableSet &vars);
  void outputInit(VariableSet &const_vars, VariableSet &vars);


};
