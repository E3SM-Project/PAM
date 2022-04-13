#pragma once

#include "common.h"

class FileIO {

public:
  
  FileIO();
  FileIO( const FileIO &fio) = delete;
  FileIO& operator=( const FileIO &fio) = delete;
  void initialize(std::string outputName, Topology &ptopo, Topology &dtopo, Parallel &par, const VariableSet<nprognostic> &progvars, const VariableSet<nconstant> &const_vars, const VariableSet<ndiagnostic> &diagvars, Stats &stats);
  void output(real time);
  void outputInit(real time);
  void outputStats(const Stats &stats);
};


FileIO::FileIO()
{

}

 void FileIO::initialize(std::string outName, Topology &ptopo, Topology &dtopo, Parallel &par, const VariableSet<nprognostic> &progvars, const VariableSet<nconstant> &constvars, const VariableSet<ndiagnostic> &diagvars, Stats &stats)
{

}
   
   void FileIO::output(real time)
   {        

   }
   
  void FileIO::outputInit(real time)
   {

     
   }
   
   void FileIO::outputStats(const Stats &stats)
   {


   }

   
   
