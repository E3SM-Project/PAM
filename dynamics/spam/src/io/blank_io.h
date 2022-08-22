#pragma once

#include "common.h"
#include "field_sets.h"
#include "model.h"
#include "stats.h"

class FileIO {

public:
  FileIO();
  FileIO(const FileIO &fio) = delete;
  FileIO &operator=(const FileIO &fio) = delete;
  void initialize(std::string outputName, Topology &ptopo, Topology &dtopo,
                  Parallel &par, const FieldSet<nprognostic> &progvars,
                  const FieldSet<nconstant> &const_vars,
                  const std::vector<std::unique_ptr<Diagnostic>> &diag,
                  Stats &stats);
  void output(real time);
  void outputInit(real time);
  void outputStats(const Stats &stats);
};

FileIO::FileIO() {}

void FileIO::initialize(std::string outName, Topology &ptopo, Topology &dtopo,
                        Parallel &par, const FieldSet<nprognostic> &progvars,
                        const FieldSet<nconstant> &constvars,
                        const std::vector<std::unique_ptr<Diagnostic>> &diag,
                        Stats &stats) {}

void FileIO::output(real time) {}

void FileIO::outputInit(real time) {}

void FileIO::outputStats(const Stats &stats) {}
