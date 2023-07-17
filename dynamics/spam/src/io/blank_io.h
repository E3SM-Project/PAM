#pragma once

#include "common.h"
#include "field_sets.h"
#include "model.h"
#include "stats.h"

namespace pamc {

class FileIO {

public:
  FileIO() {}
  FileIO(const FileIO &fio) = delete;
  FileIO &operator=(const FileIO &fio) = delete;
  void initialize(std::string outputName, Topology &ptopo, Topology &dtopo,
                  Parallel &par, const FieldSet<nprognostic> &progvars,
                  const FieldSet<nconstant> &const_vars,
                  const std::vector<std::unique_ptr<Diagnostic>> &diag,
                  Stats &stats) {}
  void output(real time) {}
  void outputInit(real time, const Geometry<Straight> &primal_geometry,
                  const Geometry<Twisted> &dual_geometry,
                  const ModelParameters &params) {}
  void outputStats(const Stats &stats) {}
};
} // namespace pamc
