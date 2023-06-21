#pragma once

#include "common.h"
#include "field_sets.h"
#include "model.h"

real norm(FieldSet<nprognostic> &x) {
  // note that this function assumes that x halos have been exchanged
  real accum = 0;
  for (auto f : x.fields_arr) {
    accum = std::max(yakl::intrinsics::maxval(yakl::intrinsics::abs(f.data)),
                     accum);
  }
  return accum;
}

class TimeIntegrator {
public:
  bool is_initialized = false;
  bool is_ssp = false;
  bool is_semi_implicit = false;
  std::string tstype;

  virtual void initialize(ModelParameters &params, Tendencies &tend,
                          LinearSystem &linsys, FieldSet<nprognostic> &xvars,
                          FieldSet<nconstant> &consts,
                          FieldSet<nauxiliary> &auxiliarys) = 0;

  virtual void step_forward(real dt) = 0;

  TimeIntegrator(const std::string tstype) : tstype(tstype) {}
  virtual ~TimeIntegrator() = default;
};
