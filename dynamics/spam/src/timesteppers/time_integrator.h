#pragma once

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
