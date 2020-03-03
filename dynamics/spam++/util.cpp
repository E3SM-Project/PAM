
#include "util.h"

template<class F> void set_n_form(F &init_cond_func, Field &x, Geometry &geometry, int ndof) {

  if (ndims == 1) { geometry.set_one_form_values(init_cond_func, x, ndof); }
  if (ndims == 2) { geometry.set_two_form_values(init_cond_func, x, ndof); }
  if (ndims == 3) { geometry.set_three_form_values(init_cond_func, x, ndof); }

}

template<class F> void set_n_minus_1_form(F &init_cond_func, Field &x, Geometry &geometry, int ndof) {

  if (ndims == 1) { geometry.set_zero_form_values(init_cond_func, x, ndof) ;}
  if (ndims == 2) { geometry.set_one_form_values(init_cond_func, x, ndof); }
  if (ndims == 3) { geometry.set_two_form_values(init_cond_func, x, ndof); }

}
