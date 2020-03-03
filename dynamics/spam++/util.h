
#ifndef _UTIL_H_
#define _UTIL_H_

#include "fields.h"
#include "geometry.h"

template<class F> void set_n_form(F &init_cond_func, Field &x, Geometry &geometry);
template<class F> void set_n_minus_1_form(F &init_cond_func, Field &x, Geometry &geometry);


#endif
