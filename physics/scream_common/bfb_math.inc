
#ifndef SCREAM_BFB_MATH_INC
#define SCREAM_BFB_MATH_INC

#define bfb_square(val) ((val)*(val))
#define bfb_cube(val)   ((val)*(val)*(val))
#define bfb_quad(val)   (bfb_square(bfb_square(val)))

#define bfb_pow(base, exp) (base)**(exp)
#define bfb_cbrt(base) (base)**(1.0D0/3.0D0)
#define bfb_gamma(val) gamma(val)
#define bfb_log(val) log(val)
#define bfb_log10(val) log10(val)
#define bfb_exp(val) exp(val)
#define bfb_expm1(val) (exp(val) - 1)
#define bfb_tanh(val) tanh(val)
#define bfb_sqrt(val) sqrt(val)
#define bfb_tanh(val) tanh(val)
#define bfb_erf(val) erf(val)

#endif
