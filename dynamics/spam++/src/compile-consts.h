
#ifndef _COMPILE_CONSTS_H_
#define _COMPILE_CONSTS_H_

// Number of Dimensions
uint constexpr number_of_dims = 2;

// Spatial order of accuracy for the model
uint constexpr differential_order = 2;

// Reconstruction type
RECONSTRUCTION_TYPE reconstruction_type = RECONSTRUCTION_TYPE::CFV;
uint constexpr reconstruction_order = 2;

// Halo sizes
uint maxhalosize = mymax(reconstruction_order,differential_order)/2; // IS THIS ALWAYS CORRECT?

// initial condition quadrature pts
uint constexpr ic_quad_pts = 3;

// Time scheme
TIME_TYPE time_type = TIME_TYPE::KGRK;
uint constexpr n_time_stages = 4;

#endif
