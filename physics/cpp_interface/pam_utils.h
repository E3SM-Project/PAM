
#pragma once

#include <type_traits>
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/ekat_session.hpp"

#include "pam_coupler.h"

// for kokkos debuge only
#if defined(DEBUG)
#include "Cuda/Kokkos_Cuda_Instance.hpp"
#endif

// scream session initialize/finalize
void initialize_session (bool print_config = true) {
  ekat::initialize_ekat_session(print_config);

  // Make sure scream only has its FPEs
  ekat::disable_all_fpes();
  ekat::enable_fpes(ekat::get_default_fpes());
}

void finalize_session () {
  ekat::finalize_ekat_session();
}

// reshape an input array (nj, nk) to output array (nk, nj)
template <typename Scalar>
void reshape(const Scalar* sv, Scalar* dv, int nj, int nk) {
  parallel_for( SimpleBounds<2>(nj, nk) , YAKL_LAMBDA (int j, int k) {
    dv[k*nj+j] = sv[j*nk+k];
  });
}

// reshape and inverse an input array (nj, nk) to output array (nk, nj)
template <typename Scalar>
void reshape_and_inverse(const Scalar* sv, Scalar* dv, int nj, int nk) {
  parallel_for( SimpleBounds<2>(nj, nk) , YAKL_LAMBDA (int j, int k) {
    dv[k*nj+j] = 1./sv[j*nk+k];
  });
}

// convert YAKL multidimensional array to Kokkos view on GPU
template <typename SizeT, typename ViewT>
void array_to_view(const typename ViewT::value_type::scalar* const data,
                   const SizeT& size,
                   ViewT& view)
{
  using PackT = typename ViewT::value_type;
  EKAT_ASSERT(PackT::n >= 1);

  const size_t npack = (size + PackT::n-1)/PackT::n;

#if defined(DEBUG)
  kokkos_impl_cuda_set_serial_execution(true);
#endif
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {npack, PackT::n}), KOKKOS_LAMBDA(int k, int s) {
    const size_t scalar_offset = k*PackT::n;
    if (scalar_offset+s < size) view(k)[s] = data[scalar_offset+s];
  });
}

// 2D YAKL array to Kokkos view
template <typename SizeT, typename ViewT>
void array_to_view(const typename ViewT::value_type::scalar* const data,
                   const SizeT& dim1_size,
                   const SizeT& dim2_size,
                   ViewT& view)
{
  using PackT = typename ViewT::value_type;
  EKAT_ASSERT(PackT::n >= 1);

  const int pack_size = static_cast<int>(PackT::n);
  const int npack     = (dim2_size+pack_size-1)/pack_size;

#if defined(DEBUG)
  kokkos_impl_cuda_set_serial_execution(true);
#endif
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {dim1_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int k, int s) {
     const int num_scalars = k*pack_size;
     const int scalar_offset = i*dim2_size+num_scalars;
     if (num_scalars+s<dim2_size)  view(i,k)[s] = data[scalar_offset+s];
  });
}

// 3D YAKL to Kokkos views
template <typename SizeT, typename ViewT>
void array_to_view(const typename ViewT::value_type::scalar* const data,
                   const SizeT& dim1_size,
                   const SizeT& dim2_size,
                   const SizeT& dim3_size,
                   ViewT& view)
{
  using PackT = typename ViewT::value_type;
  EKAT_ASSERT(PackT::n >= 1);

  const int pack_size = static_cast<int>(PackT::n);
  const int npack     = (dim3_size+pack_size-1)/pack_size;

#if defined(DEBUG)
  kokkos_impl_cuda_set_serial_execution(true);
#endif
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {dim1_size, dim2_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int j, int k, int s) {
     const int num_scalars = k*pack_size;
     const int scalar_offset = (i*dim2_size+j)*dim3_size+num_scalars;
     if (num_scalars+s<dim3_size) view(i,j,k)[s] = data[scalar_offset+s];
  });
}

// 1D Kokkos view to YAKL array
template <typename SizeT, typename ViewT, typename ArrayT>
void view_to_array(const ViewT& view,
                   const SizeT& size,
                   ArrayT& array)
{
  using PackT      = typename ViewT::value_type;
  using scalarType = typename ViewT::value_type::scalar;
  using arrayType  = typename ArrayT::type;

  auto is_same_type = std::is_same<scalarType, arrayType>::value;

  EKAT_ASSERT(is_same_type);
  EKAT_ASSERT(PackT::n >= 1);

  const size_t npack = (size + PackT::n-1)/PackT::n;
    
#if defined(DEBUG)
  kokkos_impl_cuda_set_serial_execution(true);
#endif
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {npack, PackT::n}), KOKKOS_LAMBDA(int k, int s) {
      const size_t scalar_offset = k*PackT::n;
      if (scalar_offset+s < size) array.myData[scalar_offset+s] = view(k)[s];
  });
}

// 2D Kokkos view to YAKL array
template <typename SizeT, typename ViewT, typename ArrayT>
void view_to_array(const ViewT& view,
                   const SizeT& dim1_size,
                   const SizeT& dim2_size,
                   ArrayT& array)
{
  using PackT      = typename ViewT::value_type;
  using scalarType = typename ViewT::value_type::scalar;
  using arrayType  = typename ArrayT::type;
 
  auto is_same_type = std::is_same<scalarType, arrayType>::value;

  EKAT_ASSERT(is_same_type);
  EKAT_ASSERT(PackT::n >= 1);

  const int pack_size = static_cast<int>(PackT::n);
  const int npack     = (dim2_size+pack_size-1)/pack_size;

#if defined(DEBUG)
  kokkos_impl_cuda_set_serial_execution(true);
#endif
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {dim1_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int k, int s) {
    const int num_scalars = k*pack_size;
    const int scalar_offset = i*dim2_size + num_scalars;
    if (num_scalars+s<dim2_size)  array.myData[scalar_offset+s] = view(i,k)[s];
  });
}

// 3D Kokkos view to YAKL array
template <typename SizeT, typename ViewT, typename ArrayT>
void view_to_array(const ViewT& view,
                   const SizeT& dim1_size,
                   const SizeT& dim2_size,
                   const SizeT& dim3_size,
                   ArrayT& array)
{
  using PackT      = typename ViewT::value_type;
  using scalarType = typename ViewT::value_type::scalar;
  using arrayType  = typename ArrayT::type;

  auto is_same_type = std::is_same<scalarType, arrayType>::value;

  EKAT_ASSERT(is_same_type);
  EKAT_ASSERT(PackT::n >= 1);

  const int pack_size = static_cast<int>(PackT::n);
  const int npack     = (dim3_size+pack_size-1)/pack_size;

#if defined(DEBUG)
  kokkos_impl_cuda_set_serial_execution(true);
#endif
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {dim1_size, dim2_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int j, int k, int s) {
    const int num_scalars = k*pack_size;
    const int scalar_offset = (i*dim2_size+j)*dim3_size+num_scalars;
    if (num_scalars+s<dim3_size)  array.myData[scalar_offset+s] = view(i,j,k)[s];
  });
}

// validation code 
template <typename SizeT, typename ViewT>
void array_to_view_2d(typename ViewT::value_type::scalar* data,
                      const SizeT& dim1_sizes,
                      const SizeT& dim2_sizes,
                      ViewT& view)
{
  using PackT = typename ViewT::value_type;
  EKAT_ASSERT(PackT::n >= 1);

    const int dim1_size = static_cast<int>(dim1_sizes);
    const int dim2_size = static_cast<int>(dim2_sizes);
    const int pack_size = static_cast<int>(PackT::n);
    const int npack     = (dim2_size+pack_size-1)/pack_size;

#if defined(DEBUG)
    kokkos_impl_cuda_set_serial_execution(true);
#endif
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {dim1_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int k, int s) {
       const int num_scalars = k*pack_size;
       const int scalar_offset = i*dim2_size + num_scalars;
       if (num_scalars+s<dim2_size) view(i,k)[s] = data[scalar_offset+s];
   });
}

