
#pragma once

#include <type_traits>
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/ekat_session.hpp"
#include "ArrayIR.h"

// for kokkos debuge only
#if defined(DEBUG)
#include "Cuda/Kokkos_Cuda_Instance.hpp"
#endif
//
// plan data layout in the ArrayIR, take the nomenclature from dnn
// batch N(1), channels C(2), depth D(3), height H(4), width W(5)
// for example, offset_nchw(n, c, h, w) = n * CHW + c * HW + h * W + w
//
enum DataFormat {
  NCHW = 0,
  NCWH = 1,
  NHWC = 2,
  NWCH = 3,
  NWHC = 4
};

// using arrayIR instead of carray from YAKL
template<typename T, int N>
struct IRType {
  using type = typename ArrayIR::ArrayIR<T, N>;
  static constexpr bool is_irtype = true;
  // dimensions of the arrayIR
  enum { Ndim = N };
};

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

// convert YAKL multidimensional array to Kokkos view on GPU
template <typename SizeT, typename ViewT>
void array_to_view(const typename ViewT::value_type::scalar* const data, 
                   const DataFormat& format,
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
                   const DataFormat& format,
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
 switch (format)
 {
  case DataFormat::NCHW:
     Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {dim1_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int k, int s) {
       const int num_scalars = k*pack_size;
       const int scalar_offset = i*dim2_size+num_scalars;
       if (num_scalars+s<dim2_size)  view(i,k)[s] = data[scalar_offset+s];
     });
     break;
  case DataFormat::NCWH:
     Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {dim1_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int k, int s) {
       const int num_scalars   = k*pack_size+s;
       const int scalar_offset = num_scalars*dim1_size+i;
       if (num_scalars+s<dim2_size)  view(i,k)[s] = data[scalar_offset];
     });
     break;
  default:
     std::cout << "2D plan data format error!\n";
     break;
 }
}

// 3D YAKL to Kokkos views
template <typename SizeT, typename ViewT>
void array_to_view(const typename ViewT::value_type::scalar* const data,
                   const DataFormat& format,
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
 switch (format)
 {
  case DataFormat::NCHW:
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {dim1_size, dim2_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int j, int k, int s) {
       const int num_scalars = k*pack_size;
       const int scalar_offset = (i*dim2_size+j)*dim3_size+num_scalars;
       if (num_scalars+s<dim3_size) view(i,j,k)[s] = data[scalar_offset+s];
    });
    break;
  case DataFormat::NWCH:
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {dim1_size, dim2_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int j, int k, int s) {
       const int num_scalars   = k*pack_size+s;
       const int scalar_offset = (j*dim3_size+num_scalars)*dim1_size+i;
       if (num_scalars+s<dim3_size) view(i,j,k)[s] = data[scalar_offset];
    });
    break;
  default:
     std::cout << "3D plan data format error!\n";
     break;
 }
}

// 1D Kokkos view to YAKL array
template <typename SizeT, typename ViewT, typename ArrayT>
void view_to_array(const ViewT& view,
                   const SizeT& size,
                   const DataFormat& format,
                   ArrayT& array)
{
  using PackT      = typename ViewT::value_type;
  using scalarType = typename ViewT::value_type::scalar;
  using arrayType  = typename ArrayT::remove_cv_type;

  auto is_same_type = std::is_same<scalarType, arrayType>::value;

  EKAT_ASSERT(is_same_type);
  EKAT_ASSERT(PackT::n >= 1);

  const size_t npack = (size + PackT::n-1)/PackT::n;
    
#if defined(DEBUG)
  kokkos_impl_cuda_set_serial_execution(true);
#endif
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {npack, PackT::n}), KOKKOS_LAMBDA(int k, int s) {
      const size_t scalar_offset = k*PackT::n;
      if (scalar_offset+s < size) array.data()[scalar_offset+s] = view(k)[s];
  });
}

// 2D Kokkos view to YAKL array
template <typename SizeT, typename ViewT, typename ArrayT>
void view_to_array(const ViewT& view,
                   const SizeT& dim1_size,
                   const SizeT& dim2_size,
                   const DataFormat& format,
                   ArrayT& array)
{
  using PackT      = typename ViewT::value_type;
  using scalarType = typename ViewT::value_type::scalar;
  using arrayType  = typename ArrayT::remove_cv_type;
 
  auto is_same_type = std::is_same<scalarType, arrayType>::value;

  EKAT_ASSERT(is_same_type);
  EKAT_ASSERT(PackT::n >= 1);

  const int pack_size = static_cast<int>(PackT::n);
  const int npack     = (dim2_size+pack_size-1)/pack_size;

#if defined(DEBUG)
  kokkos_impl_cuda_set_serial_execution(true);
#endif
 switch (format)
 {
  case DataFormat::NCHW:
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {dim1_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int k, int s) {
      const int num_scalars = k*pack_size;
      const int scalar_offset = i*dim2_size + num_scalars;
      if (num_scalars+s<dim2_size)  array.data()[scalar_offset+s] = view(i,k)[s];
    });
    break;
  case DataFormat::NCWH:
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {dim1_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int k, int s) {
      const int num_scalars   = k*pack_size+s;
      const int scalar_offset = num_scalars*dim1_size+i;
      if (num_scalars+s<dim2_size)  array.data()[scalar_offset] = view(i,k)[s];
    });
    break;
  default:
    std::cout << "2D plain data layout error!\n";
    break;
 }
}

// 3D Kokkos view to YAKL array
template <typename SizeT, typename ViewT, typename ArrayT>
void view_to_array(const ViewT& view,
                   const SizeT& dim1_size,
                   const SizeT& dim2_size,
                   const SizeT& dim3_size,
                   const DataFormat& format,
                   ArrayT& array)
{
  using PackT      = typename ViewT::value_type;
  using scalarType = typename ViewT::value_type::scalar;
  using arrayType  = typename ArrayT::remove_cv_type;

  auto is_same_type = std::is_same<scalarType, arrayType>::value;

  EKAT_ASSERT(is_same_type);
  EKAT_ASSERT(PackT::n >= 1);

  const int pack_size = static_cast<int>(PackT::n);
  const int npack     = (dim3_size+pack_size-1)/pack_size;

#if defined(DEBUG)
  kokkos_impl_cuda_set_serial_execution(true);
#endif
 switch (format)
 {
  case DataFormat::NCHW:
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {dim1_size, dim2_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int j, int k, int s) {
      const int num_scalars = k*pack_size;
      const int scalar_offset = (i*dim2_size+j)*dim3_size+num_scalars;
      if (num_scalars+s<dim3_size)  array.data()[scalar_offset+s] = view(i,j,k)[s];
    });
    break;
  case DataFormat::NWCH:
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {dim1_size, dim2_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int j, int k, int s) {
       const int num_scalars   = k*pack_size+s;
       const int scalar_offset = (j*dim3_size+num_scalars)*dim1_size+i;
      if (num_scalars+s<dim3_size)  array.data()[scalar_offset] = view(i,j,k)[s];
    });
    break;
  default:
    std::cout << "3D plain data layout error!\n";
    break;
 }
}

// validation code 
template <typename SizeT, typename ViewT>
void array_to_view_2d(typename ViewT::value_type::scalar* data,
                      const SizeT& dim1_sizes,
                      const SizeT& dim2_sizes,
                      const DataFormat& format,
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
 switch (format)
 {
  case DataFormat::NCHW:
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {dim1_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int k, int s) {
       const int num_scalars = k*pack_size;
       const int scalar_offset = i*dim2_size + num_scalars;
       if (num_scalars+s<dim2_size) view(i,k)[s] = data[scalar_offset+s];
    });
    break;
  case DataFormat::NCWH:
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {dim1_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int k, int s) {
       const int num_scalars   = k*pack_size+s;
       const int scalar_offset = num_scalars*dim1_size+i;
       if (num_scalars+s<dim2_size) view(i,k)[s] = data[scalar_offset];
   });
   break;
  default:
    std::cout << "2D plain data layout error!\n";
    break;
 }
}

