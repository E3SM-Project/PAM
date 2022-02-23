#ifndef SCREAM_TYPES_HPP
#define SCREAM_TYPES_HPP

#include "ekat/ekat.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"

#define SCREAM_PACK_SIZE 1
#define SCREAM_SMALL_PACK_SIZE 1

namespace scream
{

// Scalar types
using ekat::Int;
// #ifdef SCREAM_DOUBLE_PRECISION
// using Real = double;
// #else
// using Real = float;
// #endif

using Real = double;

// Kokkos types
using ekat::KokkosTypes;
using ekat::DefaultDevice;
using ekat::HostDevice;
using ekat::Unmanaged;

// Miscellanea

// An enum to be used with object that have 'repository'-like behavior
enum class RepoState {
  Clean,
  Open,
  Closed
};

// We cannot expect BFB results between f90 and cxx if optimizations are on.
// Same goes for cuda-memcheck because it makes the bfb math layer prohibitively
// expensive and so must be turned off.
static constexpr bool SCREAM_BFB_TESTING = false;

/*
 * Utility function for handling floating point literals,
 * so that they match the scream precision. This is
 * especially useful for bfb tests agaisnt fortran,
 * to ensure that literals are not a source of round-off differences.
 */
template<typename T> KOKKOS_INLINE_FUNCTION
constexpr typename std::enable_if<std::is_arithmetic<T>::value,Real>::type
sp (const T val) {
  return Real(val);
}

} // namespace scream

#endif // SCREAM_TYPES_HPP
