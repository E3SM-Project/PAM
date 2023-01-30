
#pragma once

#include <type_traits>
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/ekat_session.hpp"
#include "ArrayIR.h"
#include "share/scream_types.hpp"

// for kokkos debuge only
#if defined(DEBUG)
#include "Cuda/Kokkos_Cuda_Instance.hpp"
#endif


namespace ScreamCXX {

  using array_ir::ArrayIR;
  using Kokkos::View;
  using Kokkos::LayoutRight;
  using ekat::Pack;

  inline void die(std::string msg) {
    std::cerr << "ERROR:\n" << msg << std::endl;
    throw "";
  }

  template <typename MemorySpace>
  inline void check_memory_spaces( bool valid_on_host , bool valid_on_device ) {
    if ( std::is_same<MemorySpace,Kokkos::HostSpace>::value ) {
      if (! valid_on_host  ) die("View is only valid on the host, but ArrayIR not valid on the host");
    }
    #ifdef KOKKOS_ENABLE_CUDA
      if ( std::is_same<MemorySpace,Kokkos::CudaSpace>::value ) {
        if (! valid_on_device) die("View is only valid on the device, but ArrayIR not valid on the device");
      }
    #endif
    #ifdef KOKKOS_ENABLE_HIP
      if ( std::is_same<MemorySpace,Kokkos::HIPSpace>::value ) {
        if (! valid_on_device) die("View is only valid on the device, but ArrayIR not valid on the device");
      }
    #endif
  }

  template <class T,
            typename std::enable_if<std::is_arithmetic<T>::value,bool>::type = false >
  inline View<Pack<T,1> *,LayoutRight,scream::DefaultDevice> ArrayIR_to_View_of_Packs( ArrayIR<T,1> const &arr ) {
    check_memory_spaces<scream::DefaultDevice::memory_space>( arr.data_valid_on_host() , arr.data_valid_on_device() );
    return View<Pack<T,1> *,LayoutRight,scream::DefaultDevice::memory_space>( static_cast<Pack<T,1> *>(static_cast<void *>(arr.data())) , arr.dimension(0) );
  }

  template <class T,
            typename std::enable_if<std::is_arithmetic<T>::value,bool>::type = false >
  inline View<Pack<T,1> **,LayoutRight,scream::DefaultDevice> ArrayIR_to_View_of_Packs( ArrayIR<T,2> const &arr ) {
    check_memory_spaces<scream::DefaultDevice::memory_space>( arr.data_valid_on_host() , arr.data_valid_on_device() );
    return View<Pack<T,1> **,LayoutRight,scream::DefaultDevice::memory_space>( static_cast<Pack<T,1> *>(static_cast<void *>(arr.data())) , arr.dimension(0) , arr.dimension(1) );
  }

  template <class T,
            typename std::enable_if<std::is_arithmetic<T>::value,bool>::type = false >
  inline View<Pack<T,1> ***,LayoutRight,scream::DefaultDevice> ArrayIR_to_View_of_Packs( ArrayIR<T,3> const &arr ) {
    check_memory_spaces<scream::DefaultDevice::memory_space>( arr.data_valid_on_host() , arr.data_valid_on_device() );
    return View<Pack<T,1> ***,LayoutRight,scream::DefaultDevice::memory_space>( static_cast<Pack<T,1> *>(static_cast<void *>(arr.data())) , arr.dimension(0) , arr.dimension(1) , arr.dimension(2) );
  }

  template <class T,
            typename std::enable_if<std::is_arithmetic<T>::value,bool>::type = false >
  inline View<T *,LayoutRight,scream::DefaultDevice> ArrayIR_to_View( ArrayIR<T,1> const &arr ) {
    check_memory_spaces<scream::DefaultDevice::memory_space>( arr.data_valid_on_host() , arr.data_valid_on_device() );
    return View<T *,LayoutRight,scream::DefaultDevice::memory_space>( arr.data() , arr.dimension(0) );
  }

  template <class T,
            typename std::enable_if<std::is_arithmetic<T>::value,bool>::type = false >
  inline View<T **,LayoutRight,scream::DefaultDevice> ArrayIR_to_View( ArrayIR<T,2> const &arr ) {
    check_memory_spaces<scream::DefaultDevice::memory_space>( arr.data_valid_on_host() , arr.data_valid_on_device() );
    return View<T **,LayoutRight,scream::DefaultDevice::memory_space>( arr.data() , arr.dimension(0) , arr.dimension(1) );
  }

  template <class T,
            typename std::enable_if<std::is_arithmetic<T>::value,bool>::type = false >
  inline View<T ***,LayoutRight,scream::DefaultDevice> ArrayIR_to_View( ArrayIR<T,3> const &arr ) {
    check_memory_spaces<scream::DefaultDevice::memory_space>( arr.data_valid_on_host() , arr.data_valid_on_device() );
    return View<T ***,LayoutRight,scream::DefaultDevice::memory_space>( arr.data() , arr.dimension(0) , arr.dimension(1) , arr.dimension(2) );
  }

}



