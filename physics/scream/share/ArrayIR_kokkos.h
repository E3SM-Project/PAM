
#pragma once

namespace pam {


  template <class T>
  inline Kokkos::View<T*,Kokkos::LayoutRight,ekat::HostDevice>
  arrayIR_to_kokkos_view( ArrayIR<T,1,IRMemHost> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      std::cerr << "Error: converting non-initialized ArrayIR object into Array is not allowed" << std::endl;
      throw "";
    }
    return Kokkos::View<T*,Kokkos::LayoutRight,ekat::HostDevice>( arrayIR.get_data()  ,
                                                                  arrayIR.get_dims()[0] );
  }

  template <class T>
  inline Kokkos::View<T*,Kokkos::LayoutRight,ekat::DefaultDevice>
  arrayIR_to_kokkos_view( ArrayIR<T,1,IRMemDevice> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      std::cerr << "Error: converting non-initialized ArrayIR object into Array is not allowed" << std::endl;
      throw "";
    }
    return Kokkos::View<T*,Kokkos::LayoutRight,ekat::DefaultDevice>( arrayIR.get_data()  ,
                                                                     arrayIR.get_dims()[0] );
  }

  template <class T>
  inline Kokkos::View<T**,Kokkos::LayoutRight,ekat::HostDevice>
  arrayIR_to_kokkos_view( ArrayIR<T,2,IRMemHost> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      std::cerr << "Error: converting non-initialized ArrayIR object into Array is not allowed" << std::endl;
      throw "";
    }
    return Kokkos::View<T**,Kokkos::LayoutRight,ekat::HostDevice>( arrayIR.get_data()    ,
                                                                   arrayIR.get_dims()[0] ,
                                                                   arrayIR.get_dims()[1] );
  }

  template <class T>
  inline Kokkos::View<T**,Kokkos::LayoutRight,ekat::DefaultDevice>
  arrayIR_to_kokkos_view( ArrayIR<T,2,IRMemDevice> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      std::cerr << "Error: converting non-initialized ArrayIR object into Array is not allowed" << std::endl;
      throw "";
    }
    return Kokkos::View<T**,Kokkos::LayoutRight,ekat::DefaultDevice>( arrayIR.get_data()    ,
                                                                      arrayIR.get_dims()[0] ,
                                                                      arrayIR.get_dims()[1] );
  }

  template <class T>
  inline Kokkos::View<T***,Kokkos::LayoutRight,ekat::HostDevice>
  arrayIR_to_kokkos_view( ArrayIR<T,3,IRMemHost> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      std::cerr << "Error: converting non-initialized ArrayIR object into Array is not allowed" << std::endl;
      throw "";
    }
    return Kokkos::View<T***,Kokkos::LayoutRight,ekat::HostDevice>( arrayIR.get_data()    ,
                                                                    arrayIR.get_dims()[0] ,
                                                                    arrayIR.get_dims()[1] ,
                                                                    arrayIR.get_dims()[2]  );
  }

  template <class T>
  inline Kokkos::View<T***,Kokkos::LayoutRight,ekat::DefaultDevice>
  arrayIR_to_kokkos_view( ArrayIR<T,3,IRMemDevice> const &arrayIR ) {
    if (! arrayIR.initialized()) {
      std::cerr << "Error: converting non-initialized ArrayIR object into Array is not allowed" << std::endl;
      throw "";
    }
    return Kokkos::View<T***,Kokkos::LayoutRight,ekat::DefaultDevice>( arrayIR.get_data()    ,
                                                                       arrayIR.get_dims()[0] ,
                                                                       arrayIR.get_dims()[1] ,
                                                                       arrayIR.get_dims()[2]  );
  }

}


