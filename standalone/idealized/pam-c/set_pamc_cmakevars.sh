#!/bin/bash

add_cmake_vars="  -DPAMC_MODEL=extrudedmodel "
add_cmake_vars+=" -DPAMC_HAMIL=moistcompressibleeuler_rho       "
add_cmake_vars+=" -DPAMC_THERMO=constkappavirtualpottemp      "
add_cmake_vars+=" -DPAMC_IO=serial        "
