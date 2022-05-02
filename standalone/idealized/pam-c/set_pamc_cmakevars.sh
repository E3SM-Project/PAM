#!/bin/bash

add_cmake_vars="  -DPAMC_MODEL=extrudedmodel "
add_cmake_vars+=" -DPAMC_HAMIL=mce_rho       "
add_cmake_vars+=" -DPAMC_THERMO=constkappavirpottemp      "
add_cmake_vars+=" -DPAMC_IO=serial        "



