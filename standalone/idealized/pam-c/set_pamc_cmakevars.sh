#!/bin/bash

add_cmake_vars="  -DPAMC_MODEL=extrudedmodel "
add_cmake_vars+=" -DPAMC_HAMIL=man       "
add_cmake_vars+=" -DPAMC_THERMO=constkappavirpottemp      "
add_cmake_vars+=" -DPAMC_IO=serial        "

#ce mce_rho an man
#idealgaspottemp constkappavirpottemp
