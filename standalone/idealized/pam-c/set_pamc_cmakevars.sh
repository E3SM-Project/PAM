#!/bin/bash

add_cmake_vars="  -DPAMC_MODEL=extrudedmodel "
add_cmake_vars+=" -DPAMC_HAMIL=ce       "
add_cmake_vars+=" -DPAMC_THERMO=idealgaspottemp      "
add_cmake_vars+=" -DPAMC_IO=serial        "
