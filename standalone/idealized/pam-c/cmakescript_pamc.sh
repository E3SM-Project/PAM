#!/bin/bash

add_cmake_vars="  -DPAMC_MODEL=layermodel "
add_cmake_vars+=" -DPAMC_HAMIL=tswe       "
add_cmake_vars+=" -DPAMC_THERMO=none      "
add_cmake_vars+=" -DPAMC_IO=serial        "
add_cmake_vars+=" -DPAMC_INNERMPI=true    "



