#pragma once

#include "common.h"

#ifdef PAMC_PARALLELIO
#include "yakl_parallel_io.h"
#endif

#ifdef PAMC_SERIALIO
#include "yakl_serial_io.h"
#endif

#ifdef PAMC_NOIO
#include "blank_io.h"
#endif
