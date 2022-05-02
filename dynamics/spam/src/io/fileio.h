#pragma once

#include "common.h"

#ifdef _PARALLELIO
#include "yakl_parallel_io.h"
#endif

#ifdef _SERIALIO
#include "yakl_serial_io.h"
#endif

#ifdef _NOIO
#include "blank_io.h"
#endif
