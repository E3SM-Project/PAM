
#include "Spatial_expl.h"
#include "Temporal_ader.h"

// Define the Spatial operator based on constants from the Temporal operatora header file
typedef Spatial_operator Spatial;

// Define the Temporal operator based on the Spatial operator
typedef Temporal_operator<Spatial> Dycore;

