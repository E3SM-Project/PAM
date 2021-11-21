
#include "Spatial_cons_anelastic_fv_Agrid.h"
#include "Temporal_ader.h"

// Define the Spatial operator based on constants from the Temporal operatora header file
typedef Spatial_operator<nTimeDerivs,timeAvg,nAder> Spatial;

// Define the Temporal operator based on the Spatial operator
typedef Temporal_operator<Spatial> Dycore;

