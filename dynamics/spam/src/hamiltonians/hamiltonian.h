#pragma once

#include "anelastic.h"
#include "common.h"
#include "compressible_euler.h"
#include "functionals.h"
#include "kinetic_energy.h"
#include "layer_models.h"

#ifdef _SWE
using Hamiltonian = Hamiltonian_SWE_Hs;
#elif _TSWE
using Hamiltonian = Hamiltonian_TSWE_Hs;
#elif _CE
using Hamiltonian = Hamiltonian_CE_Hs;
#elif _AN
using Hamiltonian = Hamiltonian_AN_Hs;
#elif _MAN
using Hamiltonian = Hamiltonian_MAN_Hs;
#elif _MCErho
using Hamiltonian = Hamiltonian_MCE_Hs;
#elif _MCErhod
using Hamiltonian = Hamiltonian_MCE_Hs;
#elif _CEp
using Hamiltonian = Hamiltonian_CE_p_Hs;
#elif _MCErhop
using Hamiltonian = Hamiltonian_MCE_p_Hs;
#elif _MCErhodp
using Hamiltonian = Hamiltonian_MCE_p_Hs;
#endif
