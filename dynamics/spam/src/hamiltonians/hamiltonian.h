#pragma once

#include "anelastic.h"
#include "common.h"

#include "type_traits"

namespace pamc {
template <class HamilT, class ThermoT>
struct two_point_discrete_gradient_implemented : std::false_type {};
template <class HamilT, class ThermoT>
constexpr bool two_point_discrete_gradient_implemented_v =
    two_point_discrete_gradient_implemented<HamilT, ThermoT>::value;
} // namespace pamc

#include "compressible_euler.h"
#include "functionals.h"
#include "kinetic_energy.h"
#include "layer_models.h"

namespace pamc {
#ifdef PAMC_SWE
using Hamiltonian = Hamiltonian_SWE_Hs;
#elif PAMC_TSWE
using Hamiltonian = Hamiltonian_TSWE_Hs;
#elif PAMC_CE
using Hamiltonian = Hamiltonian_CE_Hs;
#elif PAMC_AN
using Hamiltonian = Hamiltonian_AN_Hs;
#elif PAMC_MAN
using Hamiltonian = Hamiltonian_MAN_Hs;
#elif PAMC_MCErho
using Hamiltonian = Hamiltonian_MCE_Hs;
#elif PAMC_MCErhod
using Hamiltonian = Hamiltonian_MCE_Hs;
#elif PAMC_CEp
using Hamiltonian = Hamiltonian_CE_p_Hs;
#elif PAMC_MCErhop
using Hamiltonian = Hamiltonian_MCE_p_Hs;
#elif PAMC_MCErhodp
using Hamiltonian = Hamiltonian_MCE_p_Hs;
#endif
} // namespace pamc
