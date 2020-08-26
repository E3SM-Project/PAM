
#include "const.h"
#include "PhysicsSaturationAdjustment.h"

int main() {
  PhysicsSaturationAdjustment physics;

  real rho_d = 1.16;
  real theta_d = 300;
  real p_v_sat = physics.saturationVaporPressure(theta_d);
  real rho_v_sat = p_v_sat / physics.R_v / theta_d;
  std::cout << "Saturation Vapor Pressure at 300 K: " << p_v_sat << "\n";
  std::cout << "vapor density @ SVP: " << rho_v_sat << "\n\n";

  {
    std::cout << "****** CONDENSATION TEST ******\n";
    real satRatio = 1.1;
    real rho_v = rho_v_sat * satRatio;
    real rho_c = 0;
    real rho = rho_d + rho_v + rho_c;
    real temp = theta_d;
    real theta = physics.thetaFromTemp(rho, rho_v, rho_c, temp);
    real rho_theta = rho*theta;

    std::cout << "Pre-adjusted (rho, rho_v, rho_c, rho_theta): \n" << std::setprecision(20)
                                                                   << rho       << ", \n"
                                                                   << rho_v     << ", \n"
                                                                   << rho_c     << ", \n"
                                                                   << rho_theta << "\n\n";
    real mass1 = rho_d + rho_v + rho_c;

    physics.computeAdjustedState(rho_d, rho_v, rho_c, rho_theta);

    std::cout << "Post-adjusted (rho, rho_v, rho_c, rho_theta): \n" << std::setprecision(20)
                                                                    << rho       << ", \n"
                                                                    << rho_v     << ", \n"
                                                                    << rho_c     << ", \n"
                                                                    << rho_theta << "\n\n";
    real mass2 = rho_d + rho_v + rho_c;
    std::cout << "Relative Mass Change: " << std::setprecision(20) << (mass2-mass1)/mass1 << "\n";
  }

  {
    std::cout << "****** EVAPORATION TEST ******\n";
    real satRatio = 0.9;
    real rho_v = rho_v_sat * satRatio;
    real rho_c = 0.001;
    real rho = rho_d + rho_v + rho_c;
    real temp = theta_d;
    real theta = physics.thetaFromTemp(rho, rho_v, rho_c, temp);
    real rho_theta = rho*theta;

    std::cout << "Pre-adjusted (rho, rho_v, rho_c, rho_theta): \n" << std::setprecision(20)
                                                                   << rho       << ", \n"
                                                                   << rho_v     << ", \n"
                                                                   << rho_c     << ", \n"
                                                                   << rho_theta << "\n\n";
    real mass1 = rho_d + rho_v + rho_c;

    physics.computeAdjustedState(rho_d, rho_v, rho_c, rho_theta);

    std::cout << "Post-adjusted (rho, rho_v, rho_c, rho_theta): \n" << std::setprecision(20)
                                                                    << rho       << ", \n"
                                                                    << rho_v     << ", \n"
                                                                    << rho_c     << ", \n"
                                                                    << rho_theta << "\n\n";
    real mass2 = rho_d + rho_v + rho_c;
    std::cout << "Relative Mass Change: " << std::setprecision(20) << (mass2-mass1)/mass1 << "\n";
  }

}


