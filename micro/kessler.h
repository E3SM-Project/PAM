//-----------------------------------------------------------------------
//
//  Version:  2.0
//
//  Date:  January 22nd, 2015
//
//  Change log:
//  v2 - Added sub-cycling of rain sedimentation so as not to violate
//       CFL condition.
//
//  The KESSLER subroutine implements the Kessler (1969) microphysics
//  parameterization as described by Soong and Ogura (1973) and Klemp
//  and Wilhelmson (1978, KW). KESSLER is called at the end of each
//  time step and makes the final adjustments to the potential
//  temperature and moisture variables due to microphysical processes
//  occurring during that time step. KESSLER is called once for each
//  vertical column of grid cells. Increments are computed and added
//  into the respective variables. The Kessler scheme contains three
//  moisture categories: water vapor, cloud water (liquid water that
//  moves with the flow), and rain water (liquid water that falls
//  relative to the surrounding air). There  are no ice categories.
//  Variables in the column are ordered from the surface to the top.
//
//  SUBROUTINE KESSLER(theta, qv, qc, qr, rho, pk, dt, z, nz, rainnc)
//
//  Input variables:
//     theta  - potential temperature (K)
//     qv     - water vapor mixing ratio (gm/gm)
//     qc     - cloud water mixing ratio (gm/gm)
//     qr     - rain  water mixing ratio (gm/gm)
//     rho    - dry air density (not mean state as in KW) (kg/m^3)
//     pk     - Exner function  (not mean state as in KW) (p/p0)**(R/cp)
//     dt     - time step (s)
//     z      - heights of thermodynamic levels in the grid column (m)
//     nz     - number of thermodynamic levels in the column
//     precl  - Precipitation rate (m_water/s)
//
// Output variables:
//     Increments are added into t, qv, qc, qr, and rainnc which are
//     returned to the routine from which KESSLER was called. To obtain
//     the total precip qt, after calling the KESSLER routine, compute:
//
//       qt = sum over surface grid cells of (rainnc * cell area)  (kg)
//       [here, the conversion to kg uses (10^3 kg/m^3)*(10^-3 m/mm) = 1]
//
//
//  Authors: Paul Ullrich
//           University of California, Davis
//           Email: paullrich@ucdavis.edu
//
//           Based on a code by Joseph Klemp
//           (National Center for Atmospheric Research)
//
//  Reference:
//
//    Klemp, J. B., W. C. Skamarock, W. C., and S.-H. Park, 2015:
//    Idealized Global Nonhydrostatic Atmospheric Test Cases on a Reduced
//    Radius Sphere. Journal of Advances in Modeling Earth Systems. 
//    doi:10.1002/2015MS000435
//
//=======================================================================

void kessler(real1d &theta, real1d &qv, real1d &qc, real1d &qr, real1d const &rho,
             real &precl, real1d const &z, real1d const &pk, real dt, int nz) {
  real1d r    ("r"    ,nz);
  real1d rhalf("rhalf",nz);
  real1d velqr("velqr",nz);
  real1d sed  ("sed"  ,nz);
  real1d pc   ("pc"   ,nz);

  real f2x = 17.27_fp
  real f5 = 237.3_fp * f2x * 2500000._fp / 1003._fp
  real xk = .2875_fp      //  kappa (r/cp)
  real psl    = 1000._fp  //  pressure at sea level (mb)
  real rhoqr  = 1000._fp  //  density of liquid water (kg/m^3)

  for (int k=0; k < nz; k++) {
    r(k)     = 0.001_fp * rho(k);
    rhalf(k) = sqrt( rho(0) / rho(k) );
    pc(k)    = 3.8_fp / ( pow( pk(k) , (1._fp/xk) ) * psl );

    // Liquid water terminal velocity (m/s) following KW eq. 2.15
    velqr(k)  = 36.34_fp * pow( qr(k)*r(k) , 0.1364_fp ) * rhalf(k);
  }

  // Maximum time step size in accordance with CFL condition
  if (dt <= 0) { stoprun("kessler.f90 called with nonpositive dt"); }

  real dt_max = dt;
  for (k=0; k < nz-1; k++) {
    if (velqr(k) != 0) {
      dt_max = min(dt_max, 0.8_fp*(z(k+1)-z(k))/velqr(k));
    }
  }

  // Number of subcycles
  int rainsplit = ceil(dt / dt_max);
  real dt0 = dt / static_cast<real>(rainsplit);

  // Subcycle through rain process
  real precl = 0;

  for (int nt=0; nt < rainsplit; nt++) {
    // Precipitation rate (m/s)
    precl = precl + rho(0) * qr(0) * velqr(0) / rhoqr;

    // Sedimentation term using upstream differencing
    for (int k=0; k < nz-1; k++) {
      sed(k) = dt0*(r(k+1)*qr(k+1)*velqr(k+1)-r(k)*qr(k)*velqr(k))/(r(k)*(z(k+1)-z(k)));
    }
    sed(nz-1) = -dt0*qr(nz-1)*velqr(nz-1)/(.5*(z(nz-1)-z(nz-2)));

    // Adjustment terms
    for (int k=0; k < nz; k++) {
      // Autoconversion and accretion rates following KW eq. 2.13a,b
      real qrprod = qc(k) - ( qc(k)-dt0*max(0.001_fp*(qc(k)-0.001_fp) , 0._fp) ) /
                            ( 1 + dt0 * 2.2_fp * pow( qr(k) , 0.875_fp ) );
      qc(k) = max( qc(k)-qrprod , 0._fp );
      qr(k) = max( qr(k)+qrprod+sed(k) , 0._fp );

      // Saturation vapor mixing ratio (gm/gm) following KW eq. 2.11
      real qvs = pc(k)*exp(f2x*(pk(k)*theta(k)-273._fp) / (pk(k)*theta(k)- 36._fp));
      real tmp = pk(k)*theta(k)-36._fp;
      real prod = (qv(k)-qvs)/(1._fp+qvs*f5/(tmp*tmp));

      // Evaporation rate following KW eq. 2.14a,b
      real tmp1 = dt0*( ( ( 1.6_fp + 124.9_fp * pow( r(k)*qr(k) , 0.2046_fp ) ) *
                             pow( r(k)*qr(k) , 0.525_fp ) ) /
                           ( 2550000._fp*pc(k) / (3.8_fp *qvs)+540000._fp) ) * 
                         ( max(qvs-qv(k),0._fp) / (r(k)*qvs) );
      real tmp2 = max( -prod-qc(k) , 0._fp );
      real tmp3 = qr(k);
      real ern = min( tmp1 , min( tmp2 , tmp3 ) );

      // Saturation adjustment following KW eq. 3.10
      theta(k)= theta(k) + 2500000._fp / (1003._fp*pk(k)) * ( max( prod , -qc(k) )-ern );
      qv(k) = max( qv(k) - max( prod , -qc(k) ) + ern , 0. );
      qc(k) = qc(k) + max( prod , -qc(k) );
      qr(k) = qr(k) - ern;
    }

    // Recalculate liquid water terminal velocity
    if (nt != rainsplit-1) {
      for (int k=0; k < nz; k++) {
        velqr(k)  = 36.34_fp * pow( qr(k)*r(k) , 0.1364_fp ) * rhalf(k);
      }
    }
  }

  precl = precl / static_cast<real>(rainsplit);
}



