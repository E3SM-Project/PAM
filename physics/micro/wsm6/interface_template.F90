
module blah

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! rhoair0, rhowater, rhosnow, cliq, and cpv are all defined in module_model_constants.F90

! Set to 1 for hail-like graupel parameters, otherwise, use graupel-like graupel parameters
integer, parameter :: hail_like = 0

! allowed_to_read is not used by wsm6
logical, parameter ::allowed_to_read = .false. 

call wsm6init(rhoair0,rhowater,rhosnow,cliq,cpv, config_flags%hail_opt,allowed_to_read )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!======================================================================
! Grid structure in physics part of WRF
!----------------------------------------------------------------------
! The horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. All predicted
! variables are carried at half levels except w, which is at full
! levels. Some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------
! In WRF, kms (smallest number) is the bottom level and kme (largest
! number) is the top level.  In your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
!
!======================================================================
! Definitions
!-----------
! Rho_d      dry density (kg/m^3)
! Theta_m    moist potential temperature (K)
! Qv         water vapor    mixing ratio (kg/kg)
! Qc         cloud water    mixing ratio (kg/kg)
! Qr         rain water     mixing ratio (kg/kg)
! Qi         cloud ice      mixing ratio (kg/kg)
! Qs         snow           mixing ratio (kg/kg)
! Qg         graupel        mixing ratio (kg/kg)
! Qh         hail           mixing ratio (kg/kg)
! Qndrop     droplet number mixing ratio (#/kg)
! Qni        cloud ice number concentration (#/kg)
! Qns        snow      number concentration (#/kg)
! Qnr        rain      number concentration (#/kg)
! Qng        graupel   number concentration (#/kg)
! Qnh        hail      number concentration (#/kg)

! Qzr        rain             reflectivity (m6/kg)
! Qzi        ice              reflectivity (m6/kg)
! Qzs        snow             reflectivity (m6/kg)
! Qzg        graupel          reflectivity (m6/kg)
! Qzh        hail             reflectivity (m6/kg)

! Qvolg        graupel   particle volume (m3/kg)

!----------------------------------------------------------------------
!-- th        potential temperature    (K)
!-- moist_new     updated moisture array   (kg/kg)
!-- moist_old     Old moisture array       (kg/kg)
!-- rho           density of air           (kg/m^3)
!-- pi_phy        exner function           (dimensionless)
!-- p             pressure                 (Pa)
!-- RAINNC        grid scale precipitation (mm)
!-- RAINNCV       one time step grid scale precipitation (mm/step)
!-- SNOWNC        grid scale snow and ice (mm)
!-- SNOWNCV       one time step grid scale snow and ice (mm/step)
!-- GRAUPELNC     grid scale graupel (mm)
!-- GRAUPELNCV    one time step grid scale graupel (mm/step)
!-- HAILNC        grid scale hail (mm)
!-- HAILNCV       one time step grid scale hail (mm/step)
!-- SR            one time step mass ratio of snow to total precip
!-- z             Height above sea level   (m)
!-- dt            Time step              (s)
!-- G             acceleration due to gravity  (m/s^2)
!-- CP            heat capacity at constant pressure for dry air (J/kg/K)
!-- R_d           gas constant for dry air (J/kg/K)
!-- R_v           gas constant for water vapor (J/kg/K)
!-- XLS           latent heat of sublimation   (J/kg)
!-- XLV           latent heat of vaporization  (J/kg)
!-- XLF           latent heat of melting       (J/kg)
!-- rhowater      water density                      (kg/m^3)
!-- rhosnow       snow density               (kg/m^3)
!-- F_ICE_PHY     Fraction of ice.
!-- F_RAIN_PHY    Fraction of rain.
!-- F_RIMEF_PHY   Mass ratio of rimed ice (rime factor)
!-- t8w           temperature at layer interfaces
!-- cldfra, cldfra_old, current, previous cloud fraction
!-- exch_h        vertical diffusivity (m2/s)
!-- qlsink        Fractional cloud water sink (/s)
!-- precr         rain precipitation rate at all levels (kg/m2/s)
!-- preci         ice precipitation rate at all levels (kg/m2/s)
!-- precs         snow precipitation rate at all levels (kg/m2/s)
!-- precg         graupel precipitation rate at all levels (kg/m2/s)                             &
!-- P_QV          species index for water vapor
!-- P_QC          species index for cloud water
!-- P_QR          species index for rain water
!-- P_QI          species index for cloud ice
!-- P_QS          species index for snow
!-- P_QG          species index for graupel
!-- P_QH          species index for hail
!-- P_QNDROP      species index for cloud drop mixing ratio
!-- P_QNR         species index for rain number concentration,
!-- P_QNI         species index for cloud ice number concentration
!-- P_QNS         species index for snow number concentration,
!-- P_QNG         species index for graupel number concentration,
!-- P_QNH         species index for hail number concentration,
!-- P_QZR         species index for rain    reflectivity
!-- P_QZI         species index for ice     reflectivity
!-- P_QZS         species index for snow    reflectivity
!-- P_QZG         species index for graupel reflectivity
!-- P_QZH         species index for hail    reflectivity
!-- P_QVOLG       species index for graupel particle volume,
!-- id            grid id number
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- its           start index for i in tile
!-- ite           end index for i in tile
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- num_tiles     number of tiles
!-- diagflag      Logical to tell us when to produce diagnostics for history or restart


  SUBROUTINE wsm6(th, q, qc, qr, qi, qs, qg                        &
                 ,den, pii, p, delz                                &
                 ,delt,g, cpd, cpv, rd, rv, t0c                    &
                 ,ep1, ep2, qmin                                   &
                 ,XLS, XLV0, XLF0, den0, denr                      &
                 ,cliq,cice,psat                                   &
                 ,rain, rainncv                                    &
                 ,snow, snowncv                                    &
                 ,sr                                               &
                 ,graupel, graupelncv                              &
                 ,ids,ide, jds,jde, kds,kde                        &
                 ,ims,ime, jms,jme, kms,kme                        &
                 ,its,ite, jts,jte, kts,kte                        &
                                                                   )
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
  INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(INOUT) ::                                          &
                                                             th,  &
                                                              q,  &
                                                              qc, &
                                                              qi, &
                                                              qr, &
                                                              qs, &
                                                              qg
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(IN   ) ::                                          &
                                                             den, &
                                                             pii, &
                                                               p, &
                                                            delz
  REAL, INTENT(IN   ) ::                                    delt, &
                                                               g, &
                                                              rd, &
                                                              rv, &
                                                             t0c, &
                                                            den0, &
                                                             cpd, &
                                                             cpv, &
                                                             ep1, &
                                                             ep2, &
                                                            qmin, &
                                                             XLS, &
                                                            XLV0, &
                                                            XLF0, &
                                                            cliq, &
                                                            cice, &
                                                            psat, &
                                                            denr
  REAL, DIMENSION( ims:ime , jms:jme ),                           &
        INTENT(INOUT) ::                                    rain, &
                                                         rainncv, &
                                                              sr
!+---+-----------------------------------------------------------------+

  REAL, DIMENSION( ims:ime , jms:jme ), OPTIONAL,                 &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv
  REAL, DIMENSION( ims:ime , jms:jme ), OPTIONAL,                 &
        INTENT(INOUT) ::                                 graupel, &
                                                      graupelncv



call wsm6(TH=th                     ! INOUT          : ikj : Dry potential temperature
         ,Q=qv_curr                 ! INOUT          : ikj : Water vapor mixing ratio (kg/kg) TODO: DRY or MOIST?
         ,QC=qc_curr                ! INOUT          : ikj : Cloud water mixing ratio (kg/kg) TODO: DRY or MOIST?
         ,QR=qr_curr                ! INOUT          : ikj : Rain water mixing ratio (kg/kg) TODO: DRY or MOIST?
         ,QI=qi_curr                ! INOUT          : ikj : Cloud ice mixing ratio (kg/kg) TODO: DRY or MOIST?
         ,QS=qs_curr                ! INOUT          : ikj : Snow mixing ratio (kg/kg) TODO: DRY or MOIST?
         ,QG=qg_curr                ! INOUT          : ikj : Graupel mixing ratio (kg/kg) TODO: DRY or MOIST?
         ,DEN=rho                   ! IN             : ikj : Density (kg/m^3) TODO: DRY or FULL?
         ,PII=pi_phy                ! IN             : ikj : Dry Exner function (dimensionless)
         ,P=p                       ! IN             : ikj : Pressure (Pa) TODO: DRY or FULL?
         ,DELZ=dz8w                 ! IN             : ikj : Layer thickness in (m)
         ,DELT=dt                   ! IN             :     : Time step
         ,G=g                       ! IN             :     : Gravity
         ,CPD=cp                    ! IN             :     : Heat capacity at constant pressure for dry air (J/kg/K)
         ,CPV=cpv                   ! IN             :     : Heat capacity at constant pressure for water vapor (J/kg/K)
         ,RD=r_d                    ! IN             :     : Gas constant for dry air (J/kg/K)
         ,RV=r_v                    ! IN             :     : Gas constant for water vapor (J/kg/K)
         ,T0C=svpt0                 ! IN             :     : Constant for saturation vapor pressure calculation (K)
         ,EP1=ep_1                  ! IN             :     : Constant for virtual temperature (r_v/r_d - 1) (dimensionless)
         ,EP2=ep_2                  ! IN             :     : Constant for specific humidity calculation (dimensionless)
         ,QMIN=epsilon              ! IN             :     : Machine precision (e.g., 1.e-15)
         ,XLS=xls                   ! IN             :     : Latent heat of sublimation   (J/kg)
         ,XLV0=xlv                  ! IN             :     : Latent heat of vaporization  (J/kg)
         ,XLF0=xlf                  ! IN             :     : Latent heat of melting       (J/kg)
         ,DEN0=rhoair0              ! IN             :     : Density of dry air at 0^oC and 1000mb pressure (kg m^-3)
         ,DENR=rhowater             ! IN             :     : Density of liquid water at 0^oC (kg m^-3)
         ,CLIQ=cliq                 ! IN             :     : Specific heat of liquid water at 0^oC
         ,CICE=cice                 ! IN             :     : Specific heat of ice at 0^oC
         ,PSAT=psat                 ! IN             :     : 610.78 (No documentation)
         ,RAIN=rainnc               ! INOUT          : ij  : Grid scale precipitation (mm)
         ,RAINNCV=rainncv           ! INOUT          : ij  : One time step grid scale precipitation (mm/step)
         ,SNOW=snownc               ! INOUT, OPTIONAL: ij  : Grid scale snow and ice (mm)
         ,SNOWNCV=snowncv           ! INOUT, OPTIONAL: ij  : One time step grid scale snow and ice (mm/step)
         ,SR=sr                     ! INOUT          : ij  : One time step mass ratio of snow to total precip
         ,GRAUPEL=graupelnc         ! INOUT, OPTIONAL: ij  : Grid scale graupel (mm)
         ,GRAUPELNCV=graupelncv     ! INOUT, OPTIONAL: ij  : One time step grid scale graupel (mm/step)
         ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde ! IN   :     : Starting and ending indices for the "domain" (SEEMS TO NOT BE USED)
         ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme ! IN   :     : Starting and ending indices for the "memory"
         ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte ! IN   :     : Starting and ending indices for the "tile"   (CAN PROBABLY BE SET TO MEMORY BOUNDS)





end module blah
