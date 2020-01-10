
module atmosphere_mod
	implicit none

contains

  subroutine atmosphere(alt, sigma, delta, theta)
    !   -------------------------------------------------------------------------
    ! purpose - compute the properties of the 1976 standard atmosphere to 86 km.
    ! author - ralph carmichael, public domain aeronautical software
    ! note - if alt > 86, the values returned will not be correct, but they will
    !   not be too far removed from the correct values for density.
    !   the reference document does not use the terms pressure and temperature
    !   above 86 km.
    use params
    implicit none
    !============================================================================
    !     a r g u m e n t s                                                     |
    !============================================================================
    real(crm_rknd),intent(in)::  alt  ! geometric altitude, km.
    real(crm_rknd),intent(out):: sigma! density/sea-level standard density
    real(crm_rknd),intent(out):: delta! pressure/sea-level standard pressure
    real(crm_rknd),intent(out):: theta! temperature/sea-level standard temperature
    !============================================================================
    !     l o c a l   c o n s t a n t s                                         |
    !============================================================================
    real(crm_rknd),parameter:: rearth = 6369.0 ! radius of the earth (km)
    real(crm_rknd),parameter:: gmr = 34.163195 ! gas constant
    integer(crm_iknd),parameter:: ntab=8! number of entries in the defining tables
    !============================================================================
    !     l o c a l   v a r i a b l e s                                         |
    !============================================================================
    integer(crm_iknd):: i,j,k    ! counters
    real(crm_rknd):: h           ! geopotential altitude (km)
    real(crm_rknd):: tgrad, tbase! temperature gradient and base temp of this layer
    real(crm_rknd):: tlocal      ! local temperature
    real(crm_rknd):: deltah      ! height above base of this layer
    !============================================================================
    !     l o c a l   a r r a y s   ( 1 9 7 6   s t d.  a t m o s p h e r e )   |
    !============================================================================
    real(crm_rknd),dimension(ntab),parameter:: htab= (/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0,84.852/)
    real(crm_rknd),dimension(ntab),parameter:: ttab=  (/288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946/)
    real(crm_rknd),dimension(ntab),parameter:: ptab=  (/1.0, 2.233611e-1, &
    5.403295e-2, 8.5666784e-3, 1.0945601e-3, 6.6063531e-4, 3.9046834e-5, 3.68501e-6/)
    real(crm_rknd),dimension(ntab),parameter:: gtab=  (/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0/)
    !----------------------------------------------------------------------------
    h=alt*rearth/(alt+rearth)! convert geometric to geopotential altitude

    i=1
    j=ntab                     ! setting up for=binary search
    do
      k=(i+j)/2
      if (h < htab(k)) then
        j=k
      else
        i=k
      end if
      if (j <= i+1) exit
    end do

    tgrad=gtab(i)              ! i will be in 1...ntab-1
    tbase=ttab(i)
    deltah=h-htab(i)
    tlocal=tbase+tgrad*deltah
    theta=tlocal/ttab(1)            ! temperature ratio

    if (tgrad == 0.0) then        ! pressure ratio
      delta=ptab(i)*exp(-gmr*deltah/tbase)
    else
      delta=ptab(i)*(tbase/tlocal)**(gmr/tgrad)
    end if

    sigma=delta/theta            ! density ratio
    return
  end subroutine atmosphere

end module atmosphere_mod
