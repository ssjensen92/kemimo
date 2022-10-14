module kemimo_dust_rates
    real*8:: kph_factor, k_H2O_ph_0
    real*8:: Nsites, Ffuva
    real*8 :: Ffuva_CO, Ffuva_H2, Ffuva_HD, Ffuva_N2
contains

  !*******************
  !compute rates and store into commons kall(:)
  subroutine computeDustRates(n, ngas, variable_Tgas, variable_crflux, variable_Av, &
       Tdust, a_grain, Gnot)
    use kemimo_commons
    use kemimo_sticking
    use kemimo_swappingrates
    implicit none
    real*8,intent(in):: n(nmols)
    real*8,intent(in):: ngas, variable_Tgas
    real*8,intent(in):: variable_crflux, variable_Av
    real*8,intent(in):: a_grain, Gnot, Tdust
    real*8:: prefreezeout, prefreezeout_H, indns, invTd
    real*8:: ktmp(nrea)
    real*8:: p3, p4, pexp
    real*8:: n_ice_tot
    real*8, parameter :: beta = 2.5
    real*8, parameter :: T0 = 87.0 !87.0 np-ASW ice  - 56 silicate
    real*8, parameter :: S0 = 0.76 !0.76 np-ASW ice  - 0.95 silicate
    real*8, parameter :: gamma_H2O = 2.2d0
    real*8, parameter :: gamma_CO2 = 2.03d0
    real*8, parameter :: P_H2O_isrf = 5.4d-3
    real*8, parameter :: P_H2O_cr = 4.7d-3
    integer::i


    pexp = -3.5 !size distribution exponent
    p3 = pexp + 3d0
    p4 = pexp + 4d0
    
    invTd = 1d0/variable_Tgas
    ! -------------------------------------------------
    ! Update swapping rates:
    call computeSwapRates(invTd)

    ! -------------------------------------------------
    ! Update sticking rates:
    call computeSticking(variable_Tgas)
    ! -------------------------------------------------    
    !freezeout prefactor
    prefreezeout =  sqrt(variable_Tgas) * a_grain**2 * pi * xdust * ngas

    !number of binding sites per unit of volume
    Nsites = 4.0 * pi * a_grain**2 * site_density ! Number of sites pr dust grain
    ndns = xdust * ngas * Nsites

    !inverse of ndns
    indns = 1d0/ndns

    ! -------------------------------------------------
    !FUV photons per bound molecule per unit of time, 1/s
    kph_factor = a_grain**2.0 * pi * xdust * ngas
    ! Adjust/scale F_cr (CR-UV flux) based on the current CR flux
    F_cr = variable_crflux * F_cr_not / 1.3d-17
    ! Av attenuation from Hollenbach+2009
    Ffuva = kph_factor * (F_cr + Gnot*Fnot*exp(-1.8d0*variable_Av))

    ! Total gas-phase unattenuated dissociation rate from KIDA
    k_H2O_ph_0 = 8.01d-10 
    !k_H2O_ph_0 = 8.310d-10 

    ! Specials to account for self-shielding
    Ffuva_CO = ss_CO * Ffuva
    Ffuva_H2 = ss_H2 * Ffuva
    Ffuva_HD = ss_HD * Ffuva
    Ffuva_N2 = ss_N2 * Ffuva

    !!BEGIN_RATES
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-22 21:21:48
    ! CHANGESET: xxxxxxx
    ! BY: unknown@unknown


    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_RATES


  end subroutine computeDustRates


end module kemimo_dust_rates
