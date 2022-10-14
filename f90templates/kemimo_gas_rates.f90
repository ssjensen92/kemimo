module kemimo_gas_rates
contains

  !*******************
  !compute rates and store into commons kall(:)
  subroutine computeGasRates(n, ngas, variable_Tgas, variable_crflux, variable_Av, &
       Tdust, a_grain, Gnot)
    use kemimo_commons
    implicit none
    real*8,intent(in):: n(nmols)
    real*8,intent(in):: ngas, variable_Tgas
    real*8,intent(in):: variable_crflux, variable_Av
    real*8,intent(in):: Gnot, Tdust, a_grain
    real*8:: invTd
    real*8:: ktmp(nrea)
    real*8:: k_H2O_ph_0
    real*8, parameter :: beta = 2.5
    real*8, parameter :: T0 = 87.0 !87.0 np-ASW ice  - 56 silicate
    real*8, parameter :: S0 = 0.76 !0.76 np-ASW ice  - 0.95 silicate
    real*8, parameter :: gamma_H2O = 2.2
    real*8, parameter :: gamma_CO2 = 3.0
    real*8, parameter :: P_H2O_isrf = 5.4d-3
    real*8, parameter :: P_H2O_cr = 4.7d-3
    real*8 :: cross_section_rate, tau, Jtilde_ion, Jtilde_neutral, kgr_ion, kgr_neutral
    integer::i
    character(len=255) :: fname_CO, fname_N2

    H2_coverage = theta_H2(variable_Tgas, ngas*5d-1)


    ! -------------------------------------------------
    ! calculate grain-ion recombination rates (Draine & Sutil 1987):
    cross_section_rate = sqrt(8d0 * kb * variable_Tgas / (pi * pmass)) * a_grain**2 * pi
    ! 
    tau = a_grain*kb*variable_Tgas/qe**2
    Jtilde_ion = (1d0 + 1d0/tau)*(1d0 + sqrt(2d0/(2d0+tau)))
    Jtilde_neutral = 1d0 + sqrt(pi/(2d0*tau))
    kgr_ion = cross_section_rate * Jtilde_ion
    kgr_neutral = cross_section_rate * Jtilde_neutral

    ! Other constants and stuff:
    invTd = 1d0/variable_Tgas

    k_H2O_ph_0 = 7.700d-10   ! Total gas-phase unattenuated dissociation rate from KIDA (Heays+2017)


    !!BEGIN_RATES
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-22 21:21:49
    ! CHANGESET: xxxxxxx
    ! BY: unknown@unknown

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_RATES


  end subroutine computeGasRates


end module kemimo_gas_rates
