module kemimo
  use iso_c_binding
contains
  !*****************
  !print fluxes to stdout
  subroutine kemimo_printFluxes(n, nflux, idxList, onlyDust)
    use kemimo_commons
    use kemimo_flux
    implicit none
    real*8,intent(in)::n(nmols)
    integer,intent(in)::nflux,idxList(:)
    logical,optional,intent(in)::onlyDust
    logical::onlyDustValue

    onlyDustValue = .false.
    if(present(onlyDust)) then
       onlyDustValue = onlyDust
    end if

    call printFluxes(n(:), nflux, idxList(:), onlyDustValue)

  end subroutine kemimo_printFluxes


  !*********************
  !save fluxes to file
  subroutine kemimo_dumpFluxes(n, unit, xvar) bind(C,name='kemimo_dumpFluxes')
    use kemimo_commons
    use kemimo_flux
    implicit none
    real(kind=c_double),intent(in)::n(nmols), xvar
    integer(kind=c_int),intent(in)::unit

    call dumpFluxes(n(:), unit, xvar)

  end subroutine kemimo_dumpFluxes

  !*****************
  !return index of ascending sorted array v
  ! bubble sorting, faster algorithms are welcomed
  function kemimo_sortedIndexs(v) result(idx)
    use kemimo_commons
    use kemimo_flux
    implicit none
    real*8,intent(in)::v(:)
    integer::idx(size(v))

    idx(:) = sortedIndexs(v(:))

  end function kemimo_sortedIndexs

  !***************
  !load verbatim reaction
  subroutine kemimo_loadverbatim() bind(C, name='kemimo_loadVerbatim')
    use kemimo_flux

    call loadVerbatim()

  end subroutine kemimo_loadverbatim

  !*******************
  !compute rates and store into commons kall(:)
  subroutine kemimo_computeRates(n, ngas, variable_Tgas, variable_crflux, &
              variable_Av, Td, size_min, size_max, Ghabing, onlyGas, onlyDust)
    use kemimo_commons
    use kemimo_rates
    use kemimo_reactionarray
    implicit none
    real*8,intent(in)::n(nmols)
    real*8,intent(in)::ngas, variable_Tgas
    real*8,intent(in)::variable_crflux, variable_Av
    real*8,intent(in),optional::size_min, size_max, Ghabing, Td
    logical,optional,intent(in)::onlyGas, onlyDust
    real*8::amin, amax, Gnot, Tdust

    call init_reactionarray()

    !default distribution min/max size, cm
    !amin = 5d-7
    !amax = 2.5d-5
    amin = 1d-5
    amax = 1d-5
    Gnot = 1d0 !Draine flux in Habing units
    Tdust = variable_Tgas

    !replace optinal values if present
    if(present(size_min)) amin = size_min
    if(present(size_max)) amax = size_max
    if(present(Ghabing)) Gnot = Ghabing
    if(present(Td)) Tdust = Td

    call computeRates(n, ngas, variable_Tgas, variable_crflux, variable_Av, &
        Tdust, amin, amax, Gnot)

    if(present(onlyGas) .and. onlyGas .eqv. .true.) then
      call switchOffDustRates()
    end if

    if(present(onlyDust) .and. onlyDust .eqv. .true.) then
      call switchOffGasRates()
    end if



  end subroutine kemimo_computeRates
    
  !************************
  ! a wrapper is used here because the original has optional arguments
  subroutine kemimo_computeRates_c(n, ngas, variable_Tgas, variable_crflux, &
       variable_Av) bind(C,name='kemimo_computeRates')
    use kemimo_commons
    implicit none
    real(kind=c_double), intent(inout) :: n(nmols)
    real(kind=c_double), intent(inout)::ngas, variable_Tgas
    real(kind=c_double), intent(inout)::variable_crflux, variable_Av

    call kemimo_computeRates(n, ngas, variable_Tgas, variable_crflux, variable_Av)

  end subroutine kemimo_computeRates_c

  !************************
  !evolve chemistry for a time-step dt (s)
  ! n(:) are species number densities
  subroutine kemimo_dochem(n, dt) bind(C, name='kemimo_dochem')
    use kemimo_commons
    use kemimo_ode
    use kemimo_dust_rates
    implicit none
    real(kind=c_double),intent(inout)::n(nmols)
    real(kind=c_double),intent(inout)::dt
    real(kind=c_double) :: OPR, H2_ice, delta_H2_ortho_ice, delta_H2_para_ice, R
    real(kind=c_double) :: Nsurface, Nmantle
    integer :: i
    integer,parameter:: offset = mantle_start - surface_start
    ! ----------------------------------------------------------------
    ! Update OPR
    OPR = n(idx_o_H2_gas)/(n(idx_p_H2_gas)+ n(idx_o_H2_gas))

    ! --------------------------------------------------
    ! Update H2_ice for current ndns, H2_coverage (calculated in computeRates call!)
    H2_ice = H2_coverage * ndns * layerThickness

    ! Change:
    delta_H2_para_ice = H2_ice*(1d0 - OPR) - n(idx_p_H2_0001)
    delta_H2_ortho_ice = H2_ice*OPR - n(idx_o_H2_0001)
    
    ! Update surface mask based on H2 update:
    R = (delta_H2_para_ice + delta_H2_ortho_ice)
    n(idx_surface_mask) = n(idx_surface_mask) + R/ndns
    
    n(idx_p_H2_0001) = H2_ice*(1d0 - OPR)
    n(idx_o_H2_0001) = H2_ice*OPR

    call dochem(n(:), dt)

  end subroutine kemimo_dochem


    !************************
  !get H2 coverage:
  ! n(:) are species number densities
  function kemimo_theta_H2(Tgas, n_H2_gas) bind(C, name='kemimo_theta_H2')
    use kemimo_commons
    implicit none
    real(kind=c_double), intent(in):: Tgas, n_H2_gas
    real(kind=c_double) :: kemimo_theta_H2

    kemimo_theta_H2 = theta_H2(Tgas, n_H2_gas)

  end function kemimo_theta_H2

  !*****************
  !differential equations, only formation
  function kemimo_fexForm(n) result(dnf)
    use kemimo_commons
    use kemimo_ode
    implicit none
    real*8,intent(in)::n(nmols)
    real*8::dnf(nmols)

    dnf(:) = fexForm(n(:))

  end function kemimo_fexForm

end module kemimo
