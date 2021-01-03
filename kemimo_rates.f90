module kemimo_rates
contains

  !*******************
  !compute rates and store into commons kall(:)
  subroutine computeRates(n, ngas, variable_Tgas, variable_crflux, variable_Av, &
       Tdust, a_grain, Gnot)
    use kemimo_commons
    use kemimo_gas_rates
    use kemimo_dust_rates
    implicit none
    real*8,intent(in)::n(nmols)
    real*8,intent(in)::ngas, variable_Tgas
    real*8,intent(in)::variable_crflux, variable_Av
    real*8,intent(in)::a_grain, Gnot, Tdust
    real*8 :: a, x_H2, x_HD, omega, N_H2_crit
    character(len=255) :: fname_CO, fname_N2

    ! -------------------------------------------------
    ! This is the H2 coverage factor (theta_H2) used when H2 is fixed.
    H2_coverage = theta_H2(variable_Tgas, ngas*5d-1)

    ! -------------------------------------------------
    ! Estimate N_H2, N_CO etc from Av:
    N_H2 = 0.5 * 1.59d21 * variable_Av
    N_CO = 1d-4 * N_H2
    N_HD = N_H2 * 1.5d-5
    N_HI = N_H2 * 5d-5
    N_N2 = N_H2 * 1d-5


    ! --------------------------------------------------------
    ! self-shielding functions:
    ! ---------------------
    ! CO
    fname_CO = "visser_CO_shield.dat"
    ss_CO = calc_ss(N_H2, N_CO, fname_CO, 47, 42, 8)
    ! ---------------------
    ! H2:
    omega = 1.3d-2 * (1d0 + (variable_Tgas/2.7d3)**(1.3d0))**(1d0/1.3d0) * exp(-(variable_Tgas/3.9d3)**1.46d1)
    N_H2_crit = 1.3d14 * (1d0 + (variable_Tgas/6d2)**0.8)
    a = 1.4d0
    x_H2 = N_H2 / N_H2_crit
    ! Grassi+2020:
    ss_H2 = (1d0 - omega)/(1d0 + x_H2)**a * exp(-5.0d-7 * (1d0 + x_H2)) + omega/sqrt(1d0 + x_H2) * exp(-8.5d-4 * sqrt(1d0 + x_H2))
    ! Kamp & Bertoldi:
    !ss_H2 = 9.65d-1/(1d0 + x_H2)**2.0 + 3.5d-2/sqrt(1d0 + x_H2) * exp(-8.5d-4 * sqrt(1d0 + x_H2))
    ! HD shielding of HI:
    ss_H2 = ss_H2 * 1d0/(1d0 + N_HI/2.85d23)**(1.62) * exp(-1.49d-1 * N_HI/2.85d23)
    ! ---------------------
    ! HD
    x_HD = N_HD / N_H2_crit
    ss_HD = (1d0 - omega)/(1d0 + x_HD)**a * exp(-5d-7 * (1d0 + x_HD)) + omega/sqrt(1d0 + x_HD) * exp(-8.5d-4 * sqrt(1d0 + x_HD))
    ! H2 and HI shielding of HD (Wolcot-Green + Haiman):
    ss_HD = ss_HD * 1d0/(1d0 + N_H2/2.34d19)**(2.38d-1) * exp(-5.2d-3 * N_H2/2.34d19)
    ss_HD = ss_HD * 1d0/(1d0 + N_HI/2.85d23)**(1.62) * exp(-1.49d-1 * N_HI/2.85d23)
    ! ---------------------
    ! N2
    fname_N2 = "N2_shield_50K.dat"
    ss_N2 = calc_ss(N_H2, N_N2, fname_N2, 46, 46, 7)


    ! -------------------------------------------------
    kall(:) = 0d0
    call computeGasRates(n, ngas, variable_Tgas, variable_crflux, variable_Av, &
       Tdust, a_grain, Gnot)

    call computeDustRates(n, ngas, variable_Tgas, variable_crflux, variable_Av, &
       Tdust, a_grain, Gnot)

    !small bias to favour solver convergence
    kall(:) = kall(:) + 1d-40

    !this part is intended to test reduced networks
    !ktmp(1:nreadust) = 0d0

    !kall(1:nreadust) = ktmp(1:nreadust)

  end subroutine computeRates

  !*************
  !switch off gas->gas rates
  subroutine switchOffGasRates()
    use kemimo_commons
    implicit none

    kall(nreadust+1:nrea-1) = 0d0

  end subroutine switchOffGasRates

  !*************
  !switch off rates involving dust
  ! i.e. gas->dust, dust->dust, dust->gas
  subroutine switchOffDustRates()
    use kemimo_commons
    implicit none

    kall(1:nreadust) = 0d0

  end subroutine switchOffDustRates

  !********************
  !dump rate values to a file
  subroutine dumpRates(fname)
    use kemimo_commons
    implicit none
    character(len=*),intent(in)::fname
    integer::unit,i

    open(newunit=unit, file=trim(fname), status="replace")
    do i=1,nrea
       !index, rate coefficient, verbatim
       write(unit,'(I5,E17.8e3,2a)') i, kall(i), "  ", trim(verbatim(i))
    end do
    close(unit)

  end subroutine dumpRates

  !************************
  !1d liner fit
  ! x: point where to evaluate f(x)
  ! ndata: number of points in data set
  ! xmin: x minimum value
  ! invdx: inverse of x spacing
  ! dx: x spacing
  ! ydata(ndata): set of y values
  function fit1d(x, ndata, xmin, invdx, dx, ydata) result(f)
    implicit none
    integer,intent(in)::ndata
    real*8,intent(in)::xmin, invdx, ydata(:), x, dx
    integer::idx
    real*8::x0, f, f0, f1

    !index of data point corresponding to x
    idx = int((x-xmin)*invdx*(ndata-1))

    !function at idx, idx+1
    f0 = ydata(idx)
    f1 = ydata(idx+1)

    !get x variable corresponding to idx
    x0 = xmin + dx*idx

    !compute linearly interpolated function
    f = (x-x0)*invdx*(f1-f0)+f0

  end function fit1d

end module kemimo_rates
