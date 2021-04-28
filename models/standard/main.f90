! ****************
! main program: evolve a single model for several time-steps
program main
  use kemimo_commons
  use kemimo_rates
  use kemimo
  implicit none
  real*8::n(nmols), ntot, dt, t, Gnot, nout(nmols)
  real*8::Tgas, crflux, Av, tend
  integer :: unit, tunit, i, ii, iii, j, n_loop
  real*8 :: old_mask, mask


  print *, 'Running kemimo'
  ! ------------------------------------------------
  ! mandatory: load reaction names
  call kemimo_loadVerbatim()

  ! open file to write
  open(newunit=unit, file="output.dat", status="replace")

  ! Define physical initial conditions:
  Tgas = 1.0d1 ! gas temperature, K
  ntot = 2.0d4 ! gas density, cm-3
  Gnot = 1.0d0 ! radiation G0
  tend = 1d7*spy ! ending time, s
  Av = 5d0 ! visual extinction
  crflux = 1.3d-17 !CR ionization rate, 1/s

  ! ----------------------------------------------------------------
  ! Chemical initial conditions
  ! init chemical species, cm-3
  n(:) = 0d0
  ! Have to set dummy species to one.
  n(idx_dummy) = 1d0
  ! species
  n(idx_H_gas) = 5d-5*ntot
  n(idx_p_H2_gas) = ntot*4.5d-1
  n(idx_o_H2_gas) = ntot*4.5d-2
  n(idx_He_gas) = 9.75d-2*ntot
  n(idx_Cj_gas) = 7.86d-5*ntot
  n(idx_N_gas) = 2.47d-5*ntot
  n(idx_O_gas) = 1.8d-4*ntot
  n(idx_Sj_gas) = 9.14d-8*ntot
  n(idx_Sij_gas) = 9.74d-9*ntot
  n(idx_Naj_gas) = 2.25d-9*ntot
  n(idx_Fej_gas) = 2.74d-9*ntot
  n(idx_Mgj_gas) = 1.09d-9*ntot
  n(idx_Clj_gas) = 2.16d-10*ntot
  n(idx_Pj_gas) = 1d-9*ntot
  n(idx_E_gas) = n(idx_Cj_gas) + n(idx_Clj_gas) + n(idx_Sj_gas) + n(idx_Sij_gas) + n(idx_Naj_gas) + n(idx_Fej_gas) + n(idx_Mgj_gas) + n(idx_Pj_gas)

  n(idx_HD_gas) = 1.5d-5*ntot
  ! Set GRAIN:
  n(idx_GRAIN0_gas) = xdust * ntot
  
  ! -----------------------

  ! ---------------------------------------------------------------
  ! ##### Compute rates with the given parameters
  call kemimo_computeRates(n(:), ntot, Tgas, crflux, Av, Ghabing=Gnot)

  print *, ndns, xdust
  ! Loop init:
  t = 0d0 ! total time, s
  dt = 1d-3*spy ! initial timd-step, s
  n_loop = 0

  
  old_mask = n(idx_surface_mask) 
  ! ---------------------------------------------------------------
  ! ##### RUN loop:
  do
    ! run chemistry for one timestep
    call kemimo_dochem(n(:), dt)

    ! --------------------------------------------------
    ! update dt. Relying on the mask for dt updating, which is not strictly necessary.
    if (abs(mask-old_mask) > 1.3) then
      dt = dt * 0.9
    elseif (abs(mask-old_mask) > 0.7) then
      dt = dt
    else
      dt = dt * 1.1
    endif

    t = t + dt

    ! --------------------------------------------------
    ! write GAS+DUST species, 1:time, 2:Tgas, 3:ntot, 4:Av, 5->:species
    if (mod(n_loop, 2) == 0) then
      do i=1,nmols
        nout(i) = n(i)/ntot
      end do
      write(unit,'(59999E17.8e3)') t/spy, Tgas, ntot, Av, nout(:)

    endif
    
    old_mask = mask
    ! --------------------------------------------------
    ! Optional output once in a while:
    if (mod(n_loop, 5) == 0) then
      print *, 'mask: ', n(idx_surface_mask)
      call kemimo_printFluxes(n(:), 8, (/idx_H2O_0001/), .false.)
    endif
    ! increase step count
    n_loop = n_loop + 1
    ! break when time exceed ending time
    if(t>=tend) exit

  end do


  ! ----------------------------------------------------------------
  ! close unit
  close(unit)

end program main


