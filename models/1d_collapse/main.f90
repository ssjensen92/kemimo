! ****************
! main program: evolve a single model for several time-steps
program loop
  use kemimo_commons
  use kemimo_rates
  use kemimo
  implicit none
  real*8::n(nmols), ntot, dt, t, Gnot, nout(nmols), ntot_current
  real*8::Tgas, crflux, Av, tend
  integer :: unit, tunit, i, ii, iii, j, n_loop
  integer, parameter :: model_length = 139 ! 254 ! number of lines in physical model
  real*8 :: rhod, dust_mass
  real*8 :: ndust_factor
  real*8 :: dust_sites
  real*8, parameter :: a = 1d-5 ! grain radius
  real*8 :: old_mask, mask
  real*8 :: time_arr(model_length), ntot_arr(model_length), Av_arr(model_length), temperature_arr(model_length)
  real*8 :: Av_current, Tgas_current, ntot_ratio

  print *, 'Running randice'
  ! ------------------------------------------------
  ! Read tracer:
  !open(newunit=tunit, file="physical_model_NAUTILUS.dat", status="old")
  open(newunit=tunit, file="smoothed_physical_model.dat", status="old")
  ! ! Skip header:
  ! do i=1,70
  !   read(tunit, *)
  ! enddo
  do i=1, model_length
     read(tunit, *) time_arr(i), Av_arr(i), ntot_arr(i), temperature_arr(i)
  end do
  close(tunit)
  ! ------------------------------------------------
  print *, 'Loaded physical model'
  ! mandatory: load reaction names
  call kemimo_loadVerbatim()

  ! open file to write
  open(newunit=unit, file="output_collapse.dat", status="replace")

  rhod = pmass*mu*d2g !dust mass density, g/cm3
  dust_mass = rho0 * pi43 * a**3.0d0
  ndust_factor = rhod/dust_mass
  ! Physical initial conditions:
  Tgas = 1.0d1 ! gas temperature, K
  ntot = 2.0d4 ! gas density, cm-3
  Gnot = 1.0d0 ! radiation G0
  tend = 1d7*spy ! ending time, s
  Av = 4.5d0 ! visual extinction
  crflux = 1.3d-17 !CR ionization rate, 1/s

  ! ----------------------------------------------------------------
  ! Chemical initial conditions
  ! init chemical species, cm-3
  n(:) = 0d0
  ! Have to set dummy to one!
  n(idx_dummy) = 1d0
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

  ! Set GRAIN:
  n(idx_GRAIN0_gas) = ndust_factor * ntot

  ! Add D?
  n(idx_HD_gas) = 1.5d-5*ntot

  ! -----------------------
  ! compute rates with the given parameters
  ntot = ntot_arr(1)
  Tgas = temperature_arr(1)
  Av =  Av_arr(1)
  ntot_current = ntot
  Tgas_current =  Tgas
  Av_current =  Av

  print *, 'Calculating new rates for: ', ntot_current, Tgas_current, Av_current
  call kemimo_computeRates(n(:), ntot_current, Tgas, crflux, Av, Ghabing=Gnot)
  
  ! -----------------------------------------
  n_loop = 0
  t = 0d0 ! total time, s
  dt = 1d-3*spy ! initial time-step, s
  ! Initial mask
  old_mask = n(idx_surface_mask) 

  ! RUN loop:
  do i=1, model_length
    tend = time_arr(i)*spy
    print*,'tend: ', time_arr(i), 'yr, iteration: ', i

    ii = 0
    ! ----------------------------------------------------------------
    ! Compute the new dust rates
    if ((ntot_arr(i) .ne. ntot_current) .or. (Av_current .ne. Av_arr(i)) .or. (Tgas_current .ne. temperature_arr(i))) then
      print *, '------------------------------------------------------'
      print *, 'Calculating new rates for: ', ntot_arr(i), temperature_arr(i), Av_arr(i)
      call kemimo_computeRates(n(:), ntot_arr(i), temperature_arr(i), crflux, Av_arr(i), Ghabing=Gnot)
      ntot_ratio = ntot_arr(i)/ntot_current
      ntot_current = ntot_arr(i)
      Tgas_current = temperature_arr(i)
      Av_current = Av_arr(i)
      n(:) = n(:) * ntot_ratio
      ! Maintain mask
      print*, 'updating mask: ', n(idx_surface_mask), n(idx_mantle_mask)
      n(idx_surface_mask) = n(idx_surface_mask) / ntot_ratio
      n(idx_mantle_mask) = n(idx_mantle_mask) / ntot_ratio
      print*, 'new mask: ', n(idx_surface_mask), n(idx_mantle_mask)
      
      ! Set conservative dt if physical conditions changed ( experiment with this... ):
      dt = min(500d0*spy, dt)
      if (i > 1 .and. Tgas_current > temperature_arr(i-1)) dt = min(spy*50, dt)
      if (i > 1 .and. Tgas_current > 90d0) dt = min(spy*2.5, min(dt, tend - t))
    endif


    do while (t < tend)
      ! ----------------------------------------------------------------
      ! Check that while loop does not exceed the current tend:
      if ((t+(dt)) > tend) then
        print *, 'Close to the end, adjusting dt', (t+dt)/spy, tend/spy
        dt = (tend - t)*1.001d0
        ! if (dt < 1d-3) dt = dt * 1d-3
        ! print *, 'dt: ', dt/spy
      endif

      ! ----------------------------------------------------------------
      ! run chemistry
      call kemimo_dochem(n(:), dt)
      mask = n(idx_surface_mask) + n(idx_mantle_mask)
      t = t + dt
      

      ! --------------------------------------------------
      ! update dt
      if ((t+(2.0*dt)) < tend) then
        if (abs(mask)-old_mask > 1.3) then
          dt = dt * 0.96
        elseif (abs(mask)-old_mask > 0.7) then
          dt = dt
        else
          dt = dt * 1.04
        endif
      endif

      if (mod(ii, 50) == 0) print *, 'sub-step loop: ', ii, 'dt: ', (dt/spy), 'time until end: ', (tend-t)/spy

      ii = ii+1
      old_mask = mask
      ! ----------------------------------------------------------------      
      ! We also write output from subloop, but not every time.
      if (mod(ii, 4) == 0) then
        do iii=1,nmols
          nout(iii) = n(iii)/ntot_current
        end do
        write(unit,'(59999E17.8e3)') t/spy, Tgas_current, ntot_current, Av_current, nout(:)

      endif
  
    enddo

    ! ----------------------------------------------------------------
    ! write GAS+DUST species, 1:time, 2:Tgas, 3:ntot, 4:Av, 5->:species
    if (mod(n_loop, 1) == 0) then
      do iii=1,nmols
        nout(iii) = n(iii)/ntot_current
      end do
      write(unit,'(59999E17.8e3)') t/spy, Tgas_current, ntot_current, Av_current, nout(:)
    endif
    
    old_mask = n(idx_surface_mask)
  
    ! ----------------------------------------------------------------
    if (Tgas_current > 25.0d0) then
      print *, 'mask: ', mask
      !call kemimo_printFluxes(n(:), 10, (/idx_CO_0001/), .false.)
    endif

    ! ----------------------------------------------------------------
    ! increase step count
    n_loop = n_loop + 1

  end do

  ! ----------------------------------------------------------------
  ntot_current = 0.d0
  do j=91,nmols
    if (j == (nmolsu-1)) continue
    if (j == nmolsu) continue
    ntot_current = ntot_current + n(j)
  enddo

  ! ----------------------------------------------------------------
  print *, ' Final ice abundances: ', ntot_current
  print *, ' What this should be in Mly: ', (ntot_current / dust_sites)
  print *, theta_H2(Tgas_current, ntot_current)
  ! close unit
  close(unit)

end program loop


