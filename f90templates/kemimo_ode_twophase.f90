module kemimo_ode
contains

  !************************
  !evolve chemistry for a time-step dt (s)
  ! n(:) are species number densities
  subroutine dochem(n,dt)
    use kemimo_commons
    use kemimo_rates
    implicit none
    real*8,intent(inout)::n(nmols)
    real*8 ::dt
    real*8 :: ni(nmols)
    integer::i, ncount

    !DLSODES variables
    integer,parameter::meth=2 !1=adam, 2=BDF
    integer,parameter::lrw=20+16*nmols+3*nmols**2
    integer,parameter::liw=30
    !tolerances
    real*8,parameter::rtol(nmols) = 1d-4
    real*8,parameter::atol(nmols) = 1d-20
    integer::neqa(1),itol,itask,istate,iopt,mf
    integer::iwork(liw)
    real*8::rwork(lrw),tloc

    iwork(:) = 0
    rwork(:) = 0d0
    itol = 4 !both tolerances are arrays

    !number of equations is a single sized array
    neqa(:) = nmols


    !  FROM DLSODES manual
    !  Name    Location   Meaning and default value
    !  ------  ---------  -----------------------------------------------
    !  H0      RWORK(5)   Step size to be attempted on the first step.
    !                     The default value is determined by the solver.
    !  HMAX    RWORK(6)   Maximum absolute step size allowed.  The
    !                     default value is infinite.
    !  HMIN    RWORK(7)   Minimum absolute step size allowed.  The
    !                     default value is 0.  (This lower bound is not
    !                     enforced on the final step before reaching
    !                     TCRIT when ITASK = 4 or 5.)
    !  MAXORD  IWORK(5)   Maximum order to be allowed.  The default value
    !                     is 12 if METH = 1, and 5 if METH = 2. (See the
    !                     MF description above for METH.)  If MAXORD
    !                     exceeds the default value, it will be reduced
    !                     to the default value.  If MAXORD is changed
    !                     during the problem, it may cause the current
    !                     order to be reduced.
    !  MXSTEP  IWORK(6)   Maximum number of (internally defined) steps
    !                     allowed during one call to the solver.  The
    !                     default value is 500.
    !  MXHNIL  IWORK(7)   Maximum number of messages printed (per
    !                     problem) warning that T + H = T on a step
    !                     (H = step size).  This must be positive to
    !                     result in a nondefault value.  The default
    !                     value is 10.


    itask = 1
    iopt = 1
    iwork(6) = int(1e5) !maximum number of iteration before warning
    iwork(5) = 2 !maximum integration order
    MF = 121 ! hardcode jacobian, internally generated sparsity

    istate = 1
    tloc = 0d0

    ! store initial abundances:
    ni(:) = n(:)
    ! count attempts
    ncount = 0
    do
       !call solver, see DLSODES documentation
       CALL DLSODES(fex, neqa(:), n(:), tloc, dt, &
            ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, &
            LIW, JES, MF, dewset)
       !check solver output state
       if(istate==-1) then
          !maximum number of iteration reached, continue
          istate = 1
          cycle
       elseif(istate==2) then
          !success integration
          exit
       elseif(istate==-5) then
          !problem with sparsity, need to recompute
          istate = 3
          cycle
       elseif(istate==-3) then
          ncount = ncount + 1
          !problem with input.
          if (ncount > 10) then 
            print *, 'Attempted 10 times with state -3 from solver. Giving up.'
            stop
          endif
          print *, 'Solver returned state -3. Returning abundances to the initial value and reducing dt'
          n(:) = ni(:)
          dt = dt / 3d0
          istate = 1
       else
          !unknonw problem stop program
          print *,istate
          stop
       end if
    end do

  end subroutine dochem

  !*********************
  !differential equations, returns dn(:)
  ! see DLSODES documentation
  subroutine fex(neq, tt, n, dn)
    use kemimo_commons
    use kemimo_reactionarray
    use kemimo_rates
    implicit none
    integer::neq,i
    real*8::n(neq),tt,dn(neq)
    real*8:: flux
    integer:: layer, rtype
    real*8 :: Nsurface, Nmantle
    real*8 :: alpha, dnsdt, R
    integer :: offset

    dn(:) = 0d0
    n(idx_dummy) = 1d0
    ! ----------------------------------------------------------------------
    ! Chemistry part:

    dn(:) = 0d0
    n(idx_dummy) = 1d0
    do i=1, nrea-1
      ! flux = rate * n1 * n2 * layer_selection
      ! Reaction array:
      !   1: reaction_idx
      !   2: reaction layer
      !   3,4: reactants
      !   5,6,7,8: products
      !   9: reaction type (see utils.py for getReactionType function)
        flux = kall(reactionArray(i,1)) * n(reactionArray(i,3)) * n(reactionArray(i,4))
        layer = reactionArray(i,2)
        rtype = reactionArray(i,9)
        ! Ignore if the layer is zero (gas-phase)
        if (layer > 0) then
          ! two-phase:
          if (n(idx_surface_mask)*layerThickness > 1d0) then		
            ! limit thermal desorption, CR desorption and photoprocesses to *layerthickness* of mly:
            if (rtype == 1 .or. rtype == 2 .or. rtype == 3 .or. rtype == 4) then
              flux = flux * min(1d0, layerThickness/n(idx_surface_mask))
            ! 2body reactions:
            elseif (rtype == 5) then 
              if (n(idx_surface_mask) > 1d0) then
                ! From sect. 7 of Cuppen+2017, with minor adjustments:
                flux = flux / (n(idx_surface_mask))
              endif
            endif
          endif
          
        endif

        dn(reactionArray(i,3)) = dn(reactionArray(i,3)) - flux
        dn(reactionArray(i,4)) = dn(reactionArray(i,4)) - flux
        dn(reactionArray(i,5)) = dn(reactionArray(i,5)) + flux
        dn(reactionArray(i,6)) = dn(reactionArray(i,6)) + flux
        dn(reactionArray(i,7)) = dn(reactionArray(i,7)) + flux
        dn(reactionArray(i,8)) = dn(reactionArray(i,8)) + flux

    end do
    ! "Reaction rate" for mask layer
    !!! BEGIN MASK DN
    do i=surface_start, surface_end
      dn(idx_surface_mask) = dn(idx_surface_mask) + dn(i)
    enddo
    dn(idx_surface_mask) = dn(idx_surface_mask)*kall(nrea)
    !!! END MASK DN

    end subroutine fex

  !****************
  !differential equations, only formation
  function fexForm(n) result(dnf)
    use kemimo_commons
    implicit none
    real*8,intent(in)::n(nmols)
    real*8::dnf(nmols)

    !!BEGIN_RHSFORM
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_RHSFORM

  end function fexForm

  !***************************
  !Jacobian, pd(i,j)=df(i)/dx(j), see DLSODES documentation
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use kemimo_commons
    use kemimo_reactionarray
    use kemimo_rates
    implicit none
    integer::neq, j, ian, jan, i, ii, offset, layer, rtype
    real*8::tt, n(neq), pdj(neq)
    real*8:: flux
    real*8 :: Nsurface, Nmantle
    real*8 :: alpha, dnsdt, R

    !!BEGIN_JACOBIAN
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-06-06 11:13:08
    ! CHANGESET: 342b06597f6aa77e3ec690d11d948c1af26894dc
    ! BY: unknown@unknown

    n(idx_dummy) = 1d0
    pdj(:) = 0d0
    if (j == idx_dummy) return

    ! Chemistry part:
    do ii=1, nrea-1
      ! Get reaction i from jacArray for species j, entry ii      
      i = jacArray(j, ii)
      ! If mantle species or dummy, exit loop.
      if (real(i) < 0d0 .or. j > surface_end) exit

      ! ######################
      ! flux = rate * n1 * n2
      ! Reaction array:
      !   1: reaction_idx
      !   2: reaction layer
      !   3,4: reactants
      !   5,6,7,8: products
      !   9: reaction type (see utils.py for getReactionType function)

      ! if j in in reaction, this should not happen!
      !if (all(reactionArray(i,3:4) .ne. j)) print*, j, i

      ! We ignore higher layers in the chemical part:
      layer = reactionArray(i,2)
      rtype = reactionArray(i,9)
      if (layer > 1) cycle

      ! Find flux without j:
      if (reactionArray(i,3) == j) then
        flux = kall(reactionArray(i,1)) * n(reactionArray(i,4))
      elseif (reactionArray(i,4) == j) then
        flux = kall(reactionArray(i,1)) * n(reactionArray(i,3))
      else
        flux = 0d0
      endif
      


      if (layer == 1) then
        if (rtype == 5) flux = flux * 1d0/max(1d0, min(n(idx_surface_mask), real(layerThickness)))
        if (rtype == 4) flux = flux * min(n(idx_surface_mask), 4d0)
        if (rtype == 3) flux = flux * min(n(idx_surface_mask), 1d0)
      endif


      pdj(reactionArray(i,3)) = pdj(reactionArray(i,3)) - flux
      pdj(reactionArray(i,4)) = pdj(reactionArray(i,4)) - flux
      pdj(reactionArray(i,5)) = pdj(reactionArray(i,5)) + flux
      pdj(reactionArray(i,6)) = pdj(reactionArray(i,6)) + flux
      pdj(reactionArray(i,7)) = pdj(reactionArray(i,7)) + flux
      pdj(reactionArray(i,8)) = pdj(reactionArray(i,8)) + flux

    end do

    pdj(idx_dummy) = 0d0

    do i=surface_start, surface_end
      pdj(idx_surface_mask) = pdj(idx_surface_mask) + pdj(i)
    enddo
    pdj(idx_surface_mask) = pdj(idx_surface_mask) * kall(nrea)

    ewt_fac(:) = 1d0
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_JACOBIAN

  end subroutine jes

!########################################################################
  subroutine dewset(n, itol, rtol_arr, atol_arr, ycur, ewt)
    use kemimo_commons
    implicit none

    integer :: n, itol
    real*8, dimension(nmols) :: ewt
    real*8 :: rtol_arr(*), atol_arr(*), ycur(*)
    integer :: i

    select case(itol)
    case(1)
        do i = 1, nmols
          ewt(i) = rtol_arr(1)*abs(ycur(i))*ewt_fac(i) + atol_arr(1)
        enddo
    case(2)
        do i = 1, nmols
          ewt(i) = rtol_arr(1)*abs(ycur(i))*ewt_fac(i) + atol_arr(i)
        enddo
    case(3)
        do i = 1, nmols
          ewt(i) = rtol_arr(i)*abs(ycur(i))*ewt_fac(i) + atol_arr(1)
        enddo
    case(4)
        do i = 1, nmols
          ewt(i) = rtol_arr(i)*abs(ycur(i))*ewt_fac(i) + atol_arr(i)
        enddo
    end select

    return

  end subroutine dewset


end module kemimo_ode
