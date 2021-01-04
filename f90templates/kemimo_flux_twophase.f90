module kemimo_flux
contains

  !*****************
  !print fluxes to stdout
  subroutine printFluxes(n, nflux, idxList, onlyDust)
    use kemimo_commons
    use kemimo_reactionarray
    implicit none
    real*8,intent(in)::n(nmols)
    integer,intent(in)::nflux,idxList(:)
    real*8::flux(nrea-1)
    integer::idx(nrea-1),i
    logical,intent(in)::onlyDust

    !get fluxes
    if(onlyDust) then
       flux(:) = fluxesDust(n(:), idxList)
    else
       flux(:) = fluxesAll(n(:), idxList)
    end if

    !sort indexes
    idx(:) = sortedIndexs(flux(:))

    print *,"************************"
    do i=1,min(nflux,nrea-1)
       print '(2I5,2E17.8e3,2a)', i, reactionArray(idx(i),1), kall(reactionArray(idx(i),1)), &
            flux(idx(i)), "  ", verbatim(idx(i))
    end do
  end subroutine printFluxes


  !*********************
  !save fluxes to file
  subroutine dumpFluxes(n, unit, xvar)
    use kemimo_commons
    implicit none
    real*8,intent(in)::n(nmols), xvar
    real*8::flux(nrea-1)
    integer,intent(in)::unit
    integer::i

    !load fluxes, -1 means no species filter
    flux(:) = fluxesAll(n(:), (/-1/))

    !loop on reactions
    do i=1,nrea-1
       write(unit,*) xvar, i, flux(i)
    end do
    write(unit,*)

  end subroutine dumpFluxes

  !***************
  !load verbatim reaction
  subroutine loadverbatim()
    use kemimo_commons
    implicit none
    integer::unit,ios,i
    character*50::fname

    !file with verbatim
    fname = "verbatim.dat"

    !open to read
    open(newunit=unit, file=trim(fname), status="old", iostat=ios)
    !check if file exists
    if(ios/=0) then
       print *,"ERROR: problem while loading", trim(fname)
       stop
    end if

    !loop on reactions
    do i=1,nrea
       !store into common verbatim array
       read(unit,'(a)',iostat=ios) verbatim(i)
       !check if reading after EOF
       if(ios/=0) then
          print *,"ERROR: problem while reading", trim(fname)
          stop
       end if
    end do

    !close file
    close(unit)

  end subroutine loadverbatim

  !****************
  function getkMult(idxList) result(fmul)
    use kemimo_commons
    use kemimo_reactionarray
    implicit none
    integer,intent(in)::idxList(:)
    integer::i,j,k
    real*8::fmul(nrea-1)
    logical::hasSpec

    !if species index required
    !NOTE: array access here is not well written
    if(minval(idxList)>0) then
       !loop on reactions
       do i=1,nrea-1
         !has species flag
         hasSpec = .false.
          !loop on reactants+products indexes
          do j=3,maxRP2-1
             !loop on required species list
             do k=1,size(idxList)
                !check if species is there
                if(idxList(k)==reactionArray(i,j)) then
                   hasSpec = .true.
                   exit
                end if
             end do
             if(hasSpec) exit
          end do
          !set multiplicator to 1 if has species
          if(hasSpec) then
             fmul(i) = 1d0
          else
             fmul(i) = 0d0
          end if
       end do
    end if

  end function getkMult

  !****************
  !get reaction fluxes, cm-3/s
  function fluxesDust(n, idxList) result(flux)
      use kemimo_commons
      use kemimo_reactionarray
      implicit none
      real*8::n(nmols)
      real*8::flux(nrea-1),fmul(nrea-1)
      integer,intent(in)::idxList(:)
      integer:: layer, i, rtype

      fmul(:) = getkMult(idxList)


      flux(:) = 0d0
      n(idx_dummy) = 1d0
      do i=1, nrea-1
         ! flux = rate * n1 * n2
         ! Reaction array:
         !   1: reaction_idx
         !   2: reaction layer
         !   3,4: reactants
         !   5,6,7,8: products
         !   9: mask integer (indicates change in mask for the specific reaction, 0 = no change in ice thickness)
         flux(i) = kall(reactionArray(i,1)) * n(reactionArray(i,3)) * n(reactionArray(i,4))
         layer = reactionArray(i,2)
         rtype = reactionArray(i,9)
         ! Ignore if the layer is zero (gas-phase)
         if (layer < 1) flux(i) = 0d0

         if (layer > 0) then
            ! two-phase:
            if (n(idx_surface_mask)*layerThickness > 1d0) then		
               ! limit thermal desorption, CR desorption and photoprocesses to *layerthickness* of mly:
               if (rtype == 1 .or. rtype == 2 .or. rtype == 3 .or. rtype == 4) then
                  flux(i) = flux(i) * min(1d0, layerThickness/n(idx_surface_mask))
               ! 2body reactions:
               elseif (rtype == 5) then 
                  if (n(idx_surface_mask) > 1d0) then
                     ! From sect. 7 of Cuppen+2017, with minor adjustments:
                     flux(i) = flux(i) / (n(idx_surface_mask))
                  endif
               endif
            endif
         endif

      end do
      !if species filter required use multiplicator
      if(minval(idxList)>0) then
         flux(:) = flux(:)*fmul(:)
      end if

  end function fluxesDust

   !****************
  !get reaction fluxes, cm-3/s
  function fluxesAll(n, idxList) result(flux)
      use kemimo_commons
      use kemimo_reactionarray
      implicit none
      real*8::n(nmols)
      real*8::flux(nrea-1),fmul(nrea-1)
      integer,intent(in)::idxList(:)
      integer:: layer, i, rtype

      fmul(:) = getkMult(idxList)
      flux(:) = 0d0


      flux(:) = 0d0
      n(idx_dummy) = 1d0
      do i=1, nrea-1
         ! flux = rate * n1 * n2 * layer_selection
         ! Reaction array:
         !   1: reaction_idx
         !   2: reaction layer
         !   3,4: reactants
         !   5,6,7,8: products
         !   9: mask integer (indicates change in mask for the specific reaction, 0 = no change in ice thickness)
        layer = reactionArray(i,2)
        rtype = reactionArray(i,9)
        if (layer > 1) then
          cycle
        endif
        
        flux(i) = kall(reactionArray(i,1)) * n(reactionArray(i,3)) * n(reactionArray(i,4))
        if (layer > 0) then
          ! two-phase:
          if (n(idx_surface_mask)*layerThickness > 1d0) then		
            ! limit thermal desorption, CR desorption and photoprocesses to *layerthickness* of mly:
            if (rtype == 1 .or. rtype == 2 .or. rtype == 3 .or. rtype == 4) then
              flux(i) = flux(i) * min(1d0, layerThickness/n(idx_surface_mask))
            ! 2body reactions:
            elseif (rtype == 5) then 
              if (n(idx_surface_mask) > 1d0) then
                ! From sect. 7 of Cuppen+2017, with minor adjustments:
                flux(i) = flux(i) / (n(idx_surface_mask))
              endif
            endif
          endif
         endif
         
      end do
      !if species filter required use multiplicator
      if(minval(idxList)>0) then
         flux(:) = flux(:)*fmul(:)
      end if

  end function fluxesAll

  !**************************
  !return index of ascending sorted array v
  ! bubble sorting, faster algorithms are welcomed
  function sortedIndexs(v) result(idx)
    implicit none
    real*8,intent(in)::v(:)
    real*8::vloc(size(v)),tmp
    integer::idx(size(v)),i,itmp
    logical::noSwap

    !store non-ordered indexes
    do i=1,size(v)
       idx(i) = i
    end do

    !local copy to swap
    vloc(:) = v(:)

    !loop until no swap
    do
       !logical to check if swapping
       noSwap = .true.
       !loop on items
       do i=2,size(v)
          !if descending swap
          if(vloc(i-1)<vloc(i)) then
             noSwap = .false.
             !swap value array elements
             tmp = vloc(i)
             vloc(i) = vloc(i-1)
             vloc(i-1) = tmp
             !swap index array elements
             itmp = idx(i)
             idx(i) = idx(i-1)
             idx(i-1) = itmp
          end if
       end do
       !if no swaps break loop
       if(noSwap) exit
    end do

  end function sortedIndexs

end module kemimo_flux
