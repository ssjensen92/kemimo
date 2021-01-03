module kemimo_reactionarray
  use kemimo_commons

  integer,parameter::maxRP2=9
  integer, dimension(nrea-1, maxRP2):: reactionArray
  integer, dimension(nmols-3, nrea-1):: jacArray

  contains
  
  subroutine init_reactionarray
    implicit none
    integer::unit, i, j
    open(newunit=unit, file="reactionarray.dat", status="old")
    read(unit, *) ((reactionArray(i,j), j=1,maxRP2), i=1,nrea-1)
    close(unit)

    open(newunit=unit, file="jacarray.dat", status="old")
    read(unit, *) ((jacArray(i,j), j=1,nrea-1), i=1,nmols-3)
    close(unit)
  end subroutine 

end module kemimo_reactionarray
