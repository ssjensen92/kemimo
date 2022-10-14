module kemimo_commons

  !constants
  real*8,parameter::spy = 365d0*24d0*3600d0
  real*8,parameter::pi = acos(-1d0)
  real*8,parameter::qe = 4.803204d-10
  real*8,parameter::pi43 = 4d0*pi/3d0
  real*8,parameter::pmass=1.6726219d-24 !g
  real*8,parameter::gravity = 6.67259d-8 !cgs
  real*8,parameter::au2cm = 1.49597871d13 !AU->cm
  real*8,parameter::kb = 1.38064852d-16 !erg/K

  !!BEGIN_ARRAYSIZE
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-22 21:21:33
    ! CHANGESET: xxxxxxx
    ! BY: unknown@unknown

    !number of species
    integer,parameter::nmols=1865
    !number of unique species
    integer,parameter::nmolsu=1412
    !number of dust species
    integer,parameter::nmols_dust=906
    !number of reactions (nlayer*dust+gas)
    integer,parameter::nrea=52129
    !number of unique reactions (dust+gas)
    !number of dust-phase reactions
    integer,parameter::nreadust=3400
    !number of gas-phase reactions
    integer,parameter::nreagas=nrea-nreadust

    !monolayer thickness of each layer in model:
    integer,parameter::layerThickness=4
    !do swapping?
    logical, parameter:: doSwap=.true.
    !idx for mantle and surface species:
    integer,parameter::surface_start=957
    integer,parameter::surface_end=1409
    integer,parameter::mantle_start=1413
    integer,parameter::mantle_end=1865
    integer,parameter:: CO_desorption_idx=2201

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_ARRAYSIZE
  ! GRAIN PROPERTIES:
  real*8, parameter :: d2g = 1d-2 !dust/gas mass ratio
  real*8, parameter :: rho0 = 3d0 !bulk density, g/cm3
  real*8, parameter :: mu = 1.43 !mean molecular weight
  real*8, parameter :: app = 3d-8 !binding sites separation, cm
  real*8, parameter :: agrain = 1d-5
  real*8, parameter :: xdust = d2g * pmass / (4.0/3.0 * rho0 * agrain**3.0 * pi) !1.33d-12 ! dust density relative to n_H

  ! column densities for self-shielding:
  real*8 :: N_H2, N_CO, N_N2, N_HD, N_HI
  ! self-shielding factors/coefficients
  real*8 :: ss_CO, ss_H2, ss_HD, ss_N2

  ! UV RADIATION:
  real*8, parameter :: Fnot = 2d8 ! Draine ISRF
  real*8, parameter :: F_cr_not = 3d3 ! Cosmic-ray induced UV photons (never attenuated)
  real*8 :: F_cr

  ! Photodesorption yield for CO:
  real*8:: Y_CO

  ! dust sites pr cm^3
  real*8 :: ndns
  ! site density on dust grains
  real*8, parameter :: site_density = 1.5d15 ! grain sites pr cm^2
  ! H2 coverage factor
  real*8 :: H2_coverage
  !rate coefficients
  real*8::kall(nrea)
  ! swap rates & stick. coeff.
  real*8::kswap(nmols_dust)
  real*8::kstick(nmols_dust)
  real*8::ebind(nmols_dust)

  ! DEWSET VARS:
  integer,parameter:: ewt_flag=1
  real*8 :: ewt_fac(nmols)
  real*8 :: dn_surface
  integer:: reduce_dt

  !$omp threadprivate(kall)

  !prepare verbatim reaction names arrays
  integer,parameter::maxVerbatimSize=50
  character(len=maxVerbatimSize)::verbatim(nrea)
  !$omp threadprivate(verbatim)

  !!BEGIN_SPECIESNAMES
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-22 21:21:33
    ! CHANGESET: xxxxxxx
    ! BY: unknown@unknown

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_SPECIESNAMES

  !!BEGIN_IDXLIST
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-22 21:21:33
    ! CHANGESET: xxxxxxx
    ! BY: unknown@unknown

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_IDXLIST


  !!BEGIN_RPINDEX
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-05-11 12:33:15
    ! CHANGESET: e02cc77ec7ee5d47e175d29cf394a48b3b538665
    ! BY: unknown@unknown

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_RPINDEX


  contains

  function theta_H2(Tgas, n_H2_gas)
    implicit none
    real*8, intent(in) :: Tgas, n_H2_gas
    !real*8 :: H2_ice_density
    real*8 :: theta, E_critical, stick
    real*8 :: Tdust ! For now assumed equal to Tgas
    real*8 :: v_th
    ! Data array
    real*8, parameter :: data(121,2) = reshape( (/291.81700, 0.00003, &
                                                292.33683, 0.00022, &
                                                293.35351, 0.00041, &
                                                294.36680, 0.00058, &
                                                295.37667, 0.00071, &
                                                296.38740, 0.00085, &
                                                297.39643, 0.00098, &
                                                297.90605, 0.00109, &
                                                299.24972, 0.00125, &
                                                300.42692, 0.00140, &
                                                301.43424, 0.00151, &
                                                302.60861, 0.00164, &
                                                303.95341, 0.00180, &
                                                305.46270, 0.00196, &
                                                307.46799, 0.00211, &
                                                309.47669, 0.00229, &
                                                311.48113, 0.00243, &
                                                313.48387, 0.00256, &
                                                315.48661, 0.00269, &
                                                317.98621, 0.00281, &
                                                320.98493, 0.00296, &
                                                323.98139, 0.00309, &
                                                326.97727, 0.00321, &
                                                329.47347, 0.00331, &
                                                331.47281, 0.00341, &
                                                333.47299, 0.00352, &
                                                335.47403, 0.00363, &
                                                336.98077, 0.00377, &
                                                338.49006, 0.00393, &
                                                340.49875, 0.00411, &
                                                341.00667, 0.00420, &
                                                342.35204, 0.00437, &
                                                343.69798, 0.00455, &
                                                345.37402, 0.00471, &
                                                346.71995, 0.00489, &
                                                348.39826, 0.00507, &
                                                349.57745, 0.00524, &
                                                351.58529, 0.00541, &
                                                353.42582, 0.00557, &
                                                355.60013, 0.00575, &
                                                357.60288, 0.00587, &
                                                359.60647, 0.00601, &
                                                362.10493, 0.00613, &
                                                365.09968, 0.00624, &
                                                368.58860, 0.00633, &
                                                372.57111, 0.00640, &
                                                376.54852, 0.00642, &
                                                380.52039, 0.00639, &
                                                384.48844, 0.00634, &
                                                388.45393, 0.00626, &
                                                392.41773, 0.00617, &
                                                395.88552, 0.00608, &
                                                399.35075, 0.00597, &
                                                403.30944, 0.00584, &
                                                406.77439, 0.00573, &
                                                410.23821, 0.00561, &
                                                413.70146, 0.00549, &
                                                417.16443, 0.00536, &
                                                421.12141, 0.00521, &
                                                425.07925, 0.00507, &
                                                428.54193, 0.00494, &
                                                431.50890, 0.00482, &
                                                433.98269, 0.00473, &
                                                436.94894, 0.00461, &
                                                440.40978, 0.00446, &
                                                443.37674, 0.00435, &
                                                446.34256, 0.00422, &
                                                449.31066, 0.00411, &
                                                452.27592, 0.00397, &
                                                455.24401, 0.00387, &
                                                458.21041, 0.00374, &
                                                461.17850, 0.00363, &
                                                464.64062, 0.00350, &
                                                468.10103, 0.00335, &
                                                471.07083, 0.00326, &
                                                474.03836, 0.00314, &
                                                477.00759, 0.00304, &
                                                479.97511, 0.00293, &
                                                482.94491, 0.00283, &
                                                485.91244, 0.00272, &
                                                488.88223, 0.00262, &
                                                491.85033, 0.00252, &
                                                494.81842, 0.00241, &
                                                497.78765, 0.00231, &
                                                500.75575, 0.00220, &
                                                503.72498, 0.00210, &
                                                506.69421, 0.00200, &
                                                509.66400, 0.00190, &
                                                513.12839, 0.00179, &
                                                516.59391, 0.00168, &
                                                520.05929, 0.00158, &
                                                524.01968, 0.00146, &
                                                527.98262, 0.00136, &
                                                531.94641, 0.00126, &
                                                535.91489, 0.00121, &
                                                539.88506, 0.00117, &
                                                543.85438, 0.00112, &
                                                547.32132, 0.00103, &
                                                550.28885, 0.00092, &
                                                553.75451, 0.00081, &
                                                557.72638, 0.00079, &
                                                561.69911, 0.00077, &
                                                565.67014, 0.00074, &
                                                569.63989, 0.00070, &
                                                573.61049, 0.00066, &
                                                577.58151, 0.00063, &
                                                581.55211, 0.00059, &
                                                585.52186, 0.00055, &
                                                589.49246, 0.00051, &
                                                593.46434, 0.00049, &
                                                597.43621, 0.00046, &
                                                601.40809, 0.00044, &
                                                605.37954, 0.00041, &
                                                609.35056, 0.00038, &
                                                613.32159, 0.00035, &
                                                617.29262, 0.00031, &
                                                621.26194, 0.00027, &
                                                625.22999, 0.00021, &
                                                629.19888, 0.00016, &
                                                632.66462, 0.00010, &
                                                635.14670, 0.00004/), (/121,2/), order=(/2,1/))
    ! Fixed parameters:
    real*8, parameter :: nu = 1.0d12
    real*8, parameter :: n_sites = 1.5d15
    real*8, parameter :: mass_H2 = 2.0 * pmass ! grams
    real*8, parameter :: beta = 2.5
    real*8, parameter :: T0 = 87.0 !87.0 np-ASW ice  - 56 silicate
    real*8, parameter :: S0 = 0.76 !0.76 np-ASW ice  - 0.95 silicate
    integer :: i
    real*8 :: delta_Eb, coverage_avg, Eb_avg
    real*8 :: theta_H2

    Tdust = Tgas
    stick = S0 * (1.d0+beta*Tgas/T0)/(1.d0+Tgas/T0)**beta
    v_th = sqrt(8.d0 * kb * Tgas / (mass_H2 * pi))
    E_critical = Tgas*log(4.0d0*nu*n_sites/(stick * v_th*n_H2_gas))

    ! integrate
    theta_H2 = 0.d0
    do i = 1,120
        coverage_avg = (data(i,2) + data(i+1,2))/2.d0
        Eb_avg = (data(i,1) + data(i+1,1))/2.d0
        delta_Eb = (data(i+1,1)-data(i,1))

        theta = (1.d0 + exp(-(Eb_avg-E_critical)/Tgas))**(-1.d0)
        !write(*,*) theta
        theta_H2 = theta_H2 + theta * delta_Eb * coverage_avg
    end do
    
    ! Integral should be surfae coverage factor (between 0 and 1)
    !print *, 'theta_H2: ', theta_H2
    !H2_ice_density = 4.44d-6*n_H2_gas*theta_H2 ! factor from Hocuk et al 2015
    return

  end function 

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


  !************************
  !2d bilinear interp
  ! x,y: point where to evaluate f(x,y)

  function fit2d(x, y, xarr, yarr, xydata) result(f)
    implicit none
    real*8,intent(inout) :: x, y
    real*8, intent(in), dimension(:) :: xarr, yarr
    real*8, intent(in), dimension(:,:) :: xydata
    integer :: idx, idy
    real*8 :: dx, dy, invdx, invdy
    real*8 :: x0, y0, f, f00, f01, f10, f11

    if (x > maxval(xarr)) x = maxval(xarr)
    if (y > maxval(yarr)) y = maxval(yarr)

    !index of data point corresponding to x, y
    idx = minloc(abs(xarr - x), dim=1)
    idy = minloc(abs(yarr - y), dim=1)

    !function at idx, idx+1
    f00 = xydata(idx, idy)
    f01 = xydata(idx, idy+1)
    f10 = xydata(idx+1, idy)
    f11 = xydata(idx+1, idy+1)

    ! inverse dx,dy:
    if (size(xarr) > idx+1) then
      dx = abs(xarr(idx+1) - xarr(idx))
    else
      dx = abs(xarr(idx) - xarr(idx-1))
    endif
    if (size(yarr) > idy+1) then
      dy = abs(yarr(idy+1) - yarr(idy))
    else
      dy = abs(yarr(idy) - yarr(idy-1))
    endif
    invdx = 1d0 / dx
    invdy = 1d0 / dy

    !get x variable corresponding to idx
    x0 = xarr(idx)
    y0 = yarr(idy)

    !compute linearly interpolated function
    f = (x-x0)*invdx*(y-y0)*invdy*f11 + &
      (y-y0)*invdy*(1d0-(x-x0)*invdx) * f01 + &
      (x-x0)*invdx*(1d0-(y-y0)*invdy) * f10 + &
      (1d0 - (x-x0)*invdx - (y-y0)*invdy + (x-x0)*invdx*(y-y0)*invdy) * f00

  end function fit2d



  function calc_ss(N_H2, N_X, fname, nrows, ncols, nheader) result(f)
    implicit none
    real*8, intent(inout) :: N_H2, N_X
    character(len = *), intent(inout):: fname
    integer, intent(in) :: nrows, ncols, nheader
    real*8 :: f
    integer io, i
    real*8, dimension(nrows) :: arr_X
    real*8, dimension(ncols) :: arr_H2
    real*8, dimension(ncols,nrows) :: table
    character(len=1000):: modelDir

    CALL get_environment_variable("KEMIMO_WORKDIR", modelDir)

    open(newunit=io, file=trim(modelDir)//trim(fname), status="old")
    ! read header
    do i=1, nheader
      read(io, *)
    enddo
    ! read N_H2, N_CO
    do i=1, nrows
      read(io, *) arr_X(i)
    end do
    ! read title
    read(io, *)
    do i=1, ncols
      read(io, *) arr_H2(i)
    end do
    ! read title
    read(io, *)
    ! read 12C16O ss table:
    do i=1, ncols
      read(io, *) table(i, :)
    enddo
    close(io)

    ! Interp table:
    f = fit2d(N_X, N_H2, arr_X, arr_H2, table)
    return

  end function calc_ss

end module kemimo_commons
