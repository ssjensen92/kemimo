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
  real*8, parameter :: xdust = 1.33d-12 ! dust density relative to n_H

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

    character(len=maxVerbatimSize),parameter::speciesNames(nmolsu) = (/"GRAIN0_gas     ", &
        "GRAIN-_gas     ", &
        "E_gas          ", &
        "H_gas          ", &
        "H+_gas         ", &
        "H-_gas         ", &
        "D_gas          ", &
        "D+_gas         ", &
        "D-_gas         ", &
        "p_H2_gas       ", &
        "o_H2_gas       ", &
        "p_H2+_gas      ", &
        "o_H2+_gas      ", &
        "p_D2_gas       ", &
        "o_D2_gas       ", &
        "p_D2+_gas      ", &
        "o_D2+_gas      ", &
        "p_D3+_gas      ", &
        "o_D3+_gas      ", &
        "m_D3+_gas      ", &
        "p_D2H+_gas     ", &
        "o_D2H+_gas     ", &
        "p_H2D+_gas     ", &
        "o_H2D+_gas     ", &
        "HD_gas         ", &
        "HD+_gas        ", &
        "He_gas         ", &
        "He+_gas        ", &
        "O_gas          ", &
        "O+_gas         ", &
        "O-_gas         ", &
        "O2_gas         ", &
        "O2+_gas        ", &
        "O3_gas         ", &
        "OH_gas         ", &
        "OH+_gas        ", &
        "OH-_gas        ", &
        "OD_gas         ", &
        "OD+_gas        ", &
        "OD-_gas        ", &
        "H2O_gas        ", &
        "H2O+_gas       ", &
        "HDO_gas        ", &
        "HDO+_gas       ", &
        "D2O_gas        ", &
        "D2O+_gas       ", &
        "O2H_gas        ", &
        "O2D_gas        ", &
        "HO2+_gas       ", &
        "DO2+_gas       ", &
        "HOOH_gas       ", &
        "HOOD_gas       ", &
        "DOOH_gas       ", &
        "DOOD_gas       ", &
        "p_H3+_gas      ", &
        "o_H3+_gas      ", &
        "H3O+_gas       ", &
        "H2DO+_gas      ", &
        "HD2O+_gas      ", &
        "D3O+_gas       ", &
        "F_gas          ", &
        "F+_gas         ", &
        "Cl_gas         ", &
        "Cl+_gas        ", &
        "CF+_gas        ", &
        "C_gas          ", &
        "C+_gas         ", &
        "C-_gas         ", &
        "C2_gas         ", &
        "C2+_gas        ", &
        "C2N+_gas       ", &
        "CO_gas         ", &
        "CO+_gas        ", &
        "CO2_gas        ", &
        "CO2+_gas       ", &
        "N_gas          ", &
        "N+_gas         ", &
        "N2_gas         ", &
        "N2+_gas        ", &
        "N2H+_gas       ", &
        "N2D+_gas       ", &
        "HCO_gas        ", &
        "DCO_gas        ", &
        "HCO+_gas       ", &
        "DCO+_gas       ", &
        "H2CO_gas       ", &
        "HDCO_gas       ", &
        "D2CO_gas       ", &
        "CH2OH_gas      ", &
        "CD2OD_gas      ", &
        "CH2OD_gas      ", &
        "CHDOH_gas      ", &
        "CHDOD_gas      ", &
        "CD2OH_gas      ", &
        "CH3O_gas       ", &
        "CHD2O_gas      ", &
        "CH2DO_gas      ", &
        "CD3O_gas       ", &
        "CH3OH_gas      ", &
        "CH3OD_gas      ", &
        "CHD2OH_gas     ", &
        "CHD2OD_gas     ", &
        "CH2DOH_gas     ", &
        "CH2DOD_gas     ", &
        "CD3OD_gas      ", &
        "CD3OH_gas      ", &
        "CH_gas         ", &
        "CD_gas         ", &
        "CH2_gas        ", &
        "CHD_gas        ", &
        "CD2_gas        ", &
        "CH3_gas        ", &
        "CH2D_gas       ", &
        "CHD2_gas       ", &
        "CD3_gas        ", &
        "CH4_gas        ", &
        "CH3D_gas       ", &
        "CH2D2_gas      ", &
        "CHD3_gas       ", &
        "CD4_gas        ", &
        "CH+_gas        ", &
        "CD+_gas        ", &
        "CH2+_gas       ", &
        "CHD+_gas       ", &
        "CD2+_gas       ", &
        "CH3+_gas       ", &
        "CH2D+_gas      ", &
        "CHD2+_gas      ", &
        "CD3+_gas       ", &
        "CH4+_gas       ", &
        "CH3D+_gas      ", &
        "CH2D2+_gas     ", &
        "CHD3+_gas      ", &
        "CD4+_gas       ", &
        "CH5+_gas       ", &
        "CH4D+_gas      ", &
        "CH3D2+_gas     ", &
        "CH2D3+_gas     ", &
        "CHD4+_gas      ", &
        "CD5+_gas       ", &
        "HCOOH_gas      ", &
        "HCOOD_gas      ", &
        "DCOOH_gas      ", &
        "DCOOD_gas      ", &
        "HOCO_gas       ", &
        "DOCO_gas       ", &
        "NH_gas         ", &
        "NH+_gas        ", &
        "ND_gas         ", &
        "ND+_gas        ", &
        "NH2_gas        ", &
        "NH2+_gas       ", &
        "NHD_gas        ", &
        "NHD+_gas       ", &
        "ND2_gas        ", &
        "ND2+_gas       ", &
        "NH3_gas        ", &
        "NH3+_gas       ", &
        "NH2D_gas       ", &
        "NH2D+_gas      ", &
        "NHD2_gas       ", &
        "NHD2+_gas      ", &
        "ND3_gas        ", &
        "ND3+_gas       ", &
        "NH4+_gas       ", &
        "NH3D+_gas      ", &
        "NH2D2+_gas     ", &
        "NHD3+_gas      ", &
        "ND4+_gas       ", &
        "Fe_gas         ", &
        "Fe+_gas        ", &
        "S_gas          ", &
        "S+_gas         ", &
        "S-_gas         ", &
        "C10_gas        ", &
        "C10+_gas       ", &
        "C10-_gas       ", &
        "C10H_gas       ", &
        "C10H+_gas      ", &
        "C10H-_gas      ", &
        "C10H2_gas      ", &
        "C10H2+_gas     ", &
        "C10H3+_gas     ", &
        "C10N_gas       ", &
        "C10N+_gas      ", &
        "C11_gas        ", &
        "C11+_gas       ", &
        "C2H+_gas       ", &
        "C2D+_gas       ", &
        "C2H2_gas       ", &
        "C2H2+_gas      ", &
        "C2HD_gas       ", &
        "C2HD+_gas      ", &
        "C2D2_gas       ", &
        "C2D2+_gas      ", &
        "C2DO+_gas      ", &
        "C2H3_gas       ", &
        "C2H3+_gas      ", &
        "C2H2D_gas      ", &
        "C2H2D+_gas     ", &
        "C2HD2_gas      ", &
        "C2HD2+_gas     ", &
        "C2D3_gas       ", &
        "C2D3+_gas      ", &
        "C2H4_gas       ", &
        "C2H4+_gas      ", &
        "C2H3D_gas      ", &
        "C2H3D+_gas     ", &
        "C2H2D2_gas     ", &
        "C2H2D2+_gas    ", &
        "C2HD3_gas      ", &
        "C2HD3+_gas     ", &
        "C2D4_gas       ", &
        "C2D4+_gas      ", &
        "C2H4O+_gas     ", &
        "C2H5_gas       ", &
        "C2H5+_gas      ", &
        "C2H5OH+_gas    ", &
        "C2H5OH2+_gas   ", &
        "C2H6_gas       ", &
        "C2H6+_gas      ", &
        "C2H6CO+_gas    ", &
        "C2H7+_gas      ", &
        "C2HO+_gas      ", &
        "C2N2+_gas      ", &
        "C2O+_gas       ", &
        "C2S+_gas       ", &
        "C3_gas         ", &
        "C3+_gas        ", &
        "C3-_gas        ", &
        "C3H+_gas       ", &
        "C3D+_gas       ", &
        "C3H3N+_gas     ", &
        "C3H3NH+_gas    ", &
        "C3H4_gas       ", &
        "C3H4+_gas      ", &
        "C3H5+_gas      ", &
        "C3H6OH+_gas    ", &
        "C3N_gas        ", &
        "C3N+_gas       ", &
        "C3N-_gas       ", &
        "C3O_gas        ", &
        "C3O+_gas       ", &
        "C3S_gas        ", &
        "C3S+_gas       ", &
        "C4_gas         ", &
        "C4+_gas        ", &
        "C4-_gas        ", &
        "C4H_gas        ", &
        "C4H+_gas       ", &
        "C4H-_gas       ", &
        "C4D_gas        ", &
        "C4D+_gas       ", &
        "C4D-_gas       ", &
        "C4H2_gas       ", &
        "C4H2+_gas      ", &
        "C4HD_gas       ", &
        "C4HD+_gas      ", &
        "C4D2_gas       ", &
        "C4D2+_gas      ", &
        "C4H3_gas       ", &
        "C4H3+_gas      ", &
        "C4H4+_gas      ", &
        "C4H5+_gas      ", &
        "C4H7+_gas      ", &
        "C4N_gas        ", &
        "C4N+_gas       ", &
        "C4S_gas        ", &
        "C4S+_gas       ", &
        "C5_gas         ", &
        "C5+_gas        ", &
        "C5-_gas        ", &
        "C5H_gas        ", &
        "C5H+_gas       ", &
        "C5H-_gas       ", &
        "C5D_gas        ", &
        "C5D+_gas       ", &
        "C5D-_gas       ", &
        "C5H2_gas       ", &
        "C5H2+_gas      ", &
        "C5H3_gas       ", &
        "C5H3+_gas      ", &
        "C5H3N+_gas     ", &
        "C5H4_gas       ", &
        "C5H4+_gas      ", &
        "C5H4N+_gas     ", &
        "C5H5+_gas      ", &
        "C5N_gas        ", &
        "C5N+_gas       ", &
        "C5O_gas        ", &
        "C6_gas         ", &
        "C6+_gas        ", &
        "C6-_gas        ", &
        "C6H_gas        ", &
        "C6H+_gas       ", &
        "C6H-_gas       ", &
        "C6H2_gas       ", &
        "C6H2+_gas      ", &
        "C6H3_gas       ", &
        "C6H3+_gas      ", &
        "C6H4_gas       ", &
        "C6H4+_gas      ", &
        "C6H5+_gas      ", &
        "C6H6_gas       ", &
        "C6H7+_gas      ", &
        "C6N_gas        ", &
        "C6N+_gas       ", &
        "C7_gas         ", &
        "C7+_gas        ", &
        "C7-_gas        ", &
        "C7H_gas        ", &
        "C7H+_gas       ", &
        "C7H-_gas       ", &
        "C7H2_gas       ", &
        "C7H2+_gas      ", &
        "C7H2N+_gas     ", &
        "C7H3_gas       ", &
        "C7H3+_gas      ", &
        "C7H4_gas       ", &
        "C7H4+_gas      ", &
        "C7H5+_gas      ", &
        "C7N_gas        ", &
        "C7N+_gas       ", &
        "C7O_gas        ", &
        "C8_gas         ", &
        "C8+_gas        ", &
        "C8-_gas        ", &
        "C8H_gas        ", &
        "C8H+_gas       ", &
        "C8H-_gas       ", &
        "C8H2_gas       ", &
        "C8H2+_gas      ", &
        "C8H3_gas       ", &
        "C8H3+_gas      ", &
        "C8H4_gas       ", &
        "C8H4+_gas      ", &
        "C8H4N+_gas     ", &
        "C8H5+_gas      ", &
        "C8N_gas        ", &
        "C8N+_gas       ", &
        "C9_gas         ", &
        "C9+_gas        ", &
        "C9-_gas        ", &
        "C9H_gas        ", &
        "C9H+_gas       ", &
        "C9H-_gas       ", &
        "C9H2_gas       ", &
        "C9H2+_gas      ", &
        "C9H2N+_gas     ", &
        "C9H3_gas       ", &
        "C9H3+_gas      ", &
        "C9H3N+_gas     ", &
        "C9H4_gas       ", &
        "C9H4+_gas      ", &
        "C9H5+_gas      ", &
        "C9HN+_gas      ", &
        "C9N_gas        ", &
        "C9N+_gas       ", &
        "C9O_gas        ", &
        "CCH_gas        ", &
        "CCD_gas        ", &
        "CCN_gas        ", &
        "CCO_gas        ", &
        "CCP_gas        ", &
        "CCP+_gas       ", &
        "CCS_gas        ", &
        "CD2ND2_gas     ", &
        "CD2ND2+_gas    ", &
        "CD2NH2_gas     ", &
        "CD2NH2+_gas    ", &
        "CD2NHD_gas     ", &
        "CD2NHD+_gas    ", &
        "CD3CN_gas      ", &
        "CD3CN+_gas     ", &
        "CD3CO+_gas     ", &
        "CD3O2+_gas     ", &
        "CD3OH+_gas     ", &
        "CD3OD+_gas     ", &
        "CH2CCH_gas     ", &
        "CH2CCD_gas     ", &
        "CD2CCH_gas     ", &
        "CD2CCD_gas     ", &
        "CHDCCH_gas     ", &
        "CHDCCD_gas     ", &
        "CH2CHC2H_gas   ", &
        "CH2CHCHCH2_gas ", &
        "CH2CHCN_gas    ", &
        "CH2CN+_gas     ", &
        "CHDCN+_gas     ", &
        "CD2CN+_gas     ", &
        "CH2DCN_gas     ", &
        "CH2DCN+_gas    ", &
        "CH2DCO+_gas    ", &
        "CH2DO2+_gas    ", &
        "CH2DOD+_gas    ", &
        "CH2ND2_gas     ", &
        "CH2ND2+_gas    ", &
        "CH2NH_gas      ", &
        "CH2ND_gas      ", &
        "CHDNH_gas      ", &
        "CHDND_gas      ", &
        "CD2NH_gas      ", &
        "CD2ND_gas      ", &
        "CH2NH2_gas     ", &
        "CH2NH2+_gas    ", &
        "CH2NHD_gas     ", &
        "CH2NHD+_gas    ", &
        "CH3C3N_gas     ", &
        "CH3C4H_gas     ", &
        "CH3C5N_gas     ", &
        "CH3C6H_gas     ", &
        "CH3C7N_gas     ", &
        "CH3CCH_gas     ", &
        "CH3CH2OH_gas   ", &
        "CH3CHCH2_gas   ", &
        "CH3CHO_gas     ", &
        "CH3CHOH+_gas   ", &
        "CH3CN_gas      ", &
        "CH3CN+_gas     ", &
        "CH3CNH+_gas    ", &
        "CH3CO+_gas     ", &
        "CH3COCH3_gas   ", &
        "CH3NH2_gas     ", &
        "CH3NH2+_gas    ", &
        "CH3NH3+_gas    ", &
        "CH3O2+_gas     ", &
        "CH3OCH2_gas    ", &
        "CH3OCH3_gas    ", &
        "CH3OCH3+_gas   ", &
        "CH3OCH4+_gas   ", &
        "CH3OH+_gas     ", &
        "CH2DOH+_gas    ", &
        "CH3OD+_gas     ", &
        "CH3OH2+_gas    ", &
        "CHD2CN_gas     ", &
        "CHD2CN+_gas    ", &
        "CHD2CO+_gas    ", &
        "CHD2O2+_gas    ", &
        "CHD2OH+_gas    ", &
        "CHD2OD+_gas    ", &
        "CHDNH2_gas     ", &
        "CHDNH2+_gas    ", &
        "CHDNHD_gas     ", &
        "CHDNHD+_gas    ", &
        "CHDND2_gas     ", &
        "CHDND2+_gas    ", &
        "CN_gas         ", &
        "CN+_gas        ", &
        "CN-_gas        ", &
        "CNC+_gas       ", &
        "COOCH4+_gas    ", &
        "CS_gas         ", &
        "CS+_gas        ", &
        "D2C3O+_gas     ", &
        "DC2NCH+_gas    ", &
        "DC2NCD+_gas    ", &
        "DCND+_gas      ", &
        "DCNH+_gas      ", &
        "DCS_gas        ", &
        "DCS+_gas       ", &
        "DN2O+_gas      ", &
        "DNC_gas        ", &
        "DNC+_gas       ", &
        "DNCCC_gas      ", &
        "DNCO_gas       ", &
        "DNCO+_gas      ", &
        "DNO_gas        ", &
        "DNO+_gas       ", &
        "DNS+_gas       ", &
        "DOC+_gas       ", &
        "DOCO+_gas      ", &
        "DOCS+_gas      ", &
        "HS_gas         ", &
        "HS+_gas        ", &
        "DS_gas         ", &
        "DS+_gas        ", &
        "FeH_gas        ", &
        "FeD_gas        ", &
        "H2C10N+_gas    ", &
        "H2C3O+_gas     ", &
        "H2C4N+_gas     ", &
        "H2C5N+_gas     ", &
        "H2C6N+_gas     ", &
        "H2C8N+_gas     ", &
        "H2CCN_gas      ", &
        "HDCCN_gas      ", &
        "D2CCN_gas      ", &
        "H2CCO_gas      ", &
        "H2CCO+_gas     ", &
        "HDCCO_gas      ", &
        "HDCCO+_gas     ", &
        "D2CCO_gas      ", &
        "D2CCO+_gas     ", &
        "H2CN_gas       ", &
        "H2CN+_gas      ", &
        "HDCN_gas       ", &
        "D2CN_gas       ", &
        "H2CO+_gas      ", &
        "HDCO+_gas      ", &
        "D2CO+_gas      ", &
        "H2COH+_gas     ", &
        "H2COD+_gas     ", &
        "HDCOH+_gas     ", &
        "HDCOD+_gas     ", &
        "D2COH+_gas     ", &
        "D2COD+_gas     ", &
        "H2CS_gas       ", &
        "H2CS+_gas      ", &
        "HDCS_gas       ", &
        "HDCS+_gas      ", &
        "D2CS_gas       ", &
        "D2CS+_gas      ", &
        "H2NC+_gas      ", &
        "HDNC+_gas      ", &
        "D2NC+_gas      ", &
        "H2NO+_gas      ", &
        "HDNO+_gas      ", &
        "D2NO+_gas      ", &
        "H2S_gas        ", &
        "H2S+_gas       ", &
        "HDS_gas        ", &
        "HDS+_gas       ", &
        "D2S_gas        ", &
        "D2S+_gas       ", &
        "H2S2+_gas      ", &
        "HDS2+_gas      ", &
        "D2S2+_gas      ", &
        "H3C3O+_gas     ", &
        "H3C4N+_gas     ", &
        "H3C4NH+_gas    ", &
        "H3C6NH+_gas    ", &
        "H3C7N+_gas     ", &
        "H3CS+_gas      ", &
        "H2DCS+_gas     ", &
        "HD2CS+_gas     ", &
        "D3CS+_gas      ", &
        "H3S+_gas       ", &
        "H2DS+_gas      ", &
        "HD2S+_gas      ", &
        "D3S+_gas       ", &
        "H3S2+_gas      ", &
        "H2DS2+_gas     ", &
        "HD2S2+_gas     ", &
        "D3S2+_gas      ", &
        "H5C2O2+_gas    ", &
        "HC10N+_gas     ", &
        "HC2N_gas       ", &
        "HC2N+_gas      ", &
        "DC2N_gas       ", &
        "DC2N+_gas      ", &
        "HC2NCD+_gas    ", &
        "HC2NCH+_gas    ", &
        "HC2O_gas       ", &
        "DC2O_gas       ", &
        "HC2S+_gas      ", &
        "DC2S+_gas      ", &
        "HC3N_gas       ", &
        "HC3N+_gas      ", &
        "DC3N_gas       ", &
        "DC3N+_gas      ", &
        "HC3NH+_gas     ", &
        "HC3ND+_gas     ", &
        "DC3NH+_gas     ", &
        "DC3ND+_gas     ", &
        "HC3O+_gas      ", &
        "DC3O+_gas      ", &
        "HC3S+_gas      ", &
        "DC3S+_gas      ", &
        "HC4N_gas       ", &
        "HC4N+_gas      ", &
        "DC4N_gas       ", &
        "DC4N+_gas      ", &
        "HC4O+_gas      ", &
        "DC4O+_gas      ", &
        "HC4S+_gas      ", &
        "DC4S+_gas      ", &
        "HC5N_gas       ", &
        "HC5N+_gas      ", &
        "HC5O+_gas      ", &
        "HC6N_gas       ", &
        "HC6N+_gas      ", &
        "HC7N_gas       ", &
        "HC7N+_gas      ", &
        "HC7O+_gas      ", &
        "HC8N_gas       ", &
        "HC8N+_gas      ", &
        "HC9N_gas       ", &
        "HC9O+_gas      ", &
        "HCCNC_gas      ", &
        "DCCNC_gas      ", &
        "HCN_gas        ", &
        "HCN+_gas       ", &
        "DCN_gas        ", &
        "DCN+_gas       ", &
        "HCNCC_gas      ", &
        "DCNCC_gas      ", &
        "HCNH+_gas      ", &
        "HCND+_gas      ", &
        "HCOOCH3_gas    ", &
        "HCOOH+_gas     ", &
        "HCOOD+_gas     ", &
        "DCOOH+_gas     ", &
        "DCOOD+_gas     ", &
        "HCS_gas        ", &
        "HCS+_gas       ", &
        "HDC3O+_gas     ", &
        "HN2O+_gas      ", &
        "HNC_gas        ", &
        "HNC+_gas       ", &
        "HNCCC_gas      ", &
        "HNCO_gas       ", &
        "HNCO+_gas      ", &
        "HNO_gas        ", &
        "HNO+_gas       ", &
        "HNS+_gas       ", &
        "HOC+_gas       ", &
        "HOCO+_gas      ", &
        "HOCS+_gas      ", &
        "HS2+_gas       ", &
        "DS2+_gas       ", &
        "HSO+_gas       ", &
        "DSO+_gas       ", &
        "HSO2+_gas      ", &
        "DSO2+_gas      ", &
        "HSS_gas        ", &
        "DSS_gas        ", &
        "HSSH_gas       ", &
        "HSSD_gas       ", &
        "DSSH_gas       ", &
        "DSSD_gas       ", &
        "N2O_gas        ", &
        "NC4N_gas       ", &
        "NC6N_gas       ", &
        "NC8N_gas       ", &
        "NCO+_gas       ", &
        "NH2CH2O+_gas   ", &
        "NH2CHO_gas     ", &
        "NH2CDO_gas     ", &
        "NHDCHO_gas     ", &
        "NHDCDO_gas     ", &
        "ND2CHO_gas     ", &
        "ND2CDO_gas     ", &
        "NH2CN_gas      ", &
        "NHDCN_gas      ", &
        "ND2CN_gas      ", &
        "NH2CND+_gas    ", &
        "NHDCND+_gas    ", &
        "ND2CND+_gas    ", &
        "NH2CNH+_gas    ", &
        "NHDCNH+_gas    ", &
        "ND2CNH+_gas    ", &
        "NO_gas         ", &
        "NO+_gas        ", &
        "NO2_gas        ", &
        "NO2+_gas       ", &
        "NS_gas         ", &
        "NS+_gas        ", &
        "OCN_gas        ", &
        "OCS_gas        ", &
        "OCS+_gas       ", &
        "S2_gas         ", &
        "S2+_gas        ", &
        "SO_gas         ", &
        "SO+_gas        ", &
        "SO2_gas        ", &
        "SO2+_gas       ", &
        "CH3OCHO_gas    ", &
        "Si_gas         ", &
        "Si+_gas        ", &
        "SiO_gas        ", &
        "SiO+_gas       ", &
        "SiO2_gas       ", &
        "SiS_gas        ", &
        "SiS+_gas       ", &
        "SiN_gas        ", &
        "SiN+_gas       ", &
        "SiC_gas        ", &
        "SiC+_gas       ", &
        "SiC2_gas       ", &
        "SiC2+_gas      ", &
        "SiH_gas        ", &
        "SiD_gas        ", &
        "SiH2_gas       ", &
        "SiH2+_gas      ", &
        "SiHD_gas       ", &
        "SiHD+_gas      ", &
        "SiD2_gas       ", &
        "SiD2+_gas      ", &
        "SiH3_gas       ", &
        "SiH4_gas       ", &
        "P_gas          ", &
        "P+_gas         ", &
        "PH_gas         ", &
        "PH+_gas        ", &
        "PD_gas         ", &
        "PD+_gas        ", &
        "PO_gas         ", &
        "PO+_gas        ", &
        "PH2_gas        ", &
        "PH2+_gas       ", &
        "PHD_gas        ", &
        "PHD+_gas       ", &
        "PD2_gas        ", &
        "PD2+_gas       ", &
        "C3P_gas        ", &
        "C4P_gas        ", &
        "C4P+_gas       ", &
        "CH2PH_gas      ", &
        "CH2PD_gas      ", &
        "CHDPH_gas      ", &
        "CHDPD_gas      ", &
        "CD2PH_gas      ", &
        "CD2PD_gas      ", &
        "CH2Si+_gas     ", &
        "CHDSi+_gas     ", &
        "CD2Si+_gas     ", &
        "CHSi+_gas      ", &
        "CDSi+_gas      ", &
        "CP_gas         ", &
        "CP+_gas        ", &
        "D2PO+_gas      ", &
        "D3SiO+_gas     ", &
        "DCP_gas        ", &
        "DCP+_gas       ", &
        "DF_gas         ", &
        "DF+_gas        ", &
        "DPO_gas        ", &
        "DPO+_gas       ", &
        "DSiNH+_gas     ", &
        "DSiND+_gas     ", &
        "H2CSiCH_gas    ", &
        "H2CSiCD_gas    ", &
        "HDCSiCD_gas    ", &
        "HDCSiCH_gas    ", &
        "D2CSiCH_gas    ", &
        "D2CSiCD_gas    ", &
        "H2DSiO+_gas    ", &
        "H2F+_gas       ", &
        "HDF+_gas       ", &
        "D2F+_gas       ", &
        "H2PO+_gas      ", &
        "H2SiO_gas      ", &
        "H2SiO+_gas     ", &
        "HDSiO_gas      ", &
        "HDSiO+_gas     ", &
        "D2SiO_gas      ", &
        "D2SiO+_gas     ", &
        "H3SiO+_gas     ", &
        "HCCP_gas       ", &
        "DCCP_gas       ", &
        "HCCSi_gas      ", &
        "DCCSi_gas      ", &
        "HCP_gas        ", &
        "HCP+_gas       ", &
        "HCSi_gas       ", &
        "DCSi_gas       ", &
        "HD2SiO+_gas    ", &
        "HDPO+_gas      ", &
        "HF_gas         ", &
        "HF+_gas        ", &
        "HNSi_gas       ", &
        "HNSi+_gas      ", &
        "DNSi_gas       ", &
        "DNSi+_gas      ", &
        "HPN+_gas       ", &
        "DPN+_gas       ", &
        "HPO_gas        ", &
        "HPO+_gas       ", &
        "HSiNH+_gas     ", &
        "HSiND+_gas     ", &
        "HSiO+_gas      ", &
        "DSiO+_gas      ", &
        "HSiO2+_gas     ", &
        "DSiO2+_gas     ", &
        "HSiS+_gas      ", &
        "DSiS+_gas      ", &
        "HeH+_gas       ", &
        "HeD+_gas       ", &
        "PC2H+_gas      ", &
        "PC2D+_gas      ", &
        "PC2H2+_gas     ", &
        "PC2HD+_gas     ", &
        "PC2D2+_gas     ", &
        "PC2H3+_gas     ", &
        "PC2H2D+_gas    ", &
        "PC2HD2+_gas    ", &
        "PC2D3+_gas     ", &
        "PC2H4+_gas     ", &
        "PC3H+_gas      ", &
        "PC3D+_gas      ", &
        "PC4H+_gas      ", &
        "PC4D+_gas      ", &
        "PC4H2+_gas     ", &
        "PCH2+_gas      ", &
        "PCHD+_gas      ", &
        "PCD2+_gas      ", &
        "PCH3+_gas      ", &
        "PCH2D+_gas     ", &
        "PCHD2+_gas     ", &
        "PCD3+_gas      ", &
        "PCH4+_gas      ", &
        "PCH3D+_gas     ", &
        "PCH2D2+_gas    ", &
        "PCHD3+_gas     ", &
        "PCD4+_gas      ", &
        "PH3+_gas       ", &
        "PH2D+_gas      ", &
        "PHD2+_gas      ", &
        "PD3+_gas       ", &
        "PN_gas         ", &
        "PN+_gas        ", &
        "PNH2+_gas      ", &
        "PNHD+_gas      ", &
        "PND2+_gas      ", &
        "PNH3+_gas      ", &
        "PNHD2+_gas     ", &
        "PNH2D+_gas     ", &
        "PND3+_gas      ", &
        "SiC2CH3_gas    ", &
        "SiC2H+_gas     ", &
        "SiC2D+_gas     ", &
        "SiC2H2+_gas    ", &
        "SiC2HD+_gas    ", &
        "SiC2D2+_gas    ", &
        "SiC2H3+_gas    ", &
        "SiC2H2D+_gas   ", &
        "SiC2HD2+_gas   ", &
        "SiC2D3+_gas    ", &
        "SiC3D_gas      ", &
        "SiC3D+_gas     ", &
        "SiC3H_gas      ", &
        "SiC3H+_gas     ", &
        "SiC3H2+_gas    ", &
        "SiC3HD+_gas    ", &
        "SiC3D2+_gas    ", &
        "SiC3H5_gas     ", &
        "SiC4_gas       ", &
        "SiC4+_gas      ", &
        "SiC4D_gas      ", &
        "SiC4D+_gas     ", &
        "SiC4H_gas      ", &
        "SiC4H+_gas     ", &
        "SiC6H_gas      ", &
        "SiC8H_gas      ", &
        "SiCD3_gas      ", &
        "SiCD3+_gas     ", &
        "SiCH2_gas      ", &
        "SiCHD_gas      ", &
        "SiCD2_gas      ", &
        "SiCH2D_gas     ", &
        "SiCH2D+_gas    ", &
        "SiCH3_gas      ", &
        "SiCH3+_gas     ", &
        "SiCH4+_gas     ", &
        "SiCH3D+_gas    ", &
        "SiCH2D2+_gas   ", &
        "SiCHD3+_gas    ", &
        "SiCD4+_gas     ", &
        "SiCHD2_gas     ", &
        "SiCHD2+_gas    ", &
        "SiD3_gas       ", &
        "SiD3+_gas      ", &
        "SiF+_gas       ", &
        "SiH+_gas       ", &
        "SiD+_gas       ", &
        "SiH2D_gas      ", &
        "SiH2D+_gas     ", &
        "SiH2D2_gas     ", &
        "SiH2D2+_gas    ", &
        "SiH3+_gas      ", &
        "SiH3D_gas      ", &
        "SiH3D+_gas     ", &
        "SiH4+_gas      ", &
        "SiD4_gas       ", &
        "SiD4+_gas      ", &
        "SiH5+_gas      ", &
        "SiH4D+_gas     ", &
        "SiH3D2+_gas    ", &
        "SiH2D3+_gas    ", &
        "SiHD4+_gas     ", &
        "SiD5+_gas      ", &
        "SiHD2_gas      ", &
        "SiHD2+_gas     ", &
        "SiHD3_gas      ", &
        "SiHD3+_gas     ", &
        "SiNC_gas       ", &
        "SiNC+_gas      ", &
        "SiNCH+_gas     ", &
        "SiNCD+_gas     ", &
        "c_DCCHSi_gas   ", &
        "c_DCCDSi_gas   ", &
        "c_HCCHSi_gas   ", &
        "c_HCCDSi_gas   ", &
        "c_SiC2_gas     ", &
        "l_SiC3_gas     ", &
        "l_SiC3+_gas    ", &
        "l_C3H_gas      ", &
        "l_C3D_gas      ", &
        "c_C3H_gas      ", &
        "c_C3D_gas      ", &
        "l_C3H2_gas     ", &
        "l_C3HD_gas     ", &
        "l_C3D2_gas     ", &
        "c_C3H2_gas     ", &
        "c_C3HD_gas     ", &
        "c_C3D2_gas     ", &
        "l_C3H2+_gas    ", &
        "l_C3HD+_gas    ", &
        "l_C3D2+_gas    ", &
        "c_C3H2+_gas    ", &
        "c_C3HD+_gas    ", &
        "c_C3D2+_gas    ", &
        "l_C3H3+_gas    ", &
        "l_C3H2D+_gas   ", &
        "l_C3HD2+_gas   ", &
        "l_C3D3+_gas    ", &
        "c_C3H3+_gas    ", &
        "c_C3H2D+_gas   ", &
        "c_C3HD2+_gas   ", &
        "c_C3D3+_gas    ", &
        "Mg_gas         ", &
        "Mg+_gas        ", &
        "MgH_gas        ", &
        "MgD_gas        ", &
        "MgH2_gas       ", &
        "MgHD_gas       ", &
        "MgD2_gas       ", &
        "HCl_gas        ", &
        "HCl+_gas       ", &
        "DCl_gas        ", &
        "DCl+_gas       ", &
        "H2Cl_gas       ", &
        "H2Cl+_gas      ", &
        "HDCl_gas       ", &
        "HDCl+_gas      ", &
        "D2Cl_gas       ", &
        "D2Cl+_gas      ", &
        "CCl_gas        ", &
        "CCl+_gas       ", &
        "ClO_gas        ", &
        "ClO+_gas       ", &
        "H2CCl+_gas     ", &
        "HDCCl+_gas     ", &
        "D2CCl+_gas     ", &
        "Na_gas         ", &
        "Na+_gas        ", &
        "NaH_gas        ", &
        "NaD_gas        ", &
        "NaOH_gas       ", &
        "NaOD_gas       ", &
        "NaH2+_gas      ", &
        "NaHD+_gas      ", &
        "NaD2+_gas      ", &
        "NaH2O+_gas     ", &
        "NaHDO+_gas     ", &
        "NaD2O+_gas     ", &
        "H_0001         ", &
        "D_0001         ", &
        "p_H2_0001      ", &
        "p_D2_0001      ", &
        "o_H2_0001      ", &
        "o_D2_0001      ", &
        "HD_0001        ", &
        "He_0001        ", &
        "O_0001         ", &
        "O2_0001        ", &
        "O3_0001        ", &
        "OH_0001        ", &
        "OD_0001        ", &
        "H2O_0001       ", &
        "HDO_0001       ", &
        "D2O_0001       ", &
        "O2H_0001       ", &
        "O2D_0001       ", &
        "HOOH_0001      ", &
        "HOOD_0001      ", &
        "DOOD_0001      ", &
        "Fe_0001        ", &
        "FeH_0001       ", &
        "N_0001         ", &
        "S_0001         ", &
        "C_0001         ", &
        "C2_0001        ", &
        "CO_0001        ", &
        "HCO_0001       ", &
        "DCO_0001       ", &
        "H2CO_0001      ", &
        "HDCO_0001      ", &
        "D2CO_0001      ", &
        "CH2OH_0001     ", &
        "CD2OD_0001     ", &
        "CH2OD_0001     ", &
        "CHDOH_0001     ", &
        "CHDOD_0001     ", &
        "CD2OH_0001     ", &
        "CH3O_0001      ", &
        "CHD2O_0001     ", &
        "CH2DO_0001     ", &
        "CD3O_0001      ", &
        "CH3OH_0001     ", &
        "CH3OD_0001     ", &
        "CHD2OH_0001    ", &
        "CHD2OD_0001    ", &
        "CH2DOH_0001    ", &
        "CH2DOD_0001    ", &
        "CD3OD_0001     ", &
        "CD3OH_0001     ", &
        "CH_0001        ", &
        "CD_0001        ", &
        "CH2_0001       ", &
        "CHD_0001       ", &
        "CD2_0001       ", &
        "CH3_0001       ", &
        "CH2D_0001      ", &
        "CHD2_0001      ", &
        "CD3_0001       ", &
        "CH4_0001       ", &
        "CH3D_0001      ", &
        "CH2D2_0001     ", &
        "CHD3_0001      ", &
        "CD4_0001       ", &
        "CO2_0001       ", &
        "HCOOH_0001     ", &
        "HCOOD_0001     ", &
        "DCOOH_0001     ", &
        "DCOOD_0001     ", &
        "HOCO_0001      ", &
        "DOCO_0001      ", &
        "NH_0001        ", &
        "ND_0001        ", &
        "NH2_0001       ", &
        "NHD_0001       ", &
        "ND2_0001       ", &
        "NH3_0001       ", &
        "NH2D_0001      ", &
        "NHD2_0001      ", &
        "ND3_0001       ", &
        "C10_0001       ", &
        "C10H_0001      ", &
        "C10H2_0001     ", &
        "C10N_0001      ", &
        "C11_0001       ", &
        "C2H2_0001      ", &
        "C2HD_0001      ", &
        "C2D2_0001      ", &
        "C2H3_0001      ", &
        "C2H2D_0001     ", &
        "C2HD2_0001     ", &
        "C2D3_0001      ", &
        "C2H4_0001      ", &
        "C2H3D_0001     ", &
        "C2H2D2_0001    ", &
        "C2HD3_0001     ", &
        "C2D4_0001      ", &
        "C2H5_0001      ", &
        "C2H6_0001      ", &
        "C3_0001        ", &
        "C3N_0001       ", &
        "C3O_0001       ", &
        "C3S_0001       ", &
        "C4_0001        ", &
        "C4H_0001       ", &
        "C4D_0001       ", &
        "C4H2_0001      ", &
        "C4HD_0001      ", &
        "C4D2_0001      ", &
        "C4H3_0001      ", &
        "C4N_0001       ", &
        "C4S_0001       ", &
        "C5_0001        ", &
        "C5H_0001       ", &
        "C5D_0001       ", &
        "C5H2_0001      ", &
        "C5H3_0001      ", &
        "C5H4_0001      ", &
        "C5N_0001       ", &
        "C5O_0001       ", &
        "C6_0001        ", &
        "C6H_0001       ", &
        "C6H2_0001      ", &
        "C6H3_0001      ", &
        "C6H4_0001      ", &
        "C6H6_0001      ", &
        "C6N_0001       ", &
        "C7_0001        ", &
        "C7H_0001       ", &
        "C7H2_0001      ", &
        "C7H3_0001      ", &
        "C7H4_0001      ", &
        "C7N_0001       ", &
        "C7O_0001       ", &
        "C8_0001        ", &
        "C8H_0001       ", &
        "C8H2_0001      ", &
        "C8H3_0001      ", &
        "C8H4_0001      ", &
        "C8N_0001       ", &
        "C9_0001        ", &
        "C9H_0001       ", &
        "C9H2_0001      ", &
        "C9H3_0001      ", &
        "C9H4_0001      ", &
        "C9N_0001       ", &
        "C9O_0001       ", &
        "CCH_0001       ", &
        "CCD_0001       ", &
        "CCN_0001       ", &
        "CCO_0001       ", &
        "CCS_0001       ", &
        "CD3CN_0001     ", &
        "CH2CCH_0001    ", &
        "CH2CCD_0001    ", &
        "CHDCCH_0001    ", &
        "CHDCCD_0001    ", &
        "CD2CCH_0001    ", &
        "CD2CCD_0001    ", &
        "CH2CHC2H_0001  ", &
        "CH2CHCHCH2_0001", &
        "CH2CHCN_0001   ", &
        "CH2NH_0001     ", &
        "CHDNH_0001     ", &
        "CHDND_0001     ", &
        "CH2ND_0001     ", &
        "CD2NH_0001     ", &
        "CD2ND_0001     ", &
        "CH2NH2_0001    ", &
        "CH2NHD_0001    ", &
        "CHDNH2_0001    ", &
        "CHDNHD_0001    ", &
        "CHDND2_0001    ", &
        "CH2ND2_0001    ", &
        "CD2NH2_0001    ", &
        "CD2NHD_0001    ", &
        "CD2ND2_0001    ", &
        "CH3C3N_0001    ", &
        "CH3C4H_0001    ", &
        "CH3C5N_0001    ", &
        "CH3C6H_0001    ", &
        "CH3C7N_0001    ", &
        "CH3CCH_0001    ", &
        "CH3CH2OH_0001  ", &
        "CH3CHCH2_0001  ", &
        "CH3CHO_0001    ", &
        "CH3CN_0001     ", &
        "CH2DCN_0001    ", &
        "CHD2CN_0001    ", &
        "CH3COCH3_0001  ", &
        "CH3NH2_0001    ", &
        "CH3OCH2_0001   ", &
        "CH3OCH3_0001   ", &
        "CN_0001        ", &
        "CS_0001        ", &
        "H2CCN_0001     ", &
        "HDCCN_0001     ", &
        "D2CCN_0001     ", &
        "H2CCO_0001     ", &
        "HDCCO_0001     ", &
        "D2CCO_0001     ", &
        "H2CN_0001      ", &
        "HDCN_0001      ", &
        "D2CN_0001      ", &
        "H2CS_0001      ", &
        "HDCS_0001      ", &
        "D2CS_0001      ", &
        "H2S_0001       ", &
        "HDS_0001       ", &
        "D2S_0001       ", &
        "HC2O_0001      ", &
        "DC2O_0001      ", &
        "HC3N_0001      ", &
        "DC3N_0001      ", &
        "HC4N_0001      ", &
        "DC4N_0001      ", &
        "HC5N_0001      ", &
        "HC6N_0001      ", &
        "HC7N_0001      ", &
        "HC8N_0001      ", &
        "HC9N_0001      ", &
        "HCCNC_0001     ", &
        "DCCNC_0001     ", &
        "HCN_0001       ", &
        "DCN_0001       ", &
        "HCNCC_0001     ", &
        "DCNCC_0001     ", &
        "HCOOCH3_0001   ", &
        "HCS_0001       ", &
        "DCS_0001       ", &
        "HNC_0001       ", &
        "DNC_0001       ", &
        "HNCCC_0001     ", &
        "DNCCC_0001     ", &
        "HNCO_0001      ", &
        "DNCO_0001      ", &
        "HNO_0001       ", &
        "DNO_0001       ", &
        "HS_0001        ", &
        "DS_0001        ", &
        "N2_0001        ", &
        "N2O_0001       ", &
        "NC4N_0001      ", &
        "NC6N_0001      ", &
        "NC8N_0001      ", &
        "NH2CHO_0001    ", &
        "NH2CDO_0001    ", &
        "NHDCHO_0001    ", &
        "NHDCDO_0001    ", &
        "ND2CHO_0001    ", &
        "ND2CDO_0001    ", &
        "NH2CN_0001     ", &
        "NHDCN_0001     ", &
        "ND2CN_0001     ", &
        "HSO_0001       ", &
        "DSO_0001       ", &
        "HSS_0001       ", &
        "DSS_0001       ", &
        "HSSH_0001      ", &
        "HSSD_0001      ", &
        "DSSH_0001      ", &
        "DSSD_0001      ", &
        "NO_0001        ", &
        "NO2_0001       ", &
        "NS_0001        ", &
        "OCN_0001       ", &
        "OCS_0001       ", &
        "S2_0001        ", &
        "SO_0001        ", &
        "SO2_0001       ", &
        "CH3OCHO_0001   ", &
        "Si_0001        ", &
        "SiS_0001       ", &
        "SiN_0001       ", &
        "SiC_0001       ", &
        "SiH_0001       ", &
        "SiH2_0001      ", &
        "SiH3_0001      ", &
        "SiH4_0001      ", &
        "SiC2CH3_0001   ", &
        "SiC3H_0001     ", &
        "SiC3H5_0001    ", &
        "SiC4_0001      ", &
        "SiC4H_0001     ", &
        "SiC6H_0001     ", &
        "SiC8H_0001     ", &
        "c_HCCHSi_0001  ", &
        "c_SiC2_0001    ", &
        "l_C3H_0001     ", &
        "l_C3D_0001     ", &
        "c_C3H_0001     ", &
        "c_C3D_0001     ", &
        "l_SiC3_0001    ", &
        "l_C3H2_0001    ", &
        "l_C3HD_0001    ", &
        "l_C3D2_0001    ", &
        "c_C3H2_0001    ", &
        "c_C3HD_0001    ", &
        "c_C3D2_0001    ", &
        "Mg_0001        ", &
        "MgH_0001       ", &
        "MgH2_0001      ", &
        "Na_0001        ", &
        "NaH_0001       ", &
        "F_0001         ", &
        "HF_0001        ", &
        "DF_0001        ", &
        "MgD_0001       ", &
        "MgHD_0001      ", &
        "MgD2_0001      ", &
        "NaD_0001       ", &
        "SiD_0001       ", &
        "Cl_0001        ", &
        "HCl_0001       ", &
        "DCl_0001       ", &
        "H2Cl_0001      ", &
        "HDCl_0001      ", &
        "D2Cl_0001      ", &
        "CCl_0001       ", &
        "ClO_0001       ", &
        "C3H3_0001      ", &
        "C3H2D_0001     ", &
        "C3HD2_0001     ", &
        "C3D3_0001      ", &
        "C3H4_0001      ", &
        "C3H3D_0001     ", &
        "C3H2D2_0001    ", &
        "HCCN_0001      ", &
        "DCCN_0001      ", &
        "C2H2N_0001     ", &
        "C2H5OH_0001    ", &
        "C2H3D3_0001    ", &
        "C2H5D_0001     ", &
        "C2H4D2_0001    ", &
        "C2H3N_0001     ", &
        "C2H4O_0001     ", &
        "CH2DCHO_0001   ", &
        "CH3CDO_0001    ", &
        "CD3CHO_0001    ", &
        "CHD2CDO_0001   ", &
        "CHD2CHO_0001   ", &
        "CH2DCDO_0001   ", &
        "C2H4D_0001     ", &
        "C2H2D3_0001    ", &
        "C2H3D2_0001    ", &
        "C3H3N_0001     ", &
        "H4C3N_0001     ", &
        "C3D3N_0001     ", &
        "HD3C3N_0001    ", &
        "C3H2DN_0001    ", &
        "H3DC3N_0001    ", &
        "C3HD2N_0001    ", &
        "H2D2C3N_0001   ", &
        "HC3O_0001      ", &
        "DC3O_0001      ", &
        "C4H4_0001      ", &
        "CH5N_0001      ", &
        "CH3NH_0001     ", &
        "FeD_0001       ", &
        "H2C3N_0001     ", &
        "H2C5N_0001     ", &
        "H3C5N_0001     ", &
        "H2C7N_0001     ", &
        "H3C7N_0001     ", &
        "H2C9N_0001     ", &
        "H3C9N_0001     ", &
        "H5C3N_0001     ", &
        "C2H2O_0001     ", &
        "C2HDO_0001     ", &
        "C2D2O_0001     ", &
        "HDC3N_0001     ", &
        "D2C3N_0001     ", &
        "H2C3O_0001     ", &
        "HDC3O_0001     ", &
        "D2C3O_0001     ", &
        "C2HDN_0001     ", &
        "C2D2N_0001     ", &
        "CH3N_0001      ", &
        "N2H2_0001      ", &
        "N2HD_0001      ", &
        "N2D2_0001      ", &
        "NH2OH_0001     ", &
        "NH2OD_0001     ", &
        "N2H_0001       ", &
        "N2D_0001       ", &
        "CNH2_0001      ", &
        "CHNH2_0001     ", &
        "HON_0001       ", &
        "DON_0001       ", &
        "NHNO_0001      ", &
        "CH2DN_0001     ", &
        "CHD2N_0001     ", &
        "CD3N_0001      ", &
        "NH2NO_0001     ", &
        "CH3CO_0001     ", &
        "CH2DOCH3_0001  ", &
        "CHD2OCH3_0001  ", &
        "CH2DOCH2D_0001 ", &
        "CD3OCH3_0001   ", &
        "CHD2OCH2D_0001 ", &
        "DCOOCH3_0001   ", &
        "HCOOCH2D_0001  ", &
        "DCOOCH2D_0001  ", &
        "HCOOCHD2_0001  ", &
        "DCOOCHD2_0001  ", &
        "HCOOCD3_0001   ", &
        "DCOOCD3_0001   ", &
        "NH2CO_0001     ", &
        "NHDCO_0001     ", &
        "ND2CO_0001     ", &
        "CH2DOCHO_0001  ", &
        "CH3OCDO_0001   ", &
        "CHD2OCHO_0001  ", &
        "CH2DOCDO_0001  ", &
        "CD3OCHO_0001   ", &
        "CHD2OCDO_0001  ", &
        "CHDOCDO_0001   ", &
        "CHDOCHDO_0001  ", &
        "CH3OCHD_0001   ", &
        "CH2DOCH2_0001  ", &
        "CH2DOCHD_0001  ", &
        "CHD2OCH2_0001  ", &
        "CH3OCD2_0001   ", &
        "CH2DOCD2_0001  ", &
        "CH3OCHD2_0001  ", &
        "CH3OCH2D_0001  ", &
        "CH2DCO_0001    ", &
        "CH3NHD_0001    ", &
        "CH2DNH2_0001   ", &
        "P_0001         ", &
        "PO_0001        ", &
        "PH_0001        ", &
        "PD_0001        ", &
        "PH2_0001       ", &
        "PHD_0001       ", &
        "PD2_0001       ", &
        "PN_0001        ", &
        "CP_0001        ", &
        "CCP_0001       ", &
        "C3P_0001       ", &
        "C4P_0001       ", &
        "CH2PH_0001     ", &
        "CHD2PH_0001    ", &
        "CH2PD_0001     ", &
        "CHDPD_0001     ", &
        "CD2P2_0001     ", &
        "HCP_0001       ", &
        "DCP_0001       ", &
        "HCCP_0001      ", &
        "DCCP_0001      ", &
        "HPO_0001       ", &
        "DPO_0001       ", &
        "surface_mask   ", &
        "mantle_mask    ", &
        "dummy          "/)
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_SPECIESNAMES

  !!BEGIN_IDXLIST
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-22 21:21:33
    ! CHANGESET: xxxxxxx
    ! BY: unknown@unknown

    integer,parameter::idx_GRAIN0_gas=1
    integer,parameter::idx_GRAINk_gas=2
    integer,parameter::idx_E_gas=3
    integer,parameter::idx_H_gas=4
    integer,parameter::idx_Hj_gas=5
    integer,parameter::idx_Hk_gas=6
    integer,parameter::idx_D_gas=7
    integer,parameter::idx_Dj_gas=8
    integer,parameter::idx_Dk_gas=9
    integer,parameter::idx_p_H2_gas=10
    integer,parameter::idx_o_H2_gas=11
    integer,parameter::idx_p_H2j_gas=12
    integer,parameter::idx_o_H2j_gas=13
    integer,parameter::idx_p_D2_gas=14
    integer,parameter::idx_o_D2_gas=15
    integer,parameter::idx_p_D2j_gas=16
    integer,parameter::idx_o_D2j_gas=17
    integer,parameter::idx_p_D3j_gas=18
    integer,parameter::idx_o_D3j_gas=19
    integer,parameter::idx_m_D3j_gas=20
    integer,parameter::idx_p_D2Hj_gas=21
    integer,parameter::idx_o_D2Hj_gas=22
    integer,parameter::idx_p_H2Dj_gas=23
    integer,parameter::idx_o_H2Dj_gas=24
    integer,parameter::idx_HD_gas=25
    integer,parameter::idx_HDj_gas=26
    integer,parameter::idx_He_gas=27
    integer,parameter::idx_Hej_gas=28
    integer,parameter::idx_O_gas=29
    integer,parameter::idx_Oj_gas=30
    integer,parameter::idx_Ok_gas=31
    integer,parameter::idx_O2_gas=32
    integer,parameter::idx_O2j_gas=33
    integer,parameter::idx_O3_gas=34
    integer,parameter::idx_OH_gas=35
    integer,parameter::idx_OHj_gas=36
    integer,parameter::idx_OHk_gas=37
    integer,parameter::idx_OD_gas=38
    integer,parameter::idx_ODj_gas=39
    integer,parameter::idx_ODk_gas=40
    integer,parameter::idx_H2O_gas=41
    integer,parameter::idx_H2Oj_gas=42
    integer,parameter::idx_HDO_gas=43
    integer,parameter::idx_HDOj_gas=44
    integer,parameter::idx_D2O_gas=45
    integer,parameter::idx_D2Oj_gas=46
    integer,parameter::idx_O2H_gas=47
    integer,parameter::idx_O2D_gas=48
    integer,parameter::idx_HO2j_gas=49
    integer,parameter::idx_DO2j_gas=50
    integer,parameter::idx_HOOH_gas=51
    integer,parameter::idx_HOOD_gas=52
    integer,parameter::idx_DOOH_gas=53
    integer,parameter::idx_DOOD_gas=54
    integer,parameter::idx_p_H3j_gas=55
    integer,parameter::idx_o_H3j_gas=56
    integer,parameter::idx_H3Oj_gas=57
    integer,parameter::idx_H2DOj_gas=58
    integer,parameter::idx_HD2Oj_gas=59
    integer,parameter::idx_D3Oj_gas=60
    integer,parameter::idx_F_gas=61
    integer,parameter::idx_Fj_gas=62
    integer,parameter::idx_Cl_gas=63
    integer,parameter::idx_Clj_gas=64
    integer,parameter::idx_CFj_gas=65
    integer,parameter::idx_C_gas=66
    integer,parameter::idx_Cj_gas=67
    integer,parameter::idx_Ck_gas=68
    integer,parameter::idx_C2_gas=69
    integer,parameter::idx_C2j_gas=70
    integer,parameter::idx_C2Nj_gas=71
    integer,parameter::idx_CO_gas=72
    integer,parameter::idx_COj_gas=73
    integer,parameter::idx_CO2_gas=74
    integer,parameter::idx_CO2j_gas=75
    integer,parameter::idx_N_gas=76
    integer,parameter::idx_Nj_gas=77
    integer,parameter::idx_N2_gas=78
    integer,parameter::idx_N2j_gas=79
    integer,parameter::idx_N2Hj_gas=80
    integer,parameter::idx_N2Dj_gas=81
    integer,parameter::idx_HCO_gas=82
    integer,parameter::idx_DCO_gas=83
    integer,parameter::idx_HCOj_gas=84
    integer,parameter::idx_DCOj_gas=85
    integer,parameter::idx_H2CO_gas=86
    integer,parameter::idx_HDCO_gas=87
    integer,parameter::idx_D2CO_gas=88
    integer,parameter::idx_CH2OH_gas=89
    integer,parameter::idx_CD2OD_gas=90
    integer,parameter::idx_CH2OD_gas=91
    integer,parameter::idx_CHDOH_gas=92
    integer,parameter::idx_CHDOD_gas=93
    integer,parameter::idx_CD2OH_gas=94
    integer,parameter::idx_CH3O_gas=95
    integer,parameter::idx_CHD2O_gas=96
    integer,parameter::idx_CH2DO_gas=97
    integer,parameter::idx_CD3O_gas=98
    integer,parameter::idx_CH3OH_gas=99
    integer,parameter::idx_CH3OD_gas=100
    integer,parameter::idx_CHD2OH_gas=101
    integer,parameter::idx_CHD2OD_gas=102
    integer,parameter::idx_CH2DOH_gas=103
    integer,parameter::idx_CH2DOD_gas=104
    integer,parameter::idx_CD3OD_gas=105
    integer,parameter::idx_CD3OH_gas=106
    integer,parameter::idx_CH_gas=107
    integer,parameter::idx_CD_gas=108
    integer,parameter::idx_CH2_gas=109
    integer,parameter::idx_CHD_gas=110
    integer,parameter::idx_CD2_gas=111
    integer,parameter::idx_CH3_gas=112
    integer,parameter::idx_CH2D_gas=113
    integer,parameter::idx_CHD2_gas=114
    integer,parameter::idx_CD3_gas=115
    integer,parameter::idx_CH4_gas=116
    integer,parameter::idx_CH3D_gas=117
    integer,parameter::idx_CH2D2_gas=118
    integer,parameter::idx_CHD3_gas=119
    integer,parameter::idx_CD4_gas=120
    integer,parameter::idx_CHj_gas=121
    integer,parameter::idx_CDj_gas=122
    integer,parameter::idx_CH2j_gas=123
    integer,parameter::idx_CHDj_gas=124
    integer,parameter::idx_CD2j_gas=125
    integer,parameter::idx_CH3j_gas=126
    integer,parameter::idx_CH2Dj_gas=127
    integer,parameter::idx_CHD2j_gas=128
    integer,parameter::idx_CD3j_gas=129
    integer,parameter::idx_CH4j_gas=130
    integer,parameter::idx_CH3Dj_gas=131
    integer,parameter::idx_CH2D2j_gas=132
    integer,parameter::idx_CHD3j_gas=133
    integer,parameter::idx_CD4j_gas=134
    integer,parameter::idx_CH5j_gas=135
    integer,parameter::idx_CH4Dj_gas=136
    integer,parameter::idx_CH3D2j_gas=137
    integer,parameter::idx_CH2D3j_gas=138
    integer,parameter::idx_CHD4j_gas=139
    integer,parameter::idx_CD5j_gas=140
    integer,parameter::idx_HCOOH_gas=141
    integer,parameter::idx_HCOOD_gas=142
    integer,parameter::idx_DCOOH_gas=143
    integer,parameter::idx_DCOOD_gas=144
    integer,parameter::idx_HOCO_gas=145
    integer,parameter::idx_DOCO_gas=146
    integer,parameter::idx_NH_gas=147
    integer,parameter::idx_NHj_gas=148
    integer,parameter::idx_ND_gas=149
    integer,parameter::idx_NDj_gas=150
    integer,parameter::idx_NH2_gas=151
    integer,parameter::idx_NH2j_gas=152
    integer,parameter::idx_NHD_gas=153
    integer,parameter::idx_NHDj_gas=154
    integer,parameter::idx_ND2_gas=155
    integer,parameter::idx_ND2j_gas=156
    integer,parameter::idx_NH3_gas=157
    integer,parameter::idx_NH3j_gas=158
    integer,parameter::idx_NH2D_gas=159
    integer,parameter::idx_NH2Dj_gas=160
    integer,parameter::idx_NHD2_gas=161
    integer,parameter::idx_NHD2j_gas=162
    integer,parameter::idx_ND3_gas=163
    integer,parameter::idx_ND3j_gas=164
    integer,parameter::idx_NH4j_gas=165
    integer,parameter::idx_NH3Dj_gas=166
    integer,parameter::idx_NH2D2j_gas=167
    integer,parameter::idx_NHD3j_gas=168
    integer,parameter::idx_ND4j_gas=169
    integer,parameter::idx_Fe_gas=170
    integer,parameter::idx_Fej_gas=171
    integer,parameter::idx_S_gas=172
    integer,parameter::idx_Sj_gas=173
    integer,parameter::idx_Sk_gas=174
    integer,parameter::idx_C10_gas=175
    integer,parameter::idx_C10j_gas=176
    integer,parameter::idx_C10k_gas=177
    integer,parameter::idx_C10H_gas=178
    integer,parameter::idx_C10Hj_gas=179
    integer,parameter::idx_C10Hk_gas=180
    integer,parameter::idx_C10H2_gas=181
    integer,parameter::idx_C10H2j_gas=182
    integer,parameter::idx_C10H3j_gas=183
    integer,parameter::idx_C10N_gas=184
    integer,parameter::idx_C10Nj_gas=185
    integer,parameter::idx_C11_gas=186
    integer,parameter::idx_C11j_gas=187
    integer,parameter::idx_C2Hj_gas=188
    integer,parameter::idx_C2Dj_gas=189
    integer,parameter::idx_C2H2_gas=190
    integer,parameter::idx_C2H2j_gas=191
    integer,parameter::idx_C2HD_gas=192
    integer,parameter::idx_C2HDj_gas=193
    integer,parameter::idx_C2D2_gas=194
    integer,parameter::idx_C2D2j_gas=195
    integer,parameter::idx_C2DOj_gas=196
    integer,parameter::idx_C2H3_gas=197
    integer,parameter::idx_C2H3j_gas=198
    integer,parameter::idx_C2H2D_gas=199
    integer,parameter::idx_C2H2Dj_gas=200
    integer,parameter::idx_C2HD2_gas=201
    integer,parameter::idx_C2HD2j_gas=202
    integer,parameter::idx_C2D3_gas=203
    integer,parameter::idx_C2D3j_gas=204
    integer,parameter::idx_C2H4_gas=205
    integer,parameter::idx_C2H4j_gas=206
    integer,parameter::idx_C2H3D_gas=207
    integer,parameter::idx_C2H3Dj_gas=208
    integer,parameter::idx_C2H2D2_gas=209
    integer,parameter::idx_C2H2D2j_gas=210
    integer,parameter::idx_C2HD3_gas=211
    integer,parameter::idx_C2HD3j_gas=212
    integer,parameter::idx_C2D4_gas=213
    integer,parameter::idx_C2D4j_gas=214
    integer,parameter::idx_C2H4Oj_gas=215
    integer,parameter::idx_C2H5_gas=216
    integer,parameter::idx_C2H5j_gas=217
    integer,parameter::idx_C2H5OHj_gas=218
    integer,parameter::idx_C2H5OH2j_gas=219
    integer,parameter::idx_C2H6_gas=220
    integer,parameter::idx_C2H6j_gas=221
    integer,parameter::idx_C2H6COj_gas=222
    integer,parameter::idx_C2H7j_gas=223
    integer,parameter::idx_C2HOj_gas=224
    integer,parameter::idx_C2N2j_gas=225
    integer,parameter::idx_C2Oj_gas=226
    integer,parameter::idx_C2Sj_gas=227
    integer,parameter::idx_C3_gas=228
    integer,parameter::idx_C3j_gas=229
    integer,parameter::idx_C3k_gas=230
    integer,parameter::idx_C3Hj_gas=231
    integer,parameter::idx_C3Dj_gas=232
    integer,parameter::idx_C3H3Nj_gas=233
    integer,parameter::idx_C3H3NHj_gas=234
    integer,parameter::idx_C3H4_gas=235
    integer,parameter::idx_C3H4j_gas=236
    integer,parameter::idx_C3H5j_gas=237
    integer,parameter::idx_C3H6OHj_gas=238
    integer,parameter::idx_C3N_gas=239
    integer,parameter::idx_C3Nj_gas=240
    integer,parameter::idx_C3Nk_gas=241
    integer,parameter::idx_C3O_gas=242
    integer,parameter::idx_C3Oj_gas=243
    integer,parameter::idx_C3S_gas=244
    integer,parameter::idx_C3Sj_gas=245
    integer,parameter::idx_C4_gas=246
    integer,parameter::idx_C4j_gas=247
    integer,parameter::idx_C4k_gas=248
    integer,parameter::idx_C4H_gas=249
    integer,parameter::idx_C4Hj_gas=250
    integer,parameter::idx_C4Hk_gas=251
    integer,parameter::idx_C4D_gas=252
    integer,parameter::idx_C4Dj_gas=253
    integer,parameter::idx_C4Dk_gas=254
    integer,parameter::idx_C4H2_gas=255
    integer,parameter::idx_C4H2j_gas=256
    integer,parameter::idx_C4HD_gas=257
    integer,parameter::idx_C4HDj_gas=258
    integer,parameter::idx_C4D2_gas=259
    integer,parameter::idx_C4D2j_gas=260
    integer,parameter::idx_C4H3_gas=261
    integer,parameter::idx_C4H3j_gas=262
    integer,parameter::idx_C4H4j_gas=263
    integer,parameter::idx_C4H5j_gas=264
    integer,parameter::idx_C4H7j_gas=265
    integer,parameter::idx_C4N_gas=266
    integer,parameter::idx_C4Nj_gas=267
    integer,parameter::idx_C4S_gas=268
    integer,parameter::idx_C4Sj_gas=269
    integer,parameter::idx_C5_gas=270
    integer,parameter::idx_C5j_gas=271
    integer,parameter::idx_C5k_gas=272
    integer,parameter::idx_C5H_gas=273
    integer,parameter::idx_C5Hj_gas=274
    integer,parameter::idx_C5Hk_gas=275
    integer,parameter::idx_C5D_gas=276
    integer,parameter::idx_C5Dj_gas=277
    integer,parameter::idx_C5Dk_gas=278
    integer,parameter::idx_C5H2_gas=279
    integer,parameter::idx_C5H2j_gas=280
    integer,parameter::idx_C5H3_gas=281
    integer,parameter::idx_C5H3j_gas=282
    integer,parameter::idx_C5H3Nj_gas=283
    integer,parameter::idx_C5H4_gas=284
    integer,parameter::idx_C5H4j_gas=285
    integer,parameter::idx_C5H4Nj_gas=286
    integer,parameter::idx_C5H5j_gas=287
    integer,parameter::idx_C5N_gas=288
    integer,parameter::idx_C5Nj_gas=289
    integer,parameter::idx_C5O_gas=290
    integer,parameter::idx_C6_gas=291
    integer,parameter::idx_C6j_gas=292
    integer,parameter::idx_C6k_gas=293
    integer,parameter::idx_C6H_gas=294
    integer,parameter::idx_C6Hj_gas=295
    integer,parameter::idx_C6Hk_gas=296
    integer,parameter::idx_C6H2_gas=297
    integer,parameter::idx_C6H2j_gas=298
    integer,parameter::idx_C6H3_gas=299
    integer,parameter::idx_C6H3j_gas=300
    integer,parameter::idx_C6H4_gas=301
    integer,parameter::idx_C6H4j_gas=302
    integer,parameter::idx_C6H5j_gas=303
    integer,parameter::idx_C6H6_gas=304
    integer,parameter::idx_C6H7j_gas=305
    integer,parameter::idx_C6N_gas=306
    integer,parameter::idx_C6Nj_gas=307
    integer,parameter::idx_C7_gas=308
    integer,parameter::idx_C7j_gas=309
    integer,parameter::idx_C7k_gas=310
    integer,parameter::idx_C7H_gas=311
    integer,parameter::idx_C7Hj_gas=312
    integer,parameter::idx_C7Hk_gas=313
    integer,parameter::idx_C7H2_gas=314
    integer,parameter::idx_C7H2j_gas=315
    integer,parameter::idx_C7H2Nj_gas=316
    integer,parameter::idx_C7H3_gas=317
    integer,parameter::idx_C7H3j_gas=318
    integer,parameter::idx_C7H4_gas=319
    integer,parameter::idx_C7H4j_gas=320
    integer,parameter::idx_C7H5j_gas=321
    integer,parameter::idx_C7N_gas=322
    integer,parameter::idx_C7Nj_gas=323
    integer,parameter::idx_C7O_gas=324
    integer,parameter::idx_C8_gas=325
    integer,parameter::idx_C8j_gas=326
    integer,parameter::idx_C8k_gas=327
    integer,parameter::idx_C8H_gas=328
    integer,parameter::idx_C8Hj_gas=329
    integer,parameter::idx_C8Hk_gas=330
    integer,parameter::idx_C8H2_gas=331
    integer,parameter::idx_C8H2j_gas=332
    integer,parameter::idx_C8H3_gas=333
    integer,parameter::idx_C8H3j_gas=334
    integer,parameter::idx_C8H4_gas=335
    integer,parameter::idx_C8H4j_gas=336
    integer,parameter::idx_C8H4Nj_gas=337
    integer,parameter::idx_C8H5j_gas=338
    integer,parameter::idx_C8N_gas=339
    integer,parameter::idx_C8Nj_gas=340
    integer,parameter::idx_C9_gas=341
    integer,parameter::idx_C9j_gas=342
    integer,parameter::idx_C9k_gas=343
    integer,parameter::idx_C9H_gas=344
    integer,parameter::idx_C9Hj_gas=345
    integer,parameter::idx_C9Hk_gas=346
    integer,parameter::idx_C9H2_gas=347
    integer,parameter::idx_C9H2j_gas=348
    integer,parameter::idx_C9H2Nj_gas=349
    integer,parameter::idx_C9H3_gas=350
    integer,parameter::idx_C9H3j_gas=351
    integer,parameter::idx_C9H3Nj_gas=352
    integer,parameter::idx_C9H4_gas=353
    integer,parameter::idx_C9H4j_gas=354
    integer,parameter::idx_C9H5j_gas=355
    integer,parameter::idx_C9HNj_gas=356
    integer,parameter::idx_C9N_gas=357
    integer,parameter::idx_C9Nj_gas=358
    integer,parameter::idx_C9O_gas=359
    integer,parameter::idx_CCH_gas=360
    integer,parameter::idx_CCD_gas=361
    integer,parameter::idx_CCN_gas=362
    integer,parameter::idx_CCO_gas=363
    integer,parameter::idx_CCP_gas=364
    integer,parameter::idx_CCPj_gas=365
    integer,parameter::idx_CCS_gas=366
    integer,parameter::idx_CD2ND2_gas=367
    integer,parameter::idx_CD2ND2j_gas=368
    integer,parameter::idx_CD2NH2_gas=369
    integer,parameter::idx_CD2NH2j_gas=370
    integer,parameter::idx_CD2NHD_gas=371
    integer,parameter::idx_CD2NHDj_gas=372
    integer,parameter::idx_CD3CN_gas=373
    integer,parameter::idx_CD3CNj_gas=374
    integer,parameter::idx_CD3COj_gas=375
    integer,parameter::idx_CD3O2j_gas=376
    integer,parameter::idx_CD3OHj_gas=377
    integer,parameter::idx_CD3ODj_gas=378
    integer,parameter::idx_CH2CCH_gas=379
    integer,parameter::idx_CH2CCD_gas=380
    integer,parameter::idx_CD2CCH_gas=381
    integer,parameter::idx_CD2CCD_gas=382
    integer,parameter::idx_CHDCCH_gas=383
    integer,parameter::idx_CHDCCD_gas=384
    integer,parameter::idx_CH2CHC2H_gas=385
    integer,parameter::idx_CH2CHCHCH2_gas=386
    integer,parameter::idx_CH2CHCN_gas=387
    integer,parameter::idx_CH2CNj_gas=388
    integer,parameter::idx_CHDCNj_gas=389
    integer,parameter::idx_CD2CNj_gas=390
    integer,parameter::idx_CH2DCN_gas=391
    integer,parameter::idx_CH2DCNj_gas=392
    integer,parameter::idx_CH2DCOj_gas=393
    integer,parameter::idx_CH2DO2j_gas=394
    integer,parameter::idx_CH2DODj_gas=395
    integer,parameter::idx_CH2ND2_gas=396
    integer,parameter::idx_CH2ND2j_gas=397
    integer,parameter::idx_CH2NH_gas=398
    integer,parameter::idx_CH2ND_gas=399
    integer,parameter::idx_CHDNH_gas=400
    integer,parameter::idx_CHDND_gas=401
    integer,parameter::idx_CD2NH_gas=402
    integer,parameter::idx_CD2ND_gas=403
    integer,parameter::idx_CH2NH2_gas=404
    integer,parameter::idx_CH2NH2j_gas=405
    integer,parameter::idx_CH2NHD_gas=406
    integer,parameter::idx_CH2NHDj_gas=407
    integer,parameter::idx_CH3C3N_gas=408
    integer,parameter::idx_CH3C4H_gas=409
    integer,parameter::idx_CH3C5N_gas=410
    integer,parameter::idx_CH3C6H_gas=411
    integer,parameter::idx_CH3C7N_gas=412
    integer,parameter::idx_CH3CCH_gas=413
    integer,parameter::idx_CH3CH2OH_gas=414
    integer,parameter::idx_CH3CHCH2_gas=415
    integer,parameter::idx_CH3CHO_gas=416
    integer,parameter::idx_CH3CHOHj_gas=417
    integer,parameter::idx_CH3CN_gas=418
    integer,parameter::idx_CH3CNj_gas=419
    integer,parameter::idx_CH3CNHj_gas=420
    integer,parameter::idx_CH3COj_gas=421
    integer,parameter::idx_CH3COCH3_gas=422
    integer,parameter::idx_CH3NH2_gas=423
    integer,parameter::idx_CH3NH2j_gas=424
    integer,parameter::idx_CH3NH3j_gas=425
    integer,parameter::idx_CH3O2j_gas=426
    integer,parameter::idx_CH3OCH2_gas=427
    integer,parameter::idx_CH3OCH3_gas=428
    integer,parameter::idx_CH3OCH3j_gas=429
    integer,parameter::idx_CH3OCH4j_gas=430
    integer,parameter::idx_CH3OHj_gas=431
    integer,parameter::idx_CH2DOHj_gas=432
    integer,parameter::idx_CH3ODj_gas=433
    integer,parameter::idx_CH3OH2j_gas=434
    integer,parameter::idx_CHD2CN_gas=435
    integer,parameter::idx_CHD2CNj_gas=436
    integer,parameter::idx_CHD2COj_gas=437
    integer,parameter::idx_CHD2O2j_gas=438
    integer,parameter::idx_CHD2OHj_gas=439
    integer,parameter::idx_CHD2ODj_gas=440
    integer,parameter::idx_CHDNH2_gas=441
    integer,parameter::idx_CHDNH2j_gas=442
    integer,parameter::idx_CHDNHD_gas=443
    integer,parameter::idx_CHDNHDj_gas=444
    integer,parameter::idx_CHDND2_gas=445
    integer,parameter::idx_CHDND2j_gas=446
    integer,parameter::idx_CN_gas=447
    integer,parameter::idx_CNj_gas=448
    integer,parameter::idx_CNk_gas=449
    integer,parameter::idx_CNCj_gas=450
    integer,parameter::idx_COOCH4j_gas=451
    integer,parameter::idx_CS_gas=452
    integer,parameter::idx_CSj_gas=453
    integer,parameter::idx_D2C3Oj_gas=454
    integer,parameter::idx_DC2NCHj_gas=455
    integer,parameter::idx_DC2NCDj_gas=456
    integer,parameter::idx_DCNDj_gas=457
    integer,parameter::idx_DCNHj_gas=458
    integer,parameter::idx_DCS_gas=459
    integer,parameter::idx_DCSj_gas=460
    integer,parameter::idx_DN2Oj_gas=461
    integer,parameter::idx_DNC_gas=462
    integer,parameter::idx_DNCj_gas=463
    integer,parameter::idx_DNCCC_gas=464
    integer,parameter::idx_DNCO_gas=465
    integer,parameter::idx_DNCOj_gas=466
    integer,parameter::idx_DNO_gas=467
    integer,parameter::idx_DNOj_gas=468
    integer,parameter::idx_DNSj_gas=469
    integer,parameter::idx_DOCj_gas=470
    integer,parameter::idx_DOCOj_gas=471
    integer,parameter::idx_DOCSj_gas=472
    integer,parameter::idx_HS_gas=473
    integer,parameter::idx_HSj_gas=474
    integer,parameter::idx_DS_gas=475
    integer,parameter::idx_DSj_gas=476
    integer,parameter::idx_FeH_gas=477
    integer,parameter::idx_FeD_gas=478
    integer,parameter::idx_H2C10Nj_gas=479
    integer,parameter::idx_H2C3Oj_gas=480
    integer,parameter::idx_H2C4Nj_gas=481
    integer,parameter::idx_H2C5Nj_gas=482
    integer,parameter::idx_H2C6Nj_gas=483
    integer,parameter::idx_H2C8Nj_gas=484
    integer,parameter::idx_H2CCN_gas=485
    integer,parameter::idx_HDCCN_gas=486
    integer,parameter::idx_D2CCN_gas=487
    integer,parameter::idx_H2CCO_gas=488
    integer,parameter::idx_H2CCOj_gas=489
    integer,parameter::idx_HDCCO_gas=490
    integer,parameter::idx_HDCCOj_gas=491
    integer,parameter::idx_D2CCO_gas=492
    integer,parameter::idx_D2CCOj_gas=493
    integer,parameter::idx_H2CN_gas=494
    integer,parameter::idx_H2CNj_gas=495
    integer,parameter::idx_HDCN_gas=496
    integer,parameter::idx_D2CN_gas=497
    integer,parameter::idx_H2COj_gas=498
    integer,parameter::idx_HDCOj_gas=499
    integer,parameter::idx_D2COj_gas=500
    integer,parameter::idx_H2COHj_gas=501
    integer,parameter::idx_H2CODj_gas=502
    integer,parameter::idx_HDCOHj_gas=503
    integer,parameter::idx_HDCODj_gas=504
    integer,parameter::idx_D2COHj_gas=505
    integer,parameter::idx_D2CODj_gas=506
    integer,parameter::idx_H2CS_gas=507
    integer,parameter::idx_H2CSj_gas=508
    integer,parameter::idx_HDCS_gas=509
    integer,parameter::idx_HDCSj_gas=510
    integer,parameter::idx_D2CS_gas=511
    integer,parameter::idx_D2CSj_gas=512
    integer,parameter::idx_H2NCj_gas=513
    integer,parameter::idx_HDNCj_gas=514
    integer,parameter::idx_D2NCj_gas=515
    integer,parameter::idx_H2NOj_gas=516
    integer,parameter::idx_HDNOj_gas=517
    integer,parameter::idx_D2NOj_gas=518
    integer,parameter::idx_H2S_gas=519
    integer,parameter::idx_H2Sj_gas=520
    integer,parameter::idx_HDS_gas=521
    integer,parameter::idx_HDSj_gas=522
    integer,parameter::idx_D2S_gas=523
    integer,parameter::idx_D2Sj_gas=524
    integer,parameter::idx_H2S2j_gas=525
    integer,parameter::idx_HDS2j_gas=526
    integer,parameter::idx_D2S2j_gas=527
    integer,parameter::idx_H3C3Oj_gas=528
    integer,parameter::idx_H3C4Nj_gas=529
    integer,parameter::idx_H3C4NHj_gas=530
    integer,parameter::idx_H3C6NHj_gas=531
    integer,parameter::idx_H3C7Nj_gas=532
    integer,parameter::idx_H3CSj_gas=533
    integer,parameter::idx_H2DCSj_gas=534
    integer,parameter::idx_HD2CSj_gas=535
    integer,parameter::idx_D3CSj_gas=536
    integer,parameter::idx_H3Sj_gas=537
    integer,parameter::idx_H2DSj_gas=538
    integer,parameter::idx_HD2Sj_gas=539
    integer,parameter::idx_D3Sj_gas=540
    integer,parameter::idx_H3S2j_gas=541
    integer,parameter::idx_H2DS2j_gas=542
    integer,parameter::idx_HD2S2j_gas=543
    integer,parameter::idx_D3S2j_gas=544
    integer,parameter::idx_H5C2O2j_gas=545
    integer,parameter::idx_HC10Nj_gas=546
    integer,parameter::idx_HC2N_gas=547
    integer,parameter::idx_HC2Nj_gas=548
    integer,parameter::idx_DC2N_gas=549
    integer,parameter::idx_DC2Nj_gas=550
    integer,parameter::idx_HC2NCDj_gas=551
    integer,parameter::idx_HC2NCHj_gas=552
    integer,parameter::idx_HC2O_gas=553
    integer,parameter::idx_DC2O_gas=554
    integer,parameter::idx_HC2Sj_gas=555
    integer,parameter::idx_DC2Sj_gas=556
    integer,parameter::idx_HC3N_gas=557
    integer,parameter::idx_HC3Nj_gas=558
    integer,parameter::idx_DC3N_gas=559
    integer,parameter::idx_DC3Nj_gas=560
    integer,parameter::idx_HC3NHj_gas=561
    integer,parameter::idx_HC3NDj_gas=562
    integer,parameter::idx_DC3NHj_gas=563
    integer,parameter::idx_DC3NDj_gas=564
    integer,parameter::idx_HC3Oj_gas=565
    integer,parameter::idx_DC3Oj_gas=566
    integer,parameter::idx_HC3Sj_gas=567
    integer,parameter::idx_DC3Sj_gas=568
    integer,parameter::idx_HC4N_gas=569
    integer,parameter::idx_HC4Nj_gas=570
    integer,parameter::idx_DC4N_gas=571
    integer,parameter::idx_DC4Nj_gas=572
    integer,parameter::idx_HC4Oj_gas=573
    integer,parameter::idx_DC4Oj_gas=574
    integer,parameter::idx_HC4Sj_gas=575
    integer,parameter::idx_DC4Sj_gas=576
    integer,parameter::idx_HC5N_gas=577
    integer,parameter::idx_HC5Nj_gas=578
    integer,parameter::idx_HC5Oj_gas=579
    integer,parameter::idx_HC6N_gas=580
    integer,parameter::idx_HC6Nj_gas=581
    integer,parameter::idx_HC7N_gas=582
    integer,parameter::idx_HC7Nj_gas=583
    integer,parameter::idx_HC7Oj_gas=584
    integer,parameter::idx_HC8N_gas=585
    integer,parameter::idx_HC8Nj_gas=586
    integer,parameter::idx_HC9N_gas=587
    integer,parameter::idx_HC9Oj_gas=588
    integer,parameter::idx_HCCNC_gas=589
    integer,parameter::idx_DCCNC_gas=590
    integer,parameter::idx_HCN_gas=591
    integer,parameter::idx_HCNj_gas=592
    integer,parameter::idx_DCN_gas=593
    integer,parameter::idx_DCNj_gas=594
    integer,parameter::idx_HCNCC_gas=595
    integer,parameter::idx_DCNCC_gas=596
    integer,parameter::idx_HCNHj_gas=597
    integer,parameter::idx_HCNDj_gas=598
    integer,parameter::idx_HCOOCH3_gas=599
    integer,parameter::idx_HCOOHj_gas=600
    integer,parameter::idx_HCOODj_gas=601
    integer,parameter::idx_DCOOHj_gas=602
    integer,parameter::idx_DCOODj_gas=603
    integer,parameter::idx_HCS_gas=604
    integer,parameter::idx_HCSj_gas=605
    integer,parameter::idx_HDC3Oj_gas=606
    integer,parameter::idx_HN2Oj_gas=607
    integer,parameter::idx_HNC_gas=608
    integer,parameter::idx_HNCj_gas=609
    integer,parameter::idx_HNCCC_gas=610
    integer,parameter::idx_HNCO_gas=611
    integer,parameter::idx_HNCOj_gas=612
    integer,parameter::idx_HNO_gas=613
    integer,parameter::idx_HNOj_gas=614
    integer,parameter::idx_HNSj_gas=615
    integer,parameter::idx_HOCj_gas=616
    integer,parameter::idx_HOCOj_gas=617
    integer,parameter::idx_HOCSj_gas=618
    integer,parameter::idx_HS2j_gas=619
    integer,parameter::idx_DS2j_gas=620
    integer,parameter::idx_HSOj_gas=621
    integer,parameter::idx_DSOj_gas=622
    integer,parameter::idx_HSO2j_gas=623
    integer,parameter::idx_DSO2j_gas=624
    integer,parameter::idx_HSS_gas=625
    integer,parameter::idx_DSS_gas=626
    integer,parameter::idx_HSSH_gas=627
    integer,parameter::idx_HSSD_gas=628
    integer,parameter::idx_DSSH_gas=629
    integer,parameter::idx_DSSD_gas=630
    integer,parameter::idx_N2O_gas=631
    integer,parameter::idx_NC4N_gas=632
    integer,parameter::idx_NC6N_gas=633
    integer,parameter::idx_NC8N_gas=634
    integer,parameter::idx_NCOj_gas=635
    integer,parameter::idx_NH2CH2Oj_gas=636
    integer,parameter::idx_NH2CHO_gas=637
    integer,parameter::idx_NH2CDO_gas=638
    integer,parameter::idx_NHDCHO_gas=639
    integer,parameter::idx_NHDCDO_gas=640
    integer,parameter::idx_ND2CHO_gas=641
    integer,parameter::idx_ND2CDO_gas=642
    integer,parameter::idx_NH2CN_gas=643
    integer,parameter::idx_NHDCN_gas=644
    integer,parameter::idx_ND2CN_gas=645
    integer,parameter::idx_NH2CNDj_gas=646
    integer,parameter::idx_NHDCNDj_gas=647
    integer,parameter::idx_ND2CNDj_gas=648
    integer,parameter::idx_NH2CNHj_gas=649
    integer,parameter::idx_NHDCNHj_gas=650
    integer,parameter::idx_ND2CNHj_gas=651
    integer,parameter::idx_NO_gas=652
    integer,parameter::idx_NOj_gas=653
    integer,parameter::idx_NO2_gas=654
    integer,parameter::idx_NO2j_gas=655
    integer,parameter::idx_NS_gas=656
    integer,parameter::idx_NSj_gas=657
    integer,parameter::idx_OCN_gas=658
    integer,parameter::idx_OCS_gas=659
    integer,parameter::idx_OCSj_gas=660
    integer,parameter::idx_S2_gas=661
    integer,parameter::idx_S2j_gas=662
    integer,parameter::idx_SO_gas=663
    integer,parameter::idx_SOj_gas=664
    integer,parameter::idx_SO2_gas=665
    integer,parameter::idx_SO2j_gas=666
    integer,parameter::idx_CH3OCHO_gas=667
    integer,parameter::idx_Si_gas=668
    integer,parameter::idx_Sij_gas=669
    integer,parameter::idx_SiO_gas=670
    integer,parameter::idx_SiOj_gas=671
    integer,parameter::idx_SiO2_gas=672
    integer,parameter::idx_SiS_gas=673
    integer,parameter::idx_SiSj_gas=674
    integer,parameter::idx_SiN_gas=675
    integer,parameter::idx_SiNj_gas=676
    integer,parameter::idx_SiC_gas=677
    integer,parameter::idx_SiCj_gas=678
    integer,parameter::idx_SiC2_gas=679
    integer,parameter::idx_SiC2j_gas=680
    integer,parameter::idx_SiH_gas=681
    integer,parameter::idx_SiD_gas=682
    integer,parameter::idx_SiH2_gas=683
    integer,parameter::idx_SiH2j_gas=684
    integer,parameter::idx_SiHD_gas=685
    integer,parameter::idx_SiHDj_gas=686
    integer,parameter::idx_SiD2_gas=687
    integer,parameter::idx_SiD2j_gas=688
    integer,parameter::idx_SiH3_gas=689
    integer,parameter::idx_SiH4_gas=690
    integer,parameter::idx_P_gas=691
    integer,parameter::idx_Pj_gas=692
    integer,parameter::idx_PH_gas=693
    integer,parameter::idx_PHj_gas=694
    integer,parameter::idx_PD_gas=695
    integer,parameter::idx_PDj_gas=696
    integer,parameter::idx_PO_gas=697
    integer,parameter::idx_POj_gas=698
    integer,parameter::idx_PH2_gas=699
    integer,parameter::idx_PH2j_gas=700
    integer,parameter::idx_PHD_gas=701
    integer,parameter::idx_PHDj_gas=702
    integer,parameter::idx_PD2_gas=703
    integer,parameter::idx_PD2j_gas=704
    integer,parameter::idx_C3P_gas=705
    integer,parameter::idx_C4P_gas=706
    integer,parameter::idx_C4Pj_gas=707
    integer,parameter::idx_CH2PH_gas=708
    integer,parameter::idx_CH2PD_gas=709
    integer,parameter::idx_CHDPH_gas=710
    integer,parameter::idx_CHDPD_gas=711
    integer,parameter::idx_CD2PH_gas=712
    integer,parameter::idx_CD2PD_gas=713
    integer,parameter::idx_CH2Sij_gas=714
    integer,parameter::idx_CHDSij_gas=715
    integer,parameter::idx_CD2Sij_gas=716
    integer,parameter::idx_CHSij_gas=717
    integer,parameter::idx_CDSij_gas=718
    integer,parameter::idx_CP_gas=719
    integer,parameter::idx_CPj_gas=720
    integer,parameter::idx_D2POj_gas=721
    integer,parameter::idx_D3SiOj_gas=722
    integer,parameter::idx_DCP_gas=723
    integer,parameter::idx_DCPj_gas=724
    integer,parameter::idx_DF_gas=725
    integer,parameter::idx_DFj_gas=726
    integer,parameter::idx_DPO_gas=727
    integer,parameter::idx_DPOj_gas=728
    integer,parameter::idx_DSiNHj_gas=729
    integer,parameter::idx_DSiNDj_gas=730
    integer,parameter::idx_H2CSiCH_gas=731
    integer,parameter::idx_H2CSiCD_gas=732
    integer,parameter::idx_HDCSiCD_gas=733
    integer,parameter::idx_HDCSiCH_gas=734
    integer,parameter::idx_D2CSiCH_gas=735
    integer,parameter::idx_D2CSiCD_gas=736
    integer,parameter::idx_H2DSiOj_gas=737
    integer,parameter::idx_H2Fj_gas=738
    integer,parameter::idx_HDFj_gas=739
    integer,parameter::idx_D2Fj_gas=740
    integer,parameter::idx_H2POj_gas=741
    integer,parameter::idx_H2SiO_gas=742
    integer,parameter::idx_H2SiOj_gas=743
    integer,parameter::idx_HDSiO_gas=744
    integer,parameter::idx_HDSiOj_gas=745
    integer,parameter::idx_D2SiO_gas=746
    integer,parameter::idx_D2SiOj_gas=747
    integer,parameter::idx_H3SiOj_gas=748
    integer,parameter::idx_HCCP_gas=749
    integer,parameter::idx_DCCP_gas=750
    integer,parameter::idx_HCCSi_gas=751
    integer,parameter::idx_DCCSi_gas=752
    integer,parameter::idx_HCP_gas=753
    integer,parameter::idx_HCPj_gas=754
    integer,parameter::idx_HCSi_gas=755
    integer,parameter::idx_DCSi_gas=756
    integer,parameter::idx_HD2SiOj_gas=757
    integer,parameter::idx_HDPOj_gas=758
    integer,parameter::idx_HF_gas=759
    integer,parameter::idx_HFj_gas=760
    integer,parameter::idx_HNSi_gas=761
    integer,parameter::idx_HNSij_gas=762
    integer,parameter::idx_DNSi_gas=763
    integer,parameter::idx_DNSij_gas=764
    integer,parameter::idx_HPNj_gas=765
    integer,parameter::idx_DPNj_gas=766
    integer,parameter::idx_HPO_gas=767
    integer,parameter::idx_HPOj_gas=768
    integer,parameter::idx_HSiNHj_gas=769
    integer,parameter::idx_HSiNDj_gas=770
    integer,parameter::idx_HSiOj_gas=771
    integer,parameter::idx_DSiOj_gas=772
    integer,parameter::idx_HSiO2j_gas=773
    integer,parameter::idx_DSiO2j_gas=774
    integer,parameter::idx_HSiSj_gas=775
    integer,parameter::idx_DSiSj_gas=776
    integer,parameter::idx_HeHj_gas=777
    integer,parameter::idx_HeDj_gas=778
    integer,parameter::idx_PC2Hj_gas=779
    integer,parameter::idx_PC2Dj_gas=780
    integer,parameter::idx_PC2H2j_gas=781
    integer,parameter::idx_PC2HDj_gas=782
    integer,parameter::idx_PC2D2j_gas=783
    integer,parameter::idx_PC2H3j_gas=784
    integer,parameter::idx_PC2H2Dj_gas=785
    integer,parameter::idx_PC2HD2j_gas=786
    integer,parameter::idx_PC2D3j_gas=787
    integer,parameter::idx_PC2H4j_gas=788
    integer,parameter::idx_PC3Hj_gas=789
    integer,parameter::idx_PC3Dj_gas=790
    integer,parameter::idx_PC4Hj_gas=791
    integer,parameter::idx_PC4Dj_gas=792
    integer,parameter::idx_PC4H2j_gas=793
    integer,parameter::idx_PCH2j_gas=794
    integer,parameter::idx_PCHDj_gas=795
    integer,parameter::idx_PCD2j_gas=796
    integer,parameter::idx_PCH3j_gas=797
    integer,parameter::idx_PCH2Dj_gas=798
    integer,parameter::idx_PCHD2j_gas=799
    integer,parameter::idx_PCD3j_gas=800
    integer,parameter::idx_PCH4j_gas=801
    integer,parameter::idx_PCH3Dj_gas=802
    integer,parameter::idx_PCH2D2j_gas=803
    integer,parameter::idx_PCHD3j_gas=804
    integer,parameter::idx_PCD4j_gas=805
    integer,parameter::idx_PH3j_gas=806
    integer,parameter::idx_PH2Dj_gas=807
    integer,parameter::idx_PHD2j_gas=808
    integer,parameter::idx_PD3j_gas=809
    integer,parameter::idx_PN_gas=810
    integer,parameter::idx_PNj_gas=811
    integer,parameter::idx_PNH2j_gas=812
    integer,parameter::idx_PNHDj_gas=813
    integer,parameter::idx_PND2j_gas=814
    integer,parameter::idx_PNH3j_gas=815
    integer,parameter::idx_PNHD2j_gas=816
    integer,parameter::idx_PNH2Dj_gas=817
    integer,parameter::idx_PND3j_gas=818
    integer,parameter::idx_SiC2CH3_gas=819
    integer,parameter::idx_SiC2Hj_gas=820
    integer,parameter::idx_SiC2Dj_gas=821
    integer,parameter::idx_SiC2H2j_gas=822
    integer,parameter::idx_SiC2HDj_gas=823
    integer,parameter::idx_SiC2D2j_gas=824
    integer,parameter::idx_SiC2H3j_gas=825
    integer,parameter::idx_SiC2H2Dj_gas=826
    integer,parameter::idx_SiC2HD2j_gas=827
    integer,parameter::idx_SiC2D3j_gas=828
    integer,parameter::idx_SiC3D_gas=829
    integer,parameter::idx_SiC3Dj_gas=830
    integer,parameter::idx_SiC3H_gas=831
    integer,parameter::idx_SiC3Hj_gas=832
    integer,parameter::idx_SiC3H2j_gas=833
    integer,parameter::idx_SiC3HDj_gas=834
    integer,parameter::idx_SiC3D2j_gas=835
    integer,parameter::idx_SiC3H5_gas=836
    integer,parameter::idx_SiC4_gas=837
    integer,parameter::idx_SiC4j_gas=838
    integer,parameter::idx_SiC4D_gas=839
    integer,parameter::idx_SiC4Dj_gas=840
    integer,parameter::idx_SiC4H_gas=841
    integer,parameter::idx_SiC4Hj_gas=842
    integer,parameter::idx_SiC6H_gas=843
    integer,parameter::idx_SiC8H_gas=844
    integer,parameter::idx_SiCD3_gas=845
    integer,parameter::idx_SiCD3j_gas=846
    integer,parameter::idx_SiCH2_gas=847
    integer,parameter::idx_SiCHD_gas=848
    integer,parameter::idx_SiCD2_gas=849
    integer,parameter::idx_SiCH2D_gas=850
    integer,parameter::idx_SiCH2Dj_gas=851
    integer,parameter::idx_SiCH3_gas=852
    integer,parameter::idx_SiCH3j_gas=853
    integer,parameter::idx_SiCH4j_gas=854
    integer,parameter::idx_SiCH3Dj_gas=855
    integer,parameter::idx_SiCH2D2j_gas=856
    integer,parameter::idx_SiCHD3j_gas=857
    integer,parameter::idx_SiCD4j_gas=858
    integer,parameter::idx_SiCHD2_gas=859
    integer,parameter::idx_SiCHD2j_gas=860
    integer,parameter::idx_SiD3_gas=861
    integer,parameter::idx_SiD3j_gas=862
    integer,parameter::idx_SiFj_gas=863
    integer,parameter::idx_SiHj_gas=864
    integer,parameter::idx_SiDj_gas=865
    integer,parameter::idx_SiH2D_gas=866
    integer,parameter::idx_SiH2Dj_gas=867
    integer,parameter::idx_SiH2D2_gas=868
    integer,parameter::idx_SiH2D2j_gas=869
    integer,parameter::idx_SiH3j_gas=870
    integer,parameter::idx_SiH3D_gas=871
    integer,parameter::idx_SiH3Dj_gas=872
    integer,parameter::idx_SiH4j_gas=873
    integer,parameter::idx_SiD4_gas=874
    integer,parameter::idx_SiD4j_gas=875
    integer,parameter::idx_SiH5j_gas=876
    integer,parameter::idx_SiH4Dj_gas=877
    integer,parameter::idx_SiH3D2j_gas=878
    integer,parameter::idx_SiH2D3j_gas=879
    integer,parameter::idx_SiHD4j_gas=880
    integer,parameter::idx_SiD5j_gas=881
    integer,parameter::idx_SiHD2_gas=882
    integer,parameter::idx_SiHD2j_gas=883
    integer,parameter::idx_SiHD3_gas=884
    integer,parameter::idx_SiHD3j_gas=885
    integer,parameter::idx_SiNC_gas=886
    integer,parameter::idx_SiNCj_gas=887
    integer,parameter::idx_SiNCHj_gas=888
    integer,parameter::idx_SiNCDj_gas=889
    integer,parameter::idx_c_DCCHSi_gas=890
    integer,parameter::idx_c_DCCDSi_gas=891
    integer,parameter::idx_c_HCCHSi_gas=892
    integer,parameter::idx_c_HCCDSi_gas=893
    integer,parameter::idx_c_SiC2_gas=894
    integer,parameter::idx_l_SiC3_gas=895
    integer,parameter::idx_l_SiC3j_gas=896
    integer,parameter::idx_l_C3H_gas=897
    integer,parameter::idx_l_C3D_gas=898
    integer,parameter::idx_c_C3H_gas=899
    integer,parameter::idx_c_C3D_gas=900
    integer,parameter::idx_l_C3H2_gas=901
    integer,parameter::idx_l_C3HD_gas=902
    integer,parameter::idx_l_C3D2_gas=903
    integer,parameter::idx_c_C3H2_gas=904
    integer,parameter::idx_c_C3HD_gas=905
    integer,parameter::idx_c_C3D2_gas=906
    integer,parameter::idx_l_C3H2j_gas=907
    integer,parameter::idx_l_C3HDj_gas=908
    integer,parameter::idx_l_C3D2j_gas=909
    integer,parameter::idx_c_C3H2j_gas=910
    integer,parameter::idx_c_C3HDj_gas=911
    integer,parameter::idx_c_C3D2j_gas=912
    integer,parameter::idx_l_C3H3j_gas=913
    integer,parameter::idx_l_C3H2Dj_gas=914
    integer,parameter::idx_l_C3HD2j_gas=915
    integer,parameter::idx_l_C3D3j_gas=916
    integer,parameter::idx_c_C3H3j_gas=917
    integer,parameter::idx_c_C3H2Dj_gas=918
    integer,parameter::idx_c_C3HD2j_gas=919
    integer,parameter::idx_c_C3D3j_gas=920
    integer,parameter::idx_Mg_gas=921
    integer,parameter::idx_Mgj_gas=922
    integer,parameter::idx_MgH_gas=923
    integer,parameter::idx_MgD_gas=924
    integer,parameter::idx_MgH2_gas=925
    integer,parameter::idx_MgHD_gas=926
    integer,parameter::idx_MgD2_gas=927
    integer,parameter::idx_HCl_gas=928
    integer,parameter::idx_HClj_gas=929
    integer,parameter::idx_DCl_gas=930
    integer,parameter::idx_DClj_gas=931
    integer,parameter::idx_H2Cl_gas=932
    integer,parameter::idx_H2Clj_gas=933
    integer,parameter::idx_HDCl_gas=934
    integer,parameter::idx_HDClj_gas=935
    integer,parameter::idx_D2Cl_gas=936
    integer,parameter::idx_D2Clj_gas=937
    integer,parameter::idx_CCl_gas=938
    integer,parameter::idx_CClj_gas=939
    integer,parameter::idx_ClO_gas=940
    integer,parameter::idx_ClOj_gas=941
    integer,parameter::idx_H2CClj_gas=942
    integer,parameter::idx_HDCClj_gas=943
    integer,parameter::idx_D2CClj_gas=944
    integer,parameter::idx_Na_gas=945
    integer,parameter::idx_Naj_gas=946
    integer,parameter::idx_NaH_gas=947
    integer,parameter::idx_NaD_gas=948
    integer,parameter::idx_NaOH_gas=949
    integer,parameter::idx_NaOD_gas=950
    integer,parameter::idx_NaH2j_gas=951
    integer,parameter::idx_NaHDj_gas=952
    integer,parameter::idx_NaD2j_gas=953
    integer,parameter::idx_NaH2Oj_gas=954
    integer,parameter::idx_NaHDOj_gas=955
    integer,parameter::idx_NaD2Oj_gas=956
    integer,parameter::idx_H_0001=957
    integer,parameter::idx_H_0002=1413
    integer,parameter::idx_D_0001=958
    integer,parameter::idx_D_0002=1414
    integer,parameter::idx_p_H2_0001=959
    integer,parameter::idx_p_H2_0002=1415
    integer,parameter::idx_p_D2_0001=960
    integer,parameter::idx_p_D2_0002=1416
    integer,parameter::idx_o_H2_0001=961
    integer,parameter::idx_o_H2_0002=1417
    integer,parameter::idx_o_D2_0001=962
    integer,parameter::idx_o_D2_0002=1418
    integer,parameter::idx_HD_0001=963
    integer,parameter::idx_HD_0002=1419
    integer,parameter::idx_He_0001=964
    integer,parameter::idx_He_0002=1420
    integer,parameter::idx_O_0001=965
    integer,parameter::idx_O_0002=1421
    integer,parameter::idx_O2_0001=966
    integer,parameter::idx_O2_0002=1422
    integer,parameter::idx_O3_0001=967
    integer,parameter::idx_O3_0002=1423
    integer,parameter::idx_OH_0001=968
    integer,parameter::idx_OH_0002=1424
    integer,parameter::idx_OD_0001=969
    integer,parameter::idx_OD_0002=1425
    integer,parameter::idx_H2O_0001=970
    integer,parameter::idx_H2O_0002=1426
    integer,parameter::idx_HDO_0001=971
    integer,parameter::idx_HDO_0002=1427
    integer,parameter::idx_D2O_0001=972
    integer,parameter::idx_D2O_0002=1428
    integer,parameter::idx_O2H_0001=973
    integer,parameter::idx_O2H_0002=1429
    integer,parameter::idx_O2D_0001=974
    integer,parameter::idx_O2D_0002=1430
    integer,parameter::idx_HOOH_0001=975
    integer,parameter::idx_HOOH_0002=1431
    integer,parameter::idx_HOOD_0001=976
    integer,parameter::idx_HOOD_0002=1432
    integer,parameter::idx_DOOD_0001=977
    integer,parameter::idx_DOOD_0002=1433
    integer,parameter::idx_Fe_0001=978
    integer,parameter::idx_Fe_0002=1434
    integer,parameter::idx_FeH_0001=979
    integer,parameter::idx_FeH_0002=1435
    integer,parameter::idx_N_0001=980
    integer,parameter::idx_N_0002=1436
    integer,parameter::idx_S_0001=981
    integer,parameter::idx_S_0002=1437
    integer,parameter::idx_C_0001=982
    integer,parameter::idx_C_0002=1438
    integer,parameter::idx_C2_0001=983
    integer,parameter::idx_C2_0002=1439
    integer,parameter::idx_CO_0001=984
    integer,parameter::idx_CO_0002=1440
    integer,parameter::idx_HCO_0001=985
    integer,parameter::idx_HCO_0002=1441
    integer,parameter::idx_DCO_0001=986
    integer,parameter::idx_DCO_0002=1442
    integer,parameter::idx_H2CO_0001=987
    integer,parameter::idx_H2CO_0002=1443
    integer,parameter::idx_HDCO_0001=988
    integer,parameter::idx_HDCO_0002=1444
    integer,parameter::idx_D2CO_0001=989
    integer,parameter::idx_D2CO_0002=1445
    integer,parameter::idx_CH2OH_0001=990
    integer,parameter::idx_CH2OH_0002=1446
    integer,parameter::idx_CD2OD_0001=991
    integer,parameter::idx_CD2OD_0002=1447
    integer,parameter::idx_CH2OD_0001=992
    integer,parameter::idx_CH2OD_0002=1448
    integer,parameter::idx_CHDOH_0001=993
    integer,parameter::idx_CHDOH_0002=1449
    integer,parameter::idx_CHDOD_0001=994
    integer,parameter::idx_CHDOD_0002=1450
    integer,parameter::idx_CD2OH_0001=995
    integer,parameter::idx_CD2OH_0002=1451
    integer,parameter::idx_CH3O_0001=996
    integer,parameter::idx_CH3O_0002=1452
    integer,parameter::idx_CHD2O_0001=997
    integer,parameter::idx_CHD2O_0002=1453
    integer,parameter::idx_CH2DO_0001=998
    integer,parameter::idx_CH2DO_0002=1454
    integer,parameter::idx_CD3O_0001=999
    integer,parameter::idx_CD3O_0002=1455
    integer,parameter::idx_CH3OH_0001=1000
    integer,parameter::idx_CH3OH_0002=1456
    integer,parameter::idx_CH3OD_0001=1001
    integer,parameter::idx_CH3OD_0002=1457
    integer,parameter::idx_CHD2OH_0001=1002
    integer,parameter::idx_CHD2OH_0002=1458
    integer,parameter::idx_CHD2OD_0001=1003
    integer,parameter::idx_CHD2OD_0002=1459
    integer,parameter::idx_CH2DOH_0001=1004
    integer,parameter::idx_CH2DOH_0002=1460
    integer,parameter::idx_CH2DOD_0001=1005
    integer,parameter::idx_CH2DOD_0002=1461
    integer,parameter::idx_CD3OD_0001=1006
    integer,parameter::idx_CD3OD_0002=1462
    integer,parameter::idx_CD3OH_0001=1007
    integer,parameter::idx_CD3OH_0002=1463
    integer,parameter::idx_CH_0001=1008
    integer,parameter::idx_CH_0002=1464
    integer,parameter::idx_CD_0001=1009
    integer,parameter::idx_CD_0002=1465
    integer,parameter::idx_CH2_0001=1010
    integer,parameter::idx_CH2_0002=1466
    integer,parameter::idx_CHD_0001=1011
    integer,parameter::idx_CHD_0002=1467
    integer,parameter::idx_CD2_0001=1012
    integer,parameter::idx_CD2_0002=1468
    integer,parameter::idx_CH3_0001=1013
    integer,parameter::idx_CH3_0002=1469
    integer,parameter::idx_CH2D_0001=1014
    integer,parameter::idx_CH2D_0002=1470
    integer,parameter::idx_CHD2_0001=1015
    integer,parameter::idx_CHD2_0002=1471
    integer,parameter::idx_CD3_0001=1016
    integer,parameter::idx_CD3_0002=1472
    integer,parameter::idx_CH4_0001=1017
    integer,parameter::idx_CH4_0002=1473
    integer,parameter::idx_CH3D_0001=1018
    integer,parameter::idx_CH3D_0002=1474
    integer,parameter::idx_CH2D2_0001=1019
    integer,parameter::idx_CH2D2_0002=1475
    integer,parameter::idx_CHD3_0001=1020
    integer,parameter::idx_CHD3_0002=1476
    integer,parameter::idx_CD4_0001=1021
    integer,parameter::idx_CD4_0002=1477
    integer,parameter::idx_CO2_0001=1022
    integer,parameter::idx_CO2_0002=1478
    integer,parameter::idx_HCOOH_0001=1023
    integer,parameter::idx_HCOOH_0002=1479
    integer,parameter::idx_HCOOD_0001=1024
    integer,parameter::idx_HCOOD_0002=1480
    integer,parameter::idx_DCOOH_0001=1025
    integer,parameter::idx_DCOOH_0002=1481
    integer,parameter::idx_DCOOD_0001=1026
    integer,parameter::idx_DCOOD_0002=1482
    integer,parameter::idx_HOCO_0001=1027
    integer,parameter::idx_HOCO_0002=1483
    integer,parameter::idx_DOCO_0001=1028
    integer,parameter::idx_DOCO_0002=1484
    integer,parameter::idx_NH_0001=1029
    integer,parameter::idx_NH_0002=1485
    integer,parameter::idx_ND_0001=1030
    integer,parameter::idx_ND_0002=1486
    integer,parameter::idx_NH2_0001=1031
    integer,parameter::idx_NH2_0002=1487
    integer,parameter::idx_NHD_0001=1032
    integer,parameter::idx_NHD_0002=1488
    integer,parameter::idx_ND2_0001=1033
    integer,parameter::idx_ND2_0002=1489
    integer,parameter::idx_NH3_0001=1034
    integer,parameter::idx_NH3_0002=1490
    integer,parameter::idx_NH2D_0001=1035
    integer,parameter::idx_NH2D_0002=1491
    integer,parameter::idx_NHD2_0001=1036
    integer,parameter::idx_NHD2_0002=1492
    integer,parameter::idx_ND3_0001=1037
    integer,parameter::idx_ND3_0002=1493
    integer,parameter::idx_C10_0001=1038
    integer,parameter::idx_C10_0002=1494
    integer,parameter::idx_C10H_0001=1039
    integer,parameter::idx_C10H_0002=1495
    integer,parameter::idx_C10H2_0001=1040
    integer,parameter::idx_C10H2_0002=1496
    integer,parameter::idx_C10N_0001=1041
    integer,parameter::idx_C10N_0002=1497
    integer,parameter::idx_C11_0001=1042
    integer,parameter::idx_C11_0002=1498
    integer,parameter::idx_C2H2_0001=1043
    integer,parameter::idx_C2H2_0002=1499
    integer,parameter::idx_C2HD_0001=1044
    integer,parameter::idx_C2HD_0002=1500
    integer,parameter::idx_C2D2_0001=1045
    integer,parameter::idx_C2D2_0002=1501
    integer,parameter::idx_C2H3_0001=1046
    integer,parameter::idx_C2H3_0002=1502
    integer,parameter::idx_C2H2D_0001=1047
    integer,parameter::idx_C2H2D_0002=1503
    integer,parameter::idx_C2HD2_0001=1048
    integer,parameter::idx_C2HD2_0002=1504
    integer,parameter::idx_C2D3_0001=1049
    integer,parameter::idx_C2D3_0002=1505
    integer,parameter::idx_C2H4_0001=1050
    integer,parameter::idx_C2H4_0002=1506
    integer,parameter::idx_C2H3D_0001=1051
    integer,parameter::idx_C2H3D_0002=1507
    integer,parameter::idx_C2H2D2_0001=1052
    integer,parameter::idx_C2H2D2_0002=1508
    integer,parameter::idx_C2HD3_0001=1053
    integer,parameter::idx_C2HD3_0002=1509
    integer,parameter::idx_C2D4_0001=1054
    integer,parameter::idx_C2D4_0002=1510
    integer,parameter::idx_C2H5_0001=1055
    integer,parameter::idx_C2H5_0002=1511
    integer,parameter::idx_C2H6_0001=1056
    integer,parameter::idx_C2H6_0002=1512
    integer,parameter::idx_C3_0001=1057
    integer,parameter::idx_C3_0002=1513
    integer,parameter::idx_C3N_0001=1058
    integer,parameter::idx_C3N_0002=1514
    integer,parameter::idx_C3O_0001=1059
    integer,parameter::idx_C3O_0002=1515
    integer,parameter::idx_C3S_0001=1060
    integer,parameter::idx_C3S_0002=1516
    integer,parameter::idx_C4_0001=1061
    integer,parameter::idx_C4_0002=1517
    integer,parameter::idx_C4H_0001=1062
    integer,parameter::idx_C4H_0002=1518
    integer,parameter::idx_C4D_0001=1063
    integer,parameter::idx_C4D_0002=1519
    integer,parameter::idx_C4H2_0001=1064
    integer,parameter::idx_C4H2_0002=1520
    integer,parameter::idx_C4HD_0001=1065
    integer,parameter::idx_C4HD_0002=1521
    integer,parameter::idx_C4D2_0001=1066
    integer,parameter::idx_C4D2_0002=1522
    integer,parameter::idx_C4H3_0001=1067
    integer,parameter::idx_C4H3_0002=1523
    integer,parameter::idx_C4N_0001=1068
    integer,parameter::idx_C4N_0002=1524
    integer,parameter::idx_C4S_0001=1069
    integer,parameter::idx_C4S_0002=1525
    integer,parameter::idx_C5_0001=1070
    integer,parameter::idx_C5_0002=1526
    integer,parameter::idx_C5H_0001=1071
    integer,parameter::idx_C5H_0002=1527
    integer,parameter::idx_C5D_0001=1072
    integer,parameter::idx_C5D_0002=1528
    integer,parameter::idx_C5H2_0001=1073
    integer,parameter::idx_C5H2_0002=1529
    integer,parameter::idx_C5H3_0001=1074
    integer,parameter::idx_C5H3_0002=1530
    integer,parameter::idx_C5H4_0001=1075
    integer,parameter::idx_C5H4_0002=1531
    integer,parameter::idx_C5N_0001=1076
    integer,parameter::idx_C5N_0002=1532
    integer,parameter::idx_C5O_0001=1077
    integer,parameter::idx_C5O_0002=1533
    integer,parameter::idx_C6_0001=1078
    integer,parameter::idx_C6_0002=1534
    integer,parameter::idx_C6H_0001=1079
    integer,parameter::idx_C6H_0002=1535
    integer,parameter::idx_C6H2_0001=1080
    integer,parameter::idx_C6H2_0002=1536
    integer,parameter::idx_C6H3_0001=1081
    integer,parameter::idx_C6H3_0002=1537
    integer,parameter::idx_C6H4_0001=1082
    integer,parameter::idx_C6H4_0002=1538
    integer,parameter::idx_C6H6_0001=1083
    integer,parameter::idx_C6H6_0002=1539
    integer,parameter::idx_C6N_0001=1084
    integer,parameter::idx_C6N_0002=1540
    integer,parameter::idx_C7_0001=1085
    integer,parameter::idx_C7_0002=1541
    integer,parameter::idx_C7H_0001=1086
    integer,parameter::idx_C7H_0002=1542
    integer,parameter::idx_C7H2_0001=1087
    integer,parameter::idx_C7H2_0002=1543
    integer,parameter::idx_C7H3_0001=1088
    integer,parameter::idx_C7H3_0002=1544
    integer,parameter::idx_C7H4_0001=1089
    integer,parameter::idx_C7H4_0002=1545
    integer,parameter::idx_C7N_0001=1090
    integer,parameter::idx_C7N_0002=1546
    integer,parameter::idx_C7O_0001=1091
    integer,parameter::idx_C7O_0002=1547
    integer,parameter::idx_C8_0001=1092
    integer,parameter::idx_C8_0002=1548
    integer,parameter::idx_C8H_0001=1093
    integer,parameter::idx_C8H_0002=1549
    integer,parameter::idx_C8H2_0001=1094
    integer,parameter::idx_C8H2_0002=1550
    integer,parameter::idx_C8H3_0001=1095
    integer,parameter::idx_C8H3_0002=1551
    integer,parameter::idx_C8H4_0001=1096
    integer,parameter::idx_C8H4_0002=1552
    integer,parameter::idx_C8N_0001=1097
    integer,parameter::idx_C8N_0002=1553
    integer,parameter::idx_C9_0001=1098
    integer,parameter::idx_C9_0002=1554
    integer,parameter::idx_C9H_0001=1099
    integer,parameter::idx_C9H_0002=1555
    integer,parameter::idx_C9H2_0001=1100
    integer,parameter::idx_C9H2_0002=1556
    integer,parameter::idx_C9H3_0001=1101
    integer,parameter::idx_C9H3_0002=1557
    integer,parameter::idx_C9H4_0001=1102
    integer,parameter::idx_C9H4_0002=1558
    integer,parameter::idx_C9N_0001=1103
    integer,parameter::idx_C9N_0002=1559
    integer,parameter::idx_C9O_0001=1104
    integer,parameter::idx_C9O_0002=1560
    integer,parameter::idx_CCH_0001=1105
    integer,parameter::idx_CCH_0002=1561
    integer,parameter::idx_CCD_0001=1106
    integer,parameter::idx_CCD_0002=1562
    integer,parameter::idx_CCN_0001=1107
    integer,parameter::idx_CCN_0002=1563
    integer,parameter::idx_CCO_0001=1108
    integer,parameter::idx_CCO_0002=1564
    integer,parameter::idx_CCS_0001=1109
    integer,parameter::idx_CCS_0002=1565
    integer,parameter::idx_CD3CN_0001=1110
    integer,parameter::idx_CD3CN_0002=1566
    integer,parameter::idx_CH2CCH_0001=1111
    integer,parameter::idx_CH2CCH_0002=1567
    integer,parameter::idx_CH2CCD_0001=1112
    integer,parameter::idx_CH2CCD_0002=1568
    integer,parameter::idx_CHDCCH_0001=1113
    integer,parameter::idx_CHDCCH_0002=1569
    integer,parameter::idx_CHDCCD_0001=1114
    integer,parameter::idx_CHDCCD_0002=1570
    integer,parameter::idx_CD2CCH_0001=1115
    integer,parameter::idx_CD2CCH_0002=1571
    integer,parameter::idx_CD2CCD_0001=1116
    integer,parameter::idx_CD2CCD_0002=1572
    integer,parameter::idx_CH2CHC2H_0001=1117
    integer,parameter::idx_CH2CHC2H_0002=1573
    integer,parameter::idx_CH2CHCHCH2_0001=1118
    integer,parameter::idx_CH2CHCHCH2_0002=1574
    integer,parameter::idx_CH2CHCN_0001=1119
    integer,parameter::idx_CH2CHCN_0002=1575
    integer,parameter::idx_CH2NH_0001=1120
    integer,parameter::idx_CH2NH_0002=1576
    integer,parameter::idx_CHDNH_0001=1121
    integer,parameter::idx_CHDNH_0002=1577
    integer,parameter::idx_CHDND_0001=1122
    integer,parameter::idx_CHDND_0002=1578
    integer,parameter::idx_CH2ND_0001=1123
    integer,parameter::idx_CH2ND_0002=1579
    integer,parameter::idx_CD2NH_0001=1124
    integer,parameter::idx_CD2NH_0002=1580
    integer,parameter::idx_CD2ND_0001=1125
    integer,parameter::idx_CD2ND_0002=1581
    integer,parameter::idx_CH2NH2_0001=1126
    integer,parameter::idx_CH2NH2_0002=1582
    integer,parameter::idx_CH2NHD_0001=1127
    integer,parameter::idx_CH2NHD_0002=1583
    integer,parameter::idx_CHDNH2_0001=1128
    integer,parameter::idx_CHDNH2_0002=1584
    integer,parameter::idx_CHDNHD_0001=1129
    integer,parameter::idx_CHDNHD_0002=1585
    integer,parameter::idx_CHDND2_0001=1130
    integer,parameter::idx_CHDND2_0002=1586
    integer,parameter::idx_CH2ND2_0001=1131
    integer,parameter::idx_CH2ND2_0002=1587
    integer,parameter::idx_CD2NH2_0001=1132
    integer,parameter::idx_CD2NH2_0002=1588
    integer,parameter::idx_CD2NHD_0001=1133
    integer,parameter::idx_CD2NHD_0002=1589
    integer,parameter::idx_CD2ND2_0001=1134
    integer,parameter::idx_CD2ND2_0002=1590
    integer,parameter::idx_CH3C3N_0001=1135
    integer,parameter::idx_CH3C3N_0002=1591
    integer,parameter::idx_CH3C4H_0001=1136
    integer,parameter::idx_CH3C4H_0002=1592
    integer,parameter::idx_CH3C5N_0001=1137
    integer,parameter::idx_CH3C5N_0002=1593
    integer,parameter::idx_CH3C6H_0001=1138
    integer,parameter::idx_CH3C6H_0002=1594
    integer,parameter::idx_CH3C7N_0001=1139
    integer,parameter::idx_CH3C7N_0002=1595
    integer,parameter::idx_CH3CCH_0001=1140
    integer,parameter::idx_CH3CCH_0002=1596
    integer,parameter::idx_CH3CH2OH_0001=1141
    integer,parameter::idx_CH3CH2OH_0002=1597
    integer,parameter::idx_CH3CHCH2_0001=1142
    integer,parameter::idx_CH3CHCH2_0002=1598
    integer,parameter::idx_CH3CHO_0001=1143
    integer,parameter::idx_CH3CHO_0002=1599
    integer,parameter::idx_CH3CN_0001=1144
    integer,parameter::idx_CH3CN_0002=1600
    integer,parameter::idx_CH2DCN_0001=1145
    integer,parameter::idx_CH2DCN_0002=1601
    integer,parameter::idx_CHD2CN_0001=1146
    integer,parameter::idx_CHD2CN_0002=1602
    integer,parameter::idx_CH3COCH3_0001=1147
    integer,parameter::idx_CH3COCH3_0002=1603
    integer,parameter::idx_CH3NH2_0001=1148
    integer,parameter::idx_CH3NH2_0002=1604
    integer,parameter::idx_CH3OCH2_0001=1149
    integer,parameter::idx_CH3OCH2_0002=1605
    integer,parameter::idx_CH3OCH3_0001=1150
    integer,parameter::idx_CH3OCH3_0002=1606
    integer,parameter::idx_CN_0001=1151
    integer,parameter::idx_CN_0002=1607
    integer,parameter::idx_CS_0001=1152
    integer,parameter::idx_CS_0002=1608
    integer,parameter::idx_H2CCN_0001=1153
    integer,parameter::idx_H2CCN_0002=1609
    integer,parameter::idx_HDCCN_0001=1154
    integer,parameter::idx_HDCCN_0002=1610
    integer,parameter::idx_D2CCN_0001=1155
    integer,parameter::idx_D2CCN_0002=1611
    integer,parameter::idx_H2CCO_0001=1156
    integer,parameter::idx_H2CCO_0002=1612
    integer,parameter::idx_HDCCO_0001=1157
    integer,parameter::idx_HDCCO_0002=1613
    integer,parameter::idx_D2CCO_0001=1158
    integer,parameter::idx_D2CCO_0002=1614
    integer,parameter::idx_H2CN_0001=1159
    integer,parameter::idx_H2CN_0002=1615
    integer,parameter::idx_HDCN_0001=1160
    integer,parameter::idx_HDCN_0002=1616
    integer,parameter::idx_D2CN_0001=1161
    integer,parameter::idx_D2CN_0002=1617
    integer,parameter::idx_H2CS_0001=1162
    integer,parameter::idx_H2CS_0002=1618
    integer,parameter::idx_HDCS_0001=1163
    integer,parameter::idx_HDCS_0002=1619
    integer,parameter::idx_D2CS_0001=1164
    integer,parameter::idx_D2CS_0002=1620
    integer,parameter::idx_H2S_0001=1165
    integer,parameter::idx_H2S_0002=1621
    integer,parameter::idx_HDS_0001=1166
    integer,parameter::idx_HDS_0002=1622
    integer,parameter::idx_D2S_0001=1167
    integer,parameter::idx_D2S_0002=1623
    integer,parameter::idx_HC2O_0001=1168
    integer,parameter::idx_HC2O_0002=1624
    integer,parameter::idx_DC2O_0001=1169
    integer,parameter::idx_DC2O_0002=1625
    integer,parameter::idx_HC3N_0001=1170
    integer,parameter::idx_HC3N_0002=1626
    integer,parameter::idx_DC3N_0001=1171
    integer,parameter::idx_DC3N_0002=1627
    integer,parameter::idx_HC4N_0001=1172
    integer,parameter::idx_HC4N_0002=1628
    integer,parameter::idx_DC4N_0001=1173
    integer,parameter::idx_DC4N_0002=1629
    integer,parameter::idx_HC5N_0001=1174
    integer,parameter::idx_HC5N_0002=1630
    integer,parameter::idx_HC6N_0001=1175
    integer,parameter::idx_HC6N_0002=1631
    integer,parameter::idx_HC7N_0001=1176
    integer,parameter::idx_HC7N_0002=1632
    integer,parameter::idx_HC8N_0001=1177
    integer,parameter::idx_HC8N_0002=1633
    integer,parameter::idx_HC9N_0001=1178
    integer,parameter::idx_HC9N_0002=1634
    integer,parameter::idx_HCCNC_0001=1179
    integer,parameter::idx_HCCNC_0002=1635
    integer,parameter::idx_DCCNC_0001=1180
    integer,parameter::idx_DCCNC_0002=1636
    integer,parameter::idx_HCN_0001=1181
    integer,parameter::idx_HCN_0002=1637
    integer,parameter::idx_DCN_0001=1182
    integer,parameter::idx_DCN_0002=1638
    integer,parameter::idx_HCNCC_0001=1183
    integer,parameter::idx_HCNCC_0002=1639
    integer,parameter::idx_DCNCC_0001=1184
    integer,parameter::idx_DCNCC_0002=1640
    integer,parameter::idx_HCOOCH3_0001=1185
    integer,parameter::idx_HCOOCH3_0002=1641
    integer,parameter::idx_HCS_0001=1186
    integer,parameter::idx_HCS_0002=1642
    integer,parameter::idx_DCS_0001=1187
    integer,parameter::idx_DCS_0002=1643
    integer,parameter::idx_HNC_0001=1188
    integer,parameter::idx_HNC_0002=1644
    integer,parameter::idx_DNC_0001=1189
    integer,parameter::idx_DNC_0002=1645
    integer,parameter::idx_HNCCC_0001=1190
    integer,parameter::idx_HNCCC_0002=1646
    integer,parameter::idx_DNCCC_0001=1191
    integer,parameter::idx_DNCCC_0002=1647
    integer,parameter::idx_HNCO_0001=1192
    integer,parameter::idx_HNCO_0002=1648
    integer,parameter::idx_DNCO_0001=1193
    integer,parameter::idx_DNCO_0002=1649
    integer,parameter::idx_HNO_0001=1194
    integer,parameter::idx_HNO_0002=1650
    integer,parameter::idx_DNO_0001=1195
    integer,parameter::idx_DNO_0002=1651
    integer,parameter::idx_HS_0001=1196
    integer,parameter::idx_HS_0002=1652
    integer,parameter::idx_DS_0001=1197
    integer,parameter::idx_DS_0002=1653
    integer,parameter::idx_N2_0001=1198
    integer,parameter::idx_N2_0002=1654
    integer,parameter::idx_N2O_0001=1199
    integer,parameter::idx_N2O_0002=1655
    integer,parameter::idx_NC4N_0001=1200
    integer,parameter::idx_NC4N_0002=1656
    integer,parameter::idx_NC6N_0001=1201
    integer,parameter::idx_NC6N_0002=1657
    integer,parameter::idx_NC8N_0001=1202
    integer,parameter::idx_NC8N_0002=1658
    integer,parameter::idx_NH2CHO_0001=1203
    integer,parameter::idx_NH2CHO_0002=1659
    integer,parameter::idx_NH2CDO_0001=1204
    integer,parameter::idx_NH2CDO_0002=1660
    integer,parameter::idx_NHDCHO_0001=1205
    integer,parameter::idx_NHDCHO_0002=1661
    integer,parameter::idx_NHDCDO_0001=1206
    integer,parameter::idx_NHDCDO_0002=1662
    integer,parameter::idx_ND2CHO_0001=1207
    integer,parameter::idx_ND2CHO_0002=1663
    integer,parameter::idx_ND2CDO_0001=1208
    integer,parameter::idx_ND2CDO_0002=1664
    integer,parameter::idx_NH2CN_0001=1209
    integer,parameter::idx_NH2CN_0002=1665
    integer,parameter::idx_NHDCN_0001=1210
    integer,parameter::idx_NHDCN_0002=1666
    integer,parameter::idx_ND2CN_0001=1211
    integer,parameter::idx_ND2CN_0002=1667
    integer,parameter::idx_HSO_0001=1212
    integer,parameter::idx_HSO_0002=1668
    integer,parameter::idx_DSO_0001=1213
    integer,parameter::idx_DSO_0002=1669
    integer,parameter::idx_HSS_0001=1214
    integer,parameter::idx_HSS_0002=1670
    integer,parameter::idx_DSS_0001=1215
    integer,parameter::idx_DSS_0002=1671
    integer,parameter::idx_HSSH_0001=1216
    integer,parameter::idx_HSSH_0002=1672
    integer,parameter::idx_HSSD_0001=1217
    integer,parameter::idx_HSSD_0002=1673
    integer,parameter::idx_DSSH_0001=1218
    integer,parameter::idx_DSSH_0002=1674
    integer,parameter::idx_DSSD_0001=1219
    integer,parameter::idx_DSSD_0002=1675
    integer,parameter::idx_NO_0001=1220
    integer,parameter::idx_NO_0002=1676
    integer,parameter::idx_NO2_0001=1221
    integer,parameter::idx_NO2_0002=1677
    integer,parameter::idx_NS_0001=1222
    integer,parameter::idx_NS_0002=1678
    integer,parameter::idx_OCN_0001=1223
    integer,parameter::idx_OCN_0002=1679
    integer,parameter::idx_OCS_0001=1224
    integer,parameter::idx_OCS_0002=1680
    integer,parameter::idx_S2_0001=1225
    integer,parameter::idx_S2_0002=1681
    integer,parameter::idx_SO_0001=1226
    integer,parameter::idx_SO_0002=1682
    integer,parameter::idx_SO2_0001=1227
    integer,parameter::idx_SO2_0002=1683
    integer,parameter::idx_CH3OCHO_0001=1228
    integer,parameter::idx_CH3OCHO_0002=1684
    integer,parameter::idx_Si_0001=1229
    integer,parameter::idx_Si_0002=1685
    integer,parameter::idx_SiS_0001=1230
    integer,parameter::idx_SiS_0002=1686
    integer,parameter::idx_SiN_0001=1231
    integer,parameter::idx_SiN_0002=1687
    integer,parameter::idx_SiC_0001=1232
    integer,parameter::idx_SiC_0002=1688
    integer,parameter::idx_SiH_0001=1233
    integer,parameter::idx_SiH_0002=1689
    integer,parameter::idx_SiH2_0001=1234
    integer,parameter::idx_SiH2_0002=1690
    integer,parameter::idx_SiH3_0001=1235
    integer,parameter::idx_SiH3_0002=1691
    integer,parameter::idx_SiH4_0001=1236
    integer,parameter::idx_SiH4_0002=1692
    integer,parameter::idx_SiC2CH3_0001=1237
    integer,parameter::idx_SiC2CH3_0002=1693
    integer,parameter::idx_SiC3H_0001=1238
    integer,parameter::idx_SiC3H_0002=1694
    integer,parameter::idx_SiC3H5_0001=1239
    integer,parameter::idx_SiC3H5_0002=1695
    integer,parameter::idx_SiC4_0001=1240
    integer,parameter::idx_SiC4_0002=1696
    integer,parameter::idx_SiC4H_0001=1241
    integer,parameter::idx_SiC4H_0002=1697
    integer,parameter::idx_SiC6H_0001=1242
    integer,parameter::idx_SiC6H_0002=1698
    integer,parameter::idx_SiC8H_0001=1243
    integer,parameter::idx_SiC8H_0002=1699
    integer,parameter::idx_c_HCCHSi_0001=1244
    integer,parameter::idx_c_HCCHSi_0002=1700
    integer,parameter::idx_c_SiC2_0001=1245
    integer,parameter::idx_c_SiC2_0002=1701
    integer,parameter::idx_l_C3H_0001=1246
    integer,parameter::idx_l_C3H_0002=1702
    integer,parameter::idx_l_C3D_0001=1247
    integer,parameter::idx_l_C3D_0002=1703
    integer,parameter::idx_c_C3H_0001=1248
    integer,parameter::idx_c_C3H_0002=1704
    integer,parameter::idx_c_C3D_0001=1249
    integer,parameter::idx_c_C3D_0002=1705
    integer,parameter::idx_l_SiC3_0001=1250
    integer,parameter::idx_l_SiC3_0002=1706
    integer,parameter::idx_l_C3H2_0001=1251
    integer,parameter::idx_l_C3H2_0002=1707
    integer,parameter::idx_l_C3HD_0001=1252
    integer,parameter::idx_l_C3HD_0002=1708
    integer,parameter::idx_l_C3D2_0001=1253
    integer,parameter::idx_l_C3D2_0002=1709
    integer,parameter::idx_c_C3H2_0001=1254
    integer,parameter::idx_c_C3H2_0002=1710
    integer,parameter::idx_c_C3HD_0001=1255
    integer,parameter::idx_c_C3HD_0002=1711
    integer,parameter::idx_c_C3D2_0001=1256
    integer,parameter::idx_c_C3D2_0002=1712
    integer,parameter::idx_Mg_0001=1257
    integer,parameter::idx_Mg_0002=1713
    integer,parameter::idx_MgH_0001=1258
    integer,parameter::idx_MgH_0002=1714
    integer,parameter::idx_MgH2_0001=1259
    integer,parameter::idx_MgH2_0002=1715
    integer,parameter::idx_Na_0001=1260
    integer,parameter::idx_Na_0002=1716
    integer,parameter::idx_NaH_0001=1261
    integer,parameter::idx_NaH_0002=1717
    integer,parameter::idx_F_0001=1262
    integer,parameter::idx_F_0002=1718
    integer,parameter::idx_HF_0001=1263
    integer,parameter::idx_HF_0002=1719
    integer,parameter::idx_DF_0001=1264
    integer,parameter::idx_DF_0002=1720
    integer,parameter::idx_MgD_0001=1265
    integer,parameter::idx_MgD_0002=1721
    integer,parameter::idx_MgHD_0001=1266
    integer,parameter::idx_MgHD_0002=1722
    integer,parameter::idx_MgD2_0001=1267
    integer,parameter::idx_MgD2_0002=1723
    integer,parameter::idx_NaD_0001=1268
    integer,parameter::idx_NaD_0002=1724
    integer,parameter::idx_SiD_0001=1269
    integer,parameter::idx_SiD_0002=1725
    integer,parameter::idx_Cl_0001=1270
    integer,parameter::idx_Cl_0002=1726
    integer,parameter::idx_HCl_0001=1271
    integer,parameter::idx_HCl_0002=1727
    integer,parameter::idx_DCl_0001=1272
    integer,parameter::idx_DCl_0002=1728
    integer,parameter::idx_H2Cl_0001=1273
    integer,parameter::idx_H2Cl_0002=1729
    integer,parameter::idx_HDCl_0001=1274
    integer,parameter::idx_HDCl_0002=1730
    integer,parameter::idx_D2Cl_0001=1275
    integer,parameter::idx_D2Cl_0002=1731
    integer,parameter::idx_CCl_0001=1276
    integer,parameter::idx_CCl_0002=1732
    integer,parameter::idx_ClO_0001=1277
    integer,parameter::idx_ClO_0002=1733
    integer,parameter::idx_C3H3_0001=1278
    integer,parameter::idx_C3H3_0002=1734
    integer,parameter::idx_C3H2D_0001=1279
    integer,parameter::idx_C3H2D_0002=1735
    integer,parameter::idx_C3HD2_0001=1280
    integer,parameter::idx_C3HD2_0002=1736
    integer,parameter::idx_C3D3_0001=1281
    integer,parameter::idx_C3D3_0002=1737
    integer,parameter::idx_C3H4_0001=1282
    integer,parameter::idx_C3H4_0002=1738
    integer,parameter::idx_C3H3D_0001=1283
    integer,parameter::idx_C3H3D_0002=1739
    integer,parameter::idx_C3H2D2_0001=1284
    integer,parameter::idx_C3H2D2_0002=1740
    integer,parameter::idx_HCCN_0001=1285
    integer,parameter::idx_HCCN_0002=1741
    integer,parameter::idx_DCCN_0001=1286
    integer,parameter::idx_DCCN_0002=1742
    integer,parameter::idx_C2H2N_0001=1287
    integer,parameter::idx_C2H2N_0002=1743
    integer,parameter::idx_C2H5OH_0001=1288
    integer,parameter::idx_C2H5OH_0002=1744
    integer,parameter::idx_C2H3D3_0001=1289
    integer,parameter::idx_C2H3D3_0002=1745
    integer,parameter::idx_C2H5D_0001=1290
    integer,parameter::idx_C2H5D_0002=1746
    integer,parameter::idx_C2H4D2_0001=1291
    integer,parameter::idx_C2H4D2_0002=1747
    integer,parameter::idx_C2H3N_0001=1292
    integer,parameter::idx_C2H3N_0002=1748
    integer,parameter::idx_C2H4O_0001=1293
    integer,parameter::idx_C2H4O_0002=1749
    integer,parameter::idx_CH2DCHO_0001=1294
    integer,parameter::idx_CH2DCHO_0002=1750
    integer,parameter::idx_CH3CDO_0001=1295
    integer,parameter::idx_CH3CDO_0002=1751
    integer,parameter::idx_CD3CHO_0001=1296
    integer,parameter::idx_CD3CHO_0002=1752
    integer,parameter::idx_CHD2CDO_0001=1297
    integer,parameter::idx_CHD2CDO_0002=1753
    integer,parameter::idx_CHD2CHO_0001=1298
    integer,parameter::idx_CHD2CHO_0002=1754
    integer,parameter::idx_CH2DCDO_0001=1299
    integer,parameter::idx_CH2DCDO_0002=1755
    integer,parameter::idx_C2H4D_0001=1300
    integer,parameter::idx_C2H4D_0002=1756
    integer,parameter::idx_C2H2D3_0001=1301
    integer,parameter::idx_C2H2D3_0002=1757
    integer,parameter::idx_C2H3D2_0001=1302
    integer,parameter::idx_C2H3D2_0002=1758
    integer,parameter::idx_C3H3N_0001=1303
    integer,parameter::idx_C3H3N_0002=1759
    integer,parameter::idx_H4C3N_0001=1304
    integer,parameter::idx_H4C3N_0002=1760
    integer,parameter::idx_C3D3N_0001=1305
    integer,parameter::idx_C3D3N_0002=1761
    integer,parameter::idx_HD3C3N_0001=1306
    integer,parameter::idx_HD3C3N_0002=1762
    integer,parameter::idx_C3H2DN_0001=1307
    integer,parameter::idx_C3H2DN_0002=1763
    integer,parameter::idx_H3DC3N_0001=1308
    integer,parameter::idx_H3DC3N_0002=1764
    integer,parameter::idx_C3HD2N_0001=1309
    integer,parameter::idx_C3HD2N_0002=1765
    integer,parameter::idx_H2D2C3N_0001=1310
    integer,parameter::idx_H2D2C3N_0002=1766
    integer,parameter::idx_HC3O_0001=1311
    integer,parameter::idx_HC3O_0002=1767
    integer,parameter::idx_DC3O_0001=1312
    integer,parameter::idx_DC3O_0002=1768
    integer,parameter::idx_C4H4_0001=1313
    integer,parameter::idx_C4H4_0002=1769
    integer,parameter::idx_CH5N_0001=1314
    integer,parameter::idx_CH5N_0002=1770
    integer,parameter::idx_CH3NH_0001=1315
    integer,parameter::idx_CH3NH_0002=1771
    integer,parameter::idx_FeD_0001=1316
    integer,parameter::idx_FeD_0002=1772
    integer,parameter::idx_H2C3N_0001=1317
    integer,parameter::idx_H2C3N_0002=1773
    integer,parameter::idx_H2C5N_0001=1318
    integer,parameter::idx_H2C5N_0002=1774
    integer,parameter::idx_H3C5N_0001=1319
    integer,parameter::idx_H3C5N_0002=1775
    integer,parameter::idx_H2C7N_0001=1320
    integer,parameter::idx_H2C7N_0002=1776
    integer,parameter::idx_H3C7N_0001=1321
    integer,parameter::idx_H3C7N_0002=1777
    integer,parameter::idx_H2C9N_0001=1322
    integer,parameter::idx_H2C9N_0002=1778
    integer,parameter::idx_H3C9N_0001=1323
    integer,parameter::idx_H3C9N_0002=1779
    integer,parameter::idx_H5C3N_0001=1324
    integer,parameter::idx_H5C3N_0002=1780
    integer,parameter::idx_C2H2O_0001=1325
    integer,parameter::idx_C2H2O_0002=1781
    integer,parameter::idx_C2HDO_0001=1326
    integer,parameter::idx_C2HDO_0002=1782
    integer,parameter::idx_C2D2O_0001=1327
    integer,parameter::idx_C2D2O_0002=1783
    integer,parameter::idx_HDC3N_0001=1328
    integer,parameter::idx_HDC3N_0002=1784
    integer,parameter::idx_D2C3N_0001=1329
    integer,parameter::idx_D2C3N_0002=1785
    integer,parameter::idx_H2C3O_0001=1330
    integer,parameter::idx_H2C3O_0002=1786
    integer,parameter::idx_HDC3O_0001=1331
    integer,parameter::idx_HDC3O_0002=1787
    integer,parameter::idx_D2C3O_0001=1332
    integer,parameter::idx_D2C3O_0002=1788
    integer,parameter::idx_C2HDN_0001=1333
    integer,parameter::idx_C2HDN_0002=1789
    integer,parameter::idx_C2D2N_0001=1334
    integer,parameter::idx_C2D2N_0002=1790
    integer,parameter::idx_CH3N_0001=1335
    integer,parameter::idx_CH3N_0002=1791
    integer,parameter::idx_N2H2_0001=1336
    integer,parameter::idx_N2H2_0002=1792
    integer,parameter::idx_N2HD_0001=1337
    integer,parameter::idx_N2HD_0002=1793
    integer,parameter::idx_N2D2_0001=1338
    integer,parameter::idx_N2D2_0002=1794
    integer,parameter::idx_NH2OH_0001=1339
    integer,parameter::idx_NH2OH_0002=1795
    integer,parameter::idx_NH2OD_0001=1340
    integer,parameter::idx_NH2OD_0002=1796
    integer,parameter::idx_N2H_0001=1341
    integer,parameter::idx_N2H_0002=1797
    integer,parameter::idx_N2D_0001=1342
    integer,parameter::idx_N2D_0002=1798
    integer,parameter::idx_CNH2_0001=1343
    integer,parameter::idx_CNH2_0002=1799
    integer,parameter::idx_CHNH2_0001=1344
    integer,parameter::idx_CHNH2_0002=1800
    integer,parameter::idx_HON_0001=1345
    integer,parameter::idx_HON_0002=1801
    integer,parameter::idx_DON_0001=1346
    integer,parameter::idx_DON_0002=1802
    integer,parameter::idx_NHNO_0001=1347
    integer,parameter::idx_NHNO_0002=1803
    integer,parameter::idx_CH2DN_0001=1348
    integer,parameter::idx_CH2DN_0002=1804
    integer,parameter::idx_CHD2N_0001=1349
    integer,parameter::idx_CHD2N_0002=1805
    integer,parameter::idx_CD3N_0001=1350
    integer,parameter::idx_CD3N_0002=1806
    integer,parameter::idx_NH2NO_0001=1351
    integer,parameter::idx_NH2NO_0002=1807
    integer,parameter::idx_CH3CO_0001=1352
    integer,parameter::idx_CH3CO_0002=1808
    integer,parameter::idx_CH2DOCH3_0001=1353
    integer,parameter::idx_CH2DOCH3_0002=1809
    integer,parameter::idx_CHD2OCH3_0001=1354
    integer,parameter::idx_CHD2OCH3_0002=1810
    integer,parameter::idx_CH2DOCH2D_0001=1355
    integer,parameter::idx_CH2DOCH2D_0002=1811
    integer,parameter::idx_CD3OCH3_0001=1356
    integer,parameter::idx_CD3OCH3_0002=1812
    integer,parameter::idx_CHD2OCH2D_0001=1357
    integer,parameter::idx_CHD2OCH2D_0002=1813
    integer,parameter::idx_DCOOCH3_0001=1358
    integer,parameter::idx_DCOOCH3_0002=1814
    integer,parameter::idx_HCOOCH2D_0001=1359
    integer,parameter::idx_HCOOCH2D_0002=1815
    integer,parameter::idx_DCOOCH2D_0001=1360
    integer,parameter::idx_DCOOCH2D_0002=1816
    integer,parameter::idx_HCOOCHD2_0001=1361
    integer,parameter::idx_HCOOCHD2_0002=1817
    integer,parameter::idx_DCOOCHD2_0001=1362
    integer,parameter::idx_DCOOCHD2_0002=1818
    integer,parameter::idx_HCOOCD3_0001=1363
    integer,parameter::idx_HCOOCD3_0002=1819
    integer,parameter::idx_DCOOCD3_0001=1364
    integer,parameter::idx_DCOOCD3_0002=1820
    integer,parameter::idx_NH2CO_0001=1365
    integer,parameter::idx_NH2CO_0002=1821
    integer,parameter::idx_NHDCO_0001=1366
    integer,parameter::idx_NHDCO_0002=1822
    integer,parameter::idx_ND2CO_0001=1367
    integer,parameter::idx_ND2CO_0002=1823
    integer,parameter::idx_CH2DOCHO_0001=1368
    integer,parameter::idx_CH2DOCHO_0002=1824
    integer,parameter::idx_CH3OCDO_0001=1369
    integer,parameter::idx_CH3OCDO_0002=1825
    integer,parameter::idx_CHD2OCHO_0001=1370
    integer,parameter::idx_CHD2OCHO_0002=1826
    integer,parameter::idx_CH2DOCDO_0001=1371
    integer,parameter::idx_CH2DOCDO_0002=1827
    integer,parameter::idx_CD3OCHO_0001=1372
    integer,parameter::idx_CD3OCHO_0002=1828
    integer,parameter::idx_CHD2OCDO_0001=1373
    integer,parameter::idx_CHD2OCDO_0002=1829
    integer,parameter::idx_CHDOCDO_0001=1374
    integer,parameter::idx_CHDOCDO_0002=1830
    integer,parameter::idx_CHDOCHDO_0001=1375
    integer,parameter::idx_CHDOCHDO_0002=1831
    integer,parameter::idx_CH3OCHD_0001=1376
    integer,parameter::idx_CH3OCHD_0002=1832
    integer,parameter::idx_CH2DOCH2_0001=1377
    integer,parameter::idx_CH2DOCH2_0002=1833
    integer,parameter::idx_CH2DOCHD_0001=1378
    integer,parameter::idx_CH2DOCHD_0002=1834
    integer,parameter::idx_CHD2OCH2_0001=1379
    integer,parameter::idx_CHD2OCH2_0002=1835
    integer,parameter::idx_CH3OCD2_0001=1380
    integer,parameter::idx_CH3OCD2_0002=1836
    integer,parameter::idx_CH2DOCD2_0001=1381
    integer,parameter::idx_CH2DOCD2_0002=1837
    integer,parameter::idx_CH3OCHD2_0001=1382
    integer,parameter::idx_CH3OCHD2_0002=1838
    integer,parameter::idx_CH3OCH2D_0001=1383
    integer,parameter::idx_CH3OCH2D_0002=1839
    integer,parameter::idx_CH2DCO_0001=1384
    integer,parameter::idx_CH2DCO_0002=1840
    integer,parameter::idx_CH3NHD_0001=1385
    integer,parameter::idx_CH3NHD_0002=1841
    integer,parameter::idx_CH2DNH2_0001=1386
    integer,parameter::idx_CH2DNH2_0002=1842
    integer,parameter::idx_P_0001=1387
    integer,parameter::idx_P_0002=1843
    integer,parameter::idx_PO_0001=1388
    integer,parameter::idx_PO_0002=1844
    integer,parameter::idx_PH_0001=1389
    integer,parameter::idx_PH_0002=1845
    integer,parameter::idx_PD_0001=1390
    integer,parameter::idx_PD_0002=1846
    integer,parameter::idx_PH2_0001=1391
    integer,parameter::idx_PH2_0002=1847
    integer,parameter::idx_PHD_0001=1392
    integer,parameter::idx_PHD_0002=1848
    integer,parameter::idx_PD2_0001=1393
    integer,parameter::idx_PD2_0002=1849
    integer,parameter::idx_PN_0001=1394
    integer,parameter::idx_PN_0002=1850
    integer,parameter::idx_CP_0001=1395
    integer,parameter::idx_CP_0002=1851
    integer,parameter::idx_CCP_0001=1396
    integer,parameter::idx_CCP_0002=1852
    integer,parameter::idx_C3P_0001=1397
    integer,parameter::idx_C3P_0002=1853
    integer,parameter::idx_C4P_0001=1398
    integer,parameter::idx_C4P_0002=1854
    integer,parameter::idx_CH2PH_0001=1399
    integer,parameter::idx_CH2PH_0002=1855
    integer,parameter::idx_CHD2PH_0001=1400
    integer,parameter::idx_CHD2PH_0002=1856
    integer,parameter::idx_CH2PD_0001=1401
    integer,parameter::idx_CH2PD_0002=1857
    integer,parameter::idx_CHDPD_0001=1402
    integer,parameter::idx_CHDPD_0002=1858
    integer,parameter::idx_CD2P2_0001=1403
    integer,parameter::idx_CD2P2_0002=1859
    integer,parameter::idx_HCP_0001=1404
    integer,parameter::idx_HCP_0002=1860
    integer,parameter::idx_DCP_0001=1405
    integer,parameter::idx_DCP_0002=1861
    integer,parameter::idx_HCCP_0001=1406
    integer,parameter::idx_HCCP_0002=1862
    integer,parameter::idx_DCCP_0001=1407
    integer,parameter::idx_DCCP_0002=1863
    integer,parameter::idx_HPO_0001=1408
    integer,parameter::idx_HPO_0002=1864
    integer,parameter::idx_DPO_0001=1409
    integer,parameter::idx_DPO_0002=1865
    integer,parameter::idx_surface_mask=1410
    integer,parameter::idx_mantle_mask=1411
    integer,parameter::idx_dummy=1412

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

  open(newunit=io, file=trim(fname), status="old")
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
