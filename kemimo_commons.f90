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
    ! WHEN: 2025-04-24 08:53:26
    ! CHANGESET: xxxxxxx
    ! BY: unknown@unknown

    !number of species
    integer,parameter::nmols=974
    !number of unique species
    integer,parameter::nmolsu=974
    !number of dust species
    integer,parameter::nmols_dust=451
    !number of reactions (nlayer*dust+gas)
    integer,parameter::nrea=9687
    !number of unique reactions (dust+gas)
    !number of dust-phase reactions
    integer,parameter::nreadust=2230
    !number of gas-phase reactions
    integer,parameter::nreagas=nrea-nreadust

    !monolayer thickness of each layer in model:
    integer,parameter::layerThickness=1
    !idx for mantle and surface species:
    integer,parameter::surface_start=521
    integer,parameter::surface_end=971
    integer,parameter::mantle_start=524
    integer,parameter::mantle_end=974
    !do swapping?
    logical,parameter::doSwap=.false.
    integer,parameter:: CO_desorption_idx=1585

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_ARRAYSIZE
  ! GRAIN PROPERTIES:
  real*8, parameter :: d2g = 1d-2 !dust/gas mass ratio
  real*8, parameter :: rho0 = 3d0 !bulk density, g/cm3
  real*8, parameter :: mu = 1.43 !mean molecular weight
  real*8, parameter :: app = 3d-8 !binding sites separation, cm
  real*8, parameter :: agrain = 1d-5
  real*8, parameter :: xdust = mu * d2g * pmass / ((4.0/3.0) * rho0 * agrain**3.0 * pi) !1.33d-12 ! dust density relative to n_H

  ! column densities for self-shielding:
  real*8 :: N_H2, N_CO, N_N2, N_HD, N_HI
  ! self-shielding factors/coefficients
  real*8 :: ss_CO, ss_H2, ss_HD, ss_N2

  ! UV RADIATION:
  real*8, parameter :: Fnot = 1d8 ! Draine ISRF
  real*8, parameter :: F_cr_not = 1d4 ! Cosmic-ray induced UV photons (never attenuated)
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
    ! WHEN: 2025-04-24 08:53:26
    ! CHANGESET: xxxxxxx
    ! BY: unknown@unknown

    character(len=maxVerbatimSize),parameter::speciesNames(nmolsu) = (/"GRAIN0_gas        ", &
        "GRAIN-_gas        ", &
        "E_gas             ", &
        "H_gas             ", &
        "H+_gas            ", &
        "H-_gas            ", &
        "H2_gas            ", &
        "H2+_gas           ", &
        "He_gas            ", &
        "He+_gas           ", &
        "O_gas             ", &
        "O+_gas            ", &
        "O-_gas            ", &
        "O2_gas            ", &
        "O2+_gas           ", &
        "O3_gas            ", &
        "OH_gas            ", &
        "OH+_gas           ", &
        "OH-_gas           ", &
        "H2O_gas           ", &
        "H2O+_gas          ", &
        "HO2+_gas          ", &
        "HOOH_gas          ", &
        "H3+_gas           ", &
        "H3O+_gas          ", &
        "F_gas             ", &
        "F+_gas            ", &
        "Cl_gas            ", &
        "Cl+_gas           ", &
        "CF+_gas           ", &
        "C_gas             ", &
        "C+_gas            ", &
        "C-_gas            ", &
        "C2_gas            ", &
        "C2+_gas           ", &
        "C2N+_gas          ", &
        "CO_gas            ", &
        "CO+_gas           ", &
        "CO2_gas           ", &
        "CO2+_gas          ", &
        "N_gas             ", &
        "N+_gas            ", &
        "N2_gas            ", &
        "N2+_gas           ", &
        "N2H+_gas          ", &
        "HCO_gas           ", &
        "HCO+_gas          ", &
        "H2CO_gas          ", &
        "CH2OH_gas         ", &
        "CH3O_gas          ", &
        "CH3OH_gas         ", &
        "CH_gas            ", &
        "CH2_gas           ", &
        "CH3_gas           ", &
        "CH4_gas           ", &
        "CH+_gas           ", &
        "CH2+_gas          ", &
        "CH3+_gas          ", &
        "CH4+_gas          ", &
        "CH5+_gas          ", &
        "HCOOH_gas         ", &
        "HOCO_gas          ", &
        "NH_gas            ", &
        "NH+_gas           ", &
        "NH2_gas           ", &
        "NH2+_gas          ", &
        "NH3_gas           ", &
        "NH3+_gas          ", &
        "NH4+_gas          ", &
        "Fe_gas            ", &
        "Fe+_gas           ", &
        "S_gas             ", &
        "S+_gas            ", &
        "S-_gas            ", &
        "C10_gas           ", &
        "C10+_gas          ", &
        "C10-_gas          ", &
        "C10H_gas          ", &
        "C10H+_gas         ", &
        "C10H-_gas         ", &
        "C10H2_gas         ", &
        "C10H2+_gas        ", &
        "C10H3+_gas        ", &
        "C10N_gas          ", &
        "C10N+_gas         ", &
        "C11_gas           ", &
        "C11+_gas          ", &
        "C2H+_gas          ", &
        "C2H2_gas          ", &
        "C2H2+_gas         ", &
        "C2H3_gas          ", &
        "C2H3+_gas         ", &
        "C2H4_gas          ", &
        "C2H4+_gas         ", &
        "C2H4O+_gas        ", &
        "C2H5_gas          ", &
        "C2H5+_gas         ", &
        "C2H5OH+_gas       ", &
        "C2H5OH2+_gas      ", &
        "C2H6_gas          ", &
        "C2H6+_gas         ", &
        "C2H6CO+_gas       ", &
        "C2H7+_gas         ", &
        "C2HO+_gas         ", &
        "C2N2+_gas         ", &
        "C2O+_gas          ", &
        "C2S+_gas          ", &
        "C3_gas            ", &
        "C3+_gas           ", &
        "C3-_gas           ", &
        "C3H+_gas          ", &
        "C3H3N+_gas        ", &
        "C3H3NH+_gas       ", &
        "C3H4_gas          ", &
        "C3H4+_gas         ", &
        "C3H5+_gas         ", &
        "C3H6OH+_gas       ", &
        "C3N_gas           ", &
        "C3N+_gas          ", &
        "C3N-_gas          ", &
        "C3O_gas           ", &
        "C3O+_gas          ", &
        "C3S_gas           ", &
        "C3S+_gas          ", &
        "C4_gas            ", &
        "C4+_gas           ", &
        "C4-_gas           ", &
        "C4H_gas           ", &
        "C4H+_gas          ", &
        "C4H-_gas          ", &
        "C4H2_gas          ", &
        "C4H2+_gas         ", &
        "C4H3_gas          ", &
        "C4H3+_gas         ", &
        "C4H4+_gas         ", &
        "C4H5+_gas         ", &
        "C4H7+_gas         ", &
        "C4N_gas           ", &
        "C4N+_gas          ", &
        "C4S_gas           ", &
        "C4S+_gas          ", &
        "C5_gas            ", &
        "C5+_gas           ", &
        "C5-_gas           ", &
        "C5H_gas           ", &
        "C5H+_gas          ", &
        "C5H-_gas          ", &
        "C5H2_gas          ", &
        "C5H2+_gas         ", &
        "C5H3_gas          ", &
        "C5H3+_gas         ", &
        "C5H3N+_gas        ", &
        "C5H4_gas          ", &
        "C5H4+_gas         ", &
        "C5H4N+_gas        ", &
        "C5H5+_gas         ", &
        "C5N_gas           ", &
        "C5N+_gas          ", &
        "C5O_gas           ", &
        "C6_gas            ", &
        "C6+_gas           ", &
        "C6-_gas           ", &
        "C6H_gas           ", &
        "C6H+_gas          ", &
        "C6H-_gas          ", &
        "C6H2_gas          ", &
        "C6H2+_gas         ", &
        "C6H3_gas          ", &
        "C6H3+_gas         ", &
        "C6H4_gas          ", &
        "C6H4+_gas         ", &
        "C6H5+_gas         ", &
        "C6H6_gas          ", &
        "C6H7+_gas         ", &
        "C6N_gas           ", &
        "C6N+_gas          ", &
        "C7_gas            ", &
        "C7+_gas           ", &
        "C7-_gas           ", &
        "C7H_gas           ", &
        "C7H+_gas          ", &
        "C7H-_gas          ", &
        "C7H2_gas          ", &
        "C7H2+_gas         ", &
        "C7H2N+_gas        ", &
        "C7H3_gas          ", &
        "C7H3+_gas         ", &
        "C7H4_gas          ", &
        "C7H4+_gas         ", &
        "C7H5+_gas         ", &
        "C7N_gas           ", &
        "C7N+_gas          ", &
        "C7O_gas           ", &
        "C8_gas            ", &
        "C8+_gas           ", &
        "C8-_gas           ", &
        "C8H_gas           ", &
        "C8H+_gas          ", &
        "C8H-_gas          ", &
        "C8H2_gas          ", &
        "C8H2+_gas         ", &
        "C8H3_gas          ", &
        "C8H3+_gas         ", &
        "C8H4_gas          ", &
        "C8H4+_gas         ", &
        "C8H4N+_gas        ", &
        "C8H5+_gas         ", &
        "C8N_gas           ", &
        "C8N+_gas          ", &
        "C9_gas            ", &
        "C9+_gas           ", &
        "C9-_gas           ", &
        "C9H_gas           ", &
        "C9H+_gas          ", &
        "C9H-_gas          ", &
        "C9H2_gas          ", &
        "C9H2+_gas         ", &
        "C9H2N+_gas        ", &
        "C9H3_gas          ", &
        "C9H3+_gas         ", &
        "C9H3N+_gas        ", &
        "C9H4_gas          ", &
        "C9H4+_gas         ", &
        "C9H5+_gas         ", &
        "C9HN+_gas         ", &
        "C9N_gas           ", &
        "C9N+_gas          ", &
        "C9O_gas           ", &
        "CCH_gas           ", &
        "CCN_gas           ", &
        "CCO_gas           ", &
        "CCP_gas           ", &
        "CCP+_gas          ", &
        "CCS_gas           ", &
        "CH2CCH_gas        ", &
        "CH2CHC2H_gas      ", &
        "CH2CHCHCH2_gas    ", &
        "CH2CHCN_gas       ", &
        "CH2CN+_gas        ", &
        "CH2NH_gas         ", &
        "CH2NH2_gas        ", &
        "CH2NH2+_gas       ", &
        "CH3C3N_gas        ", &
        "CH3C4H_gas        ", &
        "CH3C5N_gas        ", &
        "CH3C6H_gas        ", &
        "CH3C7N_gas        ", &
        "CH3CCH_gas        ", &
        "CH3CH2OH_gas      ", &
        "CH3CHCH2_gas      ", &
        "CH3CHO_gas        ", &
        "CH3CHOH+_gas      ", &
        "CH3CN_gas         ", &
        "CH3CN+_gas        ", &
        "CH3CNH+_gas       ", &
        "CH3CO+_gas        ", &
        "CH3COCH3_gas      ", &
        "CH3NH2_gas        ", &
        "CH3NH2+_gas       ", &
        "CH3NH3+_gas       ", &
        "CH3O2+_gas        ", &
        "CH3OCH2_gas       ", &
        "CH3OCH3_gas       ", &
        "CH3OCH3+_gas      ", &
        "CH3OCH4+_gas      ", &
        "CH3OH+_gas        ", &
        "CH3OH2+_gas       ", &
        "CN_gas            ", &
        "CN+_gas           ", &
        "CN-_gas           ", &
        "CNC+_gas          ", &
        "COOCH4+_gas       ", &
        "CS_gas            ", &
        "CS+_gas           ", &
        "HS_gas            ", &
        "HS+_gas           ", &
        "FeH_gas           ", &
        "H2C10N+_gas       ", &
        "H2C3O+_gas        ", &
        "H2C4N+_gas        ", &
        "H2C5N+_gas        ", &
        "H2C6N+_gas        ", &
        "H2C8N+_gas        ", &
        "H2CCN_gas         ", &
        "H2CCO_gas         ", &
        "H2CCO+_gas        ", &
        "H2CN_gas          ", &
        "H2CN+_gas         ", &
        "H2CO+_gas         ", &
        "H2COH+_gas        ", &
        "H2CS_gas          ", &
        "H2CS+_gas         ", &
        "H2NC+_gas         ", &
        "H2NO+_gas         ", &
        "H2S_gas           ", &
        "H2S+_gas          ", &
        "HDS_gas           ", &
        "HDS+_gas          ", &
        "H2S2+_gas         ", &
        "HDS2+_gas         ", &
        "D2S2+_gas         ", &
        "H3C3O+_gas        ", &
        "H3C4N+_gas        ", &
        "H3C4NH+_gas       ", &
        "H3C6NH+_gas       ", &
        "H3C7N+_gas        ", &
        "H3CS+_gas         ", &
        "H3S+_gas          ", &
        "H3S2+_gas         ", &
        "H5C2O2+_gas       ", &
        "HC10N+_gas        ", &
        "HC2N_gas          ", &
        "HC2N+_gas         ", &
        "HC2NCH+_gas       ", &
        "HC2O_gas          ", &
        "HC2S+_gas         ", &
        "HC3N_gas          ", &
        "HC3N+_gas         ", &
        "HC3NH+_gas        ", &
        "HC3ND+_gas        ", &
        "HC3O+_gas         ", &
        "HC3S+_gas         ", &
        "HC4N_gas          ", &
        "HC4N+_gas         ", &
        "HC4O+_gas         ", &
        "HC4S+_gas         ", &
        "HC5N_gas          ", &
        "HC5N+_gas         ", &
        "HC5O+_gas         ", &
        "HC6N_gas          ", &
        "HC6N+_gas         ", &
        "HC7N_gas          ", &
        "HC7N+_gas         ", &
        "HC7O+_gas         ", &
        "HC8N_gas          ", &
        "HC8N+_gas         ", &
        "HC9N_gas          ", &
        "HC9O+_gas         ", &
        "HCCNC_gas         ", &
        "HCN_gas           ", &
        "HCN+_gas          ", &
        "HCNCC_gas         ", &
        "HCNH+_gas         ", &
        "HCOOCH3_gas       ", &
        "HCOOH+_gas        ", &
        "HCS_gas           ", &
        "HCS+_gas          ", &
        "HN2O+_gas         ", &
        "HNC_gas           ", &
        "HNC+_gas          ", &
        "HNCCC_gas         ", &
        "HNCO_gas          ", &
        "HNCO+_gas         ", &
        "HNO_gas           ", &
        "HNO+_gas          ", &
        "HNS+_gas          ", &
        "HOC+_gas          ", &
        "HOCO+_gas         ", &
        "HOCS+_gas         ", &
        "HS2+_gas          ", &
        "HSO+_gas          ", &
        "HSO2+_gas         ", &
        "HSS_gas           ", &
        "HSSH_gas          ", &
        "N2O_gas           ", &
        "NC4N_gas          ", &
        "NC6N_gas          ", &
        "NC8N_gas          ", &
        "NCO+_gas          ", &
        "NH2CH2O+_gas      ", &
        "NH2CHO_gas        ", &
        "NH2CN_gas         ", &
        "NH2CND+_gas       ", &
        "NH2CNH+_gas       ", &
        "NO_gas            ", &
        "NO+_gas           ", &
        "NO2_gas           ", &
        "NO2+_gas          ", &
        "NS_gas            ", &
        "NS+_gas           ", &
        "OCN_gas           ", &
        "OCS_gas           ", &
        "OCS+_gas          ", &
        "S2_gas            ", &
        "S2+_gas           ", &
        "SO_gas            ", &
        "SO+_gas           ", &
        "SO2_gas           ", &
        "SO2+_gas          ", &
        "CH3OCHO_gas       ", &
        "Si_gas            ", &
        "Si+_gas           ", &
        "SiO_gas           ", &
        "SiO+_gas          ", &
        "SiO2_gas          ", &
        "SiS_gas           ", &
        "SiS+_gas          ", &
        "SiN_gas           ", &
        "SiN+_gas          ", &
        "SiC_gas           ", &
        "SiC+_gas          ", &
        "SiC2_gas          ", &
        "SiC2+_gas         ", &
        "SiH_gas           ", &
        "SiH2_gas          ", &
        "SiH2+_gas         ", &
        "SiH3_gas          ", &
        "SiH4_gas          ", &
        "P_gas             ", &
        "P+_gas            ", &
        "PH_gas            ", &
        "PH+_gas           ", &
        "PO_gas            ", &
        "PO+_gas           ", &
        "PH2_gas           ", &
        "PH2+_gas          ", &
        "C3P_gas           ", &
        "C4P_gas           ", &
        "C4P+_gas          ", &
        "CH2PH_gas         ", &
        "CH2Si+_gas        ", &
        "CHSi+_gas         ", &
        "CP_gas            ", &
        "CP+_gas           ", &
        "H2CSiCH_gas       ", &
        "H2F+_gas          ", &
        "H2PO+_gas         ", &
        "H2SiO_gas         ", &
        "H2SiO+_gas        ", &
        "H3SiO+_gas        ", &
        "HCCP_gas          ", &
        "HCCSi_gas         ", &
        "HCP_gas           ", &
        "HCP+_gas          ", &
        "HCSi_gas          ", &
        "HF_gas            ", &
        "HF+_gas           ", &
        "HNSi_gas          ", &
        "HNSi+_gas         ", &
        "HPN+_gas          ", &
        "HPO_gas           ", &
        "HPO+_gas          ", &
        "HSiNH+_gas        ", &
        "HSiO+_gas         ", &
        "HSiO2+_gas        ", &
        "HSiS+_gas         ", &
        "HeH+_gas          ", &
        "PC2H+_gas         ", &
        "PC2H2+_gas        ", &
        "PC2H3+_gas        ", &
        "PC2H4+_gas        ", &
        "PC3H+_gas         ", &
        "PC4H+_gas         ", &
        "PC4H2+_gas        ", &
        "PCH2+_gas         ", &
        "PCH3+_gas         ", &
        "PCH4+_gas         ", &
        "PH3+_gas          ", &
        "PN_gas            ", &
        "PN+_gas           ", &
        "PNH2+_gas         ", &
        "PNH3+_gas         ", &
        "SiC2CH3_gas       ", &
        "SiC2H+_gas        ", &
        "SiC2H2+_gas       ", &
        "SiC2H3+_gas       ", &
        "SiC3D_gas         ", &
        "SiC3H_gas         ", &
        "SiC3H+_gas        ", &
        "SiC3H2+_gas       ", &
        "SiC3H5_gas        ", &
        "SiC4_gas          ", &
        "SiC4+_gas         ", &
        "SiC4H_gas         ", &
        "SiC4H+_gas        ", &
        "SiC6H_gas         ", &
        "SiC8H_gas         ", &
        "SiCH2_gas         ", &
        "SiCH3_gas         ", &
        "SiCH3+_gas        ", &
        "SiCH4+_gas        ", &
        "SiF+_gas          ", &
        "SiH+_gas          ", &
        "SiH3+_gas         ", &
        "SiH4+_gas         ", &
        "SiH5+_gas         ", &
        "SiNC_gas          ", &
        "SiNC+_gas         ", &
        "SiNCH+_gas        ", &
        "c_HCCHSi_gas      ", &
        "c_SiC2_gas        ", &
        "l_SiC3_gas        ", &
        "l_SiC3+_gas       ", &
        "l_C3H_gas         ", &
        "c_C3H_gas         ", &
        "l_C3H2_gas        ", &
        "c_C3H2_gas        ", &
        "l_C3H2+_gas       ", &
        "c_C3H2+_gas       ", &
        "l_C3H3+_gas       ", &
        "c_C3H3+_gas       ", &
        "Mg_gas            ", &
        "Mg+_gas           ", &
        "MgH_gas           ", &
        "MgH2_gas          ", &
        "HCl_gas           ", &
        "HCl+_gas          ", &
        "H2Cl_gas          ", &
        "H2Cl+_gas         ", &
        "CCl_gas           ", &
        "CCl+_gas          ", &
        "ClO_gas           ", &
        "ClO+_gas          ", &
        "H2CCl+_gas        ", &
        "Na_gas            ", &
        "Na+_gas           ", &
        "NaH_gas           ", &
        "NaOH_gas          ", &
        "NaH2+_gas         ", &
        "NaH2O+_gas        ", &
        "H_surface         ", &
        "D_surface         ", &
        "H2_surface        ", &
        "D2_surface        ", &
        "HD_surface        ", &
        "He_surface        ", &
        "O_surface         ", &
        "O2_surface        ", &
        "O3_surface        ", &
        "OH_surface        ", &
        "OD_surface        ", &
        "H2O_surface       ", &
        "HDO_surface       ", &
        "D2O_surface       ", &
        "O2H_surface       ", &
        "O2D_surface       ", &
        "HOOH_surface      ", &
        "HOOD_surface      ", &
        "DOOD_surface      ", &
        "Fe_surface        ", &
        "FeH_surface       ", &
        "N_surface         ", &
        "S_surface         ", &
        "C_surface         ", &
        "C2_surface        ", &
        "CO_surface        ", &
        "HCO_surface       ", &
        "DCO_surface       ", &
        "H2CO_surface      ", &
        "HDCO_surface      ", &
        "D2CO_surface      ", &
        "CH2OH_surface     ", &
        "CD2OD_surface     ", &
        "CH2OD_surface     ", &
        "CHDOH_surface     ", &
        "CHDOD_surface     ", &
        "CD2OH_surface     ", &
        "CH3O_surface      ", &
        "CHD2O_surface     ", &
        "CH2DO_surface     ", &
        "CD3O_surface      ", &
        "CH3OH_surface     ", &
        "CH3OD_surface     ", &
        "CHD2OH_surface    ", &
        "CHD2OD_surface    ", &
        "CH2DOH_surface    ", &
        "CH2DOD_surface    ", &
        "CD3OD_surface     ", &
        "CD3OH_surface     ", &
        "CH_surface        ", &
        "CD_surface        ", &
        "CH2_surface       ", &
        "CHD_surface       ", &
        "CD2_surface       ", &
        "CH3_surface       ", &
        "CH2D_surface      ", &
        "CHD2_surface      ", &
        "CD3_surface       ", &
        "CH4_surface       ", &
        "CH3D_surface      ", &
        "CH2D2_surface     ", &
        "CHD3_surface      ", &
        "CD4_surface       ", &
        "CO2_surface       ", &
        "HCOOH_surface     ", &
        "HCOOD_surface     ", &
        "DCOOH_surface     ", &
        "DCOOD_surface     ", &
        "HOCO_surface      ", &
        "DOCO_surface      ", &
        "NH_surface        ", &
        "ND_surface        ", &
        "NH2_surface       ", &
        "NHD_surface       ", &
        "ND2_surface       ", &
        "NH3_surface       ", &
        "NH2D_surface      ", &
        "NHD2_surface      ", &
        "ND3_surface       ", &
        "C10_surface       ", &
        "C10H_surface      ", &
        "C10H2_surface     ", &
        "C10N_surface      ", &
        "C11_surface       ", &
        "C2H2_surface      ", &
        "C2HD_surface      ", &
        "C2D2_surface      ", &
        "C2H3_surface      ", &
        "C2H2D_surface     ", &
        "C2HD2_surface     ", &
        "C2D3_surface      ", &
        "C2H4_surface      ", &
        "C2H3D_surface     ", &
        "C2H2D2_surface    ", &
        "C2HD3_surface     ", &
        "C2D4_surface      ", &
        "C2H5_surface      ", &
        "C2H6_surface      ", &
        "C3_surface        ", &
        "C3N_surface       ", &
        "C3O_surface       ", &
        "C3S_surface       ", &
        "C4_surface        ", &
        "C4H_surface       ", &
        "C4D_surface       ", &
        "C4H2_surface      ", &
        "C4HD_surface      ", &
        "C4D2_surface      ", &
        "C4H3_surface      ", &
        "C4N_surface       ", &
        "C4S_surface       ", &
        "C5_surface        ", &
        "C5H_surface       ", &
        "C5D_surface       ", &
        "C5H2_surface      ", &
        "C5H3_surface      ", &
        "C5H4_surface      ", &
        "C5N_surface       ", &
        "C5O_surface       ", &
        "C6_surface        ", &
        "C6H_surface       ", &
        "C6H2_surface      ", &
        "C6H3_surface      ", &
        "C6H4_surface      ", &
        "C6H6_surface      ", &
        "C6N_surface       ", &
        "C7_surface        ", &
        "C7H_surface       ", &
        "C7H2_surface      ", &
        "C7H3_surface      ", &
        "C7H4_surface      ", &
        "C7N_surface       ", &
        "C7O_surface       ", &
        "C8_surface        ", &
        "C8H_surface       ", &
        "C8H2_surface      ", &
        "C8H3_surface      ", &
        "C8H4_surface      ", &
        "C8N_surface       ", &
        "C9_surface        ", &
        "C9H_surface       ", &
        "C9H2_surface      ", &
        "C9H3_surface      ", &
        "C9H4_surface      ", &
        "C9N_surface       ", &
        "C9O_surface       ", &
        "CCH_surface       ", &
        "CCD_surface       ", &
        "CCN_surface       ", &
        "CCO_surface       ", &
        "CCS_surface       ", &
        "CD3CN_surface     ", &
        "CH2CCH_surface    ", &
        "CH2CCD_surface    ", &
        "CHDCCH_surface    ", &
        "CHDCCD_surface    ", &
        "CD2CCH_surface    ", &
        "CD2CCD_surface    ", &
        "CH2CHC2H_surface  ", &
        "CH2CHCHCH2_surface", &
        "CH2CHCN_surface   ", &
        "CH2NH_surface     ", &
        "CHDNH_surface     ", &
        "CHDND_surface     ", &
        "CH2ND_surface     ", &
        "CD2NH_surface     ", &
        "CD2ND_surface     ", &
        "CH2NH2_surface    ", &
        "CH2NHD_surface    ", &
        "CHDNH2_surface    ", &
        "CHDNHD_surface    ", &
        "CHDND2_surface    ", &
        "CH2ND2_surface    ", &
        "CD2NH2_surface    ", &
        "CD2NHD_surface    ", &
        "CD2ND2_surface    ", &
        "CH3C3N_surface    ", &
        "CH3C4H_surface    ", &
        "CH3C5N_surface    ", &
        "CH3C6H_surface    ", &
        "CH3C7N_surface    ", &
        "CH3CCH_surface    ", &
        "CH3CH2OH_surface  ", &
        "CH3CHCH2_surface  ", &
        "CH3CHO_surface    ", &
        "CH3CN_surface     ", &
        "CH2DCN_surface    ", &
        "CHD2CN_surface    ", &
        "CH3COCH3_surface  ", &
        "CH3NH2_surface    ", &
        "CH3OCH2_surface   ", &
        "CH3OCH3_surface   ", &
        "CN_surface        ", &
        "CS_surface        ", &
        "H2CCN_surface     ", &
        "HDCCN_surface     ", &
        "D2CCN_surface     ", &
        "H2CCO_surface     ", &
        "HDCCO_surface     ", &
        "D2CCO_surface     ", &
        "H2CN_surface      ", &
        "HDCN_surface      ", &
        "D2CN_surface      ", &
        "H2CS_surface      ", &
        "HDCS_surface      ", &
        "D2CS_surface      ", &
        "H2S_surface       ", &
        "HDS_surface       ", &
        "D2S_surface       ", &
        "HC2O_surface      ", &
        "DC2O_surface      ", &
        "HC3N_surface      ", &
        "DC3N_surface      ", &
        "HC4N_surface      ", &
        "DC4N_surface      ", &
        "HC5N_surface      ", &
        "HC6N_surface      ", &
        "HC7N_surface      ", &
        "HC8N_surface      ", &
        "HC9N_surface      ", &
        "HCCNC_surface     ", &
        "DCCNC_surface     ", &
        "HCN_surface       ", &
        "DCN_surface       ", &
        "HCNCC_surface     ", &
        "DCNCC_surface     ", &
        "HCOOCH3_surface   ", &
        "HCS_surface       ", &
        "DCS_surface       ", &
        "HNC_surface       ", &
        "DNC_surface       ", &
        "HNCCC_surface     ", &
        "DNCCC_surface     ", &
        "HNCO_surface      ", &
        "DNCO_surface      ", &
        "HNO_surface       ", &
        "DNO_surface       ", &
        "HS_surface        ", &
        "DS_surface        ", &
        "N2_surface        ", &
        "N2O_surface       ", &
        "NC4N_surface      ", &
        "NC6N_surface      ", &
        "NC8N_surface      ", &
        "NH2CHO_surface    ", &
        "NH2CDO_surface    ", &
        "NHDCHO_surface    ", &
        "NHDCDO_surface    ", &
        "ND2CHO_surface    ", &
        "ND2CDO_surface    ", &
        "NH2CN_surface     ", &
        "NHDCN_surface     ", &
        "ND2CN_surface     ", &
        "HSO_surface       ", &
        "DSO_surface       ", &
        "HSS_surface       ", &
        "DSS_surface       ", &
        "HSSH_surface      ", &
        "HSSD_surface      ", &
        "DSSH_surface      ", &
        "DSSD_surface      ", &
        "NO_surface        ", &
        "NO2_surface       ", &
        "NS_surface        ", &
        "OCN_surface       ", &
        "OCS_surface       ", &
        "S2_surface        ", &
        "SO_surface        ", &
        "SO2_surface       ", &
        "CH3OCHO_surface   ", &
        "Si_surface        ", &
        "SiS_surface       ", &
        "SiN_surface       ", &
        "SiC_surface       ", &
        "SiH_surface       ", &
        "SiH2_surface      ", &
        "SiH3_surface      ", &
        "SiH4_surface      ", &
        "SiC2CH3_surface   ", &
        "SiC3H_surface     ", &
        "SiC3H5_surface    ", &
        "SiC4_surface      ", &
        "SiC4H_surface     ", &
        "SiC6H_surface     ", &
        "SiC8H_surface     ", &
        "c_HCCHSi_surface  ", &
        "c_SiC2_surface    ", &
        "l_C3H_surface     ", &
        "l_C3D_surface     ", &
        "c_C3H_surface     ", &
        "c_C3D_surface     ", &
        "l_SiC3_surface    ", &
        "l_C3H2_surface    ", &
        "l_C3HD_surface    ", &
        "l_C3D2_surface    ", &
        "c_C3H2_surface    ", &
        "c_C3HD_surface    ", &
        "c_C3D2_surface    ", &
        "Mg_surface        ", &
        "MgH_surface       ", &
        "MgH2_surface      ", &
        "Na_surface        ", &
        "NaH_surface       ", &
        "F_surface         ", &
        "HF_surface        ", &
        "DF_surface        ", &
        "MgD_surface       ", &
        "MgHD_surface      ", &
        "MgD2_surface      ", &
        "NaD_surface       ", &
        "SiD_surface       ", &
        "Cl_surface        ", &
        "HCl_surface       ", &
        "DCl_surface       ", &
        "H2Cl_surface      ", &
        "HDCl_surface      ", &
        "D2Cl_surface      ", &
        "CCl_surface       ", &
        "ClO_surface       ", &
        "C3H3_surface      ", &
        "C3H2D_surface     ", &
        "C3HD2_surface     ", &
        "C3D3_surface      ", &
        "C3H4_surface      ", &
        "C3H3D_surface     ", &
        "C3H2D2_surface    ", &
        "HCCN_surface      ", &
        "DCCN_surface      ", &
        "C2H2N_surface     ", &
        "C2H5OH_surface    ", &
        "C2H3D3_surface    ", &
        "C2H5D_surface     ", &
        "C2H4D2_surface    ", &
        "C2H3N_surface     ", &
        "C2H4O_surface     ", &
        "CH2DCHO_surface   ", &
        "CH3CDO_surface    ", &
        "CD3CHO_surface    ", &
        "CHD2CDO_surface   ", &
        "CHD2CHO_surface   ", &
        "CH2DCDO_surface   ", &
        "C2H4D_surface     ", &
        "C2H2D3_surface    ", &
        "C2H3D2_surface    ", &
        "C3H3N_surface     ", &
        "H4C3N_surface     ", &
        "C3D3N_surface     ", &
        "HD3C3N_surface    ", &
        "C3H2DN_surface    ", &
        "H3DC3N_surface    ", &
        "C3HD2N_surface    ", &
        "H2D2C3N_surface   ", &
        "HC3O_surface      ", &
        "DC3O_surface      ", &
        "C4H4_surface      ", &
        "CH5N_surface      ", &
        "CH3NH_surface     ", &
        "FeD_surface       ", &
        "H2C3N_surface     ", &
        "H2C5N_surface     ", &
        "H3C5N_surface     ", &
        "H2C7N_surface     ", &
        "H3C7N_surface     ", &
        "H2C9N_surface     ", &
        "H3C9N_surface     ", &
        "H5C3N_surface     ", &
        "C2H2O_surface     ", &
        "C2HDO_surface     ", &
        "C2D2O_surface     ", &
        "HDC3N_surface     ", &
        "D2C3N_surface     ", &
        "H2C3O_surface     ", &
        "HDC3O_surface     ", &
        "D2C3O_surface     ", &
        "C2HDN_surface     ", &
        "C2D2N_surface     ", &
        "CH3N_surface      ", &
        "N2H2_surface      ", &
        "N2HD_surface      ", &
        "N2D2_surface      ", &
        "NH2OH_surface     ", &
        "NH2OD_surface     ", &
        "N2H_surface       ", &
        "N2D_surface       ", &
        "CNH2_surface      ", &
        "CHNH2_surface     ", &
        "HON_surface       ", &
        "DON_surface       ", &
        "NHNO_surface      ", &
        "CH2DN_surface     ", &
        "CHD2N_surface     ", &
        "CD3N_surface      ", &
        "NH2NO_surface     ", &
        "CH3CO_surface     ", &
        "CH2DOCH3_surface  ", &
        "CHD2OCH3_surface  ", &
        "CH2DOCH2D_surface ", &
        "CD3OCH3_surface   ", &
        "CHD2OCH2D_surface ", &
        "DCOOCH3_surface   ", &
        "HCOOCH2D_surface  ", &
        "DCOOCH2D_surface  ", &
        "HCOOCHD2_surface  ", &
        "DCOOCHD2_surface  ", &
        "HCOOCD3_surface   ", &
        "DCOOCD3_surface   ", &
        "NH2CO_surface     ", &
        "NHDCO_surface     ", &
        "ND2CO_surface     ", &
        "CH2DOCHO_surface  ", &
        "CH3OCDO_surface   ", &
        "CHD2OCHO_surface  ", &
        "CH2DOCDO_surface  ", &
        "CD3OCHO_surface   ", &
        "CHD2OCDO_surface  ", &
        "CHDOCDO_surface   ", &
        "CHDOCHDO_surface  ", &
        "CH3OCHD_surface   ", &
        "CH2DOCH2_surface  ", &
        "CH2DOCHD_surface  ", &
        "CHD2OCH2_surface  ", &
        "CH3OCD2_surface   ", &
        "CH2DOCD2_surface  ", &
        "CH3OCHD2_surface  ", &
        "CH3OCH2D_surface  ", &
        "CH2DCO_surface    ", &
        "CH3NHD_surface    ", &
        "CH2DNH2_surface   ", &
        "P_surface         ", &
        "PO_surface        ", &
        "PH_surface        ", &
        "PD_surface        ", &
        "PH2_surface       ", &
        "PHD_surface       ", &
        "PD2_surface       ", &
        "PN_surface        ", &
        "CP_surface        ", &
        "CCP_surface       ", &
        "C3P_surface       ", &
        "C4P_surface       ", &
        "CH2PH_surface     ", &
        "CHD2PH_surface    ", &
        "CH2PD_surface     ", &
        "CHDPD_surface     ", &
        "CD2P2_surface     ", &
        "HCP_surface       ", &
        "DCP_surface       ", &
        "HCCP_surface      ", &
        "DCCP_surface      ", &
        "HPO_surface       ", &
        "DPO_surface       ", &
        "surface_mask      ", &
        "mantle_mask       ", &
        "dummy             "/)
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_SPECIESNAMES

  !!BEGIN_IDXLIST
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2025-04-24 08:53:26
    ! CHANGESET: xxxxxxx
    ! BY: unknown@unknown

    integer,parameter::idx_GRAIN0_gas=1
    integer,parameter::idx_GRAINk_gas=2
    integer,parameter::idx_E_gas=3
    integer,parameter::idx_H_gas=4
    integer,parameter::idx_Hj_gas=5
    integer,parameter::idx_Hk_gas=6
    integer,parameter::idx_H2_gas=7
    integer,parameter::idx_H2j_gas=8
    integer,parameter::idx_He_gas=9
    integer,parameter::idx_Hej_gas=10
    integer,parameter::idx_O_gas=11
    integer,parameter::idx_Oj_gas=12
    integer,parameter::idx_Ok_gas=13
    integer,parameter::idx_O2_gas=14
    integer,parameter::idx_O2j_gas=15
    integer,parameter::idx_O3_gas=16
    integer,parameter::idx_OH_gas=17
    integer,parameter::idx_OHj_gas=18
    integer,parameter::idx_OHk_gas=19
    integer,parameter::idx_H2O_gas=20
    integer,parameter::idx_H2Oj_gas=21
    integer,parameter::idx_HO2j_gas=22
    integer,parameter::idx_HOOH_gas=23
    integer,parameter::idx_H3j_gas=24
    integer,parameter::idx_H3Oj_gas=25
    integer,parameter::idx_F_gas=26
    integer,parameter::idx_Fj_gas=27
    integer,parameter::idx_Cl_gas=28
    integer,parameter::idx_Clj_gas=29
    integer,parameter::idx_CFj_gas=30
    integer,parameter::idx_C_gas=31
    integer,parameter::idx_Cj_gas=32
    integer,parameter::idx_Ck_gas=33
    integer,parameter::idx_C2_gas=34
    integer,parameter::idx_C2j_gas=35
    integer,parameter::idx_C2Nj_gas=36
    integer,parameter::idx_CO_gas=37
    integer,parameter::idx_COj_gas=38
    integer,parameter::idx_CO2_gas=39
    integer,parameter::idx_CO2j_gas=40
    integer,parameter::idx_N_gas=41
    integer,parameter::idx_Nj_gas=42
    integer,parameter::idx_N2_gas=43
    integer,parameter::idx_N2j_gas=44
    integer,parameter::idx_N2Hj_gas=45
    integer,parameter::idx_HCO_gas=46
    integer,parameter::idx_HCOj_gas=47
    integer,parameter::idx_H2CO_gas=48
    integer,parameter::idx_CH2OH_gas=49
    integer,parameter::idx_CH3O_gas=50
    integer,parameter::idx_CH3OH_gas=51
    integer,parameter::idx_CH_gas=52
    integer,parameter::idx_CH2_gas=53
    integer,parameter::idx_CH3_gas=54
    integer,parameter::idx_CH4_gas=55
    integer,parameter::idx_CHj_gas=56
    integer,parameter::idx_CH2j_gas=57
    integer,parameter::idx_CH3j_gas=58
    integer,parameter::idx_CH4j_gas=59
    integer,parameter::idx_CH5j_gas=60
    integer,parameter::idx_HCOOH_gas=61
    integer,parameter::idx_HOCO_gas=62
    integer,parameter::idx_NH_gas=63
    integer,parameter::idx_NHj_gas=64
    integer,parameter::idx_NH2_gas=65
    integer,parameter::idx_NH2j_gas=66
    integer,parameter::idx_NH3_gas=67
    integer,parameter::idx_NH3j_gas=68
    integer,parameter::idx_NH4j_gas=69
    integer,parameter::idx_Fe_gas=70
    integer,parameter::idx_Fej_gas=71
    integer,parameter::idx_S_gas=72
    integer,parameter::idx_Sj_gas=73
    integer,parameter::idx_Sk_gas=74
    integer,parameter::idx_C10_gas=75
    integer,parameter::idx_C10j_gas=76
    integer,parameter::idx_C10k_gas=77
    integer,parameter::idx_C10H_gas=78
    integer,parameter::idx_C10Hj_gas=79
    integer,parameter::idx_C10Hk_gas=80
    integer,parameter::idx_C10H2_gas=81
    integer,parameter::idx_C10H2j_gas=82
    integer,parameter::idx_C10H3j_gas=83
    integer,parameter::idx_C10N_gas=84
    integer,parameter::idx_C10Nj_gas=85
    integer,parameter::idx_C11_gas=86
    integer,parameter::idx_C11j_gas=87
    integer,parameter::idx_C2Hj_gas=88
    integer,parameter::idx_C2H2_gas=89
    integer,parameter::idx_C2H2j_gas=90
    integer,parameter::idx_C2H3_gas=91
    integer,parameter::idx_C2H3j_gas=92
    integer,parameter::idx_C2H4_gas=93
    integer,parameter::idx_C2H4j_gas=94
    integer,parameter::idx_C2H4Oj_gas=95
    integer,parameter::idx_C2H5_gas=96
    integer,parameter::idx_C2H5j_gas=97
    integer,parameter::idx_C2H5OHj_gas=98
    integer,parameter::idx_C2H5OH2j_gas=99
    integer,parameter::idx_C2H6_gas=100
    integer,parameter::idx_C2H6j_gas=101
    integer,parameter::idx_C2H6COj_gas=102
    integer,parameter::idx_C2H7j_gas=103
    integer,parameter::idx_C2HOj_gas=104
    integer,parameter::idx_C2N2j_gas=105
    integer,parameter::idx_C2Oj_gas=106
    integer,parameter::idx_C2Sj_gas=107
    integer,parameter::idx_C3_gas=108
    integer,parameter::idx_C3j_gas=109
    integer,parameter::idx_C3k_gas=110
    integer,parameter::idx_C3Hj_gas=111
    integer,parameter::idx_C3H3Nj_gas=112
    integer,parameter::idx_C3H3NHj_gas=113
    integer,parameter::idx_C3H4_gas=114
    integer,parameter::idx_C3H4j_gas=115
    integer,parameter::idx_C3H5j_gas=116
    integer,parameter::idx_C3H6OHj_gas=117
    integer,parameter::idx_C3N_gas=118
    integer,parameter::idx_C3Nj_gas=119
    integer,parameter::idx_C3Nk_gas=120
    integer,parameter::idx_C3O_gas=121
    integer,parameter::idx_C3Oj_gas=122
    integer,parameter::idx_C3S_gas=123
    integer,parameter::idx_C3Sj_gas=124
    integer,parameter::idx_C4_gas=125
    integer,parameter::idx_C4j_gas=126
    integer,parameter::idx_C4k_gas=127
    integer,parameter::idx_C4H_gas=128
    integer,parameter::idx_C4Hj_gas=129
    integer,parameter::idx_C4Hk_gas=130
    integer,parameter::idx_C4H2_gas=131
    integer,parameter::idx_C4H2j_gas=132
    integer,parameter::idx_C4H3_gas=133
    integer,parameter::idx_C4H3j_gas=134
    integer,parameter::idx_C4H4j_gas=135
    integer,parameter::idx_C4H5j_gas=136
    integer,parameter::idx_C4H7j_gas=137
    integer,parameter::idx_C4N_gas=138
    integer,parameter::idx_C4Nj_gas=139
    integer,parameter::idx_C4S_gas=140
    integer,parameter::idx_C4Sj_gas=141
    integer,parameter::idx_C5_gas=142
    integer,parameter::idx_C5j_gas=143
    integer,parameter::idx_C5k_gas=144
    integer,parameter::idx_C5H_gas=145
    integer,parameter::idx_C5Hj_gas=146
    integer,parameter::idx_C5Hk_gas=147
    integer,parameter::idx_C5H2_gas=148
    integer,parameter::idx_C5H2j_gas=149
    integer,parameter::idx_C5H3_gas=150
    integer,parameter::idx_C5H3j_gas=151
    integer,parameter::idx_C5H3Nj_gas=152
    integer,parameter::idx_C5H4_gas=153
    integer,parameter::idx_C5H4j_gas=154
    integer,parameter::idx_C5H4Nj_gas=155
    integer,parameter::idx_C5H5j_gas=156
    integer,parameter::idx_C5N_gas=157
    integer,parameter::idx_C5Nj_gas=158
    integer,parameter::idx_C5O_gas=159
    integer,parameter::idx_C6_gas=160
    integer,parameter::idx_C6j_gas=161
    integer,parameter::idx_C6k_gas=162
    integer,parameter::idx_C6H_gas=163
    integer,parameter::idx_C6Hj_gas=164
    integer,parameter::idx_C6Hk_gas=165
    integer,parameter::idx_C6H2_gas=166
    integer,parameter::idx_C6H2j_gas=167
    integer,parameter::idx_C6H3_gas=168
    integer,parameter::idx_C6H3j_gas=169
    integer,parameter::idx_C6H4_gas=170
    integer,parameter::idx_C6H4j_gas=171
    integer,parameter::idx_C6H5j_gas=172
    integer,parameter::idx_C6H6_gas=173
    integer,parameter::idx_C6H7j_gas=174
    integer,parameter::idx_C6N_gas=175
    integer,parameter::idx_C6Nj_gas=176
    integer,parameter::idx_C7_gas=177
    integer,parameter::idx_C7j_gas=178
    integer,parameter::idx_C7k_gas=179
    integer,parameter::idx_C7H_gas=180
    integer,parameter::idx_C7Hj_gas=181
    integer,parameter::idx_C7Hk_gas=182
    integer,parameter::idx_C7H2_gas=183
    integer,parameter::idx_C7H2j_gas=184
    integer,parameter::idx_C7H2Nj_gas=185
    integer,parameter::idx_C7H3_gas=186
    integer,parameter::idx_C7H3j_gas=187
    integer,parameter::idx_C7H4_gas=188
    integer,parameter::idx_C7H4j_gas=189
    integer,parameter::idx_C7H5j_gas=190
    integer,parameter::idx_C7N_gas=191
    integer,parameter::idx_C7Nj_gas=192
    integer,parameter::idx_C7O_gas=193
    integer,parameter::idx_C8_gas=194
    integer,parameter::idx_C8j_gas=195
    integer,parameter::idx_C8k_gas=196
    integer,parameter::idx_C8H_gas=197
    integer,parameter::idx_C8Hj_gas=198
    integer,parameter::idx_C8Hk_gas=199
    integer,parameter::idx_C8H2_gas=200
    integer,parameter::idx_C8H2j_gas=201
    integer,parameter::idx_C8H3_gas=202
    integer,parameter::idx_C8H3j_gas=203
    integer,parameter::idx_C8H4_gas=204
    integer,parameter::idx_C8H4j_gas=205
    integer,parameter::idx_C8H4Nj_gas=206
    integer,parameter::idx_C8H5j_gas=207
    integer,parameter::idx_C8N_gas=208
    integer,parameter::idx_C8Nj_gas=209
    integer,parameter::idx_C9_gas=210
    integer,parameter::idx_C9j_gas=211
    integer,parameter::idx_C9k_gas=212
    integer,parameter::idx_C9H_gas=213
    integer,parameter::idx_C9Hj_gas=214
    integer,parameter::idx_C9Hk_gas=215
    integer,parameter::idx_C9H2_gas=216
    integer,parameter::idx_C9H2j_gas=217
    integer,parameter::idx_C9H2Nj_gas=218
    integer,parameter::idx_C9H3_gas=219
    integer,parameter::idx_C9H3j_gas=220
    integer,parameter::idx_C9H3Nj_gas=221
    integer,parameter::idx_C9H4_gas=222
    integer,parameter::idx_C9H4j_gas=223
    integer,parameter::idx_C9H5j_gas=224
    integer,parameter::idx_C9HNj_gas=225
    integer,parameter::idx_C9N_gas=226
    integer,parameter::idx_C9Nj_gas=227
    integer,parameter::idx_C9O_gas=228
    integer,parameter::idx_CCH_gas=229
    integer,parameter::idx_CCN_gas=230
    integer,parameter::idx_CCO_gas=231
    integer,parameter::idx_CCP_gas=232
    integer,parameter::idx_CCPj_gas=233
    integer,parameter::idx_CCS_gas=234
    integer,parameter::idx_CH2CCH_gas=235
    integer,parameter::idx_CH2CHC2H_gas=236
    integer,parameter::idx_CH2CHCHCH2_gas=237
    integer,parameter::idx_CH2CHCN_gas=238
    integer,parameter::idx_CH2CNj_gas=239
    integer,parameter::idx_CH2NH_gas=240
    integer,parameter::idx_CH2NH2_gas=241
    integer,parameter::idx_CH2NH2j_gas=242
    integer,parameter::idx_CH3C3N_gas=243
    integer,parameter::idx_CH3C4H_gas=244
    integer,parameter::idx_CH3C5N_gas=245
    integer,parameter::idx_CH3C6H_gas=246
    integer,parameter::idx_CH3C7N_gas=247
    integer,parameter::idx_CH3CCH_gas=248
    integer,parameter::idx_CH3CH2OH_gas=249
    integer,parameter::idx_CH3CHCH2_gas=250
    integer,parameter::idx_CH3CHO_gas=251
    integer,parameter::idx_CH3CHOHj_gas=252
    integer,parameter::idx_CH3CN_gas=253
    integer,parameter::idx_CH3CNj_gas=254
    integer,parameter::idx_CH3CNHj_gas=255
    integer,parameter::idx_CH3COj_gas=256
    integer,parameter::idx_CH3COCH3_gas=257
    integer,parameter::idx_CH3NH2_gas=258
    integer,parameter::idx_CH3NH2j_gas=259
    integer,parameter::idx_CH3NH3j_gas=260
    integer,parameter::idx_CH3O2j_gas=261
    integer,parameter::idx_CH3OCH2_gas=262
    integer,parameter::idx_CH3OCH3_gas=263
    integer,parameter::idx_CH3OCH3j_gas=264
    integer,parameter::idx_CH3OCH4j_gas=265
    integer,parameter::idx_CH3OHj_gas=266
    integer,parameter::idx_CH3OH2j_gas=267
    integer,parameter::idx_CN_gas=268
    integer,parameter::idx_CNj_gas=269
    integer,parameter::idx_CNk_gas=270
    integer,parameter::idx_CNCj_gas=271
    integer,parameter::idx_COOCH4j_gas=272
    integer,parameter::idx_CS_gas=273
    integer,parameter::idx_CSj_gas=274
    integer,parameter::idx_HS_gas=275
    integer,parameter::idx_HSj_gas=276
    integer,parameter::idx_FeH_gas=277
    integer,parameter::idx_H2C10Nj_gas=278
    integer,parameter::idx_H2C3Oj_gas=279
    integer,parameter::idx_H2C4Nj_gas=280
    integer,parameter::idx_H2C5Nj_gas=281
    integer,parameter::idx_H2C6Nj_gas=282
    integer,parameter::idx_H2C8Nj_gas=283
    integer,parameter::idx_H2CCN_gas=284
    integer,parameter::idx_H2CCO_gas=285
    integer,parameter::idx_H2CCOj_gas=286
    integer,parameter::idx_H2CN_gas=287
    integer,parameter::idx_H2CNj_gas=288
    integer,parameter::idx_H2COj_gas=289
    integer,parameter::idx_H2COHj_gas=290
    integer,parameter::idx_H2CS_gas=291
    integer,parameter::idx_H2CSj_gas=292
    integer,parameter::idx_H2NCj_gas=293
    integer,parameter::idx_H2NOj_gas=294
    integer,parameter::idx_H2S_gas=295
    integer,parameter::idx_H2Sj_gas=296
    integer,parameter::idx_HDS_gas=297
    integer,parameter::idx_HDSj_gas=298
    integer,parameter::idx_H2S2j_gas=299
    integer,parameter::idx_HDS2j_gas=300
    integer,parameter::idx_D2S2j_gas=301
    integer,parameter::idx_H3C3Oj_gas=302
    integer,parameter::idx_H3C4Nj_gas=303
    integer,parameter::idx_H3C4NHj_gas=304
    integer,parameter::idx_H3C6NHj_gas=305
    integer,parameter::idx_H3C7Nj_gas=306
    integer,parameter::idx_H3CSj_gas=307
    integer,parameter::idx_H3Sj_gas=308
    integer,parameter::idx_H3S2j_gas=309
    integer,parameter::idx_H5C2O2j_gas=310
    integer,parameter::idx_HC10Nj_gas=311
    integer,parameter::idx_HC2N_gas=312
    integer,parameter::idx_HC2Nj_gas=313
    integer,parameter::idx_HC2NCHj_gas=314
    integer,parameter::idx_HC2O_gas=315
    integer,parameter::idx_HC2Sj_gas=316
    integer,parameter::idx_HC3N_gas=317
    integer,parameter::idx_HC3Nj_gas=318
    integer,parameter::idx_HC3NHj_gas=319
    integer,parameter::idx_HC3NDj_gas=320
    integer,parameter::idx_HC3Oj_gas=321
    integer,parameter::idx_HC3Sj_gas=322
    integer,parameter::idx_HC4N_gas=323
    integer,parameter::idx_HC4Nj_gas=324
    integer,parameter::idx_HC4Oj_gas=325
    integer,parameter::idx_HC4Sj_gas=326
    integer,parameter::idx_HC5N_gas=327
    integer,parameter::idx_HC5Nj_gas=328
    integer,parameter::idx_HC5Oj_gas=329
    integer,parameter::idx_HC6N_gas=330
    integer,parameter::idx_HC6Nj_gas=331
    integer,parameter::idx_HC7N_gas=332
    integer,parameter::idx_HC7Nj_gas=333
    integer,parameter::idx_HC7Oj_gas=334
    integer,parameter::idx_HC8N_gas=335
    integer,parameter::idx_HC8Nj_gas=336
    integer,parameter::idx_HC9N_gas=337
    integer,parameter::idx_HC9Oj_gas=338
    integer,parameter::idx_HCCNC_gas=339
    integer,parameter::idx_HCN_gas=340
    integer,parameter::idx_HCNj_gas=341
    integer,parameter::idx_HCNCC_gas=342
    integer,parameter::idx_HCNHj_gas=343
    integer,parameter::idx_HCOOCH3_gas=344
    integer,parameter::idx_HCOOHj_gas=345
    integer,parameter::idx_HCS_gas=346
    integer,parameter::idx_HCSj_gas=347
    integer,parameter::idx_HN2Oj_gas=348
    integer,parameter::idx_HNC_gas=349
    integer,parameter::idx_HNCj_gas=350
    integer,parameter::idx_HNCCC_gas=351
    integer,parameter::idx_HNCO_gas=352
    integer,parameter::idx_HNCOj_gas=353
    integer,parameter::idx_HNO_gas=354
    integer,parameter::idx_HNOj_gas=355
    integer,parameter::idx_HNSj_gas=356
    integer,parameter::idx_HOCj_gas=357
    integer,parameter::idx_HOCOj_gas=358
    integer,parameter::idx_HOCSj_gas=359
    integer,parameter::idx_HS2j_gas=360
    integer,parameter::idx_HSOj_gas=361
    integer,parameter::idx_HSO2j_gas=362
    integer,parameter::idx_HSS_gas=363
    integer,parameter::idx_HSSH_gas=364
    integer,parameter::idx_N2O_gas=365
    integer,parameter::idx_NC4N_gas=366
    integer,parameter::idx_NC6N_gas=367
    integer,parameter::idx_NC8N_gas=368
    integer,parameter::idx_NCOj_gas=369
    integer,parameter::idx_NH2CH2Oj_gas=370
    integer,parameter::idx_NH2CHO_gas=371
    integer,parameter::idx_NH2CN_gas=372
    integer,parameter::idx_NH2CNDj_gas=373
    integer,parameter::idx_NH2CNHj_gas=374
    integer,parameter::idx_NO_gas=375
    integer,parameter::idx_NOj_gas=376
    integer,parameter::idx_NO2_gas=377
    integer,parameter::idx_NO2j_gas=378
    integer,parameter::idx_NS_gas=379
    integer,parameter::idx_NSj_gas=380
    integer,parameter::idx_OCN_gas=381
    integer,parameter::idx_OCS_gas=382
    integer,parameter::idx_OCSj_gas=383
    integer,parameter::idx_S2_gas=384
    integer,parameter::idx_S2j_gas=385
    integer,parameter::idx_SO_gas=386
    integer,parameter::idx_SOj_gas=387
    integer,parameter::idx_SO2_gas=388
    integer,parameter::idx_SO2j_gas=389
    integer,parameter::idx_CH3OCHO_gas=390
    integer,parameter::idx_Si_gas=391
    integer,parameter::idx_Sij_gas=392
    integer,parameter::idx_SiO_gas=393
    integer,parameter::idx_SiOj_gas=394
    integer,parameter::idx_SiO2_gas=395
    integer,parameter::idx_SiS_gas=396
    integer,parameter::idx_SiSj_gas=397
    integer,parameter::idx_SiN_gas=398
    integer,parameter::idx_SiNj_gas=399
    integer,parameter::idx_SiC_gas=400
    integer,parameter::idx_SiCj_gas=401
    integer,parameter::idx_SiC2_gas=402
    integer,parameter::idx_SiC2j_gas=403
    integer,parameter::idx_SiH_gas=404
    integer,parameter::idx_SiH2_gas=405
    integer,parameter::idx_SiH2j_gas=406
    integer,parameter::idx_SiH3_gas=407
    integer,parameter::idx_SiH4_gas=408
    integer,parameter::idx_P_gas=409
    integer,parameter::idx_Pj_gas=410
    integer,parameter::idx_PH_gas=411
    integer,parameter::idx_PHj_gas=412
    integer,parameter::idx_PO_gas=413
    integer,parameter::idx_POj_gas=414
    integer,parameter::idx_PH2_gas=415
    integer,parameter::idx_PH2j_gas=416
    integer,parameter::idx_C3P_gas=417
    integer,parameter::idx_C4P_gas=418
    integer,parameter::idx_C4Pj_gas=419
    integer,parameter::idx_CH2PH_gas=420
    integer,parameter::idx_CH2Sij_gas=421
    integer,parameter::idx_CHSij_gas=422
    integer,parameter::idx_CP_gas=423
    integer,parameter::idx_CPj_gas=424
    integer,parameter::idx_H2CSiCH_gas=425
    integer,parameter::idx_H2Fj_gas=426
    integer,parameter::idx_H2POj_gas=427
    integer,parameter::idx_H2SiO_gas=428
    integer,parameter::idx_H2SiOj_gas=429
    integer,parameter::idx_H3SiOj_gas=430
    integer,parameter::idx_HCCP_gas=431
    integer,parameter::idx_HCCSi_gas=432
    integer,parameter::idx_HCP_gas=433
    integer,parameter::idx_HCPj_gas=434
    integer,parameter::idx_HCSi_gas=435
    integer,parameter::idx_HF_gas=436
    integer,parameter::idx_HFj_gas=437
    integer,parameter::idx_HNSi_gas=438
    integer,parameter::idx_HNSij_gas=439
    integer,parameter::idx_HPNj_gas=440
    integer,parameter::idx_HPO_gas=441
    integer,parameter::idx_HPOj_gas=442
    integer,parameter::idx_HSiNHj_gas=443
    integer,parameter::idx_HSiOj_gas=444
    integer,parameter::idx_HSiO2j_gas=445
    integer,parameter::idx_HSiSj_gas=446
    integer,parameter::idx_HeHj_gas=447
    integer,parameter::idx_PC2Hj_gas=448
    integer,parameter::idx_PC2H2j_gas=449
    integer,parameter::idx_PC2H3j_gas=450
    integer,parameter::idx_PC2H4j_gas=451
    integer,parameter::idx_PC3Hj_gas=452
    integer,parameter::idx_PC4Hj_gas=453
    integer,parameter::idx_PC4H2j_gas=454
    integer,parameter::idx_PCH2j_gas=455
    integer,parameter::idx_PCH3j_gas=456
    integer,parameter::idx_PCH4j_gas=457
    integer,parameter::idx_PH3j_gas=458
    integer,parameter::idx_PN_gas=459
    integer,parameter::idx_PNj_gas=460
    integer,parameter::idx_PNH2j_gas=461
    integer,parameter::idx_PNH3j_gas=462
    integer,parameter::idx_SiC2CH3_gas=463
    integer,parameter::idx_SiC2Hj_gas=464
    integer,parameter::idx_SiC2H2j_gas=465
    integer,parameter::idx_SiC2H3j_gas=466
    integer,parameter::idx_SiC3D_gas=467
    integer,parameter::idx_SiC3H_gas=468
    integer,parameter::idx_SiC3Hj_gas=469
    integer,parameter::idx_SiC3H2j_gas=470
    integer,parameter::idx_SiC3H5_gas=471
    integer,parameter::idx_SiC4_gas=472
    integer,parameter::idx_SiC4j_gas=473
    integer,parameter::idx_SiC4H_gas=474
    integer,parameter::idx_SiC4Hj_gas=475
    integer,parameter::idx_SiC6H_gas=476
    integer,parameter::idx_SiC8H_gas=477
    integer,parameter::idx_SiCH2_gas=478
    integer,parameter::idx_SiCH3_gas=479
    integer,parameter::idx_SiCH3j_gas=480
    integer,parameter::idx_SiCH4j_gas=481
    integer,parameter::idx_SiFj_gas=482
    integer,parameter::idx_SiHj_gas=483
    integer,parameter::idx_SiH3j_gas=484
    integer,parameter::idx_SiH4j_gas=485
    integer,parameter::idx_SiH5j_gas=486
    integer,parameter::idx_SiNC_gas=487
    integer,parameter::idx_SiNCj_gas=488
    integer,parameter::idx_SiNCHj_gas=489
    integer,parameter::idx_c_HCCHSi_gas=490
    integer,parameter::idx_c_SiC2_gas=491
    integer,parameter::idx_l_SiC3_gas=492
    integer,parameter::idx_l_SiC3j_gas=493
    integer,parameter::idx_l_C3H_gas=494
    integer,parameter::idx_c_C3H_gas=495
    integer,parameter::idx_l_C3H2_gas=496
    integer,parameter::idx_c_C3H2_gas=497
    integer,parameter::idx_l_C3H2j_gas=498
    integer,parameter::idx_c_C3H2j_gas=499
    integer,parameter::idx_l_C3H3j_gas=500
    integer,parameter::idx_c_C3H3j_gas=501
    integer,parameter::idx_Mg_gas=502
    integer,parameter::idx_Mgj_gas=503
    integer,parameter::idx_MgH_gas=504
    integer,parameter::idx_MgH2_gas=505
    integer,parameter::idx_HCl_gas=506
    integer,parameter::idx_HClj_gas=507
    integer,parameter::idx_H2Cl_gas=508
    integer,parameter::idx_H2Clj_gas=509
    integer,parameter::idx_CCl_gas=510
    integer,parameter::idx_CClj_gas=511
    integer,parameter::idx_ClO_gas=512
    integer,parameter::idx_ClOj_gas=513
    integer,parameter::idx_H2CClj_gas=514
    integer,parameter::idx_Na_gas=515
    integer,parameter::idx_Naj_gas=516
    integer,parameter::idx_NaH_gas=517
    integer,parameter::idx_NaOH_gas=518
    integer,parameter::idx_NaH2j_gas=519
    integer,parameter::idx_NaH2Oj_gas=520
    integer,parameter::idx_H_surface=521
    integer,parameter::idx_D_surface=522
    integer,parameter::idx_H2_surface=523
    integer,parameter::idx_D2_surface=524
    integer,parameter::idx_HD_surface=525
    integer,parameter::idx_He_surface=526
    integer,parameter::idx_O_surface=527
    integer,parameter::idx_O2_surface=528
    integer,parameter::idx_O3_surface=529
    integer,parameter::idx_OH_surface=530
    integer,parameter::idx_OD_surface=531
    integer,parameter::idx_H2O_surface=532
    integer,parameter::idx_HDO_surface=533
    integer,parameter::idx_D2O_surface=534
    integer,parameter::idx_O2H_surface=535
    integer,parameter::idx_O2D_surface=536
    integer,parameter::idx_HOOH_surface=537
    integer,parameter::idx_HOOD_surface=538
    integer,parameter::idx_DOOD_surface=539
    integer,parameter::idx_Fe_surface=540
    integer,parameter::idx_FeH_surface=541
    integer,parameter::idx_N_surface=542
    integer,parameter::idx_S_surface=543
    integer,parameter::idx_C_surface=544
    integer,parameter::idx_C2_surface=545
    integer,parameter::idx_CO_surface=546
    integer,parameter::idx_HCO_surface=547
    integer,parameter::idx_DCO_surface=548
    integer,parameter::idx_H2CO_surface=549
    integer,parameter::idx_HDCO_surface=550
    integer,parameter::idx_D2CO_surface=551
    integer,parameter::idx_CH2OH_surface=552
    integer,parameter::idx_CD2OD_surface=553
    integer,parameter::idx_CH2OD_surface=554
    integer,parameter::idx_CHDOH_surface=555
    integer,parameter::idx_CHDOD_surface=556
    integer,parameter::idx_CD2OH_surface=557
    integer,parameter::idx_CH3O_surface=558
    integer,parameter::idx_CHD2O_surface=559
    integer,parameter::idx_CH2DO_surface=560
    integer,parameter::idx_CD3O_surface=561
    integer,parameter::idx_CH3OH_surface=562
    integer,parameter::idx_CH3OD_surface=563
    integer,parameter::idx_CHD2OH_surface=564
    integer,parameter::idx_CHD2OD_surface=565
    integer,parameter::idx_CH2DOH_surface=566
    integer,parameter::idx_CH2DOD_surface=567
    integer,parameter::idx_CD3OD_surface=568
    integer,parameter::idx_CD3OH_surface=569
    integer,parameter::idx_CH_surface=570
    integer,parameter::idx_CD_surface=571
    integer,parameter::idx_CH2_surface=572
    integer,parameter::idx_CHD_surface=573
    integer,parameter::idx_CD2_surface=574
    integer,parameter::idx_CH3_surface=575
    integer,parameter::idx_CH2D_surface=576
    integer,parameter::idx_CHD2_surface=577
    integer,parameter::idx_CD3_surface=578
    integer,parameter::idx_CH4_surface=579
    integer,parameter::idx_CH3D_surface=580
    integer,parameter::idx_CH2D2_surface=581
    integer,parameter::idx_CHD3_surface=582
    integer,parameter::idx_CD4_surface=583
    integer,parameter::idx_CO2_surface=584
    integer,parameter::idx_HCOOH_surface=585
    integer,parameter::idx_HCOOD_surface=586
    integer,parameter::idx_DCOOH_surface=587
    integer,parameter::idx_DCOOD_surface=588
    integer,parameter::idx_HOCO_surface=589
    integer,parameter::idx_DOCO_surface=590
    integer,parameter::idx_NH_surface=591
    integer,parameter::idx_ND_surface=592
    integer,parameter::idx_NH2_surface=593
    integer,parameter::idx_NHD_surface=594
    integer,parameter::idx_ND2_surface=595
    integer,parameter::idx_NH3_surface=596
    integer,parameter::idx_NH2D_surface=597
    integer,parameter::idx_NHD2_surface=598
    integer,parameter::idx_ND3_surface=599
    integer,parameter::idx_C10_surface=600
    integer,parameter::idx_C10H_surface=601
    integer,parameter::idx_C10H2_surface=602
    integer,parameter::idx_C10N_surface=603
    integer,parameter::idx_C11_surface=604
    integer,parameter::idx_C2H2_surface=605
    integer,parameter::idx_C2HD_surface=606
    integer,parameter::idx_C2D2_surface=607
    integer,parameter::idx_C2H3_surface=608
    integer,parameter::idx_C2H2D_surface=609
    integer,parameter::idx_C2HD2_surface=610
    integer,parameter::idx_C2D3_surface=611
    integer,parameter::idx_C2H4_surface=612
    integer,parameter::idx_C2H3D_surface=613
    integer,parameter::idx_C2H2D2_surface=614
    integer,parameter::idx_C2HD3_surface=615
    integer,parameter::idx_C2D4_surface=616
    integer,parameter::idx_C2H5_surface=617
    integer,parameter::idx_C2H6_surface=618
    integer,parameter::idx_C3_surface=619
    integer,parameter::idx_C3N_surface=620
    integer,parameter::idx_C3O_surface=621
    integer,parameter::idx_C3S_surface=622
    integer,parameter::idx_C4_surface=623
    integer,parameter::idx_C4H_surface=624
    integer,parameter::idx_C4D_surface=625
    integer,parameter::idx_C4H2_surface=626
    integer,parameter::idx_C4HD_surface=627
    integer,parameter::idx_C4D2_surface=628
    integer,parameter::idx_C4H3_surface=629
    integer,parameter::idx_C4N_surface=630
    integer,parameter::idx_C4S_surface=631
    integer,parameter::idx_C5_surface=632
    integer,parameter::idx_C5H_surface=633
    integer,parameter::idx_C5D_surface=634
    integer,parameter::idx_C5H2_surface=635
    integer,parameter::idx_C5H3_surface=636
    integer,parameter::idx_C5H4_surface=637
    integer,parameter::idx_C5N_surface=638
    integer,parameter::idx_C5O_surface=639
    integer,parameter::idx_C6_surface=640
    integer,parameter::idx_C6H_surface=641
    integer,parameter::idx_C6H2_surface=642
    integer,parameter::idx_C6H3_surface=643
    integer,parameter::idx_C6H4_surface=644
    integer,parameter::idx_C6H6_surface=645
    integer,parameter::idx_C6N_surface=646
    integer,parameter::idx_C7_surface=647
    integer,parameter::idx_C7H_surface=648
    integer,parameter::idx_C7H2_surface=649
    integer,parameter::idx_C7H3_surface=650
    integer,parameter::idx_C7H4_surface=651
    integer,parameter::idx_C7N_surface=652
    integer,parameter::idx_C7O_surface=653
    integer,parameter::idx_C8_surface=654
    integer,parameter::idx_C8H_surface=655
    integer,parameter::idx_C8H2_surface=656
    integer,parameter::idx_C8H3_surface=657
    integer,parameter::idx_C8H4_surface=658
    integer,parameter::idx_C8N_surface=659
    integer,parameter::idx_C9_surface=660
    integer,parameter::idx_C9H_surface=661
    integer,parameter::idx_C9H2_surface=662
    integer,parameter::idx_C9H3_surface=663
    integer,parameter::idx_C9H4_surface=664
    integer,parameter::idx_C9N_surface=665
    integer,parameter::idx_C9O_surface=666
    integer,parameter::idx_CCH_surface=667
    integer,parameter::idx_CCD_surface=668
    integer,parameter::idx_CCN_surface=669
    integer,parameter::idx_CCO_surface=670
    integer,parameter::idx_CCS_surface=671
    integer,parameter::idx_CD3CN_surface=672
    integer,parameter::idx_CH2CCH_surface=673
    integer,parameter::idx_CH2CCD_surface=674
    integer,parameter::idx_CHDCCH_surface=675
    integer,parameter::idx_CHDCCD_surface=676
    integer,parameter::idx_CD2CCH_surface=677
    integer,parameter::idx_CD2CCD_surface=678
    integer,parameter::idx_CH2CHC2H_surface=679
    integer,parameter::idx_CH2CHCHCH2_surface=680
    integer,parameter::idx_CH2CHCN_surface=681
    integer,parameter::idx_CH2NH_surface=682
    integer,parameter::idx_CHDNH_surface=683
    integer,parameter::idx_CHDND_surface=684
    integer,parameter::idx_CH2ND_surface=685
    integer,parameter::idx_CD2NH_surface=686
    integer,parameter::idx_CD2ND_surface=687
    integer,parameter::idx_CH2NH2_surface=688
    integer,parameter::idx_CH2NHD_surface=689
    integer,parameter::idx_CHDNH2_surface=690
    integer,parameter::idx_CHDNHD_surface=691
    integer,parameter::idx_CHDND2_surface=692
    integer,parameter::idx_CH2ND2_surface=693
    integer,parameter::idx_CD2NH2_surface=694
    integer,parameter::idx_CD2NHD_surface=695
    integer,parameter::idx_CD2ND2_surface=696
    integer,parameter::idx_CH3C3N_surface=697
    integer,parameter::idx_CH3C4H_surface=698
    integer,parameter::idx_CH3C5N_surface=699
    integer,parameter::idx_CH3C6H_surface=700
    integer,parameter::idx_CH3C7N_surface=701
    integer,parameter::idx_CH3CCH_surface=702
    integer,parameter::idx_CH3CH2OH_surface=703
    integer,parameter::idx_CH3CHCH2_surface=704
    integer,parameter::idx_CH3CHO_surface=705
    integer,parameter::idx_CH3CN_surface=706
    integer,parameter::idx_CH2DCN_surface=707
    integer,parameter::idx_CHD2CN_surface=708
    integer,parameter::idx_CH3COCH3_surface=709
    integer,parameter::idx_CH3NH2_surface=710
    integer,parameter::idx_CH3OCH2_surface=711
    integer,parameter::idx_CH3OCH3_surface=712
    integer,parameter::idx_CN_surface=713
    integer,parameter::idx_CS_surface=714
    integer,parameter::idx_H2CCN_surface=715
    integer,parameter::idx_HDCCN_surface=716
    integer,parameter::idx_D2CCN_surface=717
    integer,parameter::idx_H2CCO_surface=718
    integer,parameter::idx_HDCCO_surface=719
    integer,parameter::idx_D2CCO_surface=720
    integer,parameter::idx_H2CN_surface=721
    integer,parameter::idx_HDCN_surface=722
    integer,parameter::idx_D2CN_surface=723
    integer,parameter::idx_H2CS_surface=724
    integer,parameter::idx_HDCS_surface=725
    integer,parameter::idx_D2CS_surface=726
    integer,parameter::idx_H2S_surface=727
    integer,parameter::idx_HDS_surface=728
    integer,parameter::idx_D2S_surface=729
    integer,parameter::idx_HC2O_surface=730
    integer,parameter::idx_DC2O_surface=731
    integer,parameter::idx_HC3N_surface=732
    integer,parameter::idx_DC3N_surface=733
    integer,parameter::idx_HC4N_surface=734
    integer,parameter::idx_DC4N_surface=735
    integer,parameter::idx_HC5N_surface=736
    integer,parameter::idx_HC6N_surface=737
    integer,parameter::idx_HC7N_surface=738
    integer,parameter::idx_HC8N_surface=739
    integer,parameter::idx_HC9N_surface=740
    integer,parameter::idx_HCCNC_surface=741
    integer,parameter::idx_DCCNC_surface=742
    integer,parameter::idx_HCN_surface=743
    integer,parameter::idx_DCN_surface=744
    integer,parameter::idx_HCNCC_surface=745
    integer,parameter::idx_DCNCC_surface=746
    integer,parameter::idx_HCOOCH3_surface=747
    integer,parameter::idx_HCS_surface=748
    integer,parameter::idx_DCS_surface=749
    integer,parameter::idx_HNC_surface=750
    integer,parameter::idx_DNC_surface=751
    integer,parameter::idx_HNCCC_surface=752
    integer,parameter::idx_DNCCC_surface=753
    integer,parameter::idx_HNCO_surface=754
    integer,parameter::idx_DNCO_surface=755
    integer,parameter::idx_HNO_surface=756
    integer,parameter::idx_DNO_surface=757
    integer,parameter::idx_HS_surface=758
    integer,parameter::idx_DS_surface=759
    integer,parameter::idx_N2_surface=760
    integer,parameter::idx_N2O_surface=761
    integer,parameter::idx_NC4N_surface=762
    integer,parameter::idx_NC6N_surface=763
    integer,parameter::idx_NC8N_surface=764
    integer,parameter::idx_NH2CHO_surface=765
    integer,parameter::idx_NH2CDO_surface=766
    integer,parameter::idx_NHDCHO_surface=767
    integer,parameter::idx_NHDCDO_surface=768
    integer,parameter::idx_ND2CHO_surface=769
    integer,parameter::idx_ND2CDO_surface=770
    integer,parameter::idx_NH2CN_surface=771
    integer,parameter::idx_NHDCN_surface=772
    integer,parameter::idx_ND2CN_surface=773
    integer,parameter::idx_HSO_surface=774
    integer,parameter::idx_DSO_surface=775
    integer,parameter::idx_HSS_surface=776
    integer,parameter::idx_DSS_surface=777
    integer,parameter::idx_HSSH_surface=778
    integer,parameter::idx_HSSD_surface=779
    integer,parameter::idx_DSSH_surface=780
    integer,parameter::idx_DSSD_surface=781
    integer,parameter::idx_NO_surface=782
    integer,parameter::idx_NO2_surface=783
    integer,parameter::idx_NS_surface=784
    integer,parameter::idx_OCN_surface=785
    integer,parameter::idx_OCS_surface=786
    integer,parameter::idx_S2_surface=787
    integer,parameter::idx_SO_surface=788
    integer,parameter::idx_SO2_surface=789
    integer,parameter::idx_CH3OCHO_surface=790
    integer,parameter::idx_Si_surface=791
    integer,parameter::idx_SiS_surface=792
    integer,parameter::idx_SiN_surface=793
    integer,parameter::idx_SiC_surface=794
    integer,parameter::idx_SiH_surface=795
    integer,parameter::idx_SiH2_surface=796
    integer,parameter::idx_SiH3_surface=797
    integer,parameter::idx_SiH4_surface=798
    integer,parameter::idx_SiC2CH3_surface=799
    integer,parameter::idx_SiC3H_surface=800
    integer,parameter::idx_SiC3H5_surface=801
    integer,parameter::idx_SiC4_surface=802
    integer,parameter::idx_SiC4H_surface=803
    integer,parameter::idx_SiC6H_surface=804
    integer,parameter::idx_SiC8H_surface=805
    integer,parameter::idx_c_HCCHSi_surface=806
    integer,parameter::idx_c_SiC2_surface=807
    integer,parameter::idx_l_C3H_surface=808
    integer,parameter::idx_l_C3D_surface=809
    integer,parameter::idx_c_C3H_surface=810
    integer,parameter::idx_c_C3D_surface=811
    integer,parameter::idx_l_SiC3_surface=812
    integer,parameter::idx_l_C3H2_surface=813
    integer,parameter::idx_l_C3HD_surface=814
    integer,parameter::idx_l_C3D2_surface=815
    integer,parameter::idx_c_C3H2_surface=816
    integer,parameter::idx_c_C3HD_surface=817
    integer,parameter::idx_c_C3D2_surface=818
    integer,parameter::idx_Mg_surface=819
    integer,parameter::idx_MgH_surface=820
    integer,parameter::idx_MgH2_surface=821
    integer,parameter::idx_Na_surface=822
    integer,parameter::idx_NaH_surface=823
    integer,parameter::idx_F_surface=824
    integer,parameter::idx_HF_surface=825
    integer,parameter::idx_DF_surface=826
    integer,parameter::idx_MgD_surface=827
    integer,parameter::idx_MgHD_surface=828
    integer,parameter::idx_MgD2_surface=829
    integer,parameter::idx_NaD_surface=830
    integer,parameter::idx_SiD_surface=831
    integer,parameter::idx_Cl_surface=832
    integer,parameter::idx_HCl_surface=833
    integer,parameter::idx_DCl_surface=834
    integer,parameter::idx_H2Cl_surface=835
    integer,parameter::idx_HDCl_surface=836
    integer,parameter::idx_D2Cl_surface=837
    integer,parameter::idx_CCl_surface=838
    integer,parameter::idx_ClO_surface=839
    integer,parameter::idx_C3H3_surface=840
    integer,parameter::idx_C3H2D_surface=841
    integer,parameter::idx_C3HD2_surface=842
    integer,parameter::idx_C3D3_surface=843
    integer,parameter::idx_C3H4_surface=844
    integer,parameter::idx_C3H3D_surface=845
    integer,parameter::idx_C3H2D2_surface=846
    integer,parameter::idx_HCCN_surface=847
    integer,parameter::idx_DCCN_surface=848
    integer,parameter::idx_C2H2N_surface=849
    integer,parameter::idx_C2H5OH_surface=850
    integer,parameter::idx_C2H3D3_surface=851
    integer,parameter::idx_C2H5D_surface=852
    integer,parameter::idx_C2H4D2_surface=853
    integer,parameter::idx_C2H3N_surface=854
    integer,parameter::idx_C2H4O_surface=855
    integer,parameter::idx_CH2DCHO_surface=856
    integer,parameter::idx_CH3CDO_surface=857
    integer,parameter::idx_CD3CHO_surface=858
    integer,parameter::idx_CHD2CDO_surface=859
    integer,parameter::idx_CHD2CHO_surface=860
    integer,parameter::idx_CH2DCDO_surface=861
    integer,parameter::idx_C2H4D_surface=862
    integer,parameter::idx_C2H2D3_surface=863
    integer,parameter::idx_C2H3D2_surface=864
    integer,parameter::idx_C3H3N_surface=865
    integer,parameter::idx_H4C3N_surface=866
    integer,parameter::idx_C3D3N_surface=867
    integer,parameter::idx_HD3C3N_surface=868
    integer,parameter::idx_C3H2DN_surface=869
    integer,parameter::idx_H3DC3N_surface=870
    integer,parameter::idx_C3HD2N_surface=871
    integer,parameter::idx_H2D2C3N_surface=872
    integer,parameter::idx_HC3O_surface=873
    integer,parameter::idx_DC3O_surface=874
    integer,parameter::idx_C4H4_surface=875
    integer,parameter::idx_CH5N_surface=876
    integer,parameter::idx_CH3NH_surface=877
    integer,parameter::idx_FeD_surface=878
    integer,parameter::idx_H2C3N_surface=879
    integer,parameter::idx_H2C5N_surface=880
    integer,parameter::idx_H3C5N_surface=881
    integer,parameter::idx_H2C7N_surface=882
    integer,parameter::idx_H3C7N_surface=883
    integer,parameter::idx_H2C9N_surface=884
    integer,parameter::idx_H3C9N_surface=885
    integer,parameter::idx_H5C3N_surface=886
    integer,parameter::idx_C2H2O_surface=887
    integer,parameter::idx_C2HDO_surface=888
    integer,parameter::idx_C2D2O_surface=889
    integer,parameter::idx_HDC3N_surface=890
    integer,parameter::idx_D2C3N_surface=891
    integer,parameter::idx_H2C3O_surface=892
    integer,parameter::idx_HDC3O_surface=893
    integer,parameter::idx_D2C3O_surface=894
    integer,parameter::idx_C2HDN_surface=895
    integer,parameter::idx_C2D2N_surface=896
    integer,parameter::idx_CH3N_surface=897
    integer,parameter::idx_N2H2_surface=898
    integer,parameter::idx_N2HD_surface=899
    integer,parameter::idx_N2D2_surface=900
    integer,parameter::idx_NH2OH_surface=901
    integer,parameter::idx_NH2OD_surface=902
    integer,parameter::idx_N2H_surface=903
    integer,parameter::idx_N2D_surface=904
    integer,parameter::idx_CNH2_surface=905
    integer,parameter::idx_CHNH2_surface=906
    integer,parameter::idx_HON_surface=907
    integer,parameter::idx_DON_surface=908
    integer,parameter::idx_NHNO_surface=909
    integer,parameter::idx_CH2DN_surface=910
    integer,parameter::idx_CHD2N_surface=911
    integer,parameter::idx_CD3N_surface=912
    integer,parameter::idx_NH2NO_surface=913
    integer,parameter::idx_CH3CO_surface=914
    integer,parameter::idx_CH2DOCH3_surface=915
    integer,parameter::idx_CHD2OCH3_surface=916
    integer,parameter::idx_CH2DOCH2D_surface=917
    integer,parameter::idx_CD3OCH3_surface=918
    integer,parameter::idx_CHD2OCH2D_surface=919
    integer,parameter::idx_DCOOCH3_surface=920
    integer,parameter::idx_HCOOCH2D_surface=921
    integer,parameter::idx_DCOOCH2D_surface=922
    integer,parameter::idx_HCOOCHD2_surface=923
    integer,parameter::idx_DCOOCHD2_surface=924
    integer,parameter::idx_HCOOCD3_surface=925
    integer,parameter::idx_DCOOCD3_surface=926
    integer,parameter::idx_NH2CO_surface=927
    integer,parameter::idx_NHDCO_surface=928
    integer,parameter::idx_ND2CO_surface=929
    integer,parameter::idx_CH2DOCHO_surface=930
    integer,parameter::idx_CH3OCDO_surface=931
    integer,parameter::idx_CHD2OCHO_surface=932
    integer,parameter::idx_CH2DOCDO_surface=933
    integer,parameter::idx_CD3OCHO_surface=934
    integer,parameter::idx_CHD2OCDO_surface=935
    integer,parameter::idx_CHDOCDO_surface=936
    integer,parameter::idx_CHDOCHDO_surface=937
    integer,parameter::idx_CH3OCHD_surface=938
    integer,parameter::idx_CH2DOCH2_surface=939
    integer,parameter::idx_CH2DOCHD_surface=940
    integer,parameter::idx_CHD2OCH2_surface=941
    integer,parameter::idx_CH3OCD2_surface=942
    integer,parameter::idx_CH2DOCD2_surface=943
    integer,parameter::idx_CH3OCHD2_surface=944
    integer,parameter::idx_CH3OCH2D_surface=945
    integer,parameter::idx_CH2DCO_surface=946
    integer,parameter::idx_CH3NHD_surface=947
    integer,parameter::idx_CH2DNH2_surface=948
    integer,parameter::idx_P_surface=949
    integer,parameter::idx_PO_surface=950
    integer,parameter::idx_PH_surface=951
    integer,parameter::idx_PD_surface=952
    integer,parameter::idx_PH2_surface=953
    integer,parameter::idx_PHD_surface=954
    integer,parameter::idx_PD2_surface=955
    integer,parameter::idx_PN_surface=956
    integer,parameter::idx_CP_surface=957
    integer,parameter::idx_CCP_surface=958
    integer,parameter::idx_C3P_surface=959
    integer,parameter::idx_C4P_surface=960
    integer,parameter::idx_CH2PH_surface=961
    integer,parameter::idx_CHD2PH_surface=962
    integer,parameter::idx_CH2PD_surface=963
    integer,parameter::idx_CHDPD_surface=964
    integer,parameter::idx_CD2P2_surface=965
    integer,parameter::idx_HCP_surface=966
    integer,parameter::idx_DCP_surface=967
    integer,parameter::idx_HCCP_surface=968
    integer,parameter::idx_DCCP_surface=969
    integer,parameter::idx_HPO_surface=970
    integer,parameter::idx_DPO_surface=971
    integer,parameter::idx_surface_mask=972
    integer,parameter::idx_mantle_mask=973
    integer,parameter::idx_dummy=974

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
