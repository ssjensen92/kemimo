module kemimo_dust_rates
    real*8:: kph_factor, k_H2O_ph_0
    real*8:: Nsites, Ffuva
    real*8 :: Ffuva_CO, Ffuva_H2, Ffuva_HD, Ffuva_N2
contains

  !*******************
  !compute rates and store into commons kall(:)
  subroutine computeDustRates(n, ngas, variable_Tgas, variable_crflux, variable_Av, &
       Tdust, a_grain, Gnot)
    use kemimo_commons
    use kemimo_sticking
    implicit none
    real*8,intent(in):: n(nmols)
    real*8,intent(in):: ngas, variable_Tgas
    real*8,intent(in):: variable_crflux, variable_Av
    real*8,intent(in):: a_grain, Gnot, Tdust
    real*8:: prefreezeout, prefreezeout_H, indns, invTd
    real*8:: ktmp(nrea)
    real*8:: p3, p4, pexp
    real*8:: n_ice_tot
    real*8, parameter :: beta = 2.5
    real*8, parameter :: T0 = 87.0 !87.0 np-ASW ice  - 56 silicate
    real*8, parameter :: S0 = 0.76 !0.76 np-ASW ice  - 0.95 silicate
    real*8, parameter :: gamma_H2O = 2.2d0
    real*8, parameter :: gamma_CO2 = 2.03d0
    real*8, parameter :: P_H2O_isrf = 5.4d-3
    real*8, parameter :: P_H2O_cr = 4.7d-3
    integer::i


    pexp = -3.5 !size distribution exponent
    p3 = pexp + 3d0
    p4 = pexp + 4d0
    
    invTd = 1d0/variable_Tgas
    ! -------------------------------------------------
    ! Update swapping rates:
    !call computeSwapRates(invTd)

    ! -------------------------------------------------
    ! Update sticking rates:
    call computeSticking(variable_Tgas)

    ! -------------------------------------------------    
    !freezeout prefactor
    prefreezeout =  sqrt(variable_Tgas) * a_grain**2 * pi * xdust * ngas

    !number of binding sites per unit of volume
    Nsites = 4.0 * pi * a_grain**2 * site_density ! Number of sites pr dust grain
    ndns = xdust * ngas * Nsites

    !inverse of ndns
    indns = 1d0/ndns

    ! -------------------------------------------------
    !FUV photons per bound molecule per unit of time, 1/s
    kph_factor = a_grain**2.0 * pi * xdust * ngas
    ! Adjust/scale F_cr (CR-UV flux) based on the current CR flux
    F_cr = variable_crflux * F_cr_not / 1.3d-17
    ! Av attenuation from Hollenbach+2009
    Ffuva = kph_factor * max(F_cr, max(Gnot, 1d0)*Fnot*exp(-1.8d0*variable_Av))

    ! Total gas-phase unattenuated dissociation rate from KIDA
    k_H2O_ph_0 = 8.01d-10 
    !k_H2O_ph_0 = 8.310d-10 

    ! Specials to account for self-shielding
    Ffuva_CO = ss_CO * Ffuva
    Ffuva_H2 = ss_H2 * Ffuva
    Ffuva_HD = ss_HD * Ffuva
    Ffuva_N2 = ss_N2 * Ffuva

    !!BEGIN_RATES
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-22 21:21:48
    ! CHANGESET: xxxxxxx
    ! BY: unknown@unknown

    !H_gas -> H_0001 (freezeout)
    kall(1) = 1.44981548d+04*prefreezeout * kstick(958 - surface_start)

    !D_gas -> D_0001 (freezeout)
    kall(2) = 1.02517436d+04*prefreezeout * kstick(959 - surface_start)

    !p_H2_gas -> p_H2_0001 (freezeout)
    kall(3) = 1.02517436d+04*prefreezeout * kstick(960 - surface_start)

    !o_H2_gas -> o_H2_0001 (freezeout)
    kall(4) = 1.02517436d+04*prefreezeout * kstick(962 - surface_start)

    !p_D2_gas -> p_D2_0001 (freezeout)
    kall(5) = 7.24907742d+03*prefreezeout * kstick(961 - surface_start)

    !o_D2_gas -> o_D2_0001 (freezeout)
    kall(6) = 7.24907742d+03*prefreezeout * kstick(963 - surface_start)

    !HD_gas -> HD_0001 (freezeout)
    kall(7) = 8.37051360d+03*prefreezeout * kstick(964 - surface_start)

    !He_gas -> He_0001 (freezeout)
    kall(8) = 7.24907742d+03*prefreezeout * kstick(965 - surface_start)

    !O_gas -> O_0001 (freezeout)
    kall(9) = 3.62453871d+03*prefreezeout * kstick(966 - surface_start)

    !O2_gas -> O2_0001 (freezeout)
    kall(10) = 2.56293590d+03*prefreezeout * kstick(967 - surface_start)

    !O3_gas -> O3_0001 (freezeout)
    kall(11) = 2.09262840d+03*prefreezeout * kstick(968 - surface_start)

    !OH_gas -> OH_0001 (freezeout)
    kall(12) = 3.51631905d+03*prefreezeout * kstick(969 - surface_start)

    !OD_gas -> OD_0001 (freezeout)
    kall(13) = 3.41724787d+03*prefreezeout * kstick(970 - surface_start)

    !H2O_gas -> H2O_0001 (freezeout)
    kall(14) = 3.41724787d+03*prefreezeout * kstick(971 - surface_start)

    !HDO_gas -> HDO_0001 (freezeout)
    kall(15) = 3.32610483d+03*prefreezeout * kstick(972 - surface_start)

    !D2O_gas -> D2O_0001 (freezeout)
    kall(16) = 3.24188598d+03*prefreezeout * kstick(973 - surface_start)

    !O2H_gas -> O2H_0001 (freezeout)
    kall(17) = 2.52380481d+03*prefreezeout * kstick(974 - surface_start)

    !O2D_gas -> O2D_0001 (freezeout)
    kall(18) = 2.48641304d+03*prefreezeout * kstick(975 - surface_start)

    !HOOH_gas -> HOOH_0001 (freezeout)
    kall(19) = 2.48641304d+03*prefreezeout * kstick(976 - surface_start)

    !HOOD_gas -> HOOD_0001 (freezeout)
    kall(20) = 2.45063545d+03*prefreezeout * kstick(977 - surface_start)

    !DOOD_gas -> DOOD_0001 (freezeout)
    kall(21) = 2.41635914d+03*prefreezeout * kstick(978 - surface_start)

    !F_gas -> F_0001 (freezeout)
    kall(22) = 3.32610483d+03*prefreezeout * kstick(1263 - surface_start)

    !Cl_gas -> Cl_0001 (freezeout)
    kall(23) = 2.45063545d+03*prefreezeout * kstick(1271 - surface_start)

    !C_gas -> C_0001 (freezeout)
    kall(24) = 4.18525680d+03*prefreezeout * kstick(983 - surface_start)

    !C2_gas -> C2_0001 (freezeout)
    kall(25) = 2.95942346d+03*prefreezeout * kstick(984 - surface_start)

    !CO_gas -> CO_0001 (freezeout)
    kall(26) = 2.73989373d+03*prefreezeout * kstick(985 - surface_start)

    !CO2_gas -> CO2_0001 (freezeout)
    kall(27) = 2.18567908d+03*prefreezeout * kstick(1023 - surface_start)

    !N_gas -> N_0001 (freezeout)
    kall(28) = 3.87479487d+03*prefreezeout * kstick(981 - surface_start)

    !N2_gas -> N2_0001 (freezeout)
    kall(29) = 2.73989373d+03*prefreezeout * kstick(1199 - surface_start)

    !HCO_gas -> HCO_0001 (freezeout)
    kall(30) = 2.69223977d+03*prefreezeout * kstick(986 - surface_start)

    !DCO_gas -> DCO_0001 (freezeout)
    kall(31) = 2.64698882d+03*prefreezeout * kstick(987 - surface_start)

    !H2CO_gas -> H2CO_0001 (freezeout)
    kall(32) = 2.64698882d+03*prefreezeout * kstick(988 - surface_start)

    !HDCO_gas -> HDCO_0001 (freezeout)
    kall(33) = 2.60394548d+03*prefreezeout * kstick(989 - surface_start)

    !D2CO_gas -> D2CO_0001 (freezeout)
    kall(34) = 2.56293590d+03*prefreezeout * kstick(990 - surface_start)

    !CH2OH_gas -> CH2OH_0001 (freezeout)
    kall(35) = 2.60394548d+03*prefreezeout * kstick(991 - surface_start)

    !CD2OD_gas -> CD2OD_0001 (freezeout)
    kall(36) = 2.48641304d+03*prefreezeout * kstick(992 - surface_start)

    !CH2OD_gas -> CH2OD_0001 (freezeout)
    kall(37) = 2.56293590d+03*prefreezeout * kstick(993 - surface_start)

    !CHDOH_gas -> CHDOH_0001 (freezeout)
    kall(38) = 2.56293590d+03*prefreezeout * kstick(994 - surface_start)

    !CHDOD_gas -> CHDOD_0001 (freezeout)
    kall(39) = 2.52380481d+03*prefreezeout * kstick(995 - surface_start)

    !CD2OH_gas -> CD2OH_0001 (freezeout)
    kall(40) = 2.52380481d+03*prefreezeout * kstick(996 - surface_start)

    !CH3O_gas -> CH3O_0001 (freezeout)
    kall(41) = 2.60394548d+03*prefreezeout * kstick(997 - surface_start)

    !CHD2O_gas -> CHD2O_0001 (freezeout)
    kall(42) = 2.52380481d+03*prefreezeout * kstick(998 - surface_start)

    !CH2DO_gas -> CH2DO_0001 (freezeout)
    kall(43) = 2.56293590d+03*prefreezeout * kstick(999 - surface_start)

    !CD3O_gas -> CD3O_0001 (freezeout)
    kall(44) = 2.48641304d+03*prefreezeout * kstick(1000 - surface_start)

    !CH3OH_gas -> CH3OH_0001 (freezeout)
    kall(45) = 2.56293590d+03*prefreezeout * kstick(1001 - surface_start)

    !CH3OD_gas -> CH3OD_0001 (freezeout)
    kall(46) = 2.52380481d+03*prefreezeout * kstick(1002 - surface_start)

    !CHD2OH_gas -> CHD2OH_0001 (freezeout)
    kall(47) = 2.48641304d+03*prefreezeout * kstick(1003 - surface_start)

    !CHD2OD_gas -> CHD2OD_0001 (freezeout)
    kall(48) = 2.45063545d+03*prefreezeout * kstick(1004 - surface_start)

    !CH2DOH_gas -> CH2DOH_0001 (freezeout)
    kall(49) = 2.52380481d+03*prefreezeout * kstick(1005 - surface_start)

    !CH2DOD_gas -> CH2DOD_0001 (freezeout)
    kall(50) = 2.48641304d+03*prefreezeout * kstick(1006 - surface_start)

    !CD3OD_gas -> CD3OD_0001 (freezeout)
    kall(51) = 2.41635914d+03*prefreezeout * kstick(1007 - surface_start)

    !CD3OH_gas -> CD3OH_0001 (freezeout)
    kall(52) = 2.45063545d+03*prefreezeout * kstick(1008 - surface_start)

    !CH_gas -> CH_0001 (freezeout)
    kall(53) = 4.02106467d+03*prefreezeout * kstick(1009 - surface_start)

    !CD_gas -> CD_0001 (freezeout)
    kall(54) = 3.87479487d+03*prefreezeout * kstick(1010 - surface_start)

    !CH2_gas -> CH2_0001 (freezeout)
    kall(55) = 3.87479487d+03*prefreezeout * kstick(1011 - surface_start)

    !CHD_gas -> CHD_0001 (freezeout)
    kall(56) = 3.74340748d+03*prefreezeout * kstick(1012 - surface_start)

    !CD2_gas -> CD2_0001 (freezeout)
    kall(57) = 3.62453871d+03*prefreezeout * kstick(1013 - surface_start)

    !CH3_gas -> CH3_0001 (freezeout)
    kall(58) = 3.74340748d+03*prefreezeout * kstick(1014 - surface_start)

    !CH2D_gas -> CH2D_0001 (freezeout)
    kall(59) = 3.62453871d+03*prefreezeout * kstick(1015 - surface_start)

    !CHD2_gas -> CHD2_0001 (freezeout)
    kall(60) = 3.51631905d+03*prefreezeout * kstick(1016 - surface_start)

    !CD3_gas -> CD3_0001 (freezeout)
    kall(61) = 3.41724787d+03*prefreezeout * kstick(1017 - surface_start)

    !CH4_gas -> CH4_0001 (freezeout)
    kall(62) = 3.62453871d+03*prefreezeout * kstick(1018 - surface_start)

    !CH3D_gas -> CH3D_0001 (freezeout)
    kall(63) = 3.51631905d+03*prefreezeout * kstick(1019 - surface_start)

    !CH2D2_gas -> CH2D2_0001 (freezeout)
    kall(64) = 3.41724787d+03*prefreezeout * kstick(1020 - surface_start)

    !CHD3_gas -> CHD3_0001 (freezeout)
    kall(65) = 3.32610483d+03*prefreezeout * kstick(1021 - surface_start)

    !CD4_gas -> CD4_0001 (freezeout)
    kall(66) = 3.24188598d+03*prefreezeout * kstick(1022 - surface_start)

    !HCOOH_gas -> HCOOH_0001 (freezeout)
    kall(67) = 2.13763631d+03*prefreezeout * kstick(1024 - surface_start)

    !HCOOD_gas -> HCOOD_0001 (freezeout)
    kall(68) = 2.11477323d+03*prefreezeout * kstick(1025 - surface_start)

    !DCOOH_gas -> DCOOH_0001 (freezeout)
    kall(69) = 2.11477323d+03*prefreezeout * kstick(1026 - surface_start)

    !DCOOD_gas -> DCOOD_0001 (freezeout)
    kall(70) = 2.09262840d+03*prefreezeout * kstick(1027 - surface_start)

    !HOCO_gas -> HOCO_0001 (freezeout)
    kall(71) = 2.16125732d+03*prefreezeout * kstick(1028 - surface_start)

    !DOCO_gas -> DOCO_0001 (freezeout)
    kall(72) = 2.13763631d+03*prefreezeout * kstick(1029 - surface_start)

    !NH_gas -> NH_0001 (freezeout)
    kall(73) = 3.74340748d+03*prefreezeout * kstick(1030 - surface_start)

    !ND_gas -> ND_0001 (freezeout)
    kall(74) = 3.62453871d+03*prefreezeout * kstick(1031 - surface_start)

    !NH2_gas -> NH2_0001 (freezeout)
    kall(75) = 3.62453871d+03*prefreezeout * kstick(1032 - surface_start)

    !NHD_gas -> NHD_0001 (freezeout)
    kall(76) = 3.51631905d+03*prefreezeout * kstick(1033 - surface_start)

    !ND2_gas -> ND2_0001 (freezeout)
    kall(77) = 3.41724787d+03*prefreezeout * kstick(1034 - surface_start)

    !NH3_gas -> NH3_0001 (freezeout)
    kall(78) = 3.51631905d+03*prefreezeout * kstick(1035 - surface_start)

    !NH2D_gas -> NH2D_0001 (freezeout)
    kall(79) = 3.41724787d+03*prefreezeout * kstick(1036 - surface_start)

    !NHD2_gas -> NHD2_0001 (freezeout)
    kall(80) = 3.32610483d+03*prefreezeout * kstick(1037 - surface_start)

    !ND3_gas -> ND3_0001 (freezeout)
    kall(81) = 3.24188598d+03*prefreezeout * kstick(1038 - surface_start)

    !Fe_gas -> Fe_0001 (freezeout)
    kall(82) = 1.93739743d+03*prefreezeout * kstick(979 - surface_start)

    !S_gas -> S_0001 (freezeout)
    kall(83) = 2.56293590d+03*prefreezeout * kstick(982 - surface_start)

    !C10_gas -> C10_0001 (freezeout)
    kall(84) = 1.32349441d+03*prefreezeout * kstick(1039 - surface_start)

    !C10H_gas -> C10H_0001 (freezeout)
    kall(85) = 1.31801408d+03*prefreezeout * kstick(1040 - surface_start)

    !C10H2_gas -> C10H2_0001 (freezeout)
    kall(86) = 1.31260126d+03*prefreezeout * kstick(1041 - surface_start)

    !C10N_gas -> C10N_0001 (freezeout)
    kall(87) = 1.25244982d+03*prefreezeout * kstick(1042 - surface_start)

    !C11_gas -> C11_0001 (freezeout)
    kall(88) = 1.26190240d+03*prefreezeout * kstick(1043 - surface_start)

    !C2H2_gas -> C2H2_0001 (freezeout)
    kall(89) = 2.84332209d+03*prefreezeout * kstick(1044 - surface_start)

    !C2HD_gas -> C2HD_0001 (freezeout)
    kall(90) = 2.79017120d+03*prefreezeout * kstick(1045 - surface_start)

    !C2D2_gas -> C2D2_0001 (freezeout)
    kall(91) = 2.73989373d+03*prefreezeout * kstick(1046 - surface_start)

    !C2H3_gas -> C2H3_0001 (freezeout)
    kall(92) = 2.79017120d+03*prefreezeout * kstick(1047 - surface_start)

    !C2H2D_gas -> C2H2D_0001 (freezeout)
    kall(93) = 2.73989373d+03*prefreezeout * kstick(1048 - surface_start)

    !C2HD2_gas -> C2HD2_0001 (freezeout)
    kall(94) = 2.69223977d+03*prefreezeout * kstick(1049 - surface_start)

    !C2D3_gas -> C2D3_0001 (freezeout)
    kall(95) = 2.64698882d+03*prefreezeout * kstick(1050 - surface_start)

    !C2H4_gas -> C2H4_0001 (freezeout)
    kall(96) = 2.73989373d+03*prefreezeout * kstick(1051 - surface_start)

    !C2H3D_gas -> C2H3D_0001 (freezeout)
    kall(97) = 2.69223977d+03*prefreezeout * kstick(1052 - surface_start)

    !C2H2D2_gas -> C2H2D2_0001 (freezeout)
    kall(98) = 2.64698882d+03*prefreezeout * kstick(1053 - surface_start)

    !C2HD3_gas -> C2HD3_0001 (freezeout)
    kall(99) = 2.60394548d+03*prefreezeout * kstick(1054 - surface_start)

    !C2D4_gas -> C2D4_0001 (freezeout)
    kall(100) = 2.56293590d+03*prefreezeout * kstick(1055 - surface_start)

    !C2H5_gas -> C2H5_0001 (freezeout)
    kall(101) = 2.69223977d+03*prefreezeout * kstick(1056 - surface_start)

    !C2H6_gas -> C2H6_0001 (freezeout)
    kall(102) = 2.64698882d+03*prefreezeout * kstick(1057 - surface_start)

    !C3_gas -> C3_0001 (freezeout)
    kall(103) = 2.41635914d+03*prefreezeout * kstick(1058 - surface_start)

    !C3H4_gas -> C3H4_0001 (freezeout)
    kall(104) = 2.29235956d+03*prefreezeout * kstick(1283 - surface_start)

    !C3N_gas -> C3N_0001 (freezeout)
    kall(105) = 2.05034872d+03*prefreezeout * kstick(1059 - surface_start)

    !C3O_gas -> C3O_0001 (freezeout)
    kall(106) = 2.01053233d+03*prefreezeout * kstick(1060 - surface_start)

    !C3S_gas -> C3S_0001 (freezeout)
    kall(107) = 1.75815952d+03*prefreezeout * kstick(1061 - surface_start)

    !C4_gas -> C4_0001 (freezeout)
    kall(108) = 2.09262840d+03*prefreezeout * kstick(1062 - surface_start)

    !C4H_gas -> C4H_0001 (freezeout)
    kall(109) = 2.07116498d+03*prefreezeout * kstick(1063 - surface_start)

    !C4D_gas -> C4D_0001 (freezeout)
    kall(110) = 2.05034872d+03*prefreezeout * kstick(1064 - surface_start)

    !C4H2_gas -> C4H2_0001 (freezeout)
    kall(111) = 2.05034872d+03*prefreezeout * kstick(1065 - surface_start)

    !C4HD_gas -> C4HD_0001 (freezeout)
    kall(112) = 2.03014775d+03*prefreezeout * kstick(1066 - surface_start)

    !C4D2_gas -> C4D2_0001 (freezeout)
    kall(113) = 2.01053233d+03*prefreezeout * kstick(1067 - surface_start)

    !C4H3_gas -> C4H3_0001 (freezeout)
    kall(114) = 2.03014775d+03*prefreezeout * kstick(1068 - surface_start)

    !C4N_gas -> C4N_0001 (freezeout)
    kall(115) = 1.84126751d+03*prefreezeout * kstick(1069 - surface_start)

    !C4S_gas -> C4S_0001 (freezeout)
    kall(116) = 1.62094299d+03*prefreezeout * kstick(1070 - surface_start)

    !C5_gas -> C5_0001 (freezeout)
    kall(117) = 1.87170374d+03*prefreezeout * kstick(1071 - surface_start)

    !C5H_gas -> C5H_0001 (freezeout)
    kall(118) = 1.85629851d+03*prefreezeout * kstick(1072 - surface_start)

    !C5D_gas -> C5D_0001 (freezeout)
    kall(119) = 1.84126751d+03*prefreezeout * kstick(1073 - surface_start)

    !C5H2_gas -> C5H2_0001 (freezeout)
    kall(120) = 1.84126751d+03*prefreezeout * kstick(1074 - surface_start)

    !C5H3_gas -> C5H3_0001 (freezeout)
    kall(121) = 1.82659582d+03*prefreezeout * kstick(1075 - surface_start)

    !C5H4_gas -> C5H4_0001 (freezeout)
    kall(122) = 1.81226935d+03*prefreezeout * kstick(1076 - surface_start)

    !C5N_gas -> C5N_0001 (freezeout)
    kall(123) = 1.68537627d+03*prefreezeout * kstick(1077 - surface_start)

    !C5O_gas -> C5O_0001 (freezeout)
    kall(124) = 1.66305242d+03*prefreezeout * kstick(1078 - surface_start)

    !C6_gas -> C6_0001 (freezeout)
    kall(125) = 1.70862393d+03*prefreezeout * kstick(1079 - surface_start)

    !C6H_gas -> C6H_0001 (freezeout)
    kall(126) = 1.69688067d+03*prefreezeout * kstick(1080 - surface_start)

    !C6H2_gas -> C6H2_0001 (freezeout)
    kall(127) = 1.68537627d+03*prefreezeout * kstick(1081 - surface_start)

    !C6H3_gas -> C6H3_0001 (freezeout)
    kall(128) = 1.67410272d+03*prefreezeout * kstick(1082 - surface_start)

    !C6H4_gas -> C6H4_0001 (freezeout)
    kall(129) = 1.66305242d+03*prefreezeout * kstick(1083 - surface_start)

    !C6H6_gas -> C6H6_0001 (freezeout)
    kall(130) = 1.64159278d+03*prefreezeout * kstick(1084 - surface_start)

    !C6N_gas -> C6N_0001 (freezeout)
    kall(131) = 1.56337624d+03*prefreezeout * kstick(1085 - surface_start)

    !C7_gas -> C7_0001 (freezeout)
    kall(132) = 1.58187838d+03*prefreezeout * kstick(1086 - surface_start)

    !C7H_gas -> C7H_0001 (freezeout)
    kall(133) = 1.57254568d+03*prefreezeout * kstick(1087 - surface_start)

    !C7H2_gas -> C7H2_0001 (freezeout)
    kall(134) = 1.56337624d+03*prefreezeout * kstick(1088 - surface_start)

    !C7H3_gas -> C7H3_0001 (freezeout)
    kall(135) = 1.55436535d+03*prefreezeout * kstick(1089 - surface_start)

    !C7H4_gas -> C7H4_0001 (freezeout)
    kall(136) = 1.54550850d+03*prefreezeout * kstick(1090 - surface_start)

    !C7N_gas -> C7N_0001 (freezeout)
    kall(137) = 1.46453480d+03*prefreezeout * kstick(1091 - surface_start)

    !C7O_gas -> C7O_0001 (freezeout)
    kall(138) = 1.44981548d+03*prefreezeout * kstick(1092 - surface_start)

    !C8_gas -> C8_0001 (freezeout)
    kall(139) = 1.47971173d+03*prefreezeout * kstick(1093 - surface_start)

    !C8H_gas -> C8H_0001 (freezeout)
    kall(140) = 1.47206459d+03*prefreezeout * kstick(1094 - surface_start)

    !C8H2_gas -> C8H2_0001 (freezeout)
    kall(141) = 1.46453480d+03*prefreezeout * kstick(1095 - surface_start)

    !C8H3_gas -> C8H3_0001 (freezeout)
    kall(142) = 1.45711939d+03*prefreezeout * kstick(1096 - surface_start)

    !C8H4_gas -> C8H4_0001 (freezeout)
    kall(143) = 1.44981548d+03*prefreezeout * kstick(1097 - surface_start)

    !C8N_gas -> C8N_0001 (freezeout)
    kall(144) = 1.38234482d+03*prefreezeout * kstick(1098 - surface_start)

    !C9_gas -> C9_0001 (freezeout)
    kall(145) = 1.39508560d+03*prefreezeout * kstick(1099 - surface_start)

    !C9H_gas -> C9H_0001 (freezeout)
    kall(146) = 1.38867138d+03*prefreezeout * kstick(1100 - surface_start)

    !C9H2_gas -> C9H2_0001 (freezeout)
    kall(147) = 1.38234482d+03*prefreezeout * kstick(1101 - surface_start)

    !C9H3_gas -> C9H3_0001 (freezeout)
    kall(148) = 1.37610396d+03*prefreezeout * kstick(1102 - surface_start)

    !C9H4_gas -> C9H4_0001 (freezeout)
    kall(149) = 1.36994686d+03*prefreezeout * kstick(1103 - surface_start)

    !C9N_gas -> C9N_0001 (freezeout)
    kall(150) = 1.31260126d+03*prefreezeout * kstick(1104 - surface_start)

    !C9O_gas -> C9O_0001 (freezeout)
    kall(151) = 1.30197274d+03*prefreezeout * kstick(1105 - surface_start)

    !CCH_gas -> CCH_0001 (freezeout)
    kall(152) = 2.89963097d+03*prefreezeout * kstick(1106 - surface_start)

    !CCD_gas -> CCD_0001 (freezeout)
    kall(153) = 2.84332209d+03*prefreezeout * kstick(1107 - surface_start)

    !CCN_gas -> CCN_0001 (freezeout)
    kall(154) = 2.35191128d+03*prefreezeout * kstick(1108 - surface_start)

    !CCO_gas -> CCO_0001 (freezeout)
    kall(155) = 2.29235956d+03*prefreezeout * kstick(1109 - surface_start)

    !CCP_gas -> CCP_0001 (freezeout)
    kall(156) = 1.95493080d+03*prefreezeout * kstick(1397 - surface_start)

    !CCS_gas -> CCS_0001 (freezeout)
    kall(157) = 1.93739743d+03*prefreezeout * kstick(1110 - surface_start)

    !CD2ND2_gas -> CD2ND2_0001 (freezeout)
    kall(158) = 2.48641304d+03*prefreezeout * kstick(1135 - surface_start)

    !CD2NH2_gas -> CD2NH2_0001 (freezeout)
    kall(159) = 2.56293590d+03*prefreezeout * kstick(1133 - surface_start)

    !CD2NHD_gas -> CD2NHD_0001 (freezeout)
    kall(160) = 2.52380481d+03*prefreezeout * kstick(1134 - surface_start)

    !CD3CN_gas -> CD3CN_0001 (freezeout)
    kall(161) = 2.18567908d+03*prefreezeout * kstick(1111 - surface_start)

    !CH2CCH_gas -> CH2CCH_0001 (freezeout)
    kall(162) = 2.32156277d+03*prefreezeout * kstick(1112 - surface_start)

    !CH2CCD_gas -> CH2CCD_0001 (freezeout)
    kall(163) = 2.29235956d+03*prefreezeout * kstick(1113 - surface_start)

    !CD2CCH_gas -> CD2CCH_0001 (freezeout)
    kall(164) = 2.26423138d+03*prefreezeout * kstick(1116 - surface_start)

    !CD2CCD_gas -> CD2CCD_0001 (freezeout)
    kall(165) = 2.23711386d+03*prefreezeout * kstick(1117 - surface_start)

    !CHDCCH_gas -> CHDCCH_0001 (freezeout)
    kall(166) = 2.29235956d+03*prefreezeout * kstick(1114 - surface_start)

    !CHDCCD_gas -> CHDCCD_0001 (freezeout)
    kall(167) = 2.26423138d+03*prefreezeout * kstick(1115 - surface_start)

    !CH2CHC2H_gas -> CH2CHC2H_0001 (freezeout)
    kall(168) = 2.01053233d+03*prefreezeout * kstick(1118 - surface_start)

    !CH2CHCHCH2_gas -> CH2CHCHCH2_0001 (freezeout)
    kall(169) = 1.97294898d+03*prefreezeout * kstick(1119 - surface_start)

    !CH2CHCN_gas -> CH2CHCN_0001 (freezeout)
    kall(170) = 1.99147472d+03*prefreezeout * kstick(1120 - surface_start)

    !CH2DCN_gas -> CH2DCN_0001 (freezeout)
    kall(171) = 2.23711386d+03*prefreezeout * kstick(1146 - surface_start)

    !CH2ND2_gas -> CH2ND2_0001 (freezeout)
    kall(172) = 2.56293590d+03*prefreezeout * kstick(1132 - surface_start)

    !CH2NH_gas -> CH2NH_0001 (freezeout)
    kall(173) = 2.69223977d+03*prefreezeout * kstick(1121 - surface_start)

    !CH2ND_gas -> CH2ND_0001 (freezeout)
    kall(174) = 2.64698882d+03*prefreezeout * kstick(1124 - surface_start)

    !CHDNH_gas -> CHDNH_0001 (freezeout)
    kall(175) = 2.64698882d+03*prefreezeout * kstick(1122 - surface_start)

    !CHDND_gas -> CHDND_0001 (freezeout)
    kall(176) = 2.60394548d+03*prefreezeout * kstick(1123 - surface_start)

    !CD2NH_gas -> CD2NH_0001 (freezeout)
    kall(177) = 2.60394548d+03*prefreezeout * kstick(1125 - surface_start)

    !CD2ND_gas -> CD2ND_0001 (freezeout)
    kall(178) = 2.56293590d+03*prefreezeout * kstick(1126 - surface_start)

    !CH2NH2_gas -> CH2NH2_0001 (freezeout)
    kall(179) = 2.64698882d+03*prefreezeout * kstick(1127 - surface_start)

    !CH2NHD_gas -> CH2NHD_0001 (freezeout)
    kall(180) = 2.60394548d+03*prefreezeout * kstick(1128 - surface_start)

    !CH3C3N_gas -> CH3C3N_0001 (freezeout)
    kall(181) = 1.79827479d+03*prefreezeout * kstick(1136 - surface_start)

    !CH3C4H_gas -> CH3C4H_0001 (freezeout)
    kall(182) = 1.81226935d+03*prefreezeout * kstick(1137 - surface_start)

    !CH3C5N_gas -> CH3C5N_0001 (freezeout)
    kall(183) = 1.53680134d+03*prefreezeout * kstick(1138 - surface_start)

    !CH3C6H_gas -> CH3C6H_0001 (freezeout)
    kall(184) = 1.54550850d+03*prefreezeout * kstick(1139 - surface_start)

    !CH3C7N_gas -> CH3C7N_0001 (freezeout)
    kall(185) = 1.36387168d+03*prefreezeout * kstick(1140 - surface_start)

    !CH3CCH_gas -> CH3CCH_0001 (freezeout)
    kall(186) = 2.29235956d+03*prefreezeout * kstick(1141 - surface_start)

    !CH3CH2OH_gas -> CH3CH2OH_0001 (freezeout)
    kall(187) = 2.13763631d+03*prefreezeout * kstick(1142 - surface_start)

    !CH3CHCH2_gas -> CH3CHCH2_0001 (freezeout)
    kall(188) = 2.23711386d+03*prefreezeout * kstick(1143 - surface_start)

    !CH3CHO_gas -> CH3CHO_0001 (freezeout)
    kall(189) = 2.18567908d+03*prefreezeout * kstick(1144 - surface_start)

    !CH3CN_gas -> CH3CN_0001 (freezeout)
    kall(190) = 2.26423138d+03*prefreezeout * kstick(1145 - surface_start)

    !CH3COCH3_gas -> CH3COCH3_0001 (freezeout)
    kall(191) = 1.90370099d+03*prefreezeout * kstick(1148 - surface_start)

    !CH3NH2_gas -> CH3NH2_0001 (freezeout)
    kall(192) = 2.60394548d+03*prefreezeout * kstick(1149 - surface_start)

    !CH3OCH2_gas -> CH3OCH2_0001 (freezeout)
    kall(193) = 2.16125732d+03*prefreezeout * kstick(1150 - surface_start)

    !CH3OCH3_gas -> CH3OCH3_0001 (freezeout)
    kall(194) = 2.13763631d+03*prefreezeout * kstick(1151 - surface_start)

    !CHD2CN_gas -> CHD2CN_0001 (freezeout)
    kall(195) = 2.21094789d+03*prefreezeout * kstick(1147 - surface_start)

    !CHDNH2_gas -> CHDNH2_0001 (freezeout)
    kall(196) = 2.60394548d+03*prefreezeout * kstick(1129 - surface_start)

    !CHDNHD_gas -> CHDNHD_0001 (freezeout)
    kall(197) = 2.56293590d+03*prefreezeout * kstick(1130 - surface_start)

    !CHDND2_gas -> CHDND2_0001 (freezeout)
    kall(198) = 2.52380481d+03*prefreezeout * kstick(1131 - surface_start)

    !CN_gas -> CN_0001 (freezeout)
    kall(199) = 2.84332209d+03*prefreezeout * kstick(1152 - surface_start)

    !CS_gas -> CS_0001 (freezeout)
    kall(200) = 2.18567908d+03*prefreezeout * kstick(1153 - surface_start)

    !DCS_gas -> DCS_0001 (freezeout)
    kall(201) = 2.13763631d+03*prefreezeout * kstick(1188 - surface_start)

    !DNC_gas -> DNC_0001 (freezeout)
    kall(202) = 2.73989373d+03*prefreezeout * kstick(1190 - surface_start)

    !DNCCC_gas -> DNCCC_0001 (freezeout)
    kall(203) = 2.01053233d+03*prefreezeout * kstick(1192 - surface_start)

    !DNCO_gas -> DNCO_0001 (freezeout)
    kall(204) = 2.18567908d+03*prefreezeout * kstick(1194 - surface_start)

    !DNO_gas -> DNO_0001 (freezeout)
    kall(205) = 2.56293590d+03*prefreezeout * kstick(1196 - surface_start)

    !HS_gas -> HS_0001 (freezeout)
    kall(206) = 2.52380481d+03*prefreezeout * kstick(1197 - surface_start)

    !DS_gas -> DS_0001 (freezeout)
    kall(207) = 2.48641304d+03*prefreezeout * kstick(1198 - surface_start)

    !FeH_gas -> FeH_0001 (freezeout)
    kall(208) = 1.92032752d+03*prefreezeout * kstick(980 - surface_start)

    !FeD_gas -> FeD_0001 (freezeout)
    kall(209) = 1.90370099d+03*prefreezeout * kstick(1317 - surface_start)

    !H2CCN_gas -> H2CCN_0001 (freezeout)
    kall(210) = 2.29235956d+03*prefreezeout * kstick(1154 - surface_start)

    !HDCCN_gas -> HDCCN_0001 (freezeout)
    kall(211) = 2.26423138d+03*prefreezeout * kstick(1155 - surface_start)

    !D2CCN_gas -> D2CCN_0001 (freezeout)
    kall(212) = 2.23711386d+03*prefreezeout * kstick(1156 - surface_start)

    !H2CCO_gas -> H2CCO_0001 (freezeout)
    kall(213) = 2.23711386d+03*prefreezeout * kstick(1157 - surface_start)

    !HDCCO_gas -> HDCCO_0001 (freezeout)
    kall(214) = 2.21094789d+03*prefreezeout * kstick(1158 - surface_start)

    !D2CCO_gas -> D2CCO_0001 (freezeout)
    kall(215) = 2.18567908d+03*prefreezeout * kstick(1159 - surface_start)

    !H2CN_gas -> H2CN_0001 (freezeout)
    kall(216) = 2.73989373d+03*prefreezeout * kstick(1160 - surface_start)

    !HDCN_gas -> HDCN_0001 (freezeout)
    kall(217) = 2.69223977d+03*prefreezeout * kstick(1161 - surface_start)

    !D2CN_gas -> D2CN_0001 (freezeout)
    kall(218) = 2.64698882d+03*prefreezeout * kstick(1162 - surface_start)

    !H2CS_gas -> H2CS_0001 (freezeout)
    kall(219) = 2.13763631d+03*prefreezeout * kstick(1163 - surface_start)

    !HDCS_gas -> HDCS_0001 (freezeout)
    kall(220) = 2.11477323d+03*prefreezeout * kstick(1164 - surface_start)

    !D2CS_gas -> D2CS_0001 (freezeout)
    kall(221) = 2.09262840d+03*prefreezeout * kstick(1165 - surface_start)

    !H2S_gas -> H2S_0001 (freezeout)
    kall(222) = 2.48641304d+03*prefreezeout * kstick(1166 - surface_start)

    !HDS_gas -> HDS_0001 (freezeout)
    kall(223) = 2.45063545d+03*prefreezeout * kstick(1167 - surface_start)

    !D2S_gas -> D2S_0001 (freezeout)
    kall(224) = 2.41635914d+03*prefreezeout * kstick(1168 - surface_start)

    !HC2O_gas -> HC2O_0001 (freezeout)
    kall(225) = 2.26423138d+03*prefreezeout * kstick(1169 - surface_start)

    !DC2O_gas -> DC2O_0001 (freezeout)
    kall(226) = 2.23711386d+03*prefreezeout * kstick(1170 - surface_start)

    !HC3N_gas -> HC3N_0001 (freezeout)
    kall(227) = 2.03014775d+03*prefreezeout * kstick(1171 - surface_start)

    !DC3N_gas -> DC3N_0001 (freezeout)
    kall(228) = 2.01053233d+03*prefreezeout * kstick(1172 - surface_start)

    !HC4N_gas -> HC4N_0001 (freezeout)
    kall(229) = 1.82659582d+03*prefreezeout * kstick(1173 - surface_start)

    !DC4N_gas -> DC4N_0001 (freezeout)
    kall(230) = 1.81226935d+03*prefreezeout * kstick(1174 - surface_start)

    !HC5N_gas -> HC5N_0001 (freezeout)
    kall(231) = 1.67410272d+03*prefreezeout * kstick(1175 - surface_start)

    !HC6N_gas -> HC6N_0001 (freezeout)
    kall(232) = 1.55436535d+03*prefreezeout * kstick(1176 - surface_start)

    !HC7N_gas -> HC7N_0001 (freezeout)
    kall(233) = 1.45711939d+03*prefreezeout * kstick(1177 - surface_start)

    !HC8N_gas -> HC8N_0001 (freezeout)
    kall(234) = 1.37610396d+03*prefreezeout * kstick(1178 - surface_start)

    !HC9N_gas -> HC9N_0001 (freezeout)
    kall(235) = 1.30725460d+03*prefreezeout * kstick(1179 - surface_start)

    !HCCNC_gas -> HCCNC_0001 (freezeout)
    kall(236) = 2.03014775d+03*prefreezeout * kstick(1180 - surface_start)

    !DCCNC_gas -> DCCNC_0001 (freezeout)
    kall(237) = 2.01053233d+03*prefreezeout * kstick(1181 - surface_start)

    !HCN_gas -> HCN_0001 (freezeout)
    kall(238) = 2.79017120d+03*prefreezeout * kstick(1182 - surface_start)

    !DCN_gas -> DCN_0001 (freezeout)
    kall(239) = 2.73989373d+03*prefreezeout * kstick(1183 - surface_start)

    !HCNCC_gas -> HCNCC_0001 (freezeout)
    kall(240) = 2.03014775d+03*prefreezeout * kstick(1184 - surface_start)

    !DCNCC_gas -> DCNCC_0001 (freezeout)
    kall(241) = 2.01053233d+03*prefreezeout * kstick(1185 - surface_start)

    !HCOOCH3_gas -> HCOOCH3_0001 (freezeout)
    kall(242) = 1.87170374d+03*prefreezeout * kstick(1186 - surface_start)

    !HCS_gas -> HCS_0001 (freezeout)
    kall(243) = 2.16125732d+03*prefreezeout * kstick(1187 - surface_start)

    !HNC_gas -> HNC_0001 (freezeout)
    kall(244) = 2.79017120d+03*prefreezeout * kstick(1189 - surface_start)

    !HNCCC_gas -> HNCCC_0001 (freezeout)
    kall(245) = 2.03014775d+03*prefreezeout * kstick(1191 - surface_start)

    !HNCO_gas -> HNCO_0001 (freezeout)
    kall(246) = 2.21094789d+03*prefreezeout * kstick(1193 - surface_start)

    !HNO_gas -> HNO_0001 (freezeout)
    kall(247) = 2.60394548d+03*prefreezeout * kstick(1195 - surface_start)

    !HSS_gas -> HSS_0001 (freezeout)
    kall(248) = 1.79827479d+03*prefreezeout * kstick(1215 - surface_start)

    !DSS_gas -> DSS_0001 (freezeout)
    kall(249) = 1.78459950d+03*prefreezeout * kstick(1216 - surface_start)

    !HSSH_gas -> HSSH_0001 (freezeout)
    kall(250) = 1.78459950d+03*prefreezeout * kstick(1217 - surface_start)

    !HSSD_gas -> HSSD_0001 (freezeout)
    kall(251) = 1.77123152d+03*prefreezeout * kstick(1218 - surface_start)

    !DSSH_gas -> DSSH_0001 (freezeout)
    kall(252) = 1.77123152d+03*prefreezeout * kstick(1219 - surface_start)

    !DSSD_gas -> DSSD_0001 (freezeout)
    kall(253) = 1.75815952d+03*prefreezeout * kstick(1220 - surface_start)

    !N2O_gas -> N2O_0001 (freezeout)
    kall(254) = 2.18567908d+03*prefreezeout * kstick(1200 - surface_start)

    !NC4N_gas -> NC4N_0001 (freezeout)
    kall(255) = 1.66305242d+03*prefreezeout * kstick(1201 - surface_start)

    !NC6N_gas -> NC6N_0001 (freezeout)
    kall(256) = 1.44981548d+03*prefreezeout * kstick(1202 - surface_start)

    !NC8N_gas -> NC8N_0001 (freezeout)
    kall(257) = 1.30197274d+03*prefreezeout * kstick(1203 - surface_start)

    !NH2CHO_gas -> NH2CHO_0001 (freezeout)
    kall(258) = 2.16125732d+03*prefreezeout * kstick(1204 - surface_start)

    !NH2CDO_gas -> NH2CDO_0001 (freezeout)
    kall(259) = 2.13763631d+03*prefreezeout * kstick(1205 - surface_start)

    !NHDCHO_gas -> NHDCHO_0001 (freezeout)
    kall(260) = 2.13763631d+03*prefreezeout * kstick(1206 - surface_start)

    !NHDCDO_gas -> NHDCDO_0001 (freezeout)
    kall(261) = 2.11477323d+03*prefreezeout * kstick(1207 - surface_start)

    !ND2CHO_gas -> ND2CHO_0001 (freezeout)
    kall(262) = 2.11477323d+03*prefreezeout * kstick(1208 - surface_start)

    !ND2CDO_gas -> ND2CDO_0001 (freezeout)
    kall(263) = 2.09262840d+03*prefreezeout * kstick(1209 - surface_start)

    !NH2CN_gas -> NH2CN_0001 (freezeout)
    kall(264) = 2.23711386d+03*prefreezeout * kstick(1210 - surface_start)

    !NHDCN_gas -> NHDCN_0001 (freezeout)
    kall(265) = 2.21094789d+03*prefreezeout * kstick(1211 - surface_start)

    !ND2CN_gas -> ND2CN_0001 (freezeout)
    kall(266) = 2.18567908d+03*prefreezeout * kstick(1212 - surface_start)

    !NO_gas -> NO_0001 (freezeout)
    kall(267) = 2.64698882d+03*prefreezeout * kstick(1221 - surface_start)

    !NO2_gas -> NO2_0001 (freezeout)
    kall(268) = 2.13763631d+03*prefreezeout * kstick(1222 - surface_start)

    !NS_gas -> NS_0001 (freezeout)
    kall(269) = 2.13763631d+03*prefreezeout * kstick(1223 - surface_start)

    !OCN_gas -> OCN_0001 (freezeout)
    kall(270) = 2.23711386d+03*prefreezeout * kstick(1224 - surface_start)

    !OCS_gas -> OCS_0001 (freezeout)
    kall(271) = 1.87170374d+03*prefreezeout * kstick(1225 - surface_start)

    !S2_gas -> S2_0001 (freezeout)
    kall(272) = 1.81226935d+03*prefreezeout * kstick(1226 - surface_start)

    !SO_gas -> SO_0001 (freezeout)
    kall(273) = 2.09262840d+03*prefreezeout * kstick(1227 - surface_start)

    !SO2_gas -> SO2_0001 (freezeout)
    kall(274) = 1.81226935d+03*prefreezeout * kstick(1228 - surface_start)

    !CH3OCHO_gas -> CH3OCHO_0001 (freezeout)
    kall(275) = 1.87170374d+03*prefreezeout * kstick(1229 - surface_start)

    !Si_gas -> Si_0001 (freezeout)
    kall(276) = 2.73989373d+03*prefreezeout * kstick(1230 - surface_start)

    !SiS_gas -> SiS_0001 (freezeout)
    kall(277) = 1.87170374d+03*prefreezeout * kstick(1231 - surface_start)

    !SiN_gas -> SiN_0001 (freezeout)
    kall(278) = 2.23711386d+03*prefreezeout * kstick(1232 - surface_start)

    !SiC_gas -> SiC_0001 (freezeout)
    kall(279) = 2.29235956d+03*prefreezeout * kstick(1233 - surface_start)

    !SiH_gas -> SiH_0001 (freezeout)
    kall(280) = 2.69223977d+03*prefreezeout * kstick(1234 - surface_start)

    !SiD_gas -> SiD_0001 (freezeout)
    kall(281) = 2.64698882d+03*prefreezeout * kstick(1270 - surface_start)

    !SiH2_gas -> SiH2_0001 (freezeout)
    kall(282) = 2.64698882d+03*prefreezeout * kstick(1235 - surface_start)

    !SiH3_gas -> SiH3_0001 (freezeout)
    kall(283) = 2.60394548d+03*prefreezeout * kstick(1236 - surface_start)

    !SiH4_gas -> SiH4_0001 (freezeout)
    kall(284) = 2.56293590d+03*prefreezeout * kstick(1237 - surface_start)

    !P_gas -> P_0001 (freezeout)
    kall(285) = 2.60394548d+03*prefreezeout * kstick(1388 - surface_start)

    !PH_gas -> PH_0001 (freezeout)
    kall(286) = 2.56293590d+03*prefreezeout * kstick(1390 - surface_start)

    !PD_gas -> PD_0001 (freezeout)
    kall(287) = 2.52380481d+03*prefreezeout * kstick(1391 - surface_start)

    !PO_gas -> PO_0001 (freezeout)
    kall(288) = 2.11477323d+03*prefreezeout * kstick(1389 - surface_start)

    !PH2_gas -> PH2_0001 (freezeout)
    kall(289) = 2.52380481d+03*prefreezeout * kstick(1392 - surface_start)

    !PHD_gas -> PHD_0001 (freezeout)
    kall(290) = 2.48641304d+03*prefreezeout * kstick(1393 - surface_start)

    !PD2_gas -> PD2_0001 (freezeout)
    kall(291) = 2.45063545d+03*prefreezeout * kstick(1394 - surface_start)

    !C3P_gas -> C3P_0001 (freezeout)
    kall(292) = 1.77123152d+03*prefreezeout * kstick(1398 - surface_start)

    !C4P_gas -> C4P_0001 (freezeout)
    kall(293) = 1.63116986d+03*prefreezeout * kstick(1399 - surface_start)

    !CH2PH_gas -> CH2PH_0001 (freezeout)
    kall(294) = 2.13763631d+03*prefreezeout * kstick(1400 - surface_start)

    !CH2PD_gas -> CH2PD_0001 (freezeout)
    kall(295) = 2.11477323d+03*prefreezeout * kstick(1402 - surface_start)

    !CHDPD_gas -> CHDPD_0001 (freezeout)
    kall(296) = 2.09262840d+03*prefreezeout * kstick(1403 - surface_start)

    !CP_gas -> CP_0001 (freezeout)
    kall(297) = 2.21094789d+03*prefreezeout * kstick(1396 - surface_start)

    !DCP_gas -> DCP_0001 (freezeout)
    kall(298) = 2.16125732d+03*prefreezeout * kstick(1406 - surface_start)

    !DF_gas -> DF_0001 (freezeout)
    kall(299) = 3.16375676d+03*prefreezeout * kstick(1265 - surface_start)

    !DPO_gas -> DPO_0001 (freezeout)
    kall(300) = 2.07116498d+03*prefreezeout * kstick(1410 - surface_start)

    !HCCP_gas -> HCCP_0001 (freezeout)
    kall(301) = 1.93739743d+03*prefreezeout * kstick(1407 - surface_start)

    !DCCP_gas -> DCCP_0001 (freezeout)
    kall(302) = 1.92032752d+03*prefreezeout * kstick(1408 - surface_start)

    !HCP_gas -> HCP_0001 (freezeout)
    kall(303) = 2.18567908d+03*prefreezeout * kstick(1405 - surface_start)

    !HF_gas -> HF_0001 (freezeout)
    kall(304) = 3.24188598d+03*prefreezeout * kstick(1264 - surface_start)

    !HPO_gas -> HPO_0001 (freezeout)
    kall(305) = 2.09262840d+03*prefreezeout * kstick(1409 - surface_start)

    !PN_gas -> PN_0001 (freezeout)
    kall(306) = 2.16125732d+03*prefreezeout * kstick(1395 - surface_start)

    !SiC2CH3_gas -> SiC2CH3_0001 (freezeout)
    kall(307) = 1.77123152d+03*prefreezeout * kstick(1238 - surface_start)

    !SiC3H_gas -> SiC3H_0001 (freezeout)
    kall(308) = 1.79827479d+03*prefreezeout * kstick(1239 - surface_start)

    !SiC3H5_gas -> SiC3H5_0001 (freezeout)
    kall(309) = 1.74537274d+03*prefreezeout * kstick(1240 - surface_start)

    !SiC4_gas -> SiC4_0001 (freezeout)
    kall(310) = 1.66305242d+03*prefreezeout * kstick(1241 - surface_start)

    !SiC4H_gas -> SiC4H_0001 (freezeout)
    kall(311) = 1.65221808d+03*prefreezeout * kstick(1242 - surface_start)

    !SiC6H_gas -> SiC6H_0001 (freezeout)
    kall(312) = 1.44262033d+03*prefreezeout * kstick(1243 - surface_start)

    !SiC8H_gas -> SiC8H_0001 (freezeout)
    kall(313) = 1.29675439d+03*prefreezeout * kstick(1244 - surface_start)

    !c_HCCHSi_gas -> c_HCCHSi_0001 (freezeout)
    kall(314) = 1.97294898d+03*prefreezeout * kstick(1245 - surface_start)

    !c_SiC2_gas -> c_SiC2_0001 (freezeout)
    kall(315) = 2.01053233d+03*prefreezeout * kstick(1246 - surface_start)

    !l_SiC3_gas -> l_SiC3_0001 (freezeout)
    kall(316) = 1.81226935d+03*prefreezeout * kstick(1251 - surface_start)

    !l_C3H_gas -> l_C3H_0001 (freezeout)
    kall(317) = 2.38348197d+03*prefreezeout * kstick(1247 - surface_start)

    !l_C3D_gas -> l_C3D_0001 (freezeout)
    kall(318) = 2.35191128d+03*prefreezeout * kstick(1248 - surface_start)

    !c_C3H_gas -> c_C3H_0001 (freezeout)
    kall(319) = 2.38348197d+03*prefreezeout * kstick(1249 - surface_start)

    !c_C3D_gas -> c_C3D_0001 (freezeout)
    kall(320) = 2.35191128d+03*prefreezeout * kstick(1250 - surface_start)

    !l_C3H2_gas -> l_C3H2_0001 (freezeout)
    kall(321) = 2.35191128d+03*prefreezeout * kstick(1252 - surface_start)

    !l_C3HD_gas -> l_C3HD_0001 (freezeout)
    kall(322) = 2.32156277d+03*prefreezeout * kstick(1253 - surface_start)

    !l_C3D2_gas -> l_C3D2_0001 (freezeout)
    kall(323) = 2.29235956d+03*prefreezeout * kstick(1254 - surface_start)

    !c_C3H2_gas -> c_C3H2_0001 (freezeout)
    kall(324) = 2.35191128d+03*prefreezeout * kstick(1255 - surface_start)

    !c_C3HD_gas -> c_C3HD_0001 (freezeout)
    kall(325) = 2.32156277d+03*prefreezeout * kstick(1256 - surface_start)

    !c_C3D2_gas -> c_C3D2_0001 (freezeout)
    kall(326) = 2.29235956d+03*prefreezeout * kstick(1257 - surface_start)

    !Mg_gas -> Mg_0001 (freezeout)
    kall(327) = 2.95942346d+03*prefreezeout * kstick(1258 - surface_start)

    !MgH_gas -> MgH_0001 (freezeout)
    kall(328) = 2.89963097d+03*prefreezeout * kstick(1259 - surface_start)

    !MgD_gas -> MgD_0001 (freezeout)
    kall(329) = 2.84332209d+03*prefreezeout * kstick(1266 - surface_start)

    !MgH2_gas -> MgH2_0001 (freezeout)
    kall(330) = 2.84332209d+03*prefreezeout * kstick(1260 - surface_start)

    !MgHD_gas -> MgHD_0001 (freezeout)
    kall(331) = 2.79017120d+03*prefreezeout * kstick(1267 - surface_start)

    !MgD2_gas -> MgD2_0001 (freezeout)
    kall(332) = 2.73989373d+03*prefreezeout * kstick(1268 - surface_start)

    !HCl_gas -> HCl_0001 (freezeout)
    kall(333) = 2.41635914d+03*prefreezeout * kstick(1272 - surface_start)

    !DCl_gas -> DCl_0001 (freezeout)
    kall(334) = 2.38348197d+03*prefreezeout * kstick(1273 - surface_start)

    !H2Cl_gas -> H2Cl_0001 (freezeout)
    kall(335) = 2.38348197d+03*prefreezeout * kstick(1274 - surface_start)

    !HDCl_gas -> HDCl_0001 (freezeout)
    kall(336) = 2.35191128d+03*prefreezeout * kstick(1275 - surface_start)

    !D2Cl_gas -> D2Cl_0001 (freezeout)
    kall(337) = 2.32156277d+03*prefreezeout * kstick(1276 - surface_start)

    !CCl_gas -> CCl_0001 (freezeout)
    kall(338) = 2.11477323d+03*prefreezeout * kstick(1277 - surface_start)

    !ClO_gas -> ClO_0001 (freezeout)
    kall(339) = 2.03014775d+03*prefreezeout * kstick(1278 - surface_start)

    !Na_gas -> Na_0001 (freezeout)
    kall(340) = 3.02307426d+03*prefreezeout * kstick(1261 - surface_start)

    !NaH_gas -> NaH_0001 (freezeout)
    kall(341) = 2.95942346d+03*prefreezeout * kstick(1262 - surface_start)

    !NaD_gas -> NaD_0001 (freezeout)
    kall(342) = 2.89963097d+03*prefreezeout * kstick(1269 - surface_start)

    !H_0001 -> H_gas (evaporation)
    kall(343) = 3.54191744d+12*exp(-500.0*invTd)

    !D_0001 -> D_gas (evaporation)
    kall(344) = 2.55656768d+12*exp(-521.0*invTd)

    !p_H2_0001 -> p_H2_gas (evaporation)
    kall(345) = 2.37599045d+12*exp(-450.0*invTd)

    !p_D2_0001 -> p_D2_gas (evaporation)
    kall(346) = 1.68007896d+12*exp(-450.0*invTd)

    !o_H2_0001 -> o_H2_gas (evaporation)
    kall(347) = 2.37599045d+12*exp(-450.0*invTd)

    !o_D2_0001 -> o_D2_gas (evaporation)
    kall(348) = 1.68007896d+12*exp(-450.0*invTd)

    !HD_0001 -> HD_gas (evaporation)
    kall(349) = 1.93998808d+12*exp(-450.0*invTd)

    !He_0001 -> He_gas (evaporation)
    kall(350) = 7.91996816d+11*exp(-100.0*invTd)

    !O_0001 -> O_gas (evaporation)
    kall(351) = 1.58399363d+12*exp(-1600.0*invTd)

    !O2_0001 -> O2_gas (evaporation)
    kall(352) = 9.69994038d+11*exp(-1200.0*invTd)

    !O3_0001 -> O3_gas (evaporation)
    kall(353) = 1.04771331d+12*exp(-2100.0*invTd)

    !OH_0001 -> OH_gas (evaporation)
    kall(354) = 2.60560084d+12*exp(-4600.0*invTd)

    !OD_0001 -> OD_gas (evaporation)
    kall(355) = 2.53218886d+12*exp(-4600.0*invTd)

    !H2O_0001 -> H2O_gas (evaporation)
    kall(356) = 2.79390215d+12*exp(-5600.0*invTd)

    !HDO_0001 -> HDO_gas (evaporation)
    kall(357) = 2.71938466d+12*exp(-5600.0*invTd)

    !D2O_0001 -> D2O_gas (evaporation)
    kall(358) = 2.65052831d+12*exp(-5600.0*invTd)

    !O2H_0001 -> O2H_gas (evaporation)
    kall(359) = 1.94976138d+12*exp(-5000.0*invTd)

    !O2D_0001 -> O2D_gas (evaporation)
    kall(360) = 1.92087443d+12*exp(-5000.0*invTd)

    !HOOH_0001 -> HOOH_gas (evaporation)
    kall(361) = 2.10421251d+12*exp(-6000.0*invTd)

    !HOOD_0001 -> HOOD_gas (evaporation)
    kall(362) = 2.07393449d+12*exp(-6000.0*invTd)

    !DOOD_0001 -> DOOD_gas (evaporation)
    kall(363) = 2.04492698d+12*exp(-6000.0*invTd)

    !Fe_0001 -> Fe_gas (evaporation)
    kall(364) = 1.37177872d+12*exp(-4200.0*invTd)

    !FeH_0001 -> FeH_gas (evaporation)
    kall(365) = 1.43067991d+12*exp(-4650.0*invTd)

    !N_0001 -> N_gas (evaporation)
    kall(366) = 1.13594070d+12*exp(-720.0*invTd)

    !S_0001 -> S_gas (evaporation)
    kall(367) = 1.42779256d+12*exp(-2600.0*invTd)

    !C_0001 -> C_gas (evaporation)
    kall(368) = 4.57259575d+12*exp(-10000.0*invTd)

    !C2_0001 -> C2_gas (evaporation)
    kall(369) = 3.23331346d+12*exp(-10000.0*invTd)

    !CO_0001 -> CO_gas (evaporation)
    kall(370) = 1.07930973d+12*exp(-1300.0*invTd)

    !HCO_0001 -> HCO_gas (evaporation)
    kall(371) = 1.44098697d+12*exp(-2400.0*invTd)

    !DCO_0001 -> DCO_gas (evaporation)
    kall(372) = 1.41676697d+12*exp(-2400.0*invTd)

    !H2CO_0001 -> H2CO_gas (evaporation)
    kall(373) = 1.93998808d+12*exp(-4500.0*invTd)

    !HDCO_0001 -> HDCO_gas (evaporation)
    kall(374) = 1.90844145d+12*exp(-4500.0*invTd)

    !D2CO_0001 -> D2CO_gas (evaporation)
    kall(375) = 1.87838538d+12*exp(-4500.0*invTd)

    !CH2OH_0001 -> CH2OH_gas (evaporation)
    kall(376) = 1.88711741d+12*exp(-4400.0*invTd)

    !CD2OD_0001 -> CD2OD_gas (evaporation)
    kall(377) = 1.80193994d+12*exp(-4400.0*invTd)

    !CH2OD_0001 -> CH2OD_gas (evaporation)
    kall(378) = 1.85739717d+12*exp(-4400.0*invTd)

    !CHDOH_0001 -> CHDOH_gas (evaporation)
    kall(379) = 1.85739717d+12*exp(-4400.0*invTd)

    !CHDOD_0001 -> CHDOD_gas (evaporation)
    kall(380) = 1.82903830d+12*exp(-4400.0*invTd)

    !CD2OH_0001 -> CD2OH_gas (evaporation)
    kall(381) = 1.82903830d+12*exp(-4400.0*invTd)

    !CH3O_0001 -> CH3O_gas (evaporation)
    kall(382) = 1.88711741d+12*exp(-4400.0*invTd)

    !CHD2O_0001 -> CHD2O_gas (evaporation)
    kall(383) = 1.82903830d+12*exp(-4400.0*invTd)

    !CH2DO_0001 -> CH2DO_gas (evaporation)
    kall(384) = 1.85739717d+12*exp(-4400.0*invTd)

    !CD3O_0001 -> CD3O_gas (evaporation)
    kall(385) = 1.80193994d+12*exp(-4400.0*invTd)

    !CH3OH_0001 -> CH3OH_gas (evaporation)
    kall(386) = 1.97999204d+12*exp(-5000.0*invTd)

    !CH3OD_0001 -> CH3OD_gas (evaporation)
    kall(387) = 1.94976138d+12*exp(-5000.0*invTd)

    !CHD2OH_0001 -> CHD2OH_gas (evaporation)
    kall(388) = 1.92087443d+12*exp(-5000.0*invTd)

    !CHD2OD_0001 -> CHD2OD_gas (evaporation)
    kall(389) = 1.89323451d+12*exp(-5000.0*invTd)

    !CH2DOH_0001 -> CH2DOH_gas (evaporation)
    kall(390) = 1.94976138d+12*exp(-5000.0*invTd)

    !CH2DOD_0001 -> CH2DOD_gas (evaporation)
    kall(391) = 1.92087443d+12*exp(-5000.0*invTd)

    !CD3OD_0001 -> CD3OD_gas (evaporation)
    kall(392) = 1.86675440d+12*exp(-5000.0*invTd)

    !CD3OH_0001 -> CD3OH_gas (evaporation)
    kall(393) = 1.89323451d+12*exp(-5000.0*invTd)

    !CH_0001 -> CH_gas (evaporation)
    kall(394) = 1.33614202d+12*exp(-925.0*invTd)

    !CD_0001 -> CD_gas (evaporation)
    kall(395) = 1.28753866d+12*exp(-925.0*invTd)

    !CH2_0001 -> CH2_gas (evaporation)
    kall(396) = 1.58399363d+12*exp(-1400.0*invTd)

    !CHD_0001 -> CHD_gas (evaporation)
    kall(397) = 1.53028323d+12*exp(-1400.0*invTd)

    !CD2_0001 -> CD2_gas (evaporation)
    kall(398) = 1.48169037d+12*exp(-1400.0*invTd)

    !CH3_0001 -> CH3_gas (evaporation)
    kall(399) = 1.63594159d+12*exp(-1600.0*invTd)

    !CH2D_0001 -> CH2D_gas (evaporation)
    kall(400) = 1.58399363d+12*exp(-1600.0*invTd)

    !CHD2_0001 -> CHD2_gas (evaporation)
    kall(401) = 1.88206488d+12*exp(-2400.0*invTd)

    !CD3_0001 -> CD3_gas (evaporation)
    kall(402) = 1.49340352d+12*exp(-1600.0*invTd)

    !CH4_0001 -> CH4_gas (evaporation)
    kall(403) = 1.22695619d+12*exp(-960.0*invTd)

    !CH3D_0001 -> CH3D_gas (evaporation)
    kall(404) = 1.19032235d+12*exp(-960.0*invTd)

    !CH2D2_0001 -> CH2D2_gas (evaporation)
    kall(405) = 1.15678539d+12*exp(-960.0*invTd)

    !CHD3_0001 -> CHD3_gas (evaporation)
    kall(406) = 1.12593222d+12*exp(-960.0*invTd)

    !CD4_0001 -> CD4_gas (evaporation)
    kall(407) = 1.09742298d+12*exp(-960.0*invTd)

    !CO2_0001 -> CO2_gas (evaporation)
    kall(408) = 1.21762559d+12*exp(-2600.0*invTd)

    !HCOOH_0001 -> HCOOH_gas (evaporation)
    kall(409) = 1.74301897d+12*exp(-5570.0*invTd)

    !HCOOD_0001 -> HCOOD_gas (evaporation)
    kall(410) = 1.72437652d+12*exp(-5570.0*invTd)

    !DCOOH_0001 -> DCOOH_gas (evaporation)
    kall(411) = 1.72437652d+12*exp(-5570.0*invTd)

    !DCOOD_0001 -> DCOOD_gas (evaporation)
    kall(412) = 1.70631972d+12*exp(-5570.0*invTd)

    !HOCO_0001 -> HOCO_gas (evaporation)
    kall(413) = 1.05599575d+12*exp(-2000.0*invTd)

    !DOCO_0001 -> DOCO_gas (evaporation)
    kall(414) = 1.04445447d+12*exp(-2000.0*invTd)

    !NH_0001 -> NH_gas (evaporation)
    kall(415) = 2.08542452d+12*exp(-2600.0*invTd)

    !ND_0001 -> ND_gas (evaporation)
    kall(416) = 2.01920361d+12*exp(-2600.0*invTd)

    !NH2_0001 -> NH2_gas (evaporation)
    kall(417) = 2.24010528d+12*exp(-3200.0*invTd)

    !NHD_0001 -> NHD_gas (evaporation)
    kall(418) = 2.17322133d+12*exp(-3200.0*invTd)

    !ND2_0001 -> ND2_gas (evaporation)
    kall(419) = 2.11199151d+12*exp(-3200.0*invTd)

    !NH3_0001 -> NH3_gas (evaporation)
    kall(420) = 2.84911720d+12*exp(-5500.0*invTd)

    !NH2D_0001 -> NH2D_gas (evaporation)
    kall(421) = 2.76884423d+12*exp(-5500.0*invTd)

    !NHD2_0001 -> NHD2_gas (evaporation)
    kall(422) = 2.69499507d+12*exp(-5500.0*invTd)

    !ND3_0001 -> ND3_gas (evaporation)
    kall(423) = 2.62675627d+12*exp(-5500.0*invTd)

    !C10_0001 -> C10_gas (evaporation)
    kall(424) = 1.29332538d+12*exp(-8000.0*invTd)

    !C10H_0001 -> C10H_gas (evaporation)
    kall(425) = 1.32369857d+12*exp(-8450.0*invTd)

    !C10H2_0001 -> C10H2_gas (evaporation)
    kall(426) = 1.32603984d+12*exp(-8550.0*invTd)

    !C10N_0001 -> C10N_gas (evaporation)
    kall(427) = 1.28363742d+12*exp(-8800.0*invTd)

    !C11_0001 -> C11_gas (evaporation)
    kall(428) = 1.35083431d+12*exp(-9600.0*invTd)

    !C2H2_0001 -> C2H2_gas (evaporation)
    kall(429) = 1.58002868d+12*exp(-2587.0*invTd)

    !C2HD_0001 -> C2HD_gas (evaporation)
    kall(430) = 1.55049283d+12*exp(-2587.0*invTd)

    !C2D2_0001 -> C2D2_gas (evaporation)
    kall(431) = 1.52255374d+12*exp(-2587.0*invTd)

    !C2H3_0001 -> C2H3_gas (evaporation)
    kall(432) = 1.61306016d+12*exp(-2800.0*invTd)

    !C2H2D_0001 -> C2H2D_gas (evaporation)
    kall(433) = 1.58399363d+12*exp(-2800.0*invTd)

    !C2HD2_0001 -> C2HD2_gas (evaporation)
    kall(434) = 1.55644381d+12*exp(-2800.0*invTd)

    !C2D3_0001 -> C2D3_gas (evaporation)
    kall(435) = 1.53028323d+12*exp(-2800.0*invTd)

    !C2H4_0001 -> C2H4_gas (evaporation)
    kall(436) = 1.49673330d+12*exp(-2500.0*invTd)

    !C2H3D_0001 -> C2H3D_gas (evaporation)
    kall(437) = 1.47070117d+12*exp(-2500.0*invTd)

    !C2H2D2_0001 -> C2H2D2_gas (evaporation)
    kall(438) = 1.61017712d+12*exp(-3100.0*invTd)

    !C2HD3_0001 -> C2HD3_gas (evaporation)
    kall(439) = 1.42246827d+12*exp(-2500.0*invTd)

    !C2D4_0001 -> C2D4_gas (evaporation)
    kall(440) = 1.40006580d+12*exp(-2500.0*invTd)

    !C2H5_0001 -> C2H5_gas (evaporation)
    kall(441) = 1.63770351d+12*exp(-3100.0*invTd)

    !C2H6_0001 -> C2H6_gas (evaporation)
    kall(442) = 1.15678539d+12*exp(-1600.0*invTd)

    !C3_0001 -> C3_gas (evaporation)
    kall(443) = 1.29332538d+12*exp(-2400.0*invTd)

    !C3N_0001 -> C3N_gas (evaporation)
    kall(444) = 1.26719491d+12*exp(-3200.0*invTd)

    !C3O_0001 -> C3O_gas (evaporation)
    kall(445) = 1.15190883d+12*exp(-2750.0*invTd)

    !C3S_0001 -> C3S_gas (evaporation)
    kall(446) = 1.13640464d+12*exp(-3500.0*invTd)

    !C4_0001 -> C4_gas (evaporation)
    kall(447) = 1.29332538d+12*exp(-3200.0*invTd)

    !C4H_0001 -> C4H_gas (evaporation)
    kall(448) = 1.38330179d+12*exp(-3737.0*invTd)

    !C4D_0001 -> C4D_gas (evaporation)
    kall(449) = 1.36939891d+12*exp(-3737.0*invTd)

    !C4H2_0001 -> C4H2_gas (evaporation)
    kall(450) = 1.44950564d+12*exp(-4187.0*invTd)

    !C4HD_0001 -> C4HD_gas (evaporation)
    kall(451) = 1.43522445d+12*exp(-4187.0*invTd)

    !C4D2_0001 -> C4D2_gas (evaporation)
    kall(452) = 1.42135722d+12*exp(-4187.0*invTd)

    !C4H3_0001 -> C4H3_gas (evaporation)
    kall(453) = 1.51038231d+12*exp(-4637.0*invTd)

    !C4N_0001 -> C4N_gas (evaporation)
    kall(454) = 1.27229430d+12*exp(-4000.0*invTd)

    !C4S_0001 -> C4S_gas (evaporation)
    kall(455) = 1.16129529d+12*exp(-4300.0*invTd)

    !C5_0001 -> C5_gas (evaporation)
    kall(456) = 1.29332538d+12*exp(-4000.0*invTd)

    !C5H_0001 -> C5H_gas (evaporation)
    kall(457) = 1.36606984d+12*exp(-4537.0*invTd)

    !C5D_0001 -> C5D_gas (evaporation)
    kall(458) = 1.35500836d+12*exp(-4537.0*invTd)

    !C5H2_0001 -> C5H2_gas (evaporation)
    kall(459) = 1.42061786d+12*exp(-4987.0*invTd)

    !C5H3_0001 -> C5H3_gas (evaporation)
    kall(460) = 1.47150866d+12*exp(-5437.0*invTd)

    !C5H4_0001 -> C5H4_gas (evaporation)
    kall(461) = 1.51918429d+12*exp(-5887.0*invTd)

    !C5N_0001 -> C5N_gas (evaporation)
    kall(462) = 1.27572830d+12*exp(-4800.0*invTd)

    !C5O_0001 -> C5O_gas (evaporation)
    kall(463) = 1.19837092d+12*exp(-4350.0*invTd)

    !C6_0001 -> C6_gas (evaporation)
    kall(464) = 1.29332538d+12*exp(-4800.0*invTd)

    !C6H_0001 -> C6H_gas (evaporation)
    kall(465) = 1.35438022d+12*exp(-5337.0*invTd)

    !C6H2_0001 -> C6H2_gas (evaporation)
    kall(466) = 1.40076187d+12*exp(-5787.0*invTd)

    !C6H3_0001 -> C6H3_gas (evaporation)
    kall(467) = 1.44447713d+12*exp(-6237.0*invTd)

    !C6H4_0001 -> C6H4_gas (evaporation)
    kall(468) = 1.48580665d+12*exp(-6687.0*invTd)

    !C6H6_0001 -> C6H6_gas (evaporation)
    kall(469) = 1.56221632d+12*exp(-7587.0*invTd)

    !C6N_0001 -> C6N_gas (evaporation)
    kall(470) = 1.27819825d+12*exp(-5600.0*invTd)

    !C7_0001 -> C7_gas (evaporation)
    kall(471) = 1.29332538d+12*exp(-5600.0*invTd)

    !C7H_0001 -> C7H_gas (evaporation)
    kall(472) = 1.34592863d+12*exp(-6137.0*invTd)

    !C7H2_0001 -> C7H2_gas (evaporation)
    kall(473) = 1.38627068d+12*exp(-6587.0*invTd)

    !C7H3_0001 -> C7H3_gas (evaporation)
    kall(474) = 1.42458242d+12*exp(-7037.0*invTd)

    !C7H4_0001 -> C7H4_gas (evaporation)
    kall(475) = 1.46105314d+12*exp(-7487.0*invTd)

    !C7N_0001 -> C7N_gas (evaporation)
    kall(476) = 1.28006016d+12*exp(-6400.0*invTd)

    !C7O_0001 -> C7O_gas (evaporation)
    kall(477) = 1.22183318d+12*exp(-5950.0*invTd)

    !C8_0001 -> C8_gas (evaporation)
    kall(478) = 1.29332538d+12*exp(-6400.0*invTd)

    !C8H_0001 -> C8H_gas (evaporation)
    kall(479) = 1.33953298d+12*exp(-6937.0*invTd)

    !C8H2_0001 -> C8H2_gas (evaporation)
    kall(480) = 1.37522717d+12*exp(-7387.0*invTd)

    !C8H3_0001 -> C8H3_gas (evaporation)
    kall(481) = 1.40932371d+12*exp(-7837.0*invTd)

    !C8H4_0001 -> C8H4_gas (evaporation)
    kall(482) = 1.44195631d+12*exp(-8287.0*invTd)

    !C8N_0001 -> C8N_gas (evaporation)
    kall(483) = 1.28151395d+12*exp(-7200.0*invTd)

    !C9_0001 -> C9_gas (evaporation)
    kall(484) = 1.29332538d+12*exp(-7200.0*invTd)

    !C9H_0001 -> C9H_gas (evaporation)
    kall(485) = 1.33452428d+12*exp(-7737.0*invTd)

    !C9H2_0001 -> C9H2_gas (evaporation)
    kall(486) = 1.36653099d+12*exp(-8187.0*invTd)

    !C9H3_0001 -> C9H3_gas (evaporation)
    kall(487) = 1.39724769d+12*exp(-8637.0*invTd)

    !C9H4_0001 -> C9H4_gas (evaporation)
    kall(488) = 1.42677235d+12*exp(-9087.0*invTd)

    !C9N_0001 -> C9N_gas (evaporation)
    kall(489) = 1.28268055d+12*exp(-8000.0*invTd)

    !C9O_0001 -> C9O_gas (evaporation)
    kall(490) = 1.23599315d+12*exp(-7550.0*invTd)

    !CCH_0001 -> CCH_gas (evaporation)
    kall(491) = 1.73517809d+12*exp(-3000.0*invTd)

    !CCD_0001 -> CCD_gas (evaporation)
    kall(492) = 1.70148210d+12*exp(-3000.0*invTd)

    !CCN_0001 -> CCN_gas (evaporation)
    kall(493) = 1.25883049d+12*exp(-2400.0*invTd)

    !CCO_0001 -> CCO_gas (evaporation)
    kall(494) = 1.10596336d+12*exp(-1950.0*invTd)

    !CCS_0001 -> CCS_gas (evaporation)
    kall(495) = 1.09986986d+12*exp(-2700.0*invTd)

    !CD3CN_0001 -> CD3CN_gas (evaporation)
    kall(496) = 1.63361615d+12*exp(-4680.0*invTd)

    !CH2CCH_0001 -> CH2CCH_gas (evaporation)
    kall(497) = 1.45706222d+12*exp(-3300.0*invTd)

    !CH2CCD_0001 -> CH2CCD_gas (evaporation)
    kall(498) = 1.43873366d+12*exp(-3300.0*invTd)

    !CHDCCH_0001 -> CHDCCH_gas (evaporation)
    kall(499) = 1.43873366d+12*exp(-3300.0*invTd)

    !CHDCCD_0001 -> CHDCCD_gas (evaporation)
    kall(500) = 1.42107982d+12*exp(-3300.0*invTd)

    !CD2CCH_0001 -> CD2CCH_gas (evaporation)
    kall(501) = 1.42107982d+12*exp(-3300.0*invTd)

    !CD2CCD_0001 -> CD2CCD_gas (evaporation)
    kall(502) = 1.40406029d+12*exp(-3300.0*invTd)

    !CH2CHC2H_0001 -> CH2CHC2H_gas (evaporation)
    kall(503) = 1.56668840d+12*exp(-5087.0*invTd)

    !CH2CHCHCH2_0001 -> CH2CHCHCH2_gas (evaporation)
    kall(504) = 1.66786610d+12*exp(-5987.0*invTd)

    !CH2CHCN_0001 -> CH2CHCN_gas (evaporation)
    kall(505) = 1.61066705d+12*exp(-5480.0*invTd)

    !CH2NH_0001 -> CH2NH_gas (evaporation)
    kall(506) = 2.18813448d+12*exp(-5534.0*invTd)

    !CHDNH_0001 -> CHDNH_gas (evaporation)
    kall(507) = 2.15135649d+12*exp(-5534.0*invTd)

    !CHDND_0001 -> CHDND_gas (evaporation)
    kall(508) = 2.11637275d+12*exp(-5534.0*invTd)

    !CH2ND_0001 -> CH2ND_gas (evaporation)
    kall(509) = 2.15135649d+12*exp(-5534.0*invTd)

    !CD2NH_0001 -> CD2NH_gas (evaporation)
    kall(510) = 2.11637275d+12*exp(-5534.0*invTd)

    !CD2ND_0001 -> CD2ND_gas (evaporation)
    kall(511) = 2.08304197d+12*exp(-5534.0*invTd)

    !CH2NH2_0001 -> CH2NH2_gas (evaporation)
    kall(512) = 2.15135649d+12*exp(-5534.0*invTd)

    !CH2NHD_0001 -> CH2NHD_gas (evaporation)
    kall(513) = 2.11637275d+12*exp(-5534.0*invTd)

    !CHDNH2_0001 -> CHDNH2_gas (evaporation)
    kall(514) = 2.11637275d+12*exp(-5534.0*invTd)

    !CHDNHD_0001 -> CHDNHD_gas (evaporation)
    kall(515) = 2.08304197d+12*exp(-5534.0*invTd)

    !CHDND2_0001 -> CHDND2_gas (evaporation)
    kall(516) = 2.05123793d+12*exp(-5534.0*invTd)

    !CH2ND2_0001 -> CH2ND2_gas (evaporation)
    kall(517) = 2.08304197d+12*exp(-5534.0*invTd)

    !CD2NH2_0001 -> CD2NH2_gas (evaporation)
    kall(518) = 2.08304197d+12*exp(-5534.0*invTd)

    !CD2NHD_0001 -> CD2NHD_gas (evaporation)
    kall(519) = 2.05123793d+12*exp(-5534.0*invTd)

    !CD2ND2_0001 -> CD2ND2_gas (evaporation)
    kall(520) = 2.02084754d+12*exp(-5534.0*invTd)

    !CH3C3N_0001 -> CH3C3N_gas (evaporation)
    kall(521) = 1.58155484d+12*exp(-6480.0*invTd)

    !CH3C4H_0001 -> CH3C4H_gas (evaporation)
    kall(522) = 1.51918429d+12*exp(-5887.0*invTd)

    !CH3C5N_0001 -> CH3C5N_gas (evaporation)
    kall(523) = 1.49046416d+12*exp(-7880.0*invTd)

    !CH3C6H_0001 -> CH3C6H_gas (evaporation)
    kall(524) = 1.46105314d+12*exp(-7487.0*invTd)

    !CH3C7N_0001 -> CH3C7N_gas (evaporation)
    kall(525) = 1.45083618d+12*exp(-9480.0*invTd)

    !CH3CCH_0001 -> CH3CCH_gas (evaporation)
    kall(526) = 1.54388602d+12*exp(-3800.0*invTd)

    !CH3CH2OH_0001 -> CH3CH2OH_gas (evaporation)
    kall(527) = 1.71621382d+12*exp(-5400.0*invTd)

    !CH3CHCH2_0001 -> CH3CHCH2_gas (evaporation)
    kall(528) = 1.36084804d+12*exp(-3100.0*invTd)

    !CH3CHO_0001 -> CH3CHO_gas (evaporation)
    kall(529) = 1.75478524d+12*exp(-5400.0*invTd)

    !CH3CN_0001 -> CH3CN_gas (evaporation)
    kall(530) = 1.69232757d+12*exp(-4680.0*invTd)

    !CH2DCN_0001 -> CH2DCN_gas (evaporation)
    kall(531) = 1.67205944d+12*exp(-4680.0*invTd)

    !CHD2CN_0001 -> CHD2CN_gas (evaporation)
    kall(532) = 1.65250252d+12*exp(-4680.0*invTd)

    !CH3COCH3_0001 -> CH3COCH3_gas (evaporation)
    kall(533) = 1.23047688d+12*exp(-3500.0*invTd)

    !CH3NH2_0001 -> CH3NH2_gas (evaporation)
    kall(534) = 2.30843418d+12*exp(-6584.0*invTd)

    !CH3OCH2_0001 -> CH3OCH2_gas (evaporation)
    kall(535) = 1.39695108d+12*exp(-3500.0*invTd)

    !CH3OCH3_0001 -> CH3OCH3_gas (evaporation)
    kall(536) = 1.31077996d+12*exp(-3150.0*invTd)

    !CN_0001 -> CN_gas (evaporation)
    kall(537) = 1.64378788d+12*exp(-2800.0*invTd)

    !CS_0001 -> CS_gas (evaporation)
    kall(538) = 1.35083431d+12*exp(-3200.0*invTd)

    !H2CCN_0001 -> H2CCN_gas (evaporation)
    kall(539) = 1.62889698d+12*exp(-4230.0*invTd)

    !HDCCN_0001 -> HDCCN_gas (evaporation)
    kall(540) = 1.60890976d+12*exp(-4230.0*invTd)

    !D2CCN_0001 -> D2CCN_gas (evaporation)
    kall(541) = 1.58964069d+12*exp(-4230.0*invTd)

    !H2CCO_0001 -> H2CCO_gas (evaporation)
    kall(542) = 1.29332538d+12*exp(-2800.0*invTd)

    !HDCCO_0001 -> HDCCO_gas (evaporation)
    kall(543) = 1.27819825d+12*exp(-2800.0*invTd)

    !D2CCO_0001 -> D2CCO_gas (evaporation)
    kall(544) = 1.26358979d+12*exp(-2800.0*invTd)

    !H2CN_0001 -> H2CN_gas (evaporation)
    kall(545) = 1.46649314d+12*exp(-2400.0*invTd)

    !HDCN_0001 -> HDCN_gas (evaporation)
    kall(546) = 1.44098697d+12*exp(-2400.0*invTd)

    !D2CN_0001 -> D2CN_gas (evaporation)
    kall(547) = 1.41676697d+12*exp(-2400.0*invTd)

    !H2CS_0001 -> H2CS_gas (evaporation)
    kall(548) = 1.54917633d+12*exp(-4400.0*invTd)

    !HDCS_0001 -> HDCS_gas (evaporation)
    kall(549) = 1.53260713d+12*exp(-4400.0*invTd)

    !D2CS_0001 -> D2CS_gas (evaporation)
    kall(550) = 1.51655844d+12*exp(-4400.0*invTd)

    !H2S_0001 -> H2S_gas (evaporation)
    kall(551) = 1.42274437d+12*exp(-2743.0*invTd)

    !HDS_0001 -> HDS_gas (evaporation)
    kall(552) = 1.40227216d+12*exp(-2743.0*invTd)

    !D2S_0001 -> D2S_gas (evaporation)
    kall(553) = 1.38265900d+12*exp(-2743.0*invTd)

    !HC2O_0001 -> HC2O_gas (evaporation)
    kall(554) = 1.21190094d+12*exp(-2400.0*invTd)

    !DC2O_0001 -> DC2O_gas (evaporation)
    kall(555) = 1.19738664d+12*exp(-2400.0*invTd)

    !HC3N_0001 -> HC3N_gas (evaporation)
    kall(556) = 1.50107047d+12*exp(-4580.0*invTd)

    !DC3N_0001 -> DC3N_gas (evaporation)
    kall(557) = 1.48656704d+12*exp(-4580.0*invTd)

    !HC4N_0001 -> HC4N_gas (evaporation)
    kall(558) = 1.46377489d+12*exp(-5380.0*invTd)

    !DC4N_0001 -> DC4N_gas (evaporation)
    kall(559) = 1.45229413d+12*exp(-5380.0*invTd)

    !HC5N_0001 -> HC5N_gas (evaporation)
    kall(560) = 1.43786144d+12*exp(-6180.0*invTd)

    !HC6N_0001 -> HC6N_gas (evaporation)
    kall(561) = 1.49790270d+12*exp(-7780.0*invTd)

    !HC7N_0001 -> HC7N_gas (evaporation)
    kall(562) = 1.40418921d+12*exp(-7780.0*invTd)

    !HC8N_0001 -> HC8N_gas (evaporation)
    kall(563) = 1.45610724d+12*exp(-9380.0*invTd)

    !HC9N_0001 -> HC9N_gas (evaporation)
    kall(564) = 1.38325515d+12*exp(-9380.0*invTd)

    !HCCNC_0001 -> HCCNC_gas (evaporation)
    kall(565) = 1.50107047d+12*exp(-4580.0*invTd)

    !DCCNC_0001 -> DCCNC_gas (evaporation)
    kall(566) = 1.48656704d+12*exp(-4580.0*invTd)

    !HCN_0001 -> HCN_gas (evaporation)
    kall(567) = 1.85426761d+12*exp(-3700.0*invTd)

    !DCN_0001 -> DCN_gas (evaporation)
    kall(568) = 1.82085464d+12*exp(-3700.0*invTd)

    !HCNCC_0001 -> HCNCC_gas (evaporation)
    kall(569) = 1.50107047d+12*exp(-4580.0*invTd)

    !DCNCC_0001 -> DCNCC_gas (evaporation)
    kall(570) = 1.48656704d+12*exp(-4580.0*invTd)

    !HCOOCH3_0001 -> HCOOCH3_gas (evaporation)
    kall(571) = 1.62246626d+12*exp(-6295.0*invTd)

    !HCS_0001 -> HCS_gas (evaporation)
    kall(572) = 1.27158727d+12*exp(-2900.0*invTd)

    !DCS_0001 -> DCS_gas (evaporation)
    kall(573) = 1.25768973d+12*exp(-2900.0*invTd)

    !HNC_0001 -> HNC_gas (evaporation)
    kall(574) = 1.87915822d+12*exp(-3800.0*invTd)

    !DNC_0001 -> DNC_gas (evaporation)
    kall(575) = 1.84529674d+12*exp(-3800.0*invTd)

    !HNCCC_0001 -> HNCCC_gas (evaporation)
    kall(576) = 1.50107047d+12*exp(-4580.0*invTd)

    !DNCCC_0001 -> DNCCC_gas (evaporation)
    kall(577) = 1.48656704d+12*exp(-4580.0*invTd)

    !HNCO_0001 -> HNCO_gas (evaporation)
    kall(578) = 1.60230630d+12*exp(-4400.0*invTd)

    !DNCO_0001 -> DNCO_gas (evaporation)
    kall(579) = 1.58399363d+12*exp(-4400.0*invTd)

    !HNO_0001 -> HNO_gas (evaporation)
    kall(580) = 1.55823592d+12*exp(-3000.0*invTd)

    !DNO_0001 -> DNO_gas (evaporation)
    kall(581) = 1.53369524d+12*exp(-3000.0*invTd)

    !HS_0001 -> HS_gas (evaporation)
    kall(582) = 1.43277615d+12*exp(-2700.0*invTd)

    !DS_0001 -> DS_gas (evaporation)
    kall(583) = 1.41154866d+12*exp(-2700.0*invTd)

    !N2_0001 -> N2_gas (evaporation)
    kall(584) = 1.02392311d+12*exp(-1170.0*invTd)

    !N2O_0001 -> N2O_gas (evaporation)
    kall(585) = 1.16985683d+12*exp(-2400.0*invTd)

    !NC4N_0001 -> NC4N_gas (evaporation)
    kall(586) = 1.17752810d+12*exp(-4200.0*invTd)

    !NC6N_0001 -> NC6N_gas (evaporation)
    kall(587) = 1.20633361d+12*exp(-5800.0*invTd)

    !NC8N_0001 -> NC8N_gas (evaporation)
    kall(588) = 1.22365348d+12*exp(-7400.0*invTd)

    !NH2CHO_0001 -> NH2CHO_gas (evaporation)
    kall(589) = 1.87420654d+12*exp(-6300.0*invTd)

    !NH2CDO_0001 -> NH2CDO_gas (evaporation)
    kall(590) = 1.85372279d+12*exp(-6300.0*invTd)

    !NHDCHO_0001 -> NHDCHO_gas (evaporation)
    kall(591) = 1.85372279d+12*exp(-6300.0*invTd)

    !NHDCDO_0001 -> NHDCDO_gas (evaporation)
    kall(592) = 1.83389631d+12*exp(-6300.0*invTd)

    !ND2CHO_0001 -> ND2CHO_gas (evaporation)
    kall(593) = 1.83389631d+12*exp(-6300.0*invTd)

    !ND2CDO_0001 -> ND2CDO_gas (evaporation)
    kall(594) = 1.81469268d+12*exp(-6300.0*invTd)

    !NH2CN_0001 -> NH2CN_gas (evaporation)
    kall(595) = 1.82183862d+12*exp(-5556.0*invTd)

    !NHDCN_0001 -> NHDCN_gas (evaporation)
    kall(596) = 1.80052983d+12*exp(-5556.0*invTd)

    !ND2CN_0001 -> ND2CN_gas (evaporation)
    kall(597) = 1.77995167d+12*exp(-5556.0*invTd)

    !HSS_0001 -> HSS_gas (evaporation)
    kall(598) = 1.01139241d+12*exp(-2650.0*invTd)

    !DSS_0001 -> DSS_gas (evaporation)
    kall(599) = 1.00370110d+12*exp(-2650.0*invTd)

    !HSSH_0001 -> HSSH_gas (evaporation)
    kall(600) = 1.08558119d+12*exp(-3100.0*invTd)

    !HSSD_0001 -> HSSD_gas (evaporation)
    kall(601) = 1.07744938d+12*exp(-3100.0*invTd)

    !DSSH_0001 -> DSSH_gas (evaporation)
    kall(602) = 1.07744938d+12*exp(-3100.0*invTd)

    !DSSD_0001 -> DSSD_gas (evaporation)
    kall(603) = 1.06949762d+12*exp(-3100.0*invTd)

    !NO_0001 -> NO_gas (evaporation)
    kall(604) = 1.15678539d+12*exp(-1600.0*invTd)

    !NO2_0001 -> NO2_gas (evaporation)
    kall(605) = 1.14414255d+12*exp(-2400.0*invTd)

    !NS_0001 -> NS_gas (evaporation)
    kall(606) = 1.01800829d+12*exp(-1900.0*invTd)

    !OCN_0001 -> OCN_gas (evaporation)
    kall(607) = 1.19738664d+12*exp(-2400.0*invTd)

    !OCS_0001 -> OCS_gas (evaporation)
    kall(608) = 1.09894612d+12*exp(-2888.0*invTd)

    !S2_0001 -> S2_gas (evaporation)
    kall(609) = 9.28698586d+11*exp(-2200.0*invTd)

    !SO_0001 -> SO_gas (evaporation)
    kall(610) = 1.20979512d+12*exp(-2800.0*invTd)

    !SO2_0001 -> SO2_gas (evaporation)
    kall(611) = 1.15537244d+12*exp(-3405.0*invTd)

    !CH3OCHO_0001 -> CH3OCHO_gas (evaporation)
    kall(612) = 1.62246626d+12*exp(-6295.0*invTd)

    !Si_0001 -> Si_gas (evaporation)
    kall(613) = 3.22406219d+12*exp(-11600.0*invTd)

    !SiS_0001 -> SiS_gas (evaporation)
    kall(614) = 1.26057765d+12*exp(-3800.0*invTd)

    !SiN_0001 -> SiN_gas (evaporation)
    kall(615) = 1.44598174d+12*exp(-3500.0*invTd)

    !SiC_0001 -> SiC_gas (evaporation)
    kall(616) = 1.48169037d+12*exp(-3500.0*invTd)

    !SiH_0001 -> SiH_gas (evaporation)
    kall(617) = 3.35371465d+12*exp(-13000.0*invTd)

    !SiH2_0001 -> SiH2_gas (evaporation)
    kall(618) = 1.73517809d+12*exp(-3600.0*invTd)

    !SiH3_0001 -> SiH3_gas (evaporation)
    kall(619) = 1.81050653d+12*exp(-4050.0*invTd)

    !SiH4_0001 -> SiH4_gas (evaporation)
    kall(620) = 1.87838538d+12*exp(-4500.0*invTd)

    !SiC2CH3_0001 -> SiC2CH3_gas (evaporation)
    kall(621) = 1.43188484d+12*exp(-5475.0*invTd)

    !SiC3H_0001 -> SiC3H_gas (evaporation)
    kall(622) = 1.46367025d+12*exp(-5550.0*invTd)

    !SiC3H5_0001 -> SiC3H5_gas (evaporation)
    kall(623) = 1.42699647d+12*exp(-5600.0*invTd)

    !SiC4_0001 -> SiC4_gas (evaporation)
    kall(624) = 1.39563753d+12*exp(-5900.0*invTd)

    !SiC4H_0001 -> SiC4H_gas (evaporation)
    kall(625) = 1.44827096d+12*exp(-6437.0*invTd)

    !SiC6H_0001 -> SiC6H_gas (evaporation)
    kall(626) = 1.41299208d+12*exp(-8037.0*invTd)

    !SiC8H_0001 -> SiC8H_gas (evaporation)
    kall(627) = 1.06021132d+12*exp(-5600.0*invTd)

    !c_HCCHSi_0001 -> c_HCCHSi_gas (evaporation)
    kall(628) = 1.55438366d+12*exp(-5200.0*invTd)

    !c_SiC2_0001 -> c_SiC2_gas (evaporation)
    kall(629) = 1.44040953d+12*exp(-4300.0*invTd)

    !l_C3H_0001 -> l_C3H_gas (evaporation)
    kall(630) = 1.64695815d+12*exp(-4000.0*invTd)

    !l_C3D_0001 -> l_C3D_gas (evaporation)
    kall(631) = 1.62514318d+12*exp(-4000.0*invTd)

    !c_C3H_0001 -> c_C3H_gas (evaporation)
    kall(632) = 1.64695815d+12*exp(-4000.0*invTd)

    !c_C3D_0001 -> c_C3D_gas (evaporation)
    kall(633) = 1.62514318d+12*exp(-4000.0*invTd)

    !l_SiC3_0001 -> l_SiC3_gas (evaporation)
    kall(634) = 1.41399714d+12*exp(-5100.0*invTd)

    !l_C3H2_0001 -> l_C3H2_gas (evaporation)
    kall(635) = 1.49544082d+12*exp(-3387.0*invTd)

    !l_C3HD_0001 -> l_C3HD_gas (evaporation)
    kall(636) = 1.47614400d+12*exp(-3387.0*invTd)

    !l_C3D2_0001 -> l_C3D2_gas (evaporation)
    kall(637) = 1.45757541d+12*exp(-3387.0*invTd)

    !c_C3H2_0001 -> c_C3H2_gas (evaporation)
    kall(638) = 1.49544082d+12*exp(-3387.0*invTd)

    !c_C3HD_0001 -> c_C3HD_gas (evaporation)
    kall(639) = 1.47614400d+12*exp(-3387.0*invTd)

    !c_C3D2_0001 -> c_C3D2_gas (evaporation)
    kall(640) = 1.45757541d+12*exp(-3387.0*invTd)

    !Mg_0001 -> Mg_gas (evaporation)
    kall(641) = 2.35388773d+12*exp(-5300.0*invTd)

    !MgH_0001 -> MgH_gas (evaporation)
    kall(642) = 2.40224528d+12*exp(-5750.0*invTd)

    !MgH2_0001 -> MgH2_gas (evaporation)
    kall(643) = 2.44603463d+12*exp(-6200.0*invTd)

    !Na_0001 -> Na_gas (evaporation)
    kall(644) = 3.58781893d+12*exp(-11800.0*invTd)

    !NaH_0001 -> NaH_gas (evaporation)
    kall(645) = 3.57862222d+12*exp(-12250.0*invTd)

    !F_0001 -> F_gas (evaporation)
    kall(646) = 1.02783079d+12*exp(-800.0*invTd)

    !HF_0001 -> HF_gas (evaporation)
    kall(647) = 3.06739048d+12*exp(-7500.0*invTd)

    !DF_0001 -> DF_gas (evaporation)
    kall(648) = 2.99346659d+12*exp(-7500.0*invTd)

    !MgD_0001 -> MgD_gas (evaporation)
    kall(649) = 2.35559530d+12*exp(-5750.0*invTd)

    !MgHD_0001 -> MgHD_gas (evaporation)
    kall(650) = 2.40031033d+12*exp(-6200.0*invTd)

    !MgD2_0001 -> MgD2_gas (evaporation)
    kall(651) = 2.35705795d+12*exp(-6200.0*invTd)

    !NaD_0001 -> NaD_gas (evaporation)
    kall(652) = 3.50631937d+12*exp(-12250.0*invTd)

    !SiD_0001 -> SiD_gas (evaporation)
    kall(653) = 3.29734569d+12*exp(-13000.0*invTd)

    !Cl_0001 -> Cl_gas (evaporation)
    kall(654) = 1.46649314d+12*exp(-3000.0*invTd)

    !HCl_0001 -> HCl_gas (evaporation)
    kall(655) = 1.89895815d+12*exp(-5174.0*invTd)

    !DCl_0001 -> DCl_gas (evaporation)
    kall(656) = 1.87312078d+12*exp(-5174.0*invTd)

    !H2Cl_0001 -> H2Cl_gas (evaporation)
    kall(657) = 1.94870717d+12*exp(-5600.0*invTd)

    !HDCl_0001 -> HDCl_gas (evaporation)
    kall(658) = 1.92289534d+12*exp(-5600.0*invTd)

    !D2Cl_0001 -> D2Cl_gas (evaporation)
    kall(659) = 1.89808274d+12*exp(-5600.0*invTd)

    !CCl_0001 -> CCl_gas (evaporation)
    kall(660) = 1.00712019d+12*exp(-1900.0*invTd)

    !ClO_0001 -> ClO_gas (evaporation)
    kall(661) = 9.66818926d+11*exp(-1900.0*invTd)

    !C3H4_0001 -> C3H4_gas (evaporation)
    kall(662) = 1.54388602d+12*exp(-3800.0*invTd)

    !FeD_0001 -> FeD_gas (evaporation)
    kall(663) = 1.41829284d+12*exp(-4650.0*invTd)

    !P_0001 -> P_gas (evaporation)
    kall(664) = 9.43558707d+11*exp(-1100.0*invTd)

    !PO_0001 -> PO_gas (evaporation)
    kall(665) = 1.00712019d+12*exp(-1900.0*invTd)

    !PH_0001 -> PH_gas (evaporation)
    kall(666) = 1.10241291d+12*exp(-1550.0*invTd)

    !PD_0001 -> PD_gas (evaporation)
    kall(667) = 1.08558119d+12*exp(-1550.0*invTd)

    !PH2_0001 -> PH2_gas (evaporation)
    kall(668) = 1.23313737d+12*exp(-2000.0*invTd)

    !PHD_0001 -> PHD_gas (evaporation)
    kall(669) = 1.21486766d+12*exp(-2000.0*invTd)

    !PD2_0001 -> PD2_gas (evaporation)
    kall(670) = 1.19738664d+12*exp(-2000.0*invTd)

    !PN_0001 -> PN_gas (evaporation)
    kall(671) = 1.02925734d+12*exp(-1900.0*invTd)

    !CP_0001 -> CP_gas (evaporation)
    kall(672) = 1.05292152d+12*exp(-1900.0*invTd)

    !CCP_0001 -> CCP_gas (evaporation)
    kall(673) = 1.40057482d+12*exp(-4300.0*invTd)

    !C3P_0001 -> C3P_gas (evaporation)
    kall(674) = 1.48642170d+12*exp(-5900.0*invTd)

    !C4P_0001 -> C4P_gas (evaporation)
    kall(675) = 1.54337165d+12*exp(-7500.0*invTd)

    !CH2PH_0001 -> CH2PH_gas (evaporation)
    kall(676) = 1.19086132d+12*exp(-2600.0*invTd)

    !CH2PD_0001 -> CH2PD_gas (evaporation)
    kall(677) = 1.17812447d+12*exp(-2600.0*invTd)

    !CHDPD_0001 -> CHDPD_gas (evaporation)
    kall(678) = 1.16578775d+12*exp(-2600.0*invTd)

    !HCP_0001 -> HCP_gas (evaporation)
    kall(679) = 1.15760668d+12*exp(-2350.0*invTd)

    !DCP_0001 -> DCP_gas (evaporation)
    kall(680) = 1.14467212d+12*exp(-2350.0*invTd)

    !HCCP_0001 -> HCCP_gas (evaporation)
    kall(681) = 1.45883516d+12*exp(-4750.0*invTd)

    !DCCP_0001 -> DCCP_gas (evaporation)
    kall(682) = 1.44598174d+12*exp(-4750.0*invTd)

    !HPO_0001 -> HPO_gas (evaporation)
    kall(683) = 1.10832401d+12*exp(-2350.0*invTd)

    !DPO_0001 -> DPO_gas (evaporation)
    kall(684) = 1.09695629d+12*exp(-2350.0*invTd)

    !CH2_0001 + N_0001 -> H_0001 + HCN_0001 (2body_dust_dust)
    kall(685) = 1.58399363d+12*max(3.93742211d-08, exp(-2.50000000d+02 * invTd))&
        *(1.58399363d+12*max(3.93742211d-08, exp(-2.50000000d+02 * invTd)) + (exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-3.60000000d+02*invTd)*1.13594070d+12))**(-1d0) &
        *((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-3.60000000d+02*invTd)*1.13594070d+12))*indns

    !N_0001 + CHD_0001 -> D_0001 + HCN_0001 (2body_dust_dust)
    kall(686) = 1.53028323d+12*max(2.94188189d-08, exp(-2.50000000d+02 * invTd))&
        *(1.53028323d+12*max(2.94188189d-08, exp(-2.50000000d+02 * invTd)) + (exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))**(-1d0) &
        *((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !N_0001 + CHD_0001 -> H_0001 + DCN_0001 (2body_dust_dust)
    kall(687) = 1.53028323d+12*max(2.94188189d-08, exp(-2.50000000d+02 * invTd))&
        *(1.53028323d+12*max(2.94188189d-08, exp(-2.50000000d+02 * invTd)) + (exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))**(-1d0) &
        *((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !N_0001 + CD2_0001 -> D_0001 + DCN_0001 (2body_dust_dust)
    kall(688) = 1.48169037d+12*max(2.25096018d-08, exp(-2.50000000d+02 * invTd))&
        *(1.48169037d+12*max(2.25096018d-08, exp(-2.50000000d+02 * invTd)) + (exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))**(-1d0) &
        *((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !OH_0001 + CH2NH_0001 -> H2O_0001 + H2CN_0001 (2body_dust_dust)
    kall(689) = 2.60560084d+12*max(8.17647525d-15, exp(-5.91000000d+02 * invTd))&
        *(2.60560084d+12*max(8.17647525d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.18813448d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.18813448d+12))*indns

    !OD_0001 + CH2NH_0001 -> HDO_0001 + H2CN_0001 (2body_dust_dust)
    kall(690) = 2.53218886d+12*max(4.56236913d-15, exp(-5.91000000d+02 * invTd))&
        *(2.53218886d+12*max(4.56236913d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.18813448d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.18813448d+12))*indns

    !OH_0001 + CHDNH_0001 -> HDO_0001 + H2CN_0001 (2body_dust_dust)
    kall(691) = 2.60560084d+12*max(6.68321908d-15, exp(-5.91000000d+02 * invTd))&
        *(2.60560084d+12*max(6.68321908d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !OH_0001 + CHDNH_0001 -> H2O_0001 + HDCN_0001 (2body_dust_dust)
    kall(692) = 2.60560084d+12*max(6.68321908d-15, exp(-5.91000000d+02 * invTd))&
        *(2.60560084d+12*max(6.68321908d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !OD_0001 + CHDNH_0001 -> H2O_0001 + D2CN_0001 (2body_dust_dust)
    kall(693) = 2.53218886d+12*max(3.68780110d-15, exp(-5.91000000d+02 * invTd))&
        *(2.53218886d+12*max(3.68780110d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !OD_0001 + CHDNH_0001 -> HDO_0001 + HDCN_0001 (2body_dust_dust)
    kall(694) = 2.53218886d+12*max(3.68780110d-15, exp(-5.91000000d+02 * invTd))&
        *(2.53218886d+12*max(3.68780110d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !OD_0001 + CHDNH_0001 -> D2O_0001 + H2CN_0001 (2body_dust_dust)
    kall(695) = 2.53218886d+12*max(3.68780110d-15, exp(-5.91000000d+02 * invTd))&
        *(2.53218886d+12*max(3.68780110d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !OH_0001 + CH2ND_0001 -> HDO_0001 + H2CN_0001 (2body_dust_dust)
    kall(696) = 2.60560084d+12*max(6.68321908d-15, exp(-5.91000000d+02 * invTd))&
        *(2.60560084d+12*max(6.68321908d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !OH_0001 + CH2ND_0001 -> H2O_0001 + HDCN_0001 (2body_dust_dust)
    kall(697) = 2.60560084d+12*max(6.68321908d-15, exp(-5.91000000d+02 * invTd))&
        *(2.60560084d+12*max(6.68321908d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !OD_0001 + CH2ND_0001 -> D2O_0001 + H2CN_0001 (2body_dust_dust)
    kall(698) = 2.53218886d+12*max(3.68780110d-15, exp(-5.91000000d+02 * invTd))&
        *(2.53218886d+12*max(3.68780110d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !OD_0001 + CH2ND_0001 -> HDO_0001 + HDCN_0001 (2body_dust_dust)
    kall(699) = 2.53218886d+12*max(3.68780110d-15, exp(-5.91000000d+02 * invTd))&
        *(2.53218886d+12*max(3.68780110d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !OD_0001 + CH2ND_0001 -> H2O_0001 + D2CN_0001 (2body_dust_dust)
    kall(700) = 2.53218886d+12*max(3.68780110d-15, exp(-5.91000000d+02 * invTd))&
        *(2.53218886d+12*max(3.68780110d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !OH_0001 + CD2NH_0001 -> H2O_0001 + D2CN_0001 (2body_dust_dust)
    kall(701) = 2.60560084d+12*max(5.51517427d-15, exp(-5.91000000d+02 * invTd))&
        *(2.60560084d+12*max(5.51517427d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !OH_0001 + CD2NH_0001 -> HDO_0001 + HDCN_0001 (2body_dust_dust)
    kall(702) = 2.60560084d+12*max(5.51517427d-15, exp(-5.91000000d+02 * invTd))&
        *(2.60560084d+12*max(5.51517427d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !OH_0001 + CD2NH_0001 -> D2O_0001 + H2CN_0001 (2body_dust_dust)
    kall(703) = 2.60560084d+12*max(5.51517427d-15, exp(-5.91000000d+02 * invTd))&
        *(2.60560084d+12*max(5.51517427d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !OD_0001 + CD2NH_0001 -> HDO_0001 + D2CN_0001 (2body_dust_dust)
    kall(704) = 2.53218886d+12*max(3.01071486d-15, exp(-5.91000000d+02 * invTd))&
        *(2.53218886d+12*max(3.01071486d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !OD_0001 + CD2NH_0001 -> D2O_0001 + HDCN_0001 (2body_dust_dust)
    kall(705) = 2.53218886d+12*max(3.01071486d-15, exp(-5.91000000d+02 * invTd))&
        *(2.53218886d+12*max(3.01071486d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !OH_0001 + CHDND_0001 -> H2O_0001 + D2CN_0001 (2body_dust_dust)
    kall(706) = 2.60560084d+12*max(5.51517427d-15, exp(-5.91000000d+02 * invTd))&
        *(2.60560084d+12*max(5.51517427d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !OH_0001 + CHDND_0001 -> HDO_0001 + HDCN_0001 (2body_dust_dust)
    kall(707) = 2.60560084d+12*max(5.51517427d-15, exp(-5.91000000d+02 * invTd))&
        *(2.60560084d+12*max(5.51517427d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !OH_0001 + CHDND_0001 -> D2O_0001 + H2CN_0001 (2body_dust_dust)
    kall(708) = 2.60560084d+12*max(5.51517427d-15, exp(-5.91000000d+02 * invTd))&
        *(2.60560084d+12*max(5.51517427d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !OD_0001 + CHDND_0001 -> HDO_0001 + D2CN_0001 (2body_dust_dust)
    kall(709) = 2.53218886d+12*max(3.01071486d-15, exp(-5.91000000d+02 * invTd))&
        *(2.53218886d+12*max(3.01071486d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !OD_0001 + CHDND_0001 -> D2O_0001 + HDCN_0001 (2body_dust_dust)
    kall(710) = 2.53218886d+12*max(3.01071486d-15, exp(-5.91000000d+02 * invTd))&
        *(2.53218886d+12*max(3.01071486d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !OH_0001 + CD2ND_0001 -> H2O_0001 + D2CN_0001 (2body_dust_dust)
    kall(711) = 2.60560084d+12*max(4.59191714d-15, exp(-5.91000000d+02 * invTd))&
        *(2.60560084d+12*max(4.59191714d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.08304197d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.08304197d+12))*indns

    !OH_0001 + CD2ND_0001 -> HDO_0001 + HDCN_0001 (2body_dust_dust)
    kall(712) = 2.60560084d+12*max(4.59191714d-15, exp(-5.91000000d+02 * invTd))&
        *(2.60560084d+12*max(4.59191714d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.08304197d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.08304197d+12))*indns

    !OH_0001 + CD2ND_0001 -> D2O_0001 + H2CN_0001 (2body_dust_dust)
    kall(713) = 2.60560084d+12*max(4.59191714d-15, exp(-5.91000000d+02 * invTd))&
        *(2.60560084d+12*max(4.59191714d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.08304197d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.76700000d+03*invTd)*2.08304197d+12))*indns

    !OD_0001 + CD2ND_0001 -> D2O_0001 + D2CN_0001 (2body_dust_dust)
    kall(714) = 2.53218886d+12*max(2.48082201d-15, exp(-5.91000000d+02 * invTd))&
        *(2.53218886d+12*max(2.48082201d-15, exp(-5.91000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.08304197d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.76700000d+03*invTd)*2.08304197d+12))*indns

    !OH_0001 + CH3OH_0001 -> H2O_0001 + CH2OH_0001 (2body_dust_dust)
    kall(715) = 2.60560084d+12*max(6.68553483d-12, exp(-3.59000000d+02 * invTd))&
        *(2.60560084d+12*max(6.68553483d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))*indns

    !OD_0001 + CH3OH_0001 -> HDO_0001 + CH2OH_0001 (2body_dust_dust)
    kall(716) = 2.53218886d+12*max(4.13741547d-12, exp(-3.59000000d+02 * invTd))&
        *(2.53218886d+12*max(4.13741547d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))*indns

    !OH_0001 + CH2DOH_0001 -> HDO_0001 + CH2OH_0001 (2body_dust_dust)
    kall(717) = 2.60560084d+12*max(5.83346439d-12, exp(-3.59000000d+02 * invTd))&
        *(2.60560084d+12*max(5.83346439d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !OH_0001 + CH2DOH_0001 -> H2O_0001 + CHDOH_0001 (2body_dust_dust)
    kall(718) = 2.60560084d+12*max(5.83346439d-12, exp(-3.59000000d+02 * invTd))&
        *(2.60560084d+12*max(5.83346439d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !OD_0001 + CH2DOH_0001 -> D2O_0001 + CH2OH_0001 (2body_dust_dust)
    kall(719) = 2.53218886d+12*max(3.58200253d-12, exp(-3.59000000d+02 * invTd))&
        *(2.53218886d+12*max(3.58200253d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !OD_0001 + CH2DOH_0001 -> HDO_0001 + CHDOH_0001 (2body_dust_dust)
    kall(720) = 2.53218886d+12*max(3.58200253d-12, exp(-3.59000000d+02 * invTd))&
        *(2.53218886d+12*max(3.58200253d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !OH_0001 + CH3OD_0001 -> H2O_0001 + CH2OD_0001 (2body_dust_dust)
    kall(721) = 2.60560084d+12*max(5.83346439d-12, exp(-3.59000000d+02 * invTd))&
        *(2.60560084d+12*max(5.83346439d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !OD_0001 + CH3OD_0001 -> HDO_0001 + CH2OD_0001 (2body_dust_dust)
    kall(722) = 2.53218886d+12*max(3.58200253d-12, exp(-3.59000000d+02 * invTd))&
        *(2.53218886d+12*max(3.58200253d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !OH_0001 + CHD2OH_0001 -> HDO_0001 + CHDOH_0001 (2body_dust_dust)
    kall(723) = 2.60560084d+12*max(5.12072395d-12, exp(-3.59000000d+02 * invTd))&
        *(2.60560084d+12*max(5.12072395d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !OH_0001 + CHD2OH_0001 -> D2O_0001 + CD2OH_0001 (2body_dust_dust)
    kall(724) = 2.60560084d+12*max(5.12072395d-12, exp(-3.59000000d+02 * invTd))&
        *(2.60560084d+12*max(5.12072395d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !OD_0001 + CHD2OH_0001 -> D2O_0001 + CHDOH_0001 (2body_dust_dust)
    kall(725) = 2.53218886d+12*max(3.12069758d-12, exp(-3.59000000d+02 * invTd))&
        *(2.53218886d+12*max(3.12069758d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !OD_0001 + CHD2OH_0001 -> HDO_0001 + CD2OH_0001 (2body_dust_dust)
    kall(726) = 2.53218886d+12*max(3.12069758d-12, exp(-3.59000000d+02 * invTd))&
        *(2.53218886d+12*max(3.12069758d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !OH_0001 + CH2DOD_0001 -> HDO_0001 + CH2OD_0001 (2body_dust_dust)
    kall(727) = 2.60560084d+12*max(5.12072395d-12, exp(-3.59000000d+02 * invTd))&
        *(2.60560084d+12*max(5.12072395d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !OH_0001 + CH2DOD_0001 -> H2O_0001 + CHDOD_0001 (2body_dust_dust)
    kall(728) = 2.60560084d+12*max(5.12072395d-12, exp(-3.59000000d+02 * invTd))&
        *(2.60560084d+12*max(5.12072395d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !OD_0001 + CH2DOD_0001 -> HDO_0001 + CHDOD_0001 (2body_dust_dust)
    kall(729) = 2.53218886d+12*max(3.12069758d-12, exp(-3.59000000d+02 * invTd))&
        *(2.53218886d+12*max(3.12069758d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !OD_0001 + CH2DOD_0001 -> D2O_0001 + CH2OD_0001 (2body_dust_dust)
    kall(730) = 2.53218886d+12*max(3.12069758d-12, exp(-3.59000000d+02 * invTd))&
        *(2.53218886d+12*max(3.12069758d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !OH_0001 + CHD2OD_0001 -> HDO_0001 + CHDOD_0001 (2body_dust_dust)
    kall(731) = 2.60560084d+12*max(4.52042394d-12, exp(-3.59000000d+02 * invTd))&
        *(2.60560084d+12*max(4.52042394d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))*indns

    !OD_0001 + CHD2OD_0001 -> D2O_0001 + CHDOD_0001 (2body_dust_dust)
    kall(732) = 2.53218886d+12*max(2.73482526d-12, exp(-3.59000000d+02 * invTd))&
        *(2.53218886d+12*max(2.73482526d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))*indns

    !OH_0001 + CHD2OD_0001 -> H2O_0001 + CD2OD_0001 (2body_dust_dust)
    kall(733) = 2.60560084d+12*max(4.52042394d-12, exp(-3.59000000d+02 * invTd))&
        *(2.60560084d+12*max(4.52042394d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))*indns

    !OD_0001 + CHD2OD_0001 -> HDO_0001 + CD2OD_0001 (2body_dust_dust)
    kall(734) = 2.53218886d+12*max(2.73482526d-12, exp(-3.59000000d+02 * invTd))&
        *(2.53218886d+12*max(2.73482526d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))*indns

    !OH_0001 + CD3OH_0001 -> HDO_0001 + CD2OH_0001 (2body_dust_dust)
    kall(735) = 2.60560084d+12*max(4.52042394d-12, exp(-3.59000000d+02 * invTd))&
        *(2.60560084d+12*max(4.52042394d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))*indns

    !OD_0001 + CD3OH_0001 -> D2O_0001 + CD2OH_0001 (2body_dust_dust)
    kall(736) = 2.53218886d+12*max(2.73482526d-12, exp(-3.59000000d+02 * invTd))&
        *(2.53218886d+12*max(2.73482526d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))*indns

    !OH_0001 + CD3OD_0001 -> HDO_0001 + CD2OD_0001 (2body_dust_dust)
    kall(737) = 2.60560084d+12*max(4.01156072d-12, exp(-3.59000000d+02 * invTd))&
        *(2.60560084d+12*max(4.01156072d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.86675440d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.86675440d+12))*indns

    !OD_0001 + CD3OD_0001 -> D2O_0001 + CD2OD_0001 (2body_dust_dust)
    kall(738) = 2.53218886d+12*max(2.40989334d-12, exp(-3.59000000d+02 * invTd))&
        *(2.53218886d+12*max(2.40989334d-12, exp(-3.59000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.86675440d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.86675440d+12))*indns

    !OH_0001 + CH3OH_0001 -> H2O_0001 + CH3O_0001 (2body_dust_dust)
    kall(739) = 2.60560084d+12*max(6.09103860d-18, exp(-8.52000000d+02 * invTd))&
        *(2.60560084d+12*max(6.09103860d-18, exp(-8.52000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))*indns

    !OD_0001 + CH3OH_0001 -> HDO_0001 + CH3O_0001 (2body_dust_dust)
    kall(740) = 2.53218886d+12*max(2.90825376d-18, exp(-8.52000000d+02 * invTd))&
        *(2.53218886d+12*max(2.90825376d-18, exp(-8.52000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))*indns

    !OH_0001 + CH3OD_0001 -> HDO_0001 + CH3O_0001 (2body_dust_dust)
    kall(741) = 2.60560084d+12*max(4.93715399d-18, exp(-8.52000000d+02 * invTd))&
        *(2.60560084d+12*max(4.93715399d-18, exp(-8.52000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !OD_0001 + CH3OD_0001 -> D2O_0001 + CH3O_0001 (2body_dust_dust)
    kall(742) = 2.53218886d+12*max(2.32910778d-18, exp(-8.52000000d+02 * invTd))&
        *(2.53218886d+12*max(2.32910778d-18, exp(-8.52000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !OH_0001 + CH2DOH_0001 -> H2O_0001 + CH2DO_0001 (2body_dust_dust)
    kall(743) = 2.60560084d+12*max(4.93715399d-18, exp(-8.52000000d+02 * invTd))&
        *(2.60560084d+12*max(4.93715399d-18, exp(-8.52000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !OD_0001 + CH2DOH_0001 -> HDO_0001 + CH2DO_0001 (2body_dust_dust)
    kall(744) = 2.53218886d+12*max(2.32910778d-18, exp(-8.52000000d+02 * invTd))&
        *(2.53218886d+12*max(2.32910778d-18, exp(-8.52000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !OD_0001 + CHD2OD_0001 -> HDO_0001 + CHD2O_0001 (2body_dust_dust)
    kall(745) = 2.53218886d+12*max(1.53689486d-18, exp(-8.52000000d+02 * invTd))&
        *(2.53218886d+12*max(1.53689486d-18, exp(-8.52000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))*indns

    !OD_0001 + CHD2OH_0001 -> HDO_0001 + CHD2O_0001 (2body_dust_dust)
    kall(746) = 2.53218886d+12*max(1.88343685d-18, exp(-8.52000000d+02 * invTd))&
        *(2.53218886d+12*max(1.88343685d-18, exp(-8.52000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !OH_0001 + CHD2OH_0001 -> H2O_0001 + CHD2O_0001 (2body_dust_dust)
    kall(747) = 2.60560084d+12*max(4.03914651d-18, exp(-8.52000000d+02 * invTd))&
        *(2.60560084d+12*max(4.03914651d-18, exp(-8.52000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !OH_0001 + CD3OH_0001 -> H2O_0001 + CD3O_0001 (2body_dust_dust)
    kall(748) = 2.60560084d+12*max(3.33323594d-18, exp(-8.52000000d+02 * invTd))&
        *(2.60560084d+12*max(3.33323594d-18, exp(-8.52000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))*indns

    !OD_0001 + CD3OH_0001 -> HDO_0001 + CD3O_0001 (2body_dust_dust)
    kall(749) = 2.53218886d+12*max(1.53689486d-18, exp(-8.52000000d+02 * invTd))&
        *(2.53218886d+12*max(1.53689486d-18, exp(-8.52000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))*indns

    !OH_0001 + CD3OD_0001 -> HDO_0001 + CD3O_0001 (2body_dust_dust)
    kall(750) = 2.60560084d+12*max(2.77309532d-18, exp(-8.52000000d+02 * invTd))&
        *(2.60560084d+12*max(2.77309532d-18, exp(-8.52000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.86675440d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.50000000d+03*invTd)*1.86675440d+12))*indns

    !OD_0001 + CD3OD_0001 -> D2O_0001 + CD3O_0001 (2body_dust_dust)
    kall(751) = 2.53218886d+12*max(1.26479361d-18, exp(-8.52000000d+02 * invTd))&
        *(2.53218886d+12*max(1.26479361d-18, exp(-8.52000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.86675440d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.50000000d+03*invTd)*1.86675440d+12))*indns

    !OH_0001 + CH3CHO_0001 -> H2O_0001 + CH3CO_0001 (2body_dust_dust)
    kall(752) = 2.60560084d+12*max(6.56276226d-16, exp(-6.00000000d+02 * invTd))&
        *(2.60560084d+12*max(6.56276226d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.70000000d+03*invTd)*1.75478524d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.70000000d+03*invTd)*1.75478524d+12))*indns

    !OD_0001 + CH3CHO_0001 -> HDO_0001 + CH3CO_0001 (2body_dust_dust)
    kall(753) = 2.53218886d+12*max(3.18724833d-16, exp(-6.00000000d+02 * invTd))&
        *(2.53218886d+12*max(3.18724833d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.70000000d+03*invTd)*1.75478524d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.70000000d+03*invTd)*1.75478524d+12))*indns

    !OH_0001 + CH3CDO_0001 -> HDO_0001 + CH3CO_0001 (2body_dust_dust)
    kall(754) = 2.60560084d+12*max(5.88643944d-16, exp(-6.00000000d+02 * invTd))&
        *(2.60560084d+12*max(5.88643944d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))*indns

    !OH_0001 + CH3CDO_0001 -> H2O_0001 + CH2DCO_0001 (2body_dust_dust)
    kall(755) = 2.60560084d+12*max(5.88643944d-16, exp(-6.00000000d+02 * invTd))&
        *(2.60560084d+12*max(5.88643944d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))*indns

    !OD_0001 + CH3CDO_0001 -> HDO_0001 + CH2DCO_0001 (2body_dust_dust)
    kall(756) = 2.53218886d+12*max(2.83911988d-16, exp(-6.00000000d+02 * invTd))&
        *(2.53218886d+12*max(2.83911988d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))*indns

    !OD_0001 + CH3CDO_0001 -> D2O_0001 + CH3CO_0001 (2body_dust_dust)
    kall(757) = 2.53218886d+12*max(2.83911988d-16, exp(-6.00000000d+02 * invTd))&
        *(2.53218886d+12*max(2.83911988d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))*indns

    !OH_0001 + CH2DCHO_0001 -> HDO_0001 + CH3CO_0001 (2body_dust_dust)
    kall(758) = 2.60560084d+12*max(5.88643944d-16, exp(-6.00000000d+02 * invTd))&
        *(2.60560084d+12*max(5.88643944d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))*indns

    !OH_0001 + CH2DCHO_0001 -> H2O_0001 + CH2DCO_0001 (2body_dust_dust)
    kall(759) = 2.60560084d+12*max(5.88643944d-16, exp(-6.00000000d+02 * invTd))&
        *(2.60560084d+12*max(5.88643944d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))*indns

    !OD_0001 + CH2DCHO_0001 -> HDO_0001 + CH2DCO_0001 (2body_dust_dust)
    kall(760) = 2.53218886d+12*max(2.83911988d-16, exp(-6.00000000d+02 * invTd))&
        *(2.53218886d+12*max(2.83911988d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))*indns

    !OD_0001 + CH2DCHO_0001 -> D2O_0001 + CH3CO_0001 (2body_dust_dust)
    kall(761) = 2.53218886d+12*max(2.83911988d-16, exp(-6.00000000d+02 * invTd))&
        *(2.53218886d+12*max(2.83911988d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.80000000d+03*invTd)*1.76701887d+12))*indns

    !OH_0001 + CH2DCDO_0001 -> HDO_0001 + CH2DCO_0001 (2body_dust_dust)
    kall(762) = 2.60560084d+12*max(5.29977399d-16, exp(-6.00000000d+02 * invTd))&
        *(2.60560084d+12*max(5.29977399d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.74770661d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.74770661d+12))*indns

    !OH_0001 + CH2DCDO_0001 -> D2O_0001 + CH3CO_0001 (2body_dust_dust)
    kall(763) = 2.60560084d+12*max(5.29977399d-16, exp(-6.00000000d+02 * invTd))&
        *(2.60560084d+12*max(5.29977399d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.74770661d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.74770661d+12))*indns

    !OD_0001 + CH2DCDO_0001 -> D2O_0001 + CH2DCO_0001 (2body_dust_dust)
    kall(764) = 2.53218886d+12*max(2.53907545d-16, exp(-6.00000000d+02 * invTd))&
        *(2.53218886d+12*max(2.53907545d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.80000000d+03*invTd)*1.74770661d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.80000000d+03*invTd)*1.74770661d+12))*indns

    !OH_0001 + CHD2CHO_0001 -> HDO_0001 + CH2DCO_0001 (2body_dust_dust)
    kall(765) = 2.60560084d+12*max(5.29977399d-16, exp(-6.00000000d+02 * invTd))&
        *(2.60560084d+12*max(5.29977399d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.74770661d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.74770661d+12))*indns

    !OH_0001 + CHD2CHO_0001 -> D2O_0001 + CH3CO_0001 (2body_dust_dust)
    kall(766) = 2.60560084d+12*max(5.29977399d-16, exp(-6.00000000d+02 * invTd))&
        *(2.60560084d+12*max(5.29977399d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.74770661d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.80000000d+03*invTd)*1.74770661d+12))*indns

    !OD_0001 + CHD2CHO_0001 -> D2O_0001 + CH2DCO_0001 (2body_dust_dust)
    kall(767) = 2.53218886d+12*max(2.53907545d-16, exp(-6.00000000d+02 * invTd))&
        *(2.53218886d+12*max(2.53907545d-16, exp(-6.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.80000000d+03*invTd)*1.74770661d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.80000000d+03*invTd)*1.74770661d+12))*indns

    !H_0001 + HSO_0001 -> OH_0001 + HS_0001 (2body_dust_dust)
    kall(768) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.69336042d+12))*indns

    !D_0001 + HSO_0001 -> OD_0001 + HS_0001 (2body_dust_dust)
    kall(769) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*1.69336042d+12))*indns

    !D_0001 + HSO_0001 -> OH_0001 + DS_0001 (2body_dust_dust)
    kall(770) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*1.69336042d+12))*indns

    !H_0001 + DSO_0001 -> OD_0001 + HS_0001 (2body_dust_dust)
    kall(771) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.67634129d+12))*indns

    !H_0001 + DSO_0001 -> OH_0001 + DS_0001 (2body_dust_dust)
    kall(772) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.67634129d+12))*indns

    !D_0001 + DSO_0001 -> OD_0001 + DS_0001 (2body_dust_dust)
    kall(773) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*1.67634129d+12))*indns

    !H_0001 + HCO_0001 -> p_H2_0001 + CO_0001 (2body_dust_dust)
    kall(774) = 2.50000000d-01 * 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !H_0001 + HCO_0001 -> o_H2_0001 + CO_0001 (2body_dust_dust)
    kall(775) = 7.50000000d-01 * 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !H_0001 + DCO_0001 -> HD_0001 + CO_0001 (2body_dust_dust)
    kall(776) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !D_0001 + HCO_0001 -> HD_0001 + CO_0001 (2body_dust_dust)
    kall(777) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !D_0001 + DCO_0001 -> p_D2_0001 + CO_0001 (2body_dust_dust)
    kall(778) = 3.33333333d-01 * 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !D_0001 + DCO_0001 -> o_D2_0001 + CO_0001 (2body_dust_dust)
    kall(779) = 6.66666667d-01 * 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !D_0001 + OH_0001 -> H_0001 + OD_0001 (2body_dust_dust)
    kall(780) = 2.60560084d+12*max(3.32580704d-14, exp(-8.10000000d+02 * invTd))&
        *(2.60560084d+12*max(3.32580704d-14, exp(-8.10000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !OH_0001 + CO_0001 -> HOCO_0001 (2body_dust)
    kall(781) = 2.60560084d+12*max(6.13353342d-13, exp(-4.50000000d+02 * invTd))&
        *(2.60560084d+12*max(6.13353342d-13, exp(-4.50000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))*indns

    !OH_0001 + CO_0001 -> HOCO_gas (2body_gas)
    kall(782) = 0d0*max(6.13353342d-13, exp(-4.50000000d+02 * invTd))&
        *(2.60560084d+12*max(6.13353342d-13, exp(-4.50000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))*indns

    !OD_0001 + CO_0001 -> DOCO_0001 (2body_dust)
    kall(783) = 2.53218886d+12*max(3.72391032d-13, exp(-4.50000000d+02 * invTd))&
        *(2.53218886d+12*max(3.72391032d-13, exp(-4.50000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))*indns

    !OD_0001 + CO_0001 -> DOCO_gas (2body_gas)
    kall(784) = 0d0*max(3.72391032d-13, exp(-4.50000000d+02 * invTd))&
        *(2.53218886d+12*max(3.72391032d-13, exp(-4.50000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))*indns

    !H_0001 + HOCO_0001 -> H2O_0001 + CO_0001 (2body_dust_dust)
    kall(785) = 2.00000000d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.00000000d+03*invTd)*1.05599575d+12))*indns

    !D_0001 + HOCO_0001 -> HDO_0001 + CO_0001 (2body_dust_dust)
    kall(786) = 2.00000000d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.00000000d+03*invTd)*1.05599575d+12))*indns

    !H_0001 + DOCO_0001 -> HDO_0001 + CO_0001 (2body_dust_dust)
    kall(787) = 2.00000000d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.00000000d+03*invTd)*1.04445447d+12))*indns

    !D_0001 + DOCO_0001 -> D2O_0001 + CO_0001 (2body_dust_dust)
    kall(788) = 2.00000000d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.00000000d+03*invTd)*1.04445447d+12))*indns

    !H_0001 + HOCO_0001 -> HCOOH_0001 (2body_dust)
    kall(789) = 9.93047902d-02*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.00000000d+03*invTd)*1.05599575d+12))*indns

    !H_0001 + HOCO_0001 -> HCOOH_gas (2body_gas)
    kall(790) = 6.95209761d-04*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.00000000d+03*invTd)*1.05599575d+12))*indns

    !D_0001 + HOCO_0001 -> HCOOD_0001 (2body_dust)
    kall(791) = 9.93046290d-02*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.00000000d+03*invTd)*1.05599575d+12))*indns

    !D_0001 + HOCO_0001 -> HCOOD_gas (2body_gas)
    kall(792) = 6.95371007d-04*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.00000000d+03*invTd)*1.05599575d+12))*indns

    !D_0001 + HOCO_0001 -> DCOOH_0001 (2body_dust)
    kall(793) = 9.93046290d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.00000000d+03*invTd)*1.05599575d+12))*indns

    !D_0001 + HOCO_0001 -> DCOOH_gas (2body_gas)
    kall(794) = 6.95371007d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.00000000d+03*invTd)*1.05599575d+12))*indns

    !H_0001 + DOCO_0001 -> HCOOD_0001 (2body_dust)
    kall(795) = 9.93069541d-02*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.00000000d+03*invTd)*1.04445447d+12))*indns

    !H_0001 + DOCO_0001 -> HCOOD_gas (2body_gas)
    kall(796) = 6.93045894d-04*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.00000000d+03*invTd)*1.04445447d+12))*indns

    !H_0001 + DOCO_0001 -> DCOOH_0001 (2body_dust)
    kall(797) = 9.93069541d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.00000000d+03*invTd)*1.04445447d+12))*indns

    !H_0001 + DOCO_0001 -> DCOOH_gas (2body_gas)
    kall(798) = 6.93045894d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.00000000d+03*invTd)*1.04445447d+12))*indns

    !D_0001 + DOCO_0001 -> DCOOD_0001 (2body_dust)
    kall(799) = 9.93067478d-02*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.00000000d+03*invTd)*1.04445447d+12))*indns

    !D_0001 + DOCO_0001 -> DCOOD_gas (2body_gas)
    kall(800) = 6.93252168d-04*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.00000000d+03*invTd)*1.04445447d+12))*indns

    !H_0001 + HOCO_0001 -> p_H2_0001 + CO2_0001 (2body_dust_dust)
    kall(801) = 2.50000000d-01 * 7.00000000d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.00000000d+03*invTd)*1.05599575d+12))*indns

    !H_0001 + HOCO_0001 -> o_H2_0001 + CO2_0001 (2body_dust_dust)
    kall(802) = 7.50000000d-01 * 7.00000000d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.00000000d+03*invTd)*1.05599575d+12))*indns

    !D_0001 + HOCO_0001 -> HD_0001 + CO2_0001 (2body_dust_dust)
    kall(803) = 7.00000000d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.00000000d+03*invTd)*1.05599575d+12))*indns

    !H_0001 + DOCO_0001 -> HD_0001 + CO2_0001 (2body_dust_dust)
    kall(804) = 7.00000000d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.00000000d+03*invTd)*1.04445447d+12))*indns

    !D_0001 + DOCO_0001 -> p_D2_0001 + CO2_0001 (2body_dust_dust)
    kall(805) = 3.33333333d-01 * 7.00000000d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.00000000d+03*invTd)*1.04445447d+12))*indns

    !D_0001 + DOCO_0001 -> o_D2_0001 + CO2_0001 (2body_dust_dust)
    kall(806) = 6.66666667d-01 * 7.00000000d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.00000000d+03*invTd)*1.04445447d+12))*indns

    !N_0001 + HOCO_0001 -> OH_0001 + OCN_0001 (2body_dust_dust)
    kall(807) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.00000000d+03*invTd)*1.05599575d+12))*indns

    !N_0001 + DOCO_0001 -> OD_0001 + OCN_0001 (2body_dust_dust)
    kall(808) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.00000000d+03*invTd)*1.04445447d+12))*indns

    !N_0001 + HOCO_0001 -> NH_0001 + CO2_0001 (2body_dust_dust)
    kall(809) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.00000000d+03*invTd)*1.05599575d+12))*indns

    !N_0001 + DOCO_0001 -> ND_0001 + CO2_0001 (2body_dust_dust)
    kall(810) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.00000000d+03*invTd)*1.04445447d+12))*indns

    !O_0001 + HOCO_0001 -> OH_0001 + CO2_0001 (2body_dust_dust)
    kall(811) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.00000000d+03*invTd)*1.05599575d+12))*indns

    !O_0001 + DOCO_0001 -> OD_0001 + CO2_0001 (2body_dust_dust)
    kall(812) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.00000000d+03*invTd)*1.04445447d+12))*indns

    !H_0001 + HCOOH_0001 -> p_H2_0001 + HOCO_0001 (2body_dust_dust)
    kall(813) = 2.50000000d-01 * 3.54191744d+12*max(1.84984352d-14, exp(-6.15000000d+03 * invTd))&
        *(3.54191744d+12*max(1.84984352d-14, exp(-6.15000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.74301897d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.74301897d+12))*indns

    !H_0001 + HCOOH_0001 -> o_H2_0001 + HOCO_0001 (2body_dust_dust)
    kall(814) = 7.50000000d-01 * 3.54191744d+12*max(1.84984352d-14, exp(-6.15000000d+03 * invTd))&
        *(3.54191744d+12*max(1.84984352d-14, exp(-6.15000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.74301897d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.74301897d+12))*indns

    !D_0001 + HCOOH_0001 -> HD_0001 + HOCO_0001 (2body_dust_dust)
    kall(815) = 2.55656768d+12*max(6.05570978d-20, exp(-6.15000000d+03 * invTd))&
        *(2.55656768d+12*max(6.05570978d-20, exp(-6.15000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.74301897d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.74301897d+12))*indns

    !D_0001 + HCOOH_0001 -> p_H2_0001 + DOCO_0001 (2body_dust_dust)
    kall(816) = 2.50000000d-01 * 2.55656768d+12*max(6.05570978d-20, exp(-6.15000000d+03 * invTd))&
        *(2.55656768d+12*max(6.05570978d-20, exp(-6.15000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.74301897d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.74301897d+12))*indns

    !D_0001 + HCOOH_0001 -> o_H2_0001 + DOCO_0001 (2body_dust_dust)
    kall(817) = 7.50000000d-01 * 2.55656768d+12*max(6.05570978d-20, exp(-6.15000000d+03 * invTd))&
        *(2.55656768d+12*max(6.05570978d-20, exp(-6.15000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.74301897d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.74301897d+12))*indns

    !H_0001 + HCOOD_0001 -> HD_0001 + HOCO_0001 (2body_dust_dust)
    kall(818) = 3.54191744d+12*max(1.83664638d-14, exp(-6.15000000d+03 * invTd))&
        *(3.54191744d+12*max(1.83664638d-14, exp(-6.15000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))*indns

    !H_0001 + HCOOD_0001 -> p_H2_0001 + DOCO_0001 (2body_dust_dust)
    kall(819) = 2.50000000d-01 * 3.54191744d+12*max(1.83664638d-14, exp(-6.15000000d+03 * invTd))&
        *(3.54191744d+12*max(1.83664638d-14, exp(-6.15000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))*indns

    !H_0001 + HCOOD_0001 -> o_H2_0001 + DOCO_0001 (2body_dust_dust)
    kall(820) = 7.50000000d-01 * 3.54191744d+12*max(1.83664638d-14, exp(-6.15000000d+03 * invTd))&
        *(3.54191744d+12*max(1.83664638d-14, exp(-6.15000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))*indns

    !H_0001 + DCOOH_0001 -> p_H2_0001 + DOCO_0001 (2body_dust_dust)
    kall(821) = 2.50000000d-01 * 3.54191744d+12*max(1.83664638d-14, exp(-6.15000000d+03 * invTd))&
        *(3.54191744d+12*max(1.83664638d-14, exp(-6.15000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))*indns

    !H_0001 + DCOOH_0001 -> o_H2_0001 + DOCO_0001 (2body_dust_dust)
    kall(822) = 7.50000000d-01 * 3.54191744d+12*max(1.83664638d-14, exp(-6.15000000d+03 * invTd))&
        *(3.54191744d+12*max(1.83664638d-14, exp(-6.15000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))*indns

    !H_0001 + DCOOH_0001 -> HD_0001 + HOCO_0001 (2body_dust_dust)
    kall(823) = 3.54191744d+12*max(1.83664638d-14, exp(-6.15000000d+03 * invTd))&
        *(3.54191744d+12*max(1.83664638d-14, exp(-6.15000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))*indns

    !D_0001 + DCOOH_0001 -> HD_0001 + DOCO_0001 (2body_dust_dust)
    kall(824) = 2.55656768d+12*max(5.93800886d-20, exp(-6.15000000d+03 * invTd))&
        *(2.55656768d+12*max(5.93800886d-20, exp(-6.15000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))*indns

    !D_0001 + DCOOH_0001 -> p_D2_0001 + HOCO_0001 (2body_dust_dust)
    kall(825) = 3.33333333d-01 * 2.55656768d+12*max(5.93800886d-20, exp(-6.15000000d+03 * invTd))&
        *(2.55656768d+12*max(5.93800886d-20, exp(-6.15000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))*indns

    !D_0001 + DCOOH_0001 -> o_D2_0001 + HOCO_0001 (2body_dust_dust)
    kall(826) = 6.66666667d-01 * 2.55656768d+12*max(5.93800886d-20, exp(-6.15000000d+03 * invTd))&
        *(2.55656768d+12*max(5.93800886d-20, exp(-6.15000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))*indns

    !D_0001 + HCOOD_0001 -> HD_0001 + DOCO_0001 (2body_dust_dust)
    kall(827) = 2.55656768d+12*max(5.93800886d-20, exp(-6.15000000d+03 * invTd))&
        *(2.55656768d+12*max(5.93800886d-20, exp(-6.15000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))*indns

    !D_0001 + HCOOD_0001 -> p_D2_0001 + HOCO_0001 (2body_dust_dust)
    kall(828) = 3.33333333d-01 * 2.55656768d+12*max(5.93800886d-20, exp(-6.15000000d+03 * invTd))&
        *(2.55656768d+12*max(5.93800886d-20, exp(-6.15000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))*indns

    !D_0001 + HCOOD_0001 -> o_D2_0001 + HOCO_0001 (2body_dust_dust)
    kall(829) = 6.66666667d-01 * 2.55656768d+12*max(5.93800886d-20, exp(-6.15000000d+03 * invTd))&
        *(2.55656768d+12*max(5.93800886d-20, exp(-6.15000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.72437652d+12))*indns

    !H_0001 + DCOOD_0001 -> HD_0001 + DOCO_0001 (2body_dust_dust)
    kall(830) = 3.54191744d+12*max(1.82407915d-14, exp(-6.15000000d+03 * invTd))&
        *(3.54191744d+12*max(1.82407915d-14, exp(-6.15000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.70631972d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.70631972d+12))*indns

    !H_0001 + DCOOD_0001 -> p_D2_0001 + HOCO_0001 (2body_dust_dust)
    kall(831) = 3.33333333d-01 * 3.54191744d+12*max(1.82407915d-14, exp(-6.15000000d+03 * invTd))&
        *(3.54191744d+12*max(1.82407915d-14, exp(-6.15000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.70631972d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.70631972d+12))*indns

    !H_0001 + DCOOD_0001 -> o_D2_0001 + HOCO_0001 (2body_dust_dust)
    kall(832) = 6.66666667d-01 * 3.54191744d+12*max(1.82407915d-14, exp(-6.15000000d+03 * invTd))&
        *(3.54191744d+12*max(1.82407915d-14, exp(-6.15000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.70631972d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.78500000d+03*invTd)*1.70631972d+12))*indns

    !D_0001 + DCOOD_0001 -> p_D2_0001 + DOCO_0001 (2body_dust_dust)
    kall(833) = 3.33333333d-01 * 2.55656768d+12*max(5.82721647d-20, exp(-6.15000000d+03 * invTd))&
        *(2.55656768d+12*max(5.82721647d-20, exp(-6.15000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.70631972d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.70631972d+12))*indns

    !D_0001 + DCOOD_0001 -> o_D2_0001 + DOCO_0001 (2body_dust_dust)
    kall(834) = 6.66666667d-01 * 2.55656768d+12*max(5.82721647d-20, exp(-6.15000000d+03 * invTd))&
        *(2.55656768d+12*max(5.82721647d-20, exp(-6.15000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.70631972d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.78500000d+03*invTd)*1.70631972d+12))*indns

    !C_0001 + C_0001 -> C2_0001 (2body_dust)
    kall(835) = 5d-1 * 9.91459283d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !C_0001 + C_0001 -> C2_gas (2body_gas)
    kall(836) = 5d-1 * 8.54071684d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !C_0001 + C2_0001 -> C3_0001 (2body_dust)
    kall(837) = 9.90918002d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !C_0001 + C2_0001 -> C3_gas (2body_gas)
    kall(838) = 9.08199773d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !C_0001 + CCH_0001 -> l_C3H_0001 (2body_dust)
    kall(839) = 9.94263356d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !C_0001 + CCH_0001 -> l_C3H_gas (2body_gas)
    kall(840) = 5.73664420d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !C_0001 + CCH_0001 -> c_C3H_0001 (2body_dust)
    kall(841) = 9.94263356d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !C_0001 + CCH_0001 -> c_C3H_gas (2body_gas)
    kall(842) = 5.73664420d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !C_0001 + CCD_0001 -> l_C3D_0001 (2body_dust)
    kall(843) = 9.94263356d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !C_0001 + CCD_0001 -> l_C3D_gas (2body_gas)
    kall(844) = 5.73664420d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !C_0001 + CCD_0001 -> c_C3D_0001 (2body_dust)
    kall(845) = 9.94263356d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !C_0001 + CCD_0001 -> c_C3D_gas (2body_gas)
    kall(846) = 5.73664420d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !C_0001 + C2H3_0001 -> C3H3_0001 (2body_dust)
    kall(847) = 1.00000000d+00*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.40000000d+03*invTd)*1.61306016d+12))*indns

    !C_0001 + C2D3_0001 -> C3D3_0001 (2body_dust)
    kall(848) = 1.00000000d+00*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.40000000d+03*invTd)*1.53028323d+12))*indns

    !C_0001 + C2H2D_0001 -> C3H2D_0001 (2body_dust)
    kall(849) = 1.00000000d+00*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.40000000d+03*invTd)*1.58399363d+12))*indns

    !C_0001 + C2HD2_0001 -> C3HD2_0001 (2body_dust)
    kall(850) = 1.00000000d+00*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.40000000d+03*invTd)*1.55644381d+12))*indns

    !C_0001 + CCN_0001 -> C3N_0001 (2body_dust)
    kall(851) = 9.92137047d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.20000000d+03*invTd)*1.25883049d+12))*indns

    !C_0001 + CCN_0001 -> C3N_gas (2body_gas)
    kall(852) = 7.86295286d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.20000000d+03*invTd)*1.25883049d+12))*indns

    !C_0001 + CCO_0001 -> C3O_0001 (2body_dust)
    kall(853) = 9.91896344d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-9.75000000d+02*invTd)*1.10596336d+12))*indns

    !C_0001 + CCO_0001 -> C3O_gas (2body_gas)
    kall(854) = 8.10365610d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-9.75000000d+02*invTd)*1.10596336d+12))*indns

    !C_0001 + CCS_0001 -> C3S_0001 (2body_dust)
    kall(855) = 9.90406046d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.35000000d+03*invTd)*1.09986986d+12))*indns

    !C_0001 + CCS_0001 -> C3S_gas (2body_gas)
    kall(856) = 9.59395438d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.35000000d+03*invTd)*1.09986986d+12))*indns

    !C_0001 + CH_0001 -> CCH_0001 (2body_dust)
    kall(857) = 9.91056451d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-4.62500000d+02*invTd)*1.33614202d+12))*indns

    !C_0001 + CH_0001 -> CCH_gas (2body_gas)
    kall(858) = 8.94354946d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-4.62500000d+02*invTd)*1.33614202d+12))*indns

    !C_0001 + CD_0001 -> CCD_0001 (2body_dust)
    kall(859) = 9.91062899d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-4.62500000d+02*invTd)*1.28753866d+12))*indns

    !C_0001 + CD_0001 -> CCD_gas (2body_gas)
    kall(860) = 8.93710054d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-4.62500000d+02*invTd)*1.28753866d+12))*indns

    !C_0001 + CH2_0001 -> C2H2_0001 (2body_dust)
    kall(861) = 9.91247875d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !C_0001 + CH2_0001 -> C2H2_gas (2body_gas)
    kall(862) = 8.75212489d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !C_0001 + CHD_0001 -> C2HD_0001 (2body_dust)
    kall(863) = 9.91476772d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !C_0001 + CHD_0001 -> C2HD_gas (2body_gas)
    kall(864) = 8.52322817d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !C_0001 + CD2_0001 -> C2D2_0001 (2body_dust)
    kall(865) = 9.91484301d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !C_0001 + CD2_0001 -> C2D2_gas (2body_gas)
    kall(866) = 8.51569908d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !C_0001 + CH3_0001 -> C2H3_0001 (2body_dust)
    kall(867) = 9.93049858d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !C_0001 + CH3_0001 -> C2H3_gas (2body_gas)
    kall(868) = 6.95014186d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !C_0001 + CD3_0001 -> C2D3_0001 (2body_dust)
    kall(869) = 9.93136187d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !C_0001 + CD3_0001 -> C2D3_gas (2body_gas)
    kall(870) = 6.86381342d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !C_0001 + CH2D_0001 -> C2H2D_0001 (2body_dust)
    kall(871) = 9.93092408d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !C_0001 + CH2D_0001 -> C2H2D_gas (2body_gas)
    kall(872) = 6.90759184d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !C_0001 + CHD2_0001 -> C2HD2_0001 (2body_dust)
    kall(873) = 9.93114141d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !C_0001 + CHD2_0001 -> C2HD2_gas (2body_gas)
    kall(874) = 6.88585942d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !C_0001 + CN_0001 -> CCN_0001 (2body_dust)
    kall(875) = 9.91302854d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !C_0001 + CN_0001 -> CCN_gas (2body_gas)
    kall(876) = 8.69714592d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !C_0001 + HS_0001 -> H_0001 + CS_0001 (2body_dust_dust)
    kall(877) = 1.00000000d+00*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.35000000d+03*invTd)*1.43277615d+12))*indns

    !C_0001 + DS_0001 -> D_0001 + CS_0001 (2body_dust_dust)
    kall(878) = 1.00000000d+00*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.35000000d+03*invTd)*1.41154866d+12))*indns

    !C_0001 + N_0001 -> CN_0001 (2body_dust)
    kall(879) = 9.90403878d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-3.60000000d+02*invTd)*1.13594070d+12))*indns

    !C_0001 + N_0001 -> CN_gas (2body_gas)
    kall(880) = 9.59612160d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-3.60000000d+02*invTd)*1.13594070d+12))*indns

    !C_0001 + NH_0001 -> HNC_0001 (2body_dust)
    kall(881) = 9.91100648d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !C_0001 + NH_0001 -> HNC_gas (2body_gas)
    kall(882) = 8.89935173d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !C_0001 + ND_0001 -> DNC_0001 (2body_dust)
    kall(883) = 9.91101999d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !C_0001 + ND_0001 -> DNC_gas (2body_gas)
    kall(884) = 8.89800136d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !C_0001 + NH2_0001 -> H_0001 + HNC_0001 (2body_dust_dust)
    kall(885) = 1.00000000d+00*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !C_0001 + NHD_0001 -> D_0001 + HNC_0001 (2body_dust_dust)
    kall(886) = 5.00000000d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !C_0001 + NHD_0001 -> H_0001 + DNC_0001 (2body_dust_dust)
    kall(887) = 5.00000000d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !C_0001 + ND2_0001 -> D_0001 + DNC_0001 (2body_dust_dust)
    kall(888) = 1.00000000d+00*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !C_0001 + NO_0001 -> O_0001 + CN_0001 (2body_dust_dust)
    kall(889) = 5.00000000d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !C_0001 + NO_0001 -> OCN_0001 (2body_dust)
    kall(890) = 4.95489467d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !C_0001 + NO_0001 -> OCN_gas (2body_gas)
    kall(891) = 4.51053349d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !C_0001 + NS_0001 -> CN_0001 + S_0001 (2body_dust_dust)
    kall(892) = 1.00000000d+00*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-9.50000000d+02*invTd)*1.01800829d+12))*indns

    !C_0001 + O_0001 -> CO_0001 (2body_dust)
    kall(893) = 9.90197658d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !C_0001 + O_0001 -> CO_gas (2body_gas)
    kall(894) = 9.80234157d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !C_0001 + O2_0001 -> O_0001 + CO_0001 (2body_dust_dust)
    kall(895) = 1.00000000d+00*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))*indns

    !C_0001 + OCN_0001 -> CN_0001 + CO_0001 (2body_dust_dust)
    kall(896) = 1.00000000d+00*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.20000000d+03*invTd)*1.19738664d+12))*indns

    !C_0001 + OH_0001 -> H_0001 + CO_0001 (2body_dust_dust)
    kall(897) = 1.00000000d+00*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !C_0001 + OD_0001 -> D_0001 + CO_0001 (2body_dust_dust)
    kall(898) = 1.00000000d+00*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !C_0001 + S_0001 -> CS_0001 (2body_dust)
    kall(899) = 9.90464994d-01*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !C_0001 + S_0001 -> CS_gas (2body_gas)
    kall(900) = 9.53500566d-03*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !C_0001 + SO_0001 -> CO_0001 + S_0001 (2body_dust_dust)
    kall(901) = 1.00000000d+00*((exp(-5.00000000d+03*invTd)*4.57259575d+12) + (exp(-1.40000000d+03*invTd)*1.20979512d+12))*indns

    !CH_0001 + C2_0001 -> l_C3H_0001 (2body_dust)
    kall(902) = 9.93578774d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !CH_0001 + C2_0001 -> l_C3H_gas (2body_gas)
    kall(903) = 6.42122578d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !CH_0001 + C2_0001 -> c_C3H_0001 (2body_dust)
    kall(904) = 9.93578774d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !CH_0001 + C2_0001 -> c_C3H_gas (2body_gas)
    kall(905) = 6.42122578d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !CD_0001 + C2_0001 -> l_C3D_0001 (2body_dust)
    kall(906) = 9.93599739d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !CD_0001 + C2_0001 -> l_C3D_gas (2body_gas)
    kall(907) = 6.40026121d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !CD_0001 + C2_0001 -> c_C3D_0001 (2body_dust)
    kall(908) = 9.93599739d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !CD_0001 + C2_0001 -> c_C3D_gas (2body_gas)
    kall(909) = 6.40026121d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !CH_0001 + CCH_0001 -> l_C3H2_0001 (2body_dust)
    kall(910) = 9.94983178d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CH_0001 + CCH_0001 -> l_C3H2_gas (2body_gas)
    kall(911) = 5.01682185d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CH_0001 + CCH_0001 -> c_C3H2_0001 (2body_dust)
    kall(912) = 9.94983178d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CH_0001 + CCH_0001 -> c_C3H2_gas (2body_gas)
    kall(913) = 5.01682185d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CH_0001 + CCD_0001 -> l_C3HD_0001 (2body_dust)
    kall(914) = 9.94983178d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CH_0001 + CCD_0001 -> l_C3HD_gas (2body_gas)
    kall(915) = 5.01682185d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CH_0001 + CCD_0001 -> c_C3HD_0001 (2body_dust)
    kall(916) = 9.94983178d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CH_0001 + CCD_0001 -> c_C3HD_gas (2body_gas)
    kall(917) = 5.01682185d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CD_0001 + CCH_0001 -> l_C3HD_0001 (2body_dust)
    kall(918) = 9.95018965d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CD_0001 + CCH_0001 -> l_C3HD_gas (2body_gas)
    kall(919) = 4.98103485d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CD_0001 + CCH_0001 -> c_C3HD_0001 (2body_dust)
    kall(920) = 9.95018965d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CD_0001 + CCH_0001 -> c_C3HD_gas (2body_gas)
    kall(921) = 4.98103485d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CD_0001 + CCD_0001 -> l_C3D2_0001 (2body_dust)
    kall(922) = 9.95018965d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CD_0001 + CCD_0001 -> l_C3D2_gas (2body_gas)
    kall(923) = 4.98103485d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CD_0001 + CCD_0001 -> c_C3D2_0001 (2body_dust)
    kall(924) = 9.95018965d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CD_0001 + CCD_0001 -> c_C3D2_gas (2body_gas)
    kall(925) = 4.98103485d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CH_0001 + C2H3_0001 -> C3H4_0001 (2body_dust)
    kall(926) = 9.95248605d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.40000000d+03*invTd)*1.61306016d+12))*indns

    !CH_0001 + C2H3_0001 -> C3H4_gas (2body_gas)
    kall(927) = 4.75139459d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.40000000d+03*invTd)*1.61306016d+12))*indns

    !CH_0001 + C2H2D_0001 -> C3H3D_0001 (2body_dust)
    kall(928) = 9.95248605d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.40000000d+03*invTd)*1.58399363d+12))*indns

    !CH_0001 + C2HD2_0001 -> C3H2D2_0001 (2body_dust)
    kall(929) = 9.95248605d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.40000000d+03*invTd)*1.55644381d+12))*indns

    !CD_0001 + C2H3_0001 -> C3H3D_0001 (2body_dust)
    kall(930) = 9.95276738d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.40000000d+03*invTd)*1.61306016d+12))*indns

    !CD_0001 + C2H2D_0001 -> C3H2D2_0001 (2body_dust)
    kall(931) = 9.95276738d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.40000000d+03*invTd)*1.58399363d+12))*indns

    !CH_0001 + CH_0001 -> C2H2_0001 (2body_dust)
    kall(932) = 5d-1 * 9.91164417d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-4.62500000d+02*invTd)*1.33614202d+12))*indns

    !CH_0001 + CH_0001 -> C2H2_gas (2body_gas)
    kall(933) = 5d-1 * 8.83558340d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-4.62500000d+02*invTd)*1.33614202d+12))*indns

    !CH_0001 + CD_0001 -> C2HD_0001 (2body_dust)
    kall(934) = 9.91352818d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-4.62500000d+02*invTd)*1.28753866d+12))*indns

    !CH_0001 + CD_0001 -> C2HD_gas (2body_gas)
    kall(935) = 8.64718246d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-4.62500000d+02*invTd)*1.28753866d+12))*indns

    !CD_0001 + CD_0001 -> C2D2_0001 (2body_dust)
    kall(936) = 5d-1 * 9.91359215d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-4.62500000d+02*invTd)*1.28753866d+12))*indns

    !CD_0001 + CD_0001 -> C2D2_gas (2body_gas)
    kall(937) = 5d-1 * 8.64078518d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-4.62500000d+02*invTd)*1.28753866d+12))*indns

    !CH_0001 + CH2_0001 -> C2H3_0001 (2body_dust)
    kall(938) = 9.92608042d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !CH_0001 + CH2_0001 -> C2H3_gas (2body_gas)
    kall(939) = 7.39195801d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !CH_0001 + CHD_0001 -> C2H2D_0001 (2body_dust)
    kall(940) = 9.92637909d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !CH_0001 + CHD_0001 -> C2H2D_gas (2body_gas)
    kall(941) = 7.36209106d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !CH_0001 + CD2_0001 -> C2HD2_0001 (2body_dust)
    kall(942) = 9.92653581d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !CH_0001 + CD2_0001 -> C2HD2_gas (2body_gas)
    kall(943) = 7.34641870d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !CD_0001 + CH2_0001 -> C2H2D_0001 (2body_dust)
    kall(944) = 9.92623750d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !CD_0001 + CH2_0001 -> C2H2D_gas (2body_gas)
    kall(945) = 7.37624975d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !CD_0001 + CHD_0001 -> C2HD2_0001 (2body_dust)
    kall(946) = 9.92653991d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !CD_0001 + CHD_0001 -> C2HD2_gas (2body_gas)
    kall(947) = 7.34600866d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !CD_0001 + CD2_0001 -> C2D3_0001 (2body_dust)
    kall(948) = 9.92669862d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !CD_0001 + CD2_0001 -> C2D3_gas (2body_gas)
    kall(949) = 7.33013825d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !CH_0001 + CH3_0001 -> C2H4_0001 (2body_dust)
    kall(950) = 9.93145807d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !CH_0001 + CH3_0001 -> C2H4_gas (2body_gas)
    kall(951) = 6.85419274d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !CH_0001 + CD3_0001 -> C2HD3_0001 (2body_dust)
    kall(952) = 9.93222747d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !CH_0001 + CD3_0001 -> C2HD3_gas (2body_gas)
    kall(953) = 6.77725256d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !CH_0001 + CH2D_0001 -> C2H3D_0001 (2body_dust)
    kall(954) = 9.93420974d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CH_0001 + CH2D_0001 -> C2H3D_gas (2body_gas)
    kall(955) = 6.57902617d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CH_0001 + CHD2_0001 -> C2H2D2_0001 (2body_dust)
    kall(956) = 9.93802070d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !CH_0001 + CHD2_0001 -> C2H2D2_gas (2body_gas)
    kall(957) = 6.19793038d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !CD_0001 + CH3_0001 -> C2H3D_0001 (2body_dust)
    kall(958) = 9.93399439d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !CD_0001 + CH3_0001 -> C2H3D_gas (2body_gas)
    kall(959) = 6.60056100d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !CD_0001 + CH2D_0001 -> C2H2D2_0001 (2body_dust)
    kall(960) = 9.93802647d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CD_0001 + CH2D_0001 -> C2H2D2_gas (2body_gas)
    kall(961) = 6.19735307d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CD_0001 + CHD2_0001 -> C2HD3_0001 (2body_dust)
    kall(962) = 9.93223260d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !CD_0001 + CHD2_0001 -> C2HD3_gas (2body_gas)
    kall(963) = 6.77673992d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !CH_0001 + CN_0001 -> HCCN_0001 (2body_dust)
    kall(964) = 1.00000000d+00*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !CD_0001 + CN_0001 -> DCCN_0001 (2body_dust)
    kall(965) = 1.00000000d+00*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !CH_0001 + HNO_0001 -> CH2_0001 + NO_0001 (2body_dust_dust)
    kall(966) = 1.00000000d+00*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))*indns

    !CH_0001 + DNO_0001 -> CHD_0001 + NO_0001 (2body_dust_dust)
    kall(967) = 1.00000000d+00*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))*indns

    !CD_0001 + HNO_0001 -> CHD_0001 + NO_0001 (2body_dust_dust)
    kall(968) = 1.00000000d+00*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))*indns

    !CD_0001 + DNO_0001 -> CD2_0001 + NO_0001 (2body_dust_dust)
    kall(969) = 1.00000000d+00*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))*indns

    !CH_0001 + NH_0001 -> H_0001 + HCN_0001 (2body_dust_dust)
    kall(970) = 3.33000000d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !CH_0001 + ND_0001 -> D_0001 + HCN_0001 (2body_dust_dust)
    kall(971) = 1.67000000d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !CH_0001 + ND_0001 -> H_0001 + DCN_0001 (2body_dust_dust)
    kall(972) = 1.67000000d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !CD_0001 + NH_0001 -> D_0001 + HCN_0001 (2body_dust_dust)
    kall(973) = 1.67000000d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !CD_0001 + NH_0001 -> H_0001 + DCN_0001 (2body_dust_dust)
    kall(974) = 1.67000000d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !CD_0001 + ND_0001 -> D_0001 + DCN_0001 (2body_dust_dust)
    kall(975) = 3.33000000d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !CH_0001 + NH_0001 -> H_0001 + HNC_0001 (2body_dust_dust)
    kall(976) = 3.33000000d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !CH_0001 + ND_0001 -> D_0001 + HNC_0001 (2body_dust_dust)
    kall(977) = 1.67000000d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !CH_0001 + ND_0001 -> H_0001 + DNC_0001 (2body_dust_dust)
    kall(978) = 1.67000000d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !CD_0001 + NH_0001 -> D_0001 + HNC_0001 (2body_dust_dust)
    kall(979) = 1.67000000d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !CD_0001 + NH_0001 -> H_0001 + DNC_0001 (2body_dust_dust)
    kall(980) = 1.67000000d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !CD_0001 + ND_0001 -> D_0001 + DNC_0001 (2body_dust_dust)
    kall(981) = 3.33000000d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !CH_0001 + NH2_0001 -> CH2NH_0001 (2body_dust)
    kall(982) = 9.94858633d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CH_0001 + NH2_0001 -> CH2NH_gas (2body_gas)
    kall(983) = 5.14136667d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CH_0001 + NHD_0001 -> CHDNH_0001 (2body_dust)
    kall(984) = 6.62585092d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CH_0001 + NHD_0001 -> CHDNH_gas (2body_gas)
    kall(985) = 3.41490817d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CH_0001 + NHD_0001 -> CH2ND_0001 (2body_dust)
    kall(986) = 3.31292546d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CH_0001 + NHD_0001 -> CH2ND_gas (2body_gas)
    kall(987) = 1.70745408d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CH_0001 + ND2_0001 -> CD2NH_0001 (2body_dust)
    kall(988) = 3.31296975d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CH_0001 + ND2_0001 -> CD2NH_gas (2body_gas)
    kall(989) = 1.70302496d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CH_0001 + ND2_0001 -> CHDND_0001 (2body_dust)
    kall(990) = 6.62593950d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CH_0001 + ND2_0001 -> CHDND_gas (2body_gas)
    kall(991) = 3.40604992d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CD_0001 + NH2_0001 -> CHDNH_0001 (2body_dust)
    kall(992) = 6.62594349d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CD_0001 + NH2_0001 -> CHDNH_gas (2body_gas)
    kall(993) = 3.40565110d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CD_0001 + NH2_0001 -> CH2ND_0001 (2body_dust)
    kall(994) = 3.31297174d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CD_0001 + NH2_0001 -> CH2ND_gas (2body_gas)
    kall(995) = 1.70282555d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CD_0001 + NHD_0001 -> CD2NH_0001 (2body_dust)
    kall(996) = 3.31301847d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CD_0001 + NHD_0001 -> CD2NH_gas (2body_gas)
    kall(997) = 1.69815314d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CD_0001 + NHD_0001 -> CHDND_0001 (2body_dust)
    kall(998) = 6.62603694d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CD_0001 + NHD_0001 -> CHDND_gas (2body_gas)
    kall(999) = 3.39630629d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CD_0001 + ND2_0001 -> CD2ND_0001 (2body_dust)
    kall(1000) = 9.94913890d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CD_0001 + ND2_0001 -> CD2ND_gas (2body_gas)
    kall(1001) = 5.08611010d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CH_0001 + NO_0001 -> O_0001 + HCN_0001 (2body_dust_dust)
    kall(1002) = 1.00000000d+00*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !CD_0001 + NO_0001 -> O_0001 + DCN_0001 (2body_dust_dust)
    kall(1003) = 1.00000000d+00*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !CH_0001 + O2_0001 -> O_0001 + HCO_0001 (2body_dust_dust)
    kall(1004) = 1.00000000d+00*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))*indns

    !CD_0001 + O2_0001 -> O_0001 + DCO_0001 (2body_dust_dust)
    kall(1005) = 1.00000000d+00*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))*indns

    !CH2_0001 + CH2_0001 -> C2H4_0001 (2body_dust)
    kall(1006) = 5d-1 * 9.93011977d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !CH2_0001 + CH2_0001 -> C2H4_gas (2body_gas)
    kall(1007) = 5d-1 * 6.98802283d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !CH2_0001 + CHD_0001 -> C2H3D_0001 (2body_dust)
    kall(1008) = 9.93261664d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !CH2_0001 + CHD_0001 -> C2H3D_gas (2body_gas)
    kall(1009) = 6.73833583d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !CH2_0001 + CD2_0001 -> C2H2D2_0001 (2body_dust)
    kall(1010) = 9.93642292d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !CH2_0001 + CD2_0001 -> C2H2D2_gas (2body_gas)
    kall(1011) = 6.35770776d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !CHD_0001 + CHD_0001 -> C2H2D2_0001 (2body_dust)
    kall(1012) = 5d-1 * 9.93661538d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !CHD_0001 + CHD_0001 -> C2H2D2_gas (2body_gas)
    kall(1013) = 5d-1 * 6.33846221d-03*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !CHD_0001 + CD2_0001 -> C2HD3_0001 (2body_dust)
    kall(1014) = 9.93098243d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !CHD_0001 + CD2_0001 -> C2HD3_gas (2body_gas)
    kall(1015) = 6.90175704d-03*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !CH2_0001 + CH3_0001 -> C2H5_0001 (2body_dust)
    kall(1016) = 9.96126903d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !CH2_0001 + CH3_0001 -> C2H5_gas (2body_gas)
    kall(1017) = 3.87309681d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !CH2_0001 + CD3_0001 -> C2H2D3_0001 (2body_dust)
    kall(1018) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !CH2_0001 + CH2D_0001 -> C2H4D_0001 (2body_dust)
    kall(1019) = 9.96217151d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CH2_0001 + CHD2_0001 -> C2H3D2_0001 (2body_dust)
    kall(1020) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !CHD_0001 + CH3_0001 -> C2H4D_0001 (2body_dust)
    kall(1021) = 9.96214409d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !CHD_0001 + CH2D_0001 -> C2H3D2_0001 (2body_dust)
    kall(1022) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CHD_0001 + CHD2_0001 -> C2H2D3_0001 (2body_dust)
    kall(1023) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !CH3_0001 + CD2_0001 -> C2H3D2_0001 (2body_dust)
    kall(1024) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !CD2_0001 + CH2D_0001 -> C2H2D3_0001 (2body_dust)
    kall(1025) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CH2_0001 + HNO_0001 -> CH3_0001 + NO_0001 (2body_dust_dust)
    kall(1026) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))*indns

    !CH2_0001 + DNO_0001 -> CH2D_0001 + NO_0001 (2body_dust_dust)
    kall(1027) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))*indns

    !CHD_0001 + HNO_0001 -> CH2D_0001 + NO_0001 (2body_dust_dust)
    kall(1028) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))*indns

    !CHD_0001 + DNO_0001 -> CHD2_0001 + NO_0001 (2body_dust_dust)
    kall(1029) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))*indns

    !CD2_0001 + HNO_0001 -> CHD2_0001 + NO_0001 (2body_dust_dust)
    kall(1030) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))*indns

    !CD2_0001 + DNO_0001 -> CD3_0001 + NO_0001 (2body_dust_dust)
    kall(1031) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))*indns

    !CH2_0001 + NH2_0001 -> CH2NH2_0001 (2body_dust)
    kall(1032) = 9.97509858d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CH2_0001 + NH2_0001 -> CH2NH2_gas (2body_gas)
    kall(1033) = 2.49014152d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CH2_0001 + NHD_0001 -> CHDNH2_0001 (2body_dust)
    kall(1034) = 4.97768565d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CH2_0001 + NHD_0001 -> CHDNH2_gas (2body_gas)
    kall(1035) = 1.23143499d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CH2_0001 + NHD_0001 -> CH2NHD_0001 (2body_dust)
    kall(1036) = 4.97768565d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CH2_0001 + NHD_0001 -> CH2NHD_gas (2body_gas)
    kall(1037) = 1.23143499d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CH2_0001 + ND2_0001 -> CD2NH2_0001 (2body_dust)
    kall(1038) = 1.65593893d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CH2_0001 + ND2_0001 -> CD2NH2_gas (2body_gas)
    kall(1039) = 4.06107300d-04*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CH2_0001 + ND2_0001 -> CH2ND2_0001 (2body_dust)
    kall(1040) = 1.65593893d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CH2_0001 + ND2_0001 -> CH2ND2_gas (2body_gas)
    kall(1041) = 4.06107300d-04*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CH2_0001 + ND2_0001 -> CHDNHD_0001 (2body_dust)
    kall(1042) = 6.65368232d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CH2_0001 + ND2_0001 -> CHDNHD_gas (2body_gas)
    kall(1043) = 1.63176849d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CHD_0001 + NH2_0001 -> CHDNH2_0001 (2body_dust)
    kall(1044) = 4.97799719d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CHD_0001 + NH2_0001 -> CHDNH2_gas (2body_gas)
    kall(1045) = 1.20028104d-03*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CHD_0001 + NH2_0001 -> CH2NHD_0001 (2body_dust)
    kall(1046) = 4.97799719d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CHD_0001 + NH2_0001 -> CH2NHD_gas (2body_gas)
    kall(1047) = 1.20028104d-03*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CHD_0001 + NHD_0001 -> CD2NH2_0001 (2body_dust)
    kall(1048) = 1.65604471d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CHD_0001 + NHD_0001 -> CD2NH2_gas (2body_gas)
    kall(1049) = 3.95528604d-04*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CHD_0001 + NHD_0001 -> CH2ND2_0001 (2body_dust)
    kall(1050) = 1.65604471d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CHD_0001 + NHD_0001 -> CH2ND2_gas (2body_gas)
    kall(1051) = 3.95528604d-04*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CHD_0001 + NHD_0001 -> CHDNHD_0001 (2body_dust)
    kall(1052) = 6.65410737d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CHD_0001 + NHD_0001 -> CHDNHD_gas (2body_gas)
    kall(1053) = 1.58926252d-03*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CHD_0001 + ND2_0001 -> CD2NHD_0001 (2body_dust)
    kall(1054) = 4.97821856d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CHD_0001 + ND2_0001 -> CD2NHD_gas (2body_gas)
    kall(1055) = 1.17814355d-03*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CHD_0001 + ND2_0001 -> CHDND2_0001 (2body_dust)
    kall(1056) = 4.97821856d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CHD_0001 + ND2_0001 -> CHDND2_gas (2body_gas)
    kall(1057) = 1.17814355d-03*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !CD2_0001 + NH2_0001 -> CD2NH2_0001 (2body_dust)
    kall(1058) = 1.65608043d-01*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CD2_0001 + NH2_0001 -> CD2NH2_gas (2body_gas)
    kall(1059) = 3.91956964d-04*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CD2_0001 + NH2_0001 -> CH2ND2_0001 (2body_dust)
    kall(1060) = 1.65608043d-01*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CD2_0001 + NH2_0001 -> CH2ND2_gas (2body_gas)
    kall(1061) = 3.91956964d-04*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CD2_0001 + NH2_0001 -> CHDNHD_0001 (2body_dust)
    kall(1062) = 6.65425089d-01*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CD2_0001 + NH2_0001 -> CHDNHD_gas (2body_gas)
    kall(1063) = 1.57491141d-03*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !CD2_0001 + NHD_0001 -> CD2NHD_0001 (2body_dust)
    kall(1064) = 4.97833165d-01*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CD2_0001 + NHD_0001 -> CD2NHD_gas (2body_gas)
    kall(1065) = 1.16683486d-03*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CD2_0001 + NHD_0001 -> CHDND2_0001 (2body_dust)
    kall(1066) = 4.97833165d-01*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CD2_0001 + NHD_0001 -> CHDND2_gas (2body_gas)
    kall(1067) = 1.16683486d-03*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !CH2_0001 + O2_0001 -> O_0001 + H2CO_0001 (2body_dust_dust)
    kall(1068) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))*indns

    !CHD_0001 + O2_0001 -> O_0001 + HDCO_0001 (2body_dust_dust)
    kall(1069) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))*indns

    !CD2_0001 + O2_0001 -> O_0001 + D2CO_0001 (2body_dust_dust)
    kall(1070) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))*indns

    !CH3_0001 + CH2OH_0001 -> C2H5OH_0001 (2body_dust)
    kall(1071) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !CH3_0001 + CH3O_0001 -> CH3OCH3_0001 (2body_dust)
    kall(1072) = 9.98103022d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !CH3_0001 + CH3O_0001 -> CH3OCH3_gas (2body_gas)
    kall(1073) = 1.89697805d-03*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !CH3_0001 + CH2DO_0001 -> CH2DOCH3_0001 (2body_dust)
    kall(1074) = 9.99690052d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !CH3_0001 + CHD2O_0001 -> CHD2OCH3_0001 (2body_dust)
    kall(1075) = 3.99847684d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !CH3_0001 + CHD2O_0001 -> CH2DOCH2D_0001 (2body_dust)
    kall(1076) = 5.99771525d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !CH3_0001 + CD3O_0001 -> CD3OCH3_0001 (2body_dust)
    kall(1077) = 1.00000000d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-2.20000000d+03*invTd)*1.80193994d+12))*indns

    !CH3_0001 + CD3O_0001 -> CHD2OCH2D_0001 (2body_dust)
    kall(1078) = 9.00000000d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-2.20000000d+03*invTd)*1.80193994d+12))*indns

    !CD3_0001 + CH3O_0001 -> CD3OCH3_0001 (2body_dust)
    kall(1079) = 1.00000000d-01*((exp(-8.00000000d+02*invTd)*1.49340352d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !CD3_0001 + CH3O_0001 -> CHD2OCH2D_0001 (2body_dust)
    kall(1080) = 9.00000000d-01*((exp(-8.00000000d+02*invTd)*1.49340352d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !CH2D_0001 + CH3O_0001 -> CH2DOCH3_0001 (2body_dust)
    kall(1081) = 9.99570615d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !CH2D_0001 + CH2DO_0001 -> CHD2OCH3_0001 (2body_dust)
    kall(1082) = 3.99890941d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !CH2D_0001 + CH2DO_0001 -> CH2DOCH2D_0001 (2body_dust)
    kall(1083) = 5.99836412d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !CH2D_0001 + CHD2O_0001 -> CD3OCH3_0001 (2body_dust)
    kall(1084) = 1.00000000d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !CH2D_0001 + CHD2O_0001 -> CHD2OCH2D_0001 (2body_dust)
    kall(1085) = 9.00000000d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !CHD2_0001 + CH3O_0001 -> CHD2OCH3_0001 (2body_dust)
    kall(1086) = 3.99836879d-01*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !CHD2_0001 + CH3O_0001 -> CH2DOCH2D_0001 (2body_dust)
    kall(1087) = 5.99755318d-01*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !CHD2_0001 + CH2DO_0001 -> CD3OCH3_0001 (2body_dust)
    kall(1088) = 1.00000000d-01*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !CHD2_0001 + CH2DO_0001 -> CHD2OCH2D_0001 (2body_dust)
    kall(1089) = 9.00000000d-01*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !CH3_0001 + CH3_0001 -> C2H6_0001 (2body_dust)
    kall(1090) = 5d-1 * 9.94870199d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !CH3_0001 + CH3_0001 -> C2H6_gas (2body_gas)
    kall(1091) = 5d-1 * 5.12980076d-03*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !CH3_0001 + CD3_0001 -> C2H3D3_0001 (2body_dust)
    kall(1092) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !CH3_0001 + CH2D_0001 -> C2H5D_0001 (2body_dust)
    kall(1093) = 9.94965489d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CH3_0001 + CHD2_0001 -> C2H4D2_0001 (2body_dust)
    kall(1094) = 9.95014495d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !CH2D_0001 + CH2D_0001 -> C2H4D2_0001 (2body_dust)
    kall(1095) = 5d-1 * 9.95064442d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CH2D_0001 + CHD2_0001 -> C2H3D3_0001 (2body_dust)
    kall(1096) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !CH3_0001 + CN_0001 -> CH3CN_0001 (2body_dust)
    kall(1097) = 9.96187180d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !CH3_0001 + CN_0001 -> CH3CN_gas (2body_gas)
    kall(1098) = 3.81282026d-03*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !CD3_0001 + CN_0001 -> CD3CN_0001 (2body_dust)
    kall(1099) = 9.96340607d-01*((exp(-8.00000000d+02*invTd)*1.49340352d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !CD3_0001 + CN_0001 -> CD3CN_gas (2body_gas)
    kall(1100) = 3.65939274d-03*((exp(-8.00000000d+02*invTd)*1.49340352d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !CH2D_0001 + CN_0001 -> CH2DCN_0001 (2body_dust)
    kall(1101) = 9.96263065d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !CH2D_0001 + CN_0001 -> CH2DCN_gas (2body_gas)
    kall(1102) = 3.73693483d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !CHD2_0001 + CN_0001 -> CHD2CN_0001 (2body_dust)
    kall(1103) = 9.96301627d-01*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !CHD2_0001 + CN_0001 -> CHD2CN_gas (2body_gas)
    kall(1104) = 3.69837321d-03*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !CH3_0001 + HCO_0001 -> CH3CHO_0001 (2body_dust)
    kall(1105) = 9.98752331d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !CH3_0001 + HCO_0001 -> CH3CHO_gas (2body_gas)
    kall(1106) = 1.24766855d-03*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !CH3_0001 + DCO_0001 -> CH2DCHO_0001 (2body_dust)
    kall(1107) = 7.49158365d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !CH3_0001 + DCO_0001 -> CH3CDO_0001 (2body_dust)
    kall(1108) = 2.49719455d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !CD3_0001 + HCO_0001 -> CD3CHO_0001 (2body_dust)
    kall(1109) = 2.50000000d-01*((exp(-8.00000000d+02*invTd)*1.49340352d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !CD3_0001 + HCO_0001 -> CHD2CDO_0001 (2body_dust)
    kall(1110) = 7.50000000d-01*((exp(-8.00000000d+02*invTd)*1.49340352d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !CH2D_0001 + HCO_0001 -> CH2DCHO_0001 (2body_dust)
    kall(1111) = 7.49196154d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !CH2D_0001 + HCO_0001 -> CH3CDO_0001 (2body_dust)
    kall(1112) = 2.49732051d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !CH2D_0001 + DCO_0001 -> CHD2CHO_0001 (2body_dust)
    kall(1113) = 4.99477210d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !CH2D_0001 + DCO_0001 -> CH2DCDO_0001 (2body_dust)
    kall(1114) = 4.99477210d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !CHD2_0001 + HCO_0001 -> CHD2CHO_0001 (2body_dust)
    kall(1115) = 4.99483197d-01*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !CHD2_0001 + HCO_0001 -> CH2DCDO_0001 (2body_dust)
    kall(1116) = 4.99483197d-01*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !CHD2_0001 + DCO_0001 -> CD3CHO_0001 (2body_dust)
    kall(1117) = 2.50000000d-01*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !CHD2_0001 + DCO_0001 -> CHD2CDO_0001 (2body_dust)
    kall(1118) = 7.50000000d-01*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !CH3_0001 + HNO_0001 -> CH4_0001 + NO_0001 (2body_dust_dust)
    kall(1119) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))*indns

    !CH3_0001 + DNO_0001 -> CH3D_0001 + NO_0001 (2body_dust_dust)
    kall(1120) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))*indns

    !CD3_0001 + HNO_0001 -> CHD3_0001 + NO_0001 (2body_dust_dust)
    kall(1121) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.49340352d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))*indns

    !CH2D_0001 + HNO_0001 -> CH3D_0001 + NO_0001 (2body_dust_dust)
    kall(1122) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))*indns

    !CH2D_0001 + DNO_0001 -> CH2D2_0001 + NO_0001 (2body_dust_dust)
    kall(1123) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))*indns

    !CHD2_0001 + HNO_0001 -> CH2D2_0001 + NO_0001 (2body_dust_dust)
    kall(1124) = 1.00000000d+00*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))*indns

    !CHD2_0001 + DNO_0001 -> CHD3_0001 + NO_0001 (2body_dust_dust)
    kall(1125) = 1.00000000d+00*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))*indns

    !CH4_0001 + CCH_0001 -> CH3_0001 + C2H2_0001 (2body_dust_dust)
    kall(1126) = 1.73517809d+12*max(2.26960086d-28, exp(-2.50000000d+03 * invTd))&
        *(1.73517809d+12*max(2.26960086d-28, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.22695619d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.22695619d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CH4_0001 + CCD_0001 -> CH2D_0001 + C2H2_0001 (2body_dust_dust)
    kall(1127) = 1.02088926d+12*max(1.39998848d-28, exp(-2.50000000d+03 * invTd))&
        *(1.70148210d+12*max(1.39998848d-28, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.22695619d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.22695619d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CH4_0001 + CCD_0001 -> CH3_0001 + C2HD_0001 (2body_dust_dust)
    kall(1128) = 6.80592840d+11*max(1.39998848d-28, exp(-2.50000000d+03 * invTd))&
        *(1.70148210d+12*max(1.39998848d-28, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.22695619d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.22695619d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CH3D_0001 + CCH_0001 -> CH2D_0001 + C2H2_0001 (2body_dust_dust)
    kall(1129) = 1.04110685d+12*max(7.02154646d-29, exp(-2.50000000d+03 * invTd))&
        *(1.73517809d+12*max(7.02154646d-29, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.19032235d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.19032235d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CH3D_0001 + CCH_0001 -> CH3_0001 + C2HD_0001 (2body_dust_dust)
    kall(1130) = 6.94071234d+11*max(7.02154646d-29, exp(-2.50000000d+03 * invTd))&
        *(1.73517809d+12*max(7.02154646d-29, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.19032235d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.19032235d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CH3D_0001 + CCD_0001 -> CHD2_0001 + C2H2_0001 (2body_dust_dust)
    kall(1131) = 5.10444630d+11*max(4.21401780d-29, exp(-2.50000000d+03 * invTd))&
        *(1.70148210d+12*max(4.21401780d-29, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.19032235d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.19032235d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CH3D_0001 + CCD_0001 -> CH2D_0001 + C2HD_0001 (2body_dust_dust)
    kall(1132) = 1.02088926d+12*max(4.21401780d-29, exp(-2.50000000d+03 * invTd))&
        *(1.70148210d+12*max(4.21401780d-29, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.19032235d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.19032235d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CH3D_0001 + CCD_0001 -> CH3_0001 + C2D2_0001 (2body_dust_dust)
    kall(1133) = 1.70148210d+11*max(4.21401780d-29, exp(-2.50000000d+03 * invTd))&
        *(1.70148210d+12*max(4.21401780d-29, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.19032235d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.19032235d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CHD3_0001 + CCH_0001 -> CD3_0001 + C2H2_0001 (2body_dust_dust)
    kall(1134) = 1.73517809d+11*max(8.33088901d-30, exp(-2.50000000d+03 * invTd))&
        *(1.73517809d+12*max(8.33088901d-30, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.12593222d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.12593222d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CHD3_0001 + CCH_0001 -> CHD2_0001 + C2HD_0001 (2body_dust_dust)
    kall(1135) = 1.04110685d+12*max(8.33088901d-30, exp(-2.50000000d+03 * invTd))&
        *(1.73517809d+12*max(8.33088901d-30, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.12593222d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.12593222d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CHD3_0001 + CCH_0001 -> CH2D_0001 + C2D2_0001 (2body_dust_dust)
    kall(1136) = 5.20553426d+11*max(8.33088901d-30, exp(-2.50000000d+03 * invTd))&
        *(1.73517809d+12*max(8.33088901d-30, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.12593222d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.12593222d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CHD3_0001 + CCD_0001 -> CD3_0001 + C2HD_0001 (2body_dust_dust)
    kall(1137) = 6.80592840d+11*max(4.74417697d-30, exp(-2.50000000d+03 * invTd))&
        *(1.70148210d+12*max(4.74417697d-30, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.12593222d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.12593222d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CHD3_0001 + CCD_0001 -> CHD2_0001 + C2D2_0001 (2body_dust_dust)
    kall(1138) = 1.02088926d+12*max(4.74417697d-30, exp(-2.50000000d+03 * invTd))&
        *(1.70148210d+12*max(4.74417697d-30, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.12593222d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.12593222d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CH2D2_0001 + CCH_0001 -> CHD2_0001 + C2H2_0001 (2body_dust_dust)
    kall(1139) = 5.20553426d+11*max(2.33914873d-29, exp(-2.50000000d+03 * invTd))&
        *(1.73517809d+12*max(2.33914873d-29, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.15678539d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.15678539d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CH2D2_0001 + CCH_0001 -> CH2D_0001 + C2HD_0001 (2body_dust_dust)
    kall(1140) = 1.04110685d+12*max(2.33914873d-29, exp(-2.50000000d+03 * invTd))&
        *(1.73517809d+12*max(2.33914873d-29, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.15678539d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.15678539d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CH2D2_0001 + CCH_0001 -> CH3_0001 + C2D2_0001 (2body_dust_dust)
    kall(1141) = 1.73517809d+11*max(2.33914873d-29, exp(-2.50000000d+03 * invTd))&
        *(1.73517809d+12*max(2.33914873d-29, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.15678539d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.15678539d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !CH2D2_0001 + CCD_0001 -> CD3_0001 + C2H2_0001 (2body_dust_dust)
    kall(1142) = 1.70148210d+11*max(1.36695815d-29, exp(-2.50000000d+03 * invTd))&
        *(1.70148210d+12*max(1.36695815d-29, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.15678539d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.15678539d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CH2D2_0001 + CCD_0001 -> CHD2_0001 + C2HD_0001 (2body_dust_dust)
    kall(1143) = 1.02088926d+12*max(1.36695815d-29, exp(-2.50000000d+03 * invTd))&
        *(1.70148210d+12*max(1.36695815d-29, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.15678539d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.15678539d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !CH2D2_0001 + CCD_0001 -> CH2D_0001 + C2D2_0001 (2body_dust_dust)
    kall(1144) = 5.10444630d+11*max(1.36695815d-29, exp(-2.50000000d+03 * invTd))&
        *(1.70148210d+12*max(1.36695815d-29, exp(-2.50000000d+03 * invTd)) + (exp(-4.80000000d+02*invTd)*1.15678539d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-4.80000000d+02*invTd)*1.15678539d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !H_0001 + C_0001 -> CH_0001 (2body_dust)
    kall(1145) = 9.90325396d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !H_0001 + C_0001 -> CH_gas (2body_gas)
    kall(1146) = 9.67460425d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !D_0001 + C_0001 -> CD_0001 (2body_dust)
    kall(1147) = 9.90318096d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !D_0001 + C_0001 -> CD_gas (2body_gas)
    kall(1148) = 9.68190392d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !H_0001 + C2_0001 -> CCH_0001 (2body_dust)
    kall(1149) = 9.91567341d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !H_0001 + C2_0001 -> CCH_gas (2body_gas)
    kall(1150) = 8.43265868d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !D_0001 + C2_0001 -> CCD_0001 (2body_dust)
    kall(1151) = 9.91550046d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !D_0001 + C2_0001 -> CCD_gas (2body_gas)
    kall(1152) = 8.44995383d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !H_0001 + CCH_0001 -> C2H2_0001 (2body_dust)
    kall(1153) = 9.91721252d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !H_0001 + CCH_0001 -> C2H2_gas (2body_gas)
    kall(1154) = 8.27874757d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !H_0001 + CCD_0001 -> C2HD_0001 (2body_dust)
    kall(1155) = 9.92183449d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !H_0001 + CCD_0001 -> C2HD_gas (2body_gas)
    kall(1156) = 7.81655108d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !D_0001 + CCH_0001 -> C2HD_0001 (2body_dust)
    kall(1157) = 9.92163250d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !D_0001 + CCH_0001 -> C2HD_gas (2body_gas)
    kall(1158) = 7.83674967d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !D_0001 + CCD_0001 -> C2D2_0001 (2body_dust)
    kall(1159) = 9.92163250d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !D_0001 + CCD_0001 -> C2D2_gas (2body_gas)
    kall(1160) = 7.83674967d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !H_0001 + C2H2_0001 -> C2H3_0001 (2body_dust)
    kall(1161) = 0d0*max(9.07519227d-07, exp(-1.21000000d+03 * invTd))&
        *(3.54191744d+12*max(9.07519227d-07, exp(-1.21000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.29350000d+03*invTd)*1.58002868d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.29350000d+03*invTd)*1.58002868d+12))*indns

    !H_0001 + C2H2_0001 -> C2H3_gas (2body_gas)
    kall(1162) = 0d0*max(9.07519227d-07, exp(-1.21000000d+03 * invTd))&
        *(3.54191744d+12*max(9.07519227d-07, exp(-1.21000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.29350000d+03*invTd)*1.58002868d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.29350000d+03*invTd)*1.58002868d+12))*indns

    !H_0001 + C2HD_0001 -> C2H2D_0001 (2body_dust)
    kall(1163) = 3.53163380d+12*max(8.98891833d-07, exp(-1.21000000d+03 * invTd))&
        *(3.54191744d+12*max(8.98891833d-07, exp(-1.21000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.29350000d+03*invTd)*1.55049283d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.29350000d+03*invTd)*1.55049283d+12))*indns

    !H_0001 + C2HD_0001 -> C2H2D_gas (2body_gas)
    kall(1164) = 1.02836357d+10*max(8.98891833d-07, exp(-1.21000000d+03 * invTd))&
        *(3.54191744d+12*max(8.98891833d-07, exp(-1.21000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.29350000d+03*invTd)*1.55049283d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.29350000d+03*invTd)*1.55049283d+12))*indns

    !H_0001 + C2D2_0001 -> C2HD2_0001 (2body_dust)
    kall(1165) = 3.53163380d+12*max(8.90938421d-07, exp(-1.21000000d+03 * invTd))&
        *(3.54191744d+12*max(8.90938421d-07, exp(-1.21000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.29350000d+03*invTd)*1.52255374d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.29350000d+03*invTd)*1.52255374d+12))*indns

    !H_0001 + C2D2_0001 -> C2HD2_gas (2body_gas)
    kall(1166) = 1.02836357d+10*max(8.90938421d-07, exp(-1.21000000d+03 * invTd))&
        *(3.54191744d+12*max(8.90938421d-07, exp(-1.21000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.29350000d+03*invTd)*1.52255374d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.29350000d+03*invTd)*1.52255374d+12))*indns

    !D_0001 + C2H2_0001 -> C2H2D_0001 (2body_dust)
    kall(1167) = 0d0*max(4.06529653d-09, exp(-1.21000000d+03 * invTd))&
        *(2.55656768d+12*max(4.06529653d-09, exp(-1.21000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.29350000d+03*invTd)*1.58002868d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.29350000d+03*invTd)*1.58002868d+12))*indns

    !D_0001 + C2H2_0001 -> C2H2D_gas (2body_gas)
    kall(1168) = 0d0*max(4.06529653d-09, exp(-1.21000000d+03 * invTd))&
        *(2.55656768d+12*max(4.06529653d-09, exp(-1.21000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.29350000d+03*invTd)*1.58002868d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.29350000d+03*invTd)*1.58002868d+12))*indns

    !D_0001 + C2HD_0001 -> C2HD2_0001 (2body_dust)
    kall(1169) = 2.54882799d+12*max(3.96251641d-09, exp(-1.21000000d+03 * invTd))&
        *(2.55656768d+12*max(3.96251641d-09, exp(-1.21000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.29350000d+03*invTd)*1.55049283d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.29350000d+03*invTd)*1.55049283d+12))*indns

    !D_0001 + C2HD_0001 -> C2HD2_gas (2body_gas)
    kall(1170) = 7.73969266d+09*max(3.96251641d-09, exp(-1.21000000d+03 * invTd))&
        *(2.55656768d+12*max(3.96251641d-09, exp(-1.21000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.29350000d+03*invTd)*1.55049283d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.29350000d+03*invTd)*1.55049283d+12))*indns

    !D_0001 + C2D2_0001 -> C2D3_0001 (2body_dust)
    kall(1171) = 2.54882799d+12*max(3.86905225d-09, exp(-1.21000000d+03 * invTd))&
        *(2.55656768d+12*max(3.86905225d-09, exp(-1.21000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.29350000d+03*invTd)*1.52255374d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.29350000d+03*invTd)*1.52255374d+12))*indns

    !D_0001 + C2D2_0001 -> C2D3_gas (2body_gas)
    kall(1172) = 7.73969266d+09*max(3.86905225d-09, exp(-1.21000000d+03 * invTd))&
        *(2.55656768d+12*max(3.86905225d-09, exp(-1.21000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.29350000d+03*invTd)*1.52255374d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.29350000d+03*invTd)*1.52255374d+12))*indns

    !H_0001 + C2H3_0001 -> C2H4_0001 (2body_dust)
    kall(1173) = 9.94612749d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.40000000d+03*invTd)*1.61306016d+12))*indns

    !H_0001 + C2H3_0001 -> C2H4_gas (2body_gas)
    kall(1174) = 5.38725064d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.40000000d+03*invTd)*1.61306016d+12))*indns

    !H_0001 + C2D3_0001 -> C2HD3_0001 (2body_dust)
    kall(1175) = 9.94612749d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.40000000d+03*invTd)*1.53028323d+12))*indns

    !H_0001 + C2D3_0001 -> C2HD3_gas (2body_gas)
    kall(1176) = 5.38725064d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.40000000d+03*invTd)*1.53028323d+12))*indns

    !H_0001 + C2H2D_0001 -> C2H3D_0001 (2body_dust)
    kall(1177) = 9.95125805d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.40000000d+03*invTd)*1.58399363d+12))*indns

    !H_0001 + C2H2D_0001 -> C2H3D_gas (2body_gas)
    kall(1178) = 4.87419462d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.40000000d+03*invTd)*1.58399363d+12))*indns

    !H_0001 + C2HD2_0001 -> C2H2D2_0001 (2body_dust)
    kall(1179) = 9.95368610d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.40000000d+03*invTd)*1.55644381d+12))*indns

    !H_0001 + C2HD2_0001 -> C2H2D2_gas (2body_gas)
    kall(1180) = 4.63139025d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.40000000d+03*invTd)*1.55644381d+12))*indns

    !D_0001 + C2H3_0001 -> C2H3D_0001 (2body_dust)
    kall(1181) = 9.95067333d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.40000000d+03*invTd)*1.61306016d+12))*indns

    !D_0001 + C2H3_0001 -> C2H3D_gas (2body_gas)
    kall(1182) = 4.93266706d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.40000000d+03*invTd)*1.61306016d+12))*indns

    !D_0001 + C2H2D_0001 -> C2H2D2_0001 (2body_dust)
    kall(1183) = 9.95317015d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.40000000d+03*invTd)*1.58399363d+12))*indns

    !D_0001 + C2H2D_0001 -> C2H2D2_gas (2body_gas)
    kall(1184) = 4.68298492d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.40000000d+03*invTd)*1.58399363d+12))*indns

    !D_0001 + C2HD2_0001 -> C2HD3_0001 (2body_dust)
    kall(1185) = 9.94565045d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.40000000d+03*invTd)*1.55644381d+12))*indns

    !D_0001 + C2HD2_0001 -> C2HD3_gas (2body_gas)
    kall(1186) = 5.43495521d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.40000000d+03*invTd)*1.55644381d+12))*indns

    !H_0001 + C2H4_0001 -> C2H5_0001 (2body_dust)
    kall(1187) = 3.53943710d+12*max(1.72480644d-05, exp(-7.50000000d+02 * invTd))&
        *(3.54191744d+12*max(1.72480644d-05, exp(-7.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.25000000d+03*invTd)*1.49673330d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.25000000d+03*invTd)*1.49673330d+12))*indns

    !H_0001 + C2H4_0001 -> C2H5_gas (2body_gas)
    kall(1188) = 2.48033757d+09*max(1.72480644d-05, exp(-7.50000000d+02 * invTd))&
        *(3.54191744d+12*max(1.72480644d-05, exp(-7.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.25000000d+03*invTd)*1.49673330d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.25000000d+03*invTd)*1.49673330d+12))*indns

    !H_0001 + C2H3D_0001 -> C2H4D_0001 (2body_dust)
    kall(1189) = 3.53665778d+12*max(1.71358611d-05, exp(-7.50000000d+02 * invTd))&
        *(3.54191744d+12*max(1.71358611d-05, exp(-7.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.25000000d+03*invTd)*1.47070117d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.25000000d+03*invTd)*1.47070117d+12))*indns

    !H_0001 + C2HD3_0001 -> C2H2D3_0001 (2body_dust)
    kall(1190) = 3.54191744d+12*max(1.69345173d-05, exp(-7.50000000d+02 * invTd))&
        *(3.54191744d+12*max(1.69345173d-05, exp(-7.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.25000000d+03*invTd)*1.42246827d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.25000000d+03*invTd)*1.42246827d+12))*indns

    !H_0001 + C2H2D2_0001 -> C2H3D2_0001 (2body_dust)
    kall(1191) = 3.54191744d+12*max(1.70316175d-05, exp(-7.50000000d+02 * invTd))&
        *(3.54191744d+12*max(1.70316175d-05, exp(-7.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.55000000d+03*invTd)*1.61017712d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.55000000d+03*invTd)*1.61017712d+12))*indns

    !D_0001 + C2H4_0001 -> C2H4D_0001 (2body_dust)
    kall(1192) = 2.55457933d+12*max(2.38207948d-07, exp(-7.50000000d+02 * invTd))&
        *(2.55656768d+12*max(2.38207948d-07, exp(-7.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.25000000d+03*invTd)*1.49673330d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.25000000d+03*invTd)*1.49673330d+12))*indns

    !D_0001 + C2H3D_0001 -> C2H3D2_0001 (2body_dust)
    kall(1193) = 2.55656768d+12*max(2.34061726d-07, exp(-7.50000000d+02 * invTd))&
        *(2.55656768d+12*max(2.34061726d-07, exp(-7.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.25000000d+03*invTd)*1.47070117d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.25000000d+03*invTd)*1.47070117d+12))*indns

    !D_0001 + C2H2D2_0001 -> C2H2D3_0001 (2body_dust)
    kall(1194) = 2.55656768d+12*max(2.30244429d-07, exp(-7.50000000d+02 * invTd))&
        *(2.55656768d+12*max(2.30244429d-07, exp(-7.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.55000000d+03*invTd)*1.61017712d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.55000000d+03*invTd)*1.61017712d+12))*indns

    !H_0001 + C2H5_0001 -> C2H6_0001 (2body_dust)
    kall(1195) = 9.94567583d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.55000000d+03*invTd)*1.63770351d+12))*indns

    !H_0001 + C2H5_0001 -> C2H6_gas (2body_gas)
    kall(1196) = 5.43241697d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.55000000d+03*invTd)*1.63770351d+12))*indns

    !H_0001 + C2H4D_0001 -> C2H5D_0001 (2body_dust)
    kall(1197) = 9.94567583d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.55000000d+03*invTd)*1.61017712d+12))*indns

    !H_0001 + C2H3D2_0001 -> C2H4D2_0001 (2body_dust)
    kall(1198) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*2.12895557d+12))*indns

    !H_0001 + C2H2D3_0001 -> C2H3D3_0001 (2body_dust)
    kall(1199) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*2.09542661d+12))*indns

    !D_0001 + C2H5_0001 -> C2H5D_0001 (2body_dust)
    kall(1200) = 9.94518861d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.55000000d+03*invTd)*1.63770351d+12))*indns

    !D_0001 + C2H4D_0001 -> C2H4D2_0001 (2body_dust)
    kall(1201) = 9.94518861d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.55000000d+03*invTd)*1.61017712d+12))*indns

    !D_0001 + C2H3D2_0001 -> C2H3D3_0001 (2body_dust)
    kall(1202) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*2.12895557d+12))*indns

    !H_0001 + C2H5D_0001 -> HD_0001 + C2H5_0001 (2body_dust_dust)
    kall(1203) = 1.01298839d+12*max(6.56160128d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.56160128d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))*indns

    !H_0001 + C2H4D2_0001 -> HD_0001 + C2H4D_0001 (2body_dust_dust)
    kall(1204) = 1.68595270d+12*max(6.47227185d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.47227185d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))*indns

    !H_0001 + C2H3D3_0001 -> HD_0001 + C2H3D2_0001 (2body_dust_dust)
    kall(1205) = 2.02243486d+12*max(6.38934738d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.38934738d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))*indns

    !D_0001 + C2H6_0001 -> HD_0001 + C2H5_0001 (2body_dust_dust)
    kall(1206) = 7.31178358d+11*max(1.12414099d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.12414099d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !D_0001 + C2H5D_0001 -> HD_0001 + C2H4D_0001 (2body_dust_dust)
    kall(1207) = 1.21692622d+12*max(1.08070967d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.08070967d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))*indns

    !D_0001 + C2H4D2_0001 -> HD_0001 + C2H3D2_0001 (2body_dust_dust)
    kall(1208) = 1.45980015d+12*max(1.04140490d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.04140490d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))*indns

    !D_0001 + C2H3D3_0001 -> HD_0001 + C2H2D3_0001 (2body_dust_dust)
    kall(1209) = 1.45980015d+12*max(1.00568863d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.00568863d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))*indns

    !H_0001 + CCN_0001 -> HCCN_0001 (2body_dust)
    kall(1210) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.25883049d+12))*indns

    !D_0001 + CCN_0001 -> DCCN_0001 (2body_dust)
    kall(1211) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.25883049d+12))*indns

    !H_0001 + CCO_0001 -> HC2O_0001 (2body_dust)
    kall(1212) = 9.92623000d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-9.75000000d+02*invTd)*1.10596336d+12))*indns

    !H_0001 + CCO_0001 -> HC2O_gas (2body_gas)
    kall(1213) = 7.37699972d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-9.75000000d+02*invTd)*1.10596336d+12))*indns

    !D_0001 + CCO_0001 -> DC2O_0001 (2body_dust)
    kall(1214) = 9.92591249d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-9.75000000d+02*invTd)*1.10596336d+12))*indns

    !D_0001 + CCO_0001 -> DC2O_gas (2body_gas)
    kall(1215) = 7.40875093d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-9.75000000d+02*invTd)*1.10596336d+12))*indns

    !H_0001 + C3_0001 -> l_C3H_0001 (2body_dust)
    kall(1216) = 9.96070170d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !H_0001 + C3_0001 -> l_C3H_gas (2body_gas)
    kall(1217) = 3.92982956d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !H_0001 + C3_0001 -> c_C3H_0001 (2body_dust)
    kall(1218) = 9.96070170d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !H_0001 + C3_0001 -> c_C3H_gas (2body_gas)
    kall(1219) = 3.92982956d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !D_0001 + C3_0001 -> l_C3D_0001 (2body_dust)
    kall(1220) = 9.96003293d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !D_0001 + C3_0001 -> l_C3D_gas (2body_gas)
    kall(1221) = 3.99670711d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !D_0001 + C3_0001 -> c_C3D_0001 (2body_dust)
    kall(1222) = 9.96003293d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !D_0001 + C3_0001 -> c_C3D_gas (2body_gas)
    kall(1223) = 3.99670711d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !H_0001 + l_C3H_0001 -> l_C3H2_0001 (2body_dust)
    kall(1224) = 9.97208228d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))*indns

    !H_0001 + l_C3H_0001 -> l_C3H2_gas (2body_gas)
    kall(1225) = 2.79177175d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))*indns

    !H_0001 + l_C3D_0001 -> l_C3HD_0001 (2body_dust)
    kall(1226) = 9.97208228d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))*indns

    !H_0001 + l_C3D_0001 -> l_C3HD_gas (2body_gas)
    kall(1227) = 2.79177175d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))*indns

    !D_0001 + l_C3H_0001 -> l_C3HD_0001 (2body_dust)
    kall(1228) = 9.97129150d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))*indns

    !D_0001 + l_C3H_0001 -> l_C3HD_gas (2body_gas)
    kall(1229) = 2.87084996d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))*indns

    !D_0001 + l_C3D_0001 -> l_C3D2_0001 (2body_dust)
    kall(1230) = 9.97129150d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))*indns

    !D_0001 + l_C3D_0001 -> l_C3D2_gas (2body_gas)
    kall(1231) = 2.87084996d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))*indns

    !H_0001 + l_C3H2_0001 -> C3H3_0001 (2body_dust)
    kall(1232) = 3.54191744d+12*max(8.36004581d-07, exp(-1.21000000d+03 * invTd))&
        *(3.54191744d+12*max(8.36004581d-07, exp(-1.21000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.69350000d+03*invTd)*1.49544082d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.69350000d+03*invTd)*1.49544082d+12))*indns

    !H_0001 + c_C3H2_0001 -> C3H3_0001 (2body_dust)
    kall(1233) = 3.54191744d+12*max(8.36004581d-07, exp(-1.21000000d+03 * invTd))&
        *(3.54191744d+12*max(8.36004581d-07, exp(-1.21000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.69350000d+03*invTd)*1.49544082d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.69350000d+03*invTd)*1.49544082d+12))*indns

    !H_0001 + l_C3HD_0001 -> C3H2D_0001 (2body_dust)
    kall(1234) = 3.54191744d+12*max(8.32165511d-07, exp(-1.21000000d+03 * invTd))&
        *(3.54191744d+12*max(8.32165511d-07, exp(-1.21000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.69350000d+03*invTd)*1.47614400d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.69350000d+03*invTd)*1.47614400d+12))*indns

    !H_0001 + c_C3HD_0001 -> C3H2D_0001 (2body_dust)
    kall(1235) = 3.54191744d+12*max(8.32165511d-07, exp(-1.21000000d+03 * invTd))&
        *(3.54191744d+12*max(8.32165511d-07, exp(-1.21000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.69350000d+03*invTd)*1.47614400d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.69350000d+03*invTd)*1.47614400d+12))*indns

    !H_0001 + l_C3D2_0001 -> C3HD2_0001 (2body_dust)
    kall(1236) = 3.54191744d+12*max(8.28531238d-07, exp(-1.21000000d+03 * invTd))&
        *(3.54191744d+12*max(8.28531238d-07, exp(-1.21000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.69350000d+03*invTd)*1.45757541d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.69350000d+03*invTd)*1.45757541d+12))*indns

    !H_0001 + c_C3D2_0001 -> C3HD2_0001 (2body_dust)
    kall(1237) = 3.54191744d+12*max(8.28531238d-07, exp(-1.21000000d+03 * invTd))&
        *(3.54191744d+12*max(8.28531238d-07, exp(-1.21000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.69350000d+03*invTd)*1.45757541d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.69350000d+03*invTd)*1.45757541d+12))*indns

    !D_0001 + l_C3H2_0001 -> C3H2D_0001 (2body_dust)
    kall(1238) = 2.55656768d+12*max(3.25706387d-09, exp(-1.21000000d+03 * invTd))&
        *(2.55656768d+12*max(3.25706387d-09, exp(-1.21000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.69350000d+03*invTd)*1.49544082d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.69350000d+03*invTd)*1.49544082d+12))*indns

    !D_0001 + c_C3H2_0001 -> C3H2D_0001 (2body_dust)
    kall(1239) = 2.55656768d+12*max(3.25706387d-09, exp(-1.21000000d+03 * invTd))&
        *(2.55656768d+12*max(3.25706387d-09, exp(-1.21000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.69350000d+03*invTd)*1.49544082d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.69350000d+03*invTd)*1.49544082d+12))*indns

    !D_0001 + l_C3HD_0001 -> C3HD2_0001 (2body_dust)
    kall(1240) = 2.55656768d+12*max(3.21647768d-09, exp(-1.21000000d+03 * invTd))&
        *(2.55656768d+12*max(3.21647768d-09, exp(-1.21000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.69350000d+03*invTd)*1.47614400d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.69350000d+03*invTd)*1.47614400d+12))*indns

    !D_0001 + c_C3HD_0001 -> C3HD2_0001 (2body_dust)
    kall(1241) = 2.55656768d+12*max(3.21647768d-09, exp(-1.21000000d+03 * invTd))&
        *(2.55656768d+12*max(3.21647768d-09, exp(-1.21000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.69350000d+03*invTd)*1.47614400d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.69350000d+03*invTd)*1.47614400d+12))*indns

    !D_0001 + l_C3D2_0001 -> C3D3_0001 (2body_dust)
    kall(1242) = 2.55656768d+12*max(3.17831820d-09, exp(-1.21000000d+03 * invTd))&
        *(2.55656768d+12*max(3.17831820d-09, exp(-1.21000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.69350000d+03*invTd)*1.45757541d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.69350000d+03*invTd)*1.45757541d+12))*indns

    !D_0001 + c_C3D2_0001 -> C3D3_0001 (2body_dust)
    kall(1243) = 2.55656768d+12*max(3.17831820d-09, exp(-1.21000000d+03 * invTd))&
        *(2.55656768d+12*max(3.17831820d-09, exp(-1.21000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.69350000d+03*invTd)*1.45757541d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.69350000d+03*invTd)*1.45757541d+12))*indns

    !H_0001 + C3H3_0001 -> C3H4_0001 (2body_dust)
    kall(1244) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.91850000d+03*invTd)*1.57114775d+12))*indns

    !H_0001 + C3H3_0001 -> C3H4_gas (2body_gas)
    kall(1245) = 0d0*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.91850000d+03*invTd)*1.57114775d+12))*indns

    !H_0001 + C3H2D_0001 -> C3H3D_0001 (2body_dust)
    kall(1246) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.91850000d+03*invTd)*1.55138410d+12))*indns

    !H_0001 + C3HD2_0001 -> C3H2D2_0001 (2body_dust)
    kall(1247) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.91850000d+03*invTd)*1.53234799d+12))*indns

    !D_0001 + C3H3_0001 -> C3H3D_0001 (2body_dust)
    kall(1248) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.91850000d+03*invTd)*1.57114775d+12))*indns

    !D_0001 + C3H2D_0001 -> C3H2D2_0001 (2body_dust)
    kall(1249) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.91850000d+03*invTd)*1.55138410d+12))*indns

    !H_0001 + C3H3N_0001 -> H4C3N_0001 (2body_dust)
    kall(1250) = 3.54191744d+12*max(1.57587951d-05, exp(-7.50000000d+02 * invTd))&
        *(3.54191744d+12*max(1.57587951d-05, exp(-7.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.62820660d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.62820660d+12))*indns

    !H_0001 + C3D3N_0001 -> HD3C3N_0001 (2body_dust)
    kall(1251) = 3.54191744d+12*max(1.56725280d-05, exp(-7.50000000d+02 * invTd))&
        *(3.54191744d+12*max(1.56725280d-05, exp(-7.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.58399363d+12))*indns

    !H_0001 + C3H2DN_0001 -> H3DC3N_0001 (2body_dust)
    kall(1252) = 3.54191744d+12*max(1.57289353d-05, exp(-7.50000000d+02 * invTd))&
        *(3.54191744d+12*max(1.57289353d-05, exp(-7.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.61306016d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.61306016d+12))*indns

    !H_0001 + C3HD2N_0001 -> H2D2C3N_0001 (2body_dust)
    kall(1253) = 3.54191744d+12*max(1.57002004d-05, exp(-7.50000000d+02 * invTd))&
        *(3.54191744d+12*max(1.57002004d-05, exp(-7.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.59832871d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.59832871d+12))*indns

    !D_0001 + C3H3N_0001 -> H3DC3N_0001 (2body_dust)
    kall(1254) = 2.55656768d+12*max(1.86336580d-07, exp(-7.50000000d+02 * invTd))&
        *(2.55656768d+12*max(1.86336580d-07, exp(-7.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*1.62820660d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*1.62820660d+12))*indns

    !D_0001 + C3H2DN_0001 -> H2D2C3N_0001 (2body_dust)
    kall(1255) = 2.55656768d+12*max(1.85366428d-07, exp(-7.50000000d+02 * invTd))&
        *(2.55656768d+12*max(1.85366428d-07, exp(-7.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*1.61306016d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*1.61306016d+12))*indns

    !D_0001 + C3HD2N_0001 -> HD3C3N_0001 (2body_dust)
    kall(1256) = 2.55656768d+12*max(1.84435413d-07, exp(-7.50000000d+02 * invTd))&
        *(2.55656768d+12*max(1.84435413d-07, exp(-7.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*1.59832871d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*1.59832871d+12))*indns

    !H_0001 + C3N_0001 -> HC3N_0001 (2body_dust)
    kall(1257) = 9.94753633d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.60000000d+03*invTd)*1.26719491d+12))*indns

    !H_0001 + C3N_0001 -> HC3N_gas (2body_gas)
    kall(1258) = 5.24636702d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.60000000d+03*invTd)*1.26719491d+12))*indns

    !D_0001 + C3N_0001 -> DC3N_0001 (2body_dust)
    kall(1259) = 9.94716711d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.60000000d+03*invTd)*1.26719491d+12))*indns

    !D_0001 + C3N_0001 -> DC3N_gas (2body_gas)
    kall(1260) = 5.28328851d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.60000000d+03*invTd)*1.26719491d+12))*indns

    !H_0001 + C3O_0001 -> HC3O_0001 (2body_dust)
    kall(1261) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37500000d+03*invTd)*1.15190883d+12))*indns

    !D_0001 + C3O_0001 -> DC3O_0001 (2body_dust)
    kall(1262) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37500000d+03*invTd)*1.15190883d+12))*indns

    !H_0001 + CH_0001 -> CH2_0001 (2body_dust)
    kall(1263) = 9.90893303d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.62500000d+02*invTd)*1.33614202d+12))*indns

    !H_0001 + CH_0001 -> CH2_gas (2body_gas)
    kall(1264) = 9.10669662d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.62500000d+02*invTd)*1.33614202d+12))*indns

    !H_0001 + CD_0001 -> CHD_0001 (2body_dust)
    kall(1265) = 9.90884974d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.62500000d+02*invTd)*1.28753866d+12))*indns

    !H_0001 + CD_0001 -> CHD_gas (2body_gas)
    kall(1266) = 9.11502649d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.62500000d+02*invTd)*1.28753866d+12))*indns

    !D_0001 + CH_0001 -> CHD_0001 (2body_dust)
    kall(1267) = 9.90865457d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.62500000d+02*invTd)*1.33614202d+12))*indns

    !D_0001 + CH_0001 -> CHD_gas (2body_gas)
    kall(1268) = 9.13454281d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.62500000d+02*invTd)*1.33614202d+12))*indns

    !D_0001 + CD_0001 -> CD2_0001 (2body_dust)
    kall(1269) = 9.90865680d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.62500000d+02*invTd)*1.28753866d+12))*indns

    !D_0001 + CD_0001 -> CD2_gas (2body_gas)
    kall(1270) = 9.13431976d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.62500000d+02*invTd)*1.28753866d+12))*indns

    !H_0001 + CH2_0001 -> CH3_0001 (2body_dust)
    kall(1271) = 9.91694716d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !H_0001 + CH2_0001 -> CH3_gas (2body_gas)
    kall(1272) = 8.30528417d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !H_0001 + CHD_0001 -> CH2D_0001 (2body_dust)
    kall(1273) = 9.91693743d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !H_0001 + CHD_0001 -> CH2D_gas (2body_gas)
    kall(1274) = 8.30625737d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !H_0001 + CD2_0001 -> CHD2_0001 (2body_dust)
    kall(1275) = 9.92407896d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !H_0001 + CD2_0001 -> CHD2_gas (2body_gas)
    kall(1276) = 7.59210406d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !D_0001 + CH2_0001 -> CH2D_0001 (2body_dust)
    kall(1277) = 9.91644517d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !D_0001 + CH2_0001 -> CH2D_gas (2body_gas)
    kall(1278) = 8.35548313d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !D_0001 + CHD_0001 -> CHD2_0001 (2body_dust)
    kall(1279) = 9.92359409d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !D_0001 + CHD_0001 -> CHD2_gas (2body_gas)
    kall(1280) = 7.64059107d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !D_0001 + CD2_0001 -> CD3_0001 (2body_dust)
    kall(1281) = 9.91658964d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !D_0001 + CD2_0001 -> CD3_gas (2body_gas)
    kall(1282) = 8.34103600d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !H_0001 + CH2NH_0001 -> CH2NH2_0001 (2body_dust)
    kall(1283) = 4.99820046d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.18813448d+12))*indns

    !H_0001 + CH2NH_0001 -> CH2NH2_gas (2body_gas)
    kall(1284) = 1.79954431d-04*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.18813448d+12))*indns

    !H_0001 + CHDNH_0001 -> CHDNH2_0001 (2body_dust)
    kall(1285) = 2.49910023d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !H_0001 + CHDNH_0001 -> CHDNH2_gas (2body_gas)
    kall(1286) = 8.99772157d-05*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !H_0001 + CHDNH_0001 -> CH2NHD_0001 (2body_dust)
    kall(1287) = 2.49910023d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !H_0001 + CHDNH_0001 -> CH2NHD_gas (2body_gas)
    kall(1288) = 8.99772157d-05*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !H_0001 + CH2ND_0001 -> CHDNH2_0001 (2body_dust)
    kall(1289) = 2.49910023d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !H_0001 + CH2ND_0001 -> CHDNH2_gas (2body_gas)
    kall(1290) = 8.99772157d-05*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !H_0001 + CH2ND_0001 -> CH2NHD_0001 (2body_dust)
    kall(1291) = 2.49910023d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !H_0001 + CH2ND_0001 -> CH2NHD_gas (2body_gas)
    kall(1292) = 8.99772157d-05*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !H_0001 + CD2NH_0001 -> CD2NH2_0001 (2body_dust)
    kall(1293) = 8.39697677d-02*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !H_0001 + CD2NH_0001 -> CD2NH2_gas (2body_gas)
    kall(1294) = 3.02323445d-05*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !H_0001 + CD2NH_0001 -> CH2ND2_0001 (2body_dust)
    kall(1295) = 8.39697677d-02*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !H_0001 + CD2NH_0001 -> CH2ND2_gas (2body_gas)
    kall(1296) = 3.02323445d-05*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !H_0001 + CD2NH_0001 -> CHDNHD_0001 (2body_dust)
    kall(1297) = 3.33879790d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !H_0001 + CD2NH_0001 -> CHDNHD_gas (2body_gas)
    kall(1298) = 1.20209560d-04*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !H_0001 + CHDND_0001 -> CD2NH2_0001 (2body_dust)
    kall(1299) = 8.39697677d-02*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !H_0001 + CHDND_0001 -> CD2NH2_gas (2body_gas)
    kall(1300) = 3.02323445d-05*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !H_0001 + CHDND_0001 -> CH2ND2_0001 (2body_dust)
    kall(1301) = 8.39697677d-02*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !H_0001 + CHDND_0001 -> CH2ND2_gas (2body_gas)
    kall(1302) = 3.02323445d-05*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !H_0001 + CHDND_0001 -> CHDNHD_0001 (2body_dust)
    kall(1303) = 3.33879790d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !H_0001 + CHDND_0001 -> CHDNHD_gas (2body_gas)
    kall(1304) = 1.20209560d-04*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !H_0001 + CD2ND_0001 -> CD2NHD_0001 (2body_dust)
    kall(1305) = 2.49910023d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.08304197d+12))*indns

    !H_0001 + CD2ND_0001 -> CD2NHD_gas (2body_gas)
    kall(1306) = 8.99772157d-05*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.08304197d+12))*indns

    !H_0001 + CD2ND_0001 -> CHDND2_0001 (2body_dust)
    kall(1307) = 2.49910023d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.08304197d+12))*indns

    !H_0001 + CD2ND_0001 -> CHDND2_gas (2body_gas)
    kall(1308) = 8.99772157d-05*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.08304197d+12))*indns

    !D_0001 + CH2NH_0001 -> CHDNH2_0001 (2body_dust)
    kall(1309) = 2.49898966d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.18813448d+12))*indns

    !D_0001 + CH2NH_0001 -> CHDNH2_gas (2body_gas)
    kall(1310) = 1.01033975d-04*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.18813448d+12))*indns

    !D_0001 + CH2NH_0001 -> CH2NHD_0001 (2body_dust)
    kall(1311) = 2.49898966d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.18813448d+12))*indns

    !D_0001 + CH2NH_0001 -> CH2NHD_gas (2body_gas)
    kall(1312) = 1.01033975d-04*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.18813448d+12))*indns

    !D_0001 + CHDNH_0001 -> CD2NH2_0001 (2body_dust)
    kall(1313) = 8.39660526d-02*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !D_0001 + CHDNH_0001 -> CD2NH2_gas (2body_gas)
    kall(1314) = 3.39474154d-05*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !D_0001 + CHDNH_0001 -> CH2ND2_0001 (2body_dust)
    kall(1315) = 8.39660526d-02*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !D_0001 + CHDNH_0001 -> CH2ND2_gas (2body_gas)
    kall(1316) = 3.39474154d-05*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !D_0001 + CHDNH_0001 -> CHDNHD_0001 (2body_dust)
    kall(1317) = 3.33865019d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !D_0001 + CHDNH_0001 -> CHDNHD_gas (2body_gas)
    kall(1318) = 1.34981390d-04*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !D_0001 + CH2ND_0001 -> CD2NH2_0001 (2body_dust)
    kall(1319) = 8.39660526d-02*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !D_0001 + CH2ND_0001 -> CD2NH2_gas (2body_gas)
    kall(1320) = 3.39474154d-05*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !D_0001 + CH2ND_0001 -> CH2ND2_0001 (2body_dust)
    kall(1321) = 8.39660526d-02*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !D_0001 + CH2ND_0001 -> CH2ND2_gas (2body_gas)
    kall(1322) = 3.39474154d-05*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !D_0001 + CH2ND_0001 -> CHDNHD_0001 (2body_dust)
    kall(1323) = 3.33865019d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !D_0001 + CH2ND_0001 -> CHDNHD_gas (2body_gas)
    kall(1324) = 1.34981390d-04*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !D_0001 + CD2NH_0001 -> CD2NHD_0001 (2body_dust)
    kall(1325) = 2.49898966d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !D_0001 + CD2NH_0001 -> CD2NHD_gas (2body_gas)
    kall(1326) = 1.01033975d-04*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !D_0001 + CD2NH_0001 -> CHDND2_0001 (2body_dust)
    kall(1327) = 2.49898966d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !D_0001 + CD2NH_0001 -> CHDND2_gas (2body_gas)
    kall(1328) = 1.01033975d-04*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !D_0001 + CHDND_0001 -> CD2NHD_0001 (2body_dust)
    kall(1329) = 2.49898966d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !D_0001 + CHDND_0001 -> CD2NHD_gas (2body_gas)
    kall(1330) = 1.01033975d-04*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !D_0001 + CHDND_0001 -> CHDND2_0001 (2body_dust)
    kall(1331) = 2.49898966d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !D_0001 + CHDND_0001 -> CHDND2_gas (2body_gas)
    kall(1332) = 1.01033975d-04*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.76700000d+03*invTd)*2.11637275d+12))*indns

    !H_0001 + CH2NH_0001 -> CH3NH_0001 (2body_dust)
    kall(1333) = 5.00000000d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.18813448d+12))*indns

    !H_0001 + CH2NH2_0001 -> CH5N_0001 (2body_dust)
    kall(1334) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !H_0001 + CH3_0001 -> CH4_0001 (2body_dust)
    kall(1335) = 9.91613950d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !H_0001 + CH3_0001 -> CH4_gas (2body_gas)
    kall(1336) = 8.38604988d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !H_0001 + CD3_0001 -> CHD3_0001 (2body_dust)
    kall(1337) = 9.91705663d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !H_0001 + CD3_0001 -> CHD3_gas (2body_gas)
    kall(1338) = 8.29433711d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !H_0001 + CH2D_0001 -> CH3D_0001 (2body_dust)
    kall(1339) = 9.91635718d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !H_0001 + CH2D_0001 -> CH3D_gas (2body_gas)
    kall(1340) = 8.36428221d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !H_0001 + CHD2_0001 -> CH2D2_0001 (2body_dust)
    kall(1341) = 9.91669914d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !H_0001 + CHD2_0001 -> CH2D2_gas (2body_gas)
    kall(1342) = 8.33008573d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !D_0001 + CH3_0001 -> CH3D_0001 (2body_dust)
    kall(1343) = 9.91584095d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !D_0001 + CH3_0001 -> CH3D_gas (2body_gas)
    kall(1344) = 8.41590505d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !D_0001 + CH2D_0001 -> CH2D2_0001 (2body_dust)
    kall(1345) = 9.91632493d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !D_0001 + CH2D_0001 -> CH2D2_gas (2body_gas)
    kall(1346) = 8.36750710d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !D_0001 + CHD2_0001 -> CHD3_0001 (2body_dust)
    kall(1347) = 9.91666545d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !D_0001 + CHD2_0001 -> CHD3_gas (2body_gas)
    kall(1348) = 8.33345520d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !H_0001 + CH3NH_0001 -> CH5N_0001 (2body_dust)
    kall(1349) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.76700000d+03*invTd)*2.15135649d+12))*indns

    !H_0001 + CH3D_0001 -> HD_0001 + CH3_0001 (2body_dust_dust)
    kall(1350) = 1.41676697d+12*max(1.69956299d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.69956299d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))*indns

    !H_0001 + CHD3_0001 -> HD_0001 + CHD2_0001 (2body_dust_dust)
    kall(1351) = 2.12515046d+12*max(1.39908687d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.39908687d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))*indns

    !H_0001 + CH2D2_0001 -> HD_0001 + CH2D_0001 (2body_dust_dust)
    kall(1352) = 2.12515046d+12*max(1.53403970d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.53403970d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))*indns

    !D_0001 + CH4_0001 -> HD_0001 + CH3_0001 (2body_dust_dust)
    kall(1353) = 1.02262707d+12*exp(-5.94000000d+03 * invTd)&
        *(2.55656768d+12*exp(-5.94000000d+03 * invTd) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.22695619d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.22695619d+12))*indns

    !D_0001 + CH3D_0001 -> HD_0001 + CH2D_0001 (2body_dust_dust)
    kall(1354) = 1.53394061d+12*exp(-5.94000000d+03 * invTd)&
        *(2.55656768d+12*exp(-5.94000000d+03 * invTd) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))*indns

    !D_0001 + CHD3_0001 -> HD_0001 + CD3_0001 (2body_dust_dust)
    kall(1355) = 1.02262707d+12*exp(-5.94000000d+03 * invTd)&
        *(2.55656768d+12*exp(-5.94000000d+03 * invTd) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))*indns

    !D_0001 + CH2D2_0001 -> HD_0001 + CHD2_0001 (2body_dust_dust)
    kall(1356) = 1.53394061d+12*exp(-5.94000000d+03 * invTd)&
        *(2.55656768d+12*exp(-5.94000000d+03 * invTd) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))*indns

    !H_0001 + CN_0001 -> HCN_0001 (2body_dust)
    kall(1357) = 9.91735428d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !H_0001 + CN_0001 -> HCN_gas (2body_gas)
    kall(1358) = 8.26457231d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !D_0001 + CN_0001 -> DCN_0001 (2body_dust)
    kall(1359) = 9.91714357d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !D_0001 + CN_0001 -> DCN_gas (2body_gas)
    kall(1360) = 8.28564286d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !H_0001 + CS_0001 -> HCS_0001 (2body_dust)
    kall(1361) = 3.51908031d+12*max(2.91776415d-06, exp(-1.00000000d+03 * invTd))&
        *(3.54191744d+12*max(2.91776415d-06, exp(-1.00000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.60000000d+03*invTd)*1.35083431d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.60000000d+03*invTd)*1.35083431d+12))*indns

    !H_0001 + CS_0001 -> HCS_gas (2body_gas)
    kall(1362) = 2.28371238d+10*max(2.91776415d-06, exp(-1.00000000d+03 * invTd))&
        *(3.54191744d+12*max(2.91776415d-06, exp(-1.00000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.60000000d+03*invTd)*1.35083431d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.60000000d+03*invTd)*1.35083431d+12))*indns

    !D_0001 + CS_0001 -> DCS_0001 (2body_dust)
    kall(1363) = 2.53983772d+12*max(1.81112750d-08, exp(-1.00000000d+03 * invTd))&
        *(2.55656768d+12*max(1.81112750d-08, exp(-1.00000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.60000000d+03*invTd)*1.35083431d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.60000000d+03*invTd)*1.35083431d+12))*indns

    !D_0001 + CS_0001 -> DCS_gas (2body_gas)
    kall(1364) = 1.67299681d+10*max(1.81112750d-08, exp(-1.00000000d+03 * invTd))&
        *(2.55656768d+12*max(1.81112750d-08, exp(-1.00000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.60000000d+03*invTd)*1.35083431d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.60000000d+03*invTd)*1.35083431d+12))*indns

    !H_0001 + Fe_0001 -> FeH_0001 (2body_dust)
    kall(1365) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.10000000d+03*invTd)*1.37177872d+12))*indns

    !H_0001 + Fe_0001 -> FeH_gas (2body_gas)
    kall(1366) = 0d0*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.10000000d+03*invTd)*1.37177872d+12))*indns

    !H_0001 + D_0001 -> HD_0001 (2body_dust)
    kall(1367) = 9.90182749d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.60500000d+02*invTd)*2.55656768d+12))*indns

    !H_0001 + D_0001 -> HD_gas (2body_gas)
    kall(1368) = 9.81725126d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.60500000d+02*invTd)*2.55656768d+12))*indns

    !H_0001 + H2C3N_0001 -> C3H3N_0001 (2body_dust)
    kall(1369) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.64378788d+12))*indns

    !H_0001 + HDC3N_0001 -> C3H2DN_0001 (2body_dust)
    kall(1370) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.62820660d+12))*indns

    !H_0001 + D2C3N_0001 -> C3HD2N_0001 (2body_dust)
    kall(1371) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.61306016d+12))*indns

    !D_0001 + H2C3N_0001 -> C3H2DN_0001 (2body_dust)
    kall(1372) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*1.64378788d+12))*indns

    !D_0001 + HDC3N_0001 -> C3HD2N_0001 (2body_dust)
    kall(1373) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*1.62820660d+12))*indns

    !D_0001 + D2C3N_0001 -> C3D3N_0001 (2body_dust)
    kall(1374) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*1.61306016d+12))*indns

    !H_0001 + H2CN_0001 -> CH2NH_0001 (2body_dust)
    kall(1375) = 9.97589141d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.46649314d+12))*indns

    !H_0001 + H2CN_0001 -> CH2NH_gas (2body_gas)
    kall(1376) = 2.41085884d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.46649314d+12))*indns

    !H_0001 + HDCN_0001 -> CHDNH_0001 (2body_dust)
    kall(1377) = 6.65391957d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !H_0001 + HDCN_0001 -> CHDNH_gas (2body_gas)
    kall(1378) = 1.60804285d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !H_0001 + HDCN_0001 -> CH2ND_0001 (2body_dust)
    kall(1379) = 3.32197184d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !H_0001 + HDCN_0001 -> CH2ND_gas (2body_gas)
    kall(1380) = 8.02815994d-04*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !H_0001 + D2CN_0001 -> CD2NH_0001 (2body_dust)
    kall(1381) = 3.32197184d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !H_0001 + D2CN_0001 -> CD2NH_gas (2body_gas)
    kall(1382) = 8.02815994d-04*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !H_0001 + D2CN_0001 -> CHDND_0001 (2body_dust)
    kall(1383) = 6.65391957d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !H_0001 + D2CN_0001 -> CHDND_gas (2body_gas)
    kall(1384) = 1.60804285d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !D_0001 + H2CN_0001 -> CHDNH_0001 (2body_dust)
    kall(1385) = 6.65345531d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.46649314d+12))*indns

    !D_0001 + H2CN_0001 -> CHDNH_gas (2body_gas)
    kall(1386) = 1.65446909d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.46649314d+12))*indns

    !D_0001 + H2CN_0001 -> CH2ND_0001 (2body_dust)
    kall(1387) = 3.32174006d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.46649314d+12))*indns

    !D_0001 + H2CN_0001 -> CH2ND_gas (2body_gas)
    kall(1388) = 8.25994313d-04*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.46649314d+12))*indns

    !D_0001 + HDCN_0001 -> CD2NH_0001 (2body_dust)
    kall(1389) = 3.32174006d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !D_0001 + HDCN_0001 -> CD2NH_gas (2body_gas)
    kall(1390) = 8.25994313d-04*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !D_0001 + HDCN_0001 -> CHDND_0001 (2body_dust)
    kall(1391) = 6.65345531d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !D_0001 + HDCN_0001 -> CHDND_gas (2body_gas)
    kall(1392) = 1.65446909d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !D_0001 + D2CN_0001 -> CD2ND_0001 (2body_dust)
    kall(1393) = 9.97519537d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !D_0001 + D2CN_0001 -> CD2ND_gas (2body_gas)
    kall(1394) = 2.48046340d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !H_0001 + CH2OH_0001 -> CH3OH_0001 (2body_dust)
    kall(1395) = 9.97369765d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !H_0001 + CH2OH_0001 -> CH3OH_gas (2body_gas)
    kall(1396) = 2.63023458d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !H_0001 + CHDOH_0001 -> CH2DOH_0001 (2body_dust)
    kall(1397) = 9.97074008d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !H_0001 + CHDOH_0001 -> CH2DOH_gas (2body_gas)
    kall(1398) = 2.92599240d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !H_0001 + CH2OD_0001 -> CH3OD_0001 (2body_dust)
    kall(1399) = 9.97328967d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !H_0001 + CH2OD_0001 -> CH3OD_gas (2body_gas)
    kall(1400) = 2.67103338d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !H_0001 + CD2OH_0001 -> CHD2OH_0001 (2body_dust)
    kall(1401) = 9.97162838d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !H_0001 + CD2OH_0001 -> CHD2OH_gas (2body_gas)
    kall(1402) = 2.83716174d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !H_0001 + CHDOD_0001 -> CH2DOD_0001 (2body_dust)
    kall(1403) = 9.97162838d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !H_0001 + CHDOD_0001 -> CH2DOD_gas (2body_gas)
    kall(1404) = 2.83716174d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !H_0001 + CD2OD_0001 -> CHD2OD_0001 (2body_dust)
    kall(1405) = 9.97253369d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.80193994d+12))*indns

    !H_0001 + CD2OD_0001 -> CHD2OD_gas (2body_gas)
    kall(1406) = 2.74663077d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.80193994d+12))*indns

    !D_0001 + CH2OH_0001 -> CH2DOH_0001 (2body_dust)
    kall(1407) = 9.97223615d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !D_0001 + CH2OH_0001 -> CH2DOH_gas (2body_gas)
    kall(1408) = 2.77638532d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !D_0001 + CHDOH_0001 -> CHD2OH_0001 (2body_dust)
    kall(1409) = 9.97065516d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !D_0001 + CHDOH_0001 -> CHD2OH_gas (2body_gas)
    kall(1410) = 2.93448357d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !D_0001 + CH2OD_0001 -> CH2DOD_0001 (2body_dust)
    kall(1411) = 9.97269240d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !D_0001 + CH2OD_0001 -> CH2DOD_gas (2body_gas)
    kall(1412) = 2.73075960d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !D_0001 + CD2OH_0001 -> CD3OH_0001 (2body_dust)
    kall(1413) = 9.97154183d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !D_0001 + CD2OH_0001 -> CD3OH_gas (2body_gas)
    kall(1414) = 2.84581692d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !D_0001 + CHDOD_0001 -> CHD2OD_0001 (2body_dust)
    kall(1415) = 9.97154183d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !D_0001 + CHDOD_0001 -> CHD2OD_gas (2body_gas)
    kall(1416) = 2.84581692d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !H_0001 + CH3O_0001 -> CH3OH_0001 (2body_dust)
    kall(1417) = 9.97017090d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !H_0001 + CH3O_0001 -> CH3OH_gas (2body_gas)
    kall(1418) = 2.98291023d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !H_0001 + CH2DO_0001 -> CH2DOH_0001 (2body_dust)
    kall(1419) = 9.97278104d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !H_0001 + CH2DO_0001 -> CH2DOH_gas (2body_gas)
    kall(1420) = 2.72189569d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !H_0001 + CHD2O_0001 -> CHD2OH_0001 (2body_dust)
    kall(1421) = 9.97162838d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !H_0001 + CHD2O_0001 -> CHD2OH_gas (2body_gas)
    kall(1422) = 2.83716174d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !H_0001 + CD3O_0001 -> CD3OH_0001 (2body_dust)
    kall(1423) = 9.97253369d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.80193994d+12))*indns

    !H_0001 + CD3O_0001 -> CD3OH_gas (2body_gas)
    kall(1424) = 2.74663077d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.80193994d+12))*indns

    !D_0001 + CH3O_0001 -> CH3OD_0001 (2body_dust)
    kall(1425) = 9.96928103d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !D_0001 + CH3O_0001 -> CH3OD_gas (2body_gas)
    kall(1426) = 3.07189715d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !D_0001 + CH2DO_0001 -> CH2DOD_0001 (2body_dust)
    kall(1427) = 9.97269240d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !D_0001 + CH2DO_0001 -> CH2DOD_gas (2body_gas)
    kall(1428) = 2.73075960d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !D_0001 + CHD2O_0001 -> CHD2OD_0001 (2body_dust)
    kall(1429) = 9.97154183d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !D_0001 + CHD2O_0001 -> CHD2OD_gas (2body_gas)
    kall(1430) = 2.84581692d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !D_0001 + CD2OD_0001 -> CD3OD_0001 (2body_dust)
    kall(1431) = 9.97244550d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.80193994d+12))*indns

    !D_0001 + CD2OD_0001 -> CD3OD_gas (2body_gas)
    kall(1432) = 2.75545032d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.80193994d+12))*indns

    !D_0001 + CD3O_0001 -> CD3OD_0001 (2body_dust)
    kall(1433) = 9.97244550d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.80193994d+12))*indns

    !D_0001 + CD3O_0001 -> CD3OD_gas (2body_gas)
    kall(1434) = 2.75545032d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.80193994d+12))*indns

    !H_0001 + HOOH_0001 -> OH_0001 + H2O_0001 (2body_dust_dust)
    kall(1435) = 1.77095872d+12*max(2.96658509d-07, exp(-1.40000000d+03 * invTd))&
        *(3.54191744d+12*max(2.96658509d-07, exp(-1.40000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))*indns

    !H_0001 + HOOD_0001 -> H2O_0001 + OD_0001 (2body_dust_dust)
    kall(1436) = 5.91500212d+11*max(2.94842980d-07, exp(-1.40000000d+03 * invTd))&
        *(3.54191744d+12*max(2.94842980d-07, exp(-1.40000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))*indns

    !H_0001 + HOOD_0001 -> OH_0001 + HDO_0001 (2body_dust_dust)
    kall(1437) = 1.17945851d+12*max(2.94842980d-07, exp(-1.40000000d+03 * invTd))&
        *(3.54191744d+12*max(2.94842980d-07, exp(-1.40000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))*indns

    !H_0001 + DOOD_0001 -> OD_0001 + HDO_0001 (2body_dust_dust)
    kall(1438) = 1.17945851d+12*max(2.93136491d-07, exp(-1.40000000d+03 * invTd))&
        *(3.54191744d+12*max(2.93136491d-07, exp(-1.40000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))*indns

    !H_0001 + DOOD_0001 -> OH_0001 + D2O_0001 (2body_dust_dust)
    kall(1439) = 5.91500212d+11*max(2.93136491d-07, exp(-1.40000000d+03 * invTd))&
        *(3.54191744d+12*max(2.93136491d-07, exp(-1.40000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))*indns

    !D_0001 + HOOH_0001 -> H2O_0001 + OD_0001 (2body_dust_dust)
    kall(1440) = 4.26946803d+11*max(7.89763421d-10, exp(-1.40000000d+03 * invTd))&
        *(2.55656768d+12*max(7.89763421d-10, exp(-1.40000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))*indns

    !D_0001 + HOOH_0001 -> OH_0001 + HDO_0001 (2body_dust_dust)
    kall(1441) = 8.51337039d+11*max(7.89763421d-10, exp(-1.40000000d+03 * invTd))&
        *(2.55656768d+12*max(7.89763421d-10, exp(-1.40000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))*indns

    !D_0001 + HOOD_0001 -> OD_0001 + HDO_0001 (2body_dust_dust)
    kall(1442) = 8.51337039d+11*max(7.76719475d-10, exp(-1.40000000d+03 * invTd))&
        *(2.55656768d+12*max(7.76719475d-10, exp(-1.40000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))*indns

    !D_0001 + HOOD_0001 -> OH_0001 + D2O_0001 (2body_dust_dust)
    kall(1443) = 4.26946803d+11*max(7.76719475d-10, exp(-1.40000000d+03 * invTd))&
        *(2.55656768d+12*max(7.76719475d-10, exp(-1.40000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))*indns

    !D_0001 + DOOD_0001 -> OD_0001 + D2O_0001 (2body_dust_dust)
    kall(1444) = 1.27828384d+12*max(7.64570156d-10, exp(-1.40000000d+03 * invTd))&
        *(2.55656768d+12*max(7.64570156d-10, exp(-1.40000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))*indns

    !H_0001 + HOOD_0001 -> HD_0001 + O2H_0001 (2body_dust_dust)
    kall(1445) = 1.17945851d+12*max(2.46776330d-08, exp(-1.90000000d+03 * invTd))&
        *(3.54191744d+12*max(2.46776330d-08, exp(-1.90000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))*indns

    !H_0001 + DOOD_0001 -> HD_0001 + O2D_0001 (2body_dust_dust)
    kall(1446) = 1.17945851d+12*max(2.45113219d-08, exp(-1.90000000d+03 * invTd))&
        *(3.54191744d+12*max(2.45113219d-08, exp(-1.90000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))*indns

    !D_0001 + HOOH_0001 -> HD_0001 + O2H_0001 (2body_dust_dust)
    kall(1447) = 8.51337039d+11*max(2.48831255d-11, exp(-1.90000000d+03 * invTd))&
        *(2.55656768d+12*max(2.48831255d-11, exp(-1.90000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))*indns

    !D_0001 + HOOD_0001 -> HD_0001 + O2D_0001 (2body_dust_dust)
    kall(1448) = 8.51337039d+11*max(2.44050077d-11, exp(-1.90000000d+03 * invTd))&
        *(2.55656768d+12*max(2.44050077d-11, exp(-1.90000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))*indns

    !H_0001 + HDS_0001 -> HD_0001 + HS_0001 (2body_dust_dust)
    kall(1449) = 2.36245893d+12*max(1.27796071d-07, exp(-1.56000000d+03 * invTd))&
        *(3.54191744d+12*max(1.27796071d-07, exp(-1.56000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.40227216d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.40227216d+12))*indns

    !H_0001 + D2S_0001 -> HD_0001 + DS_0001 (2body_dust_dust)
    kall(1450) = 2.36245893d+12*max(1.27015418d-07, exp(-1.56000000d+03 * invTd))&
        *(3.54191744d+12*max(1.27015418d-07, exp(-1.56000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.38265900d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.38265900d+12))*indns

    !D_0001 + H2S_0001 -> HD_0001 + HS_0001 (2body_dust_dust)
    kall(1451) = 1.70523065d+12*max(2.46276056d-10, exp(-1.56000000d+03 * invTd))&
        *(2.55656768d+12*max(2.46276056d-10, exp(-1.56000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.42274437d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.42274437d+12))*indns

    !D_0001 + HDS_0001 -> HD_0001 + DS_0001 (2body_dust_dust)
    kall(1452) = 1.70523065d+12*max(2.41984331d-10, exp(-1.56000000d+03 * invTd))&
        *(2.55656768d+12*max(2.41984331d-10, exp(-1.56000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.40227216d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.40227216d+12))*indns

    !H_0001 + H4C3N_0001 -> H5C3N_0001 (2body_dust)
    kall(1453) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.61306016d+12))*indns

    !H_0001 + HC3N_0001 -> H2C3N_0001 (2body_dust)
    kall(1454) = 3.54191744d+12*max(7.98451788d-07, exp(-1.21000000d+03 * invTd))&
        *(3.54191744d+12*max(7.98451788d-07, exp(-1.21000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.29000000d+03*invTd)*1.50107047d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.29000000d+03*invTd)*1.50107047d+12))*indns

    !H_0001 + DC3N_0001 -> HDC3N_0001 (2body_dust)
    kall(1455) = 3.54191744d+12*max(7.96380912d-07, exp(-1.21000000d+03 * invTd))&
        *(3.54191744d+12*max(7.96380912d-07, exp(-1.21000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.29000000d+03*invTd)*1.48656704d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.29000000d+03*invTd)*1.48656704d+12))*indns

    !D_0001 + HC3N_0001 -> HDC3N_0001 (2body_dust)
    kall(1456) = 2.55656768d+12*max(2.87222729d-09, exp(-1.21000000d+03 * invTd))&
        *(2.55656768d+12*max(2.87222729d-09, exp(-1.21000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.29000000d+03*invTd)*1.50107047d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.29000000d+03*invTd)*1.50107047d+12))*indns

    !D_0001 + DC3N_0001 -> D2C3N_0001 (2body_dust)
    kall(1457) = 2.55656768d+12*max(2.85179154d-09, exp(-1.21000000d+03 * invTd))&
        *(2.55656768d+12*max(2.85179154d-09, exp(-1.21000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.29000000d+03*invTd)*1.48656704d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.29000000d+03*invTd)*1.48656704d+12))*indns

    !H_0001 + HC3O_0001 -> H2C3O_0001 (2body_dust)
    kall(1458) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.62820660d+12))*indns

    !H_0001 + DC3O_0001 -> HDC3O_0001 (2body_dust)
    kall(1459) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.80000000d+03*invTd)*1.61306016d+12))*indns

    !D_0001 + HC3O_0001 -> HDC3O_0001 (2body_dust)
    kall(1460) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*1.62820660d+12))*indns

    !D_0001 + DC3O_0001 -> D2C3O_0001 (2body_dust)
    kall(1461) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.80000000d+03*invTd)*1.61306016d+12))*indns

    !H_0001 + HCO_0001 -> H2CO_0001 (2body_dust)
    kall(1462) = 9.94817620d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !H_0001 + HCO_0001 -> H2CO_gas (2body_gas)
    kall(1463) = 5.18237955d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !H_0001 + DCO_0001 -> HDCO_0001 (2body_dust)
    kall(1464) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !H_0001 + DCO_0001 -> HDCO_gas (2body_gas)
    kall(1465) = 0d0*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !D_0001 + HCO_0001 -> HDCO_0001 (2body_dust)
    kall(1466) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !D_0001 + HCO_0001 -> HDCO_gas (2body_gas)
    kall(1467) = 0d0*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !D_0001 + DCO_0001 -> D2CO_0001 (2body_dust)
    kall(1468) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !D_0001 + DCO_0001 -> D2CO_gas (2body_gas)
    kall(1469) = 0d0*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !H_0001 + HCS_0001 -> H2CS_0001 (2body_dust)
    kall(1470) = 9.94347123d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.45000000d+03*invTd)*1.27158727d+12))*indns

    !H_0001 + HCS_0001 -> H2CS_gas (2body_gas)
    kall(1471) = 5.65287660d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.45000000d+03*invTd)*1.27158727d+12))*indns

    !H_0001 + DCS_0001 -> HDCS_0001 (2body_dust)
    kall(1472) = 9.94347123d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.45000000d+03*invTd)*1.25768973d+12))*indns

    !H_0001 + DCS_0001 -> HDCS_gas (2body_gas)
    kall(1473) = 5.65287660d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.45000000d+03*invTd)*1.25768973d+12))*indns

    !D_0001 + HCS_0001 -> HDCS_0001 (2body_dust)
    kall(1474) = 9.94298895d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.45000000d+03*invTd)*1.27158727d+12))*indns

    !D_0001 + HCS_0001 -> HDCS_gas (2body_gas)
    kall(1475) = 5.70110517d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.45000000d+03*invTd)*1.27158727d+12))*indns

    !D_0001 + DCS_0001 -> D2CS_0001 (2body_dust)
    kall(1476) = 9.94298895d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.45000000d+03*invTd)*1.25768973d+12))*indns

    !D_0001 + DCS_0001 -> D2CS_gas (2body_gas)
    kall(1477) = 5.70110517d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.45000000d+03*invTd)*1.25768973d+12))*indns

    !H_0001 + HS_0001 -> H2S_0001 (2body_dust)
    kall(1478) = 9.91784275d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.35000000d+03*invTd)*1.43277615d+12))*indns

    !H_0001 + HS_0001 -> H2S_gas (2body_gas)
    kall(1479) = 8.21572548d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.35000000d+03*invTd)*1.43277615d+12))*indns

    !H_0001 + DS_0001 -> HDS_0001 (2body_dust)
    kall(1480) = 9.91784275d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.35000000d+03*invTd)*1.41154866d+12))*indns

    !H_0001 + DS_0001 -> HDS_gas (2body_gas)
    kall(1481) = 8.21572548d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.35000000d+03*invTd)*1.41154866d+12))*indns

    !D_0001 + HS_0001 -> HDS_0001 (2body_dust)
    kall(1482) = 9.91759448d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.35000000d+03*invTd)*1.43277615d+12))*indns

    !D_0001 + HS_0001 -> HDS_gas (2body_gas)
    kall(1483) = 8.24055222d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.35000000d+03*invTd)*1.43277615d+12))*indns

    !D_0001 + DS_0001 -> D2S_0001 (2body_dust)
    kall(1484) = 9.91759448d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.35000000d+03*invTd)*1.41154866d+12))*indns

    !D_0001 + DS_0001 -> D2S_gas (2body_gas)
    kall(1485) = 8.24055222d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.35000000d+03*invTd)*1.41154866d+12))*indns

    !H_0001 + Mg_0001 -> MgH_0001 (2body_dust)
    kall(1486) = 9.93672611d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.65000000d+03*invTd)*2.35388773d+12))*indns

    !H_0001 + Mg_0001 -> MgH_gas (2body_gas)
    kall(1487) = 6.32738921d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.65000000d+03*invTd)*2.35388773d+12))*indns

    !H_0001 + MgH_0001 -> MgH2_0001 (2body_dust)
    kall(1488) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.87500000d+03*invTd)*2.40224528d+12))*indns

    !H_0001 + MgH_0001 -> MgH2_gas (2body_gas)
    kall(1489) = 0d0*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.87500000d+03*invTd)*2.40224528d+12))*indns

    !H_0001 + N_0001 -> NH_0001 (2body_dust)
    kall(1490) = 9.90777655d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.60000000d+02*invTd)*1.13594070d+12))*indns

    !H_0001 + N_0001 -> NH_gas (2body_gas)
    kall(1491) = 9.22234528d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.60000000d+02*invTd)*1.13594070d+12))*indns

    !D_0001 + N_0001 -> ND_0001 (2body_dust)
    kall(1492) = 9.90762326d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.60000000d+02*invTd)*1.13594070d+12))*indns

    !D_0001 + N_0001 -> ND_gas (2body_gas)
    kall(1493) = 9.23767383d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.60000000d+02*invTd)*1.13594070d+12))*indns

    !H_0001 + N2HD_0001 -> H_0001 + HD_0001 + N2_0001 (2body_dust_dust)
    kall(1494) = 2.36245893d+12*max(3.61600669d-05, exp(-6.50000000d+02 * invTd))&
        *(3.54191744d+12*max(3.61600669d-05, exp(-6.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.96197517d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.96197517d+12))*indns

    !H_0001 + N2D2_0001 -> D_0001 + HD_0001 + N2_0001 (2body_dust_dust)
    kall(1495) = 2.36245893d+12*max(3.59798047d-05, exp(-6.50000000d+02 * invTd))&
        *(3.54191744d+12*max(3.59798047d-05, exp(-6.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.93107599d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.93107599d+12))*indns

    !D_0001 + N2H2_0001 -> H_0001 + HD_0001 + N2_0001 (2body_dust_dust)
    kall(1496) = 1.70523065d+12*max(6.61506897d-07, exp(-6.50000000d+02 * invTd))&
        *(2.55656768d+12*max(6.61506897d-07, exp(-6.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.99440671d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.99440671d+12))*indns

    !D_0001 + N2HD_0001 -> D_0001 + HD_0001 + N2_0001 (2body_dust_dust)
    kall(1497) = 1.70523065d+12*max(6.52072134d-07, exp(-6.50000000d+02 * invTd))&
        *(2.55656768d+12*max(6.52072134d-07, exp(-6.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.96197517d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.96197517d+12))*indns

    !H_0001 + Na_0001 -> NaH_0001 (2body_dust)
    kall(1498) = 9.95178094d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-5.90000000d+03*invTd)*3.58781893d+12))*indns

    !H_0001 + Na_0001 -> NaH_gas (2body_gas)
    kall(1499) = 4.82190604d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-5.90000000d+03*invTd)*3.58781893d+12))*indns

    !H_0001 + NH_0001 -> NH2_0001 (2body_dust)
    kall(1500) = 9.91921830d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !H_0001 + NH_0001 -> NH2_gas (2body_gas)
    kall(1501) = 8.07816992d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !H_0001 + ND_0001 -> NHD_0001 (2body_dust)
    kall(1502) = 9.91916259d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !H_0001 + ND_0001 -> NHD_gas (2body_gas)
    kall(1503) = 8.08374081d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !D_0001 + NH_0001 -> NHD_0001 (2body_dust)
    kall(1504) = 9.91886395d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !D_0001 + NH_0001 -> NHD_gas (2body_gas)
    kall(1505) = 8.11360498d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !D_0001 + ND_0001 -> ND2_0001 (2body_dust)
    kall(1506) = 9.91881523d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !D_0001 + ND_0001 -> ND2_gas (2body_gas)
    kall(1507) = 8.11847669d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !H_0001 + NH2_0001 -> NH3_0001 (2body_dust)
    kall(1508) = 9.94751893d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !H_0001 + NH2_0001 -> NH3_gas (2body_gas)
    kall(1509) = 5.24810678d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !H_0001 + NHD_0001 -> NH2D_0001 (2body_dust)
    kall(1510) = 9.94819078d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !H_0001 + NHD_0001 -> NH2D_gas (2body_gas)
    kall(1511) = 5.18092158d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !H_0001 + ND2_0001 -> NHD2_0001 (2body_dust)
    kall(1512) = 9.94879429d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !H_0001 + ND2_0001 -> NHD2_gas (2body_gas)
    kall(1513) = 5.12057111d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !D_0001 + NH2_0001 -> NH2D_0001 (2body_dust)
    kall(1514) = 9.94751691d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !D_0001 + NH2_0001 -> NH2D_gas (2body_gas)
    kall(1515) = 5.24830926d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !D_0001 + NHD_0001 -> NHD2_0001 (2body_dust)
    kall(1516) = 9.94811359d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !D_0001 + NHD_0001 -> NHD2_gas (2body_gas)
    kall(1517) = 5.18864083d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !D_0001 + ND2_0001 -> ND3_0001 (2body_dust)
    kall(1518) = 9.94871521d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !D_0001 + ND2_0001 -> ND3_gas (2body_gas)
    kall(1519) = 5.12847904d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !H_0001 + NO_0001 -> HNO_0001 (2body_dust)
    kall(1520) = 9.93385580d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !H_0001 + NO_0001 -> HNO_gas (2body_gas)
    kall(1521) = 6.61441990d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !D_0001 + NO_0001 -> DNO_0001 (2body_dust)
    kall(1522) = 9.93300999d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !D_0001 + NO_0001 -> DNO_gas (2body_gas)
    kall(1523) = 6.69900075d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !H_0001 + O_0001 -> OH_0001 (2body_dust)
    kall(1524) = 9.90975543d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !H_0001 + O_0001 -> OH_gas (2body_gas)
    kall(1525) = 9.02445722d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !D_0001 + O_0001 -> OD_0001 (2body_dust)
    kall(1526) = 9.90962490d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !D_0001 + O_0001 -> OD_gas (2body_gas)
    kall(1527) = 9.03751030d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !H_0001 + O2_0001 -> O2H_0001 (2body_dust)
    kall(1528) = 3.52434273d+12*max(9.15952449d-07, exp(-1.20000000d+03 * invTd))&
        *(3.54191744d+12*max(9.15952449d-07, exp(-1.20000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))*indns

    !H_0001 + O2_0001 -> O2H_gas (2body_gas)
    kall(1529) = 1.75747061d+10*max(9.15952449d-07, exp(-1.20000000d+03 * invTd))&
        *(3.54191744d+12*max(9.15952449d-07, exp(-1.20000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))*indns

    !D_0001 + O2_0001 -> O2D_0001 (2body_dust)
    kall(1530) = 2.54359480d+12*max(3.86645074d-09, exp(-1.20000000d+03 * invTd))&
        *(2.55656768d+12*max(3.86645074d-09, exp(-1.20000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))*indns

    !D_0001 + O2_0001 -> O2D_gas (2body_gas)
    kall(1531) = 1.29728858d+10*max(3.86645074d-09, exp(-1.20000000d+03 * invTd))&
        *(2.55656768d+12*max(3.86645074d-09, exp(-1.20000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))*indns

    !H_0001 + O2H_0001 -> HOOH_0001 (2body_dust)
    kall(1532) = 3.58530551d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !H_0001 + O2H_0001 -> HOOH_gas (2body_gas)
    kall(1533) = 1.46944905d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !H_0001 + O2D_0001 -> HOOD_0001 (2body_dust)
    kall(1534) = 3.58526241d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !H_0001 + O2D_0001 -> HOOD_gas (2body_gas)
    kall(1535) = 1.47375899d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !D_0001 + O2H_0001 -> HOOD_0001 (2body_dust)
    kall(1536) = 3.58503298d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !D_0001 + O2H_0001 -> HOOD_gas (2body_gas)
    kall(1537) = 1.49670158d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !D_0001 + O2D_0001 -> DOOD_0001 (2body_dust)
    kall(1538) = 3.58522496d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !D_0001 + O2D_0001 -> DOOD_gas (2body_gas)
    kall(1539) = 1.47750358d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !H_0001 + O3_0001 -> OH_0001 + O2_0001 (2body_dust_dust)
    kall(1540) = 3.54191744d+12*max(1.92132168d-04, exp(-4.50000000d+02 * invTd))&
        *(3.54191744d+12*max(1.92132168d-04, exp(-4.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.05000000d+03*invTd)*1.04771331d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.05000000d+03*invTd)*1.04771331d+12))*indns

    !D_0001 + O3_0001 -> OD_0001 + O2_0001 (2body_dust_dust)
    kall(1541) = 2.55656768d+12*max(6.26671747d-06, exp(-4.50000000d+02 * invTd))&
        *(2.55656768d+12*max(6.26671747d-06, exp(-4.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.05000000d+03*invTd)*1.04771331d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.05000000d+03*invTd)*1.04771331d+12))*indns

    !H_0001 + HNCO_0001 -> NH2_0001 + CO_0001 (2body_dust_dust)
    kall(1542) = 3.54191744d+12*max(4.05509643d-09, exp(-2.30000000d+03 * invTd))&
        *(3.54191744d+12*max(4.05509643d-09, exp(-2.30000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.60230630d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.60230630d+12))*indns

    !D_0001 + HNCO_0001 -> NHD_0001 + CO_0001 (2body_dust_dust)
    kall(1543) = 2.55656768d+12*max(1.83878035d-12, exp(-2.30000000d+03 * invTd))&
        *(2.55656768d+12*max(1.83878035d-12, exp(-2.30000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.60230630d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.60230630d+12))*indns

    !H_0001 + DNCO_0001 -> NHD_0001 + CO_0001 (2body_dust_dust)
    kall(1544) = 3.54191744d+12*max(4.03490200d-09, exp(-2.30000000d+03 * invTd))&
        *(3.54191744d+12*max(4.03490200d-09, exp(-2.30000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.58399363d+12))*indns

    !D_0001 + DNCO_0001 -> ND2_0001 + CO_0001 (2body_dust_dust)
    kall(1545) = 2.55656768d+12*max(1.81383741d-12, exp(-2.30000000d+03 * invTd))&
        *(2.55656768d+12*max(1.81383741d-12, exp(-2.30000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.58399363d+12))*indns

    !H_0001 + HNCO_0001 -> NH2CO_0001 (2body_dust)
    kall(1546) = 3.54191744d+12*max(2.99278803d-07, exp(-1.39000000d+03 * invTd))&
        *(3.54191744d+12*max(2.99278803d-07, exp(-1.39000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.60230630d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.60230630d+12))*indns

    !D_0001 + HNCO_0001 -> NHDCO_0001 (2body_dust)
    kall(1547) = 2.55656768d+12*max(7.53126022d-10, exp(-1.39000000d+03 * invTd))&
        *(2.55656768d+12*max(7.53126022d-10, exp(-1.39000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.60230630d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.60230630d+12))*indns

    !H_0001 + DNCO_0001 -> NHDCO_0001 (2body_dust)
    kall(1548) = 3.54191744d+12*max(2.98119516d-07, exp(-1.39000000d+03 * invTd))&
        *(3.54191744d+12*max(2.98119516d-07, exp(-1.39000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.20000000d+03*invTd)*1.58399363d+12))*indns

    !D_0001 + DNCO_0001 -> ND2CO_0001 (2body_dust)
    kall(1549) = 2.55656768d+12*max(7.45171982d-10, exp(-1.39000000d+03 * invTd))&
        *(2.55656768d+12*max(7.45171982d-10, exp(-1.39000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.20000000d+03*invTd)*1.58399363d+12))*indns

    !OH_0001 + HNCO_0001 -> H2O_0001 + OCN_0001 (2body_dust_dust)
    kall(1550) = 2.60560084d+12*max(6.44574769d-23, exp(-1.29000000d+03 * invTd))&
        *(2.60560084d+12*max(6.44574769d-23, exp(-1.29000000d+03 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.20000000d+03*invTd)*1.60230630d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.20000000d+03*invTd)*1.60230630d+12))*indns

    !OD_0001 + HNCO_0001 -> HDO_0001 + OCN_0001 (2body_dust_dust)
    kall(1551) = 2.53218886d+12*max(2.25875190d-23, exp(-1.29000000d+03 * invTd))&
        *(2.53218886d+12*max(2.25875190d-23, exp(-1.29000000d+03 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.20000000d+03*invTd)*1.60230630d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.20000000d+03*invTd)*1.60230630d+12))*indns

    !OH_0001 + DNCO_0001 -> HDO_0001 + OCN_0001 (2body_dust_dust)
    kall(1552) = 2.60560084d+12*max(5.46360087d-23, exp(-1.29000000d+03 * invTd))&
        *(2.60560084d+12*max(5.46360087d-23, exp(-1.29000000d+03 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.20000000d+03*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.20000000d+03*invTd)*1.58399363d+12))*indns

    !OD_0001 + DNCO_0001 -> D2O_0001 + OCN_0001 (2body_dust_dust)
    kall(1553) = 2.53218886d+12*max(1.89473000d-23, exp(-1.29000000d+03 * invTd))&
        *(2.53218886d+12*max(1.89473000d-23, exp(-1.29000000d+03 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.20000000d+03*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.20000000d+03*invTd)*1.58399363d+12))*indns

    !H_0001 + OCN_0001 -> HNCO_0001 (2body_dust)
    kall(1554) = 9.95812663d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.19738664d+12))*indns

    !H_0001 + OCN_0001 -> HNCO_gas (2body_gas)
    kall(1555) = 4.18733655d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.20000000d+03*invTd)*1.19738664d+12))*indns

    !D_0001 + OCN_0001 -> DNCO_0001 (2body_dust)
    kall(1556) = 9.95728593d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.19738664d+12))*indns

    !D_0001 + OCN_0001 -> DNCO_gas (2body_gas)
    kall(1557) = 4.27140734d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.20000000d+03*invTd)*1.19738664d+12))*indns

    !H_0001 + OCS_0001 -> CO_0001 + HS_0001 (2body_dust_dust)
    kall(1558) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.44400000d+03*invTd)*1.09894612d+12))*indns

    !D_0001 + OCS_0001 -> CO_0001 + DS_0001 (2body_dust_dust)
    kall(1559) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.44400000d+03*invTd)*1.09894612d+12))*indns

    !H_0001 + OH_0001 -> H2O_0001 (2body_dust)
    kall(1560) = 9.92633255d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !H_0001 + OH_0001 -> H2O_gas (2body_gas)
    kall(1561) = 7.36674467d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !H_0001 + OD_0001 -> HDO_0001 (2body_dust)
    kall(1562) = 9.92605838d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !H_0001 + OD_0001 -> HDO_gas (2body_gas)
    kall(1563) = 7.39416209d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !D_0001 + OH_0001 -> HDO_0001 (2body_dust)
    kall(1564) = 9.92576693d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !D_0001 + OH_0001 -> HDO_gas (2body_gas)
    kall(1565) = 7.42330691d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !D_0001 + OD_0001 -> D2O_0001 (2body_dust)
    kall(1566) = 9.92561944d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !D_0001 + OD_0001 -> D2O_gas (2body_gas)
    kall(1567) = 7.43805561d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !H_0001 + S_0001 -> HS_0001 (2body_dust)
    kall(1568) = 9.90728616d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !H_0001 + S_0001 -> HS_gas (2body_gas)
    kall(1569) = 9.27138409d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !D_0001 + S_0001 -> DS_0001 (2body_dust)
    kall(1570) = 9.90718050d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !D_0001 + S_0001 -> DS_gas (2body_gas)
    kall(1571) = 9.28195036d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !H_0001 + SiH_0001 -> SiH2_0001 (2body_dust)
    kall(1572) = 9.92695229d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-6.50000000d+03*invTd)*3.35371465d+12))*indns

    !H_0001 + SiH_0001 -> SiH2_gas (2body_gas)
    kall(1573) = 7.30477085d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-6.50000000d+03*invTd)*3.35371465d+12))*indns

    !H_0001 + SiH2_0001 -> SiH3_0001 (2body_dust)
    kall(1574) = 9.95293735d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.80000000d+03*invTd)*1.73517809d+12))*indns

    !H_0001 + SiH2_0001 -> SiH3_gas (2body_gas)
    kall(1575) = 4.70626510d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.80000000d+03*invTd)*1.73517809d+12))*indns

    !H_0001 + SiH3_0001 -> SiH4_0001 (2body_dust)
    kall(1576) = 9.96123709d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.02500000d+03*invTd)*1.81050653d+12))*indns

    !H_0001 + SiH3_0001 -> SiH4_gas (2body_gas)
    kall(1577) = 3.87629116d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.02500000d+03*invTd)*1.81050653d+12))*indns

    !H_0001 + SO2_0001 -> O2_0001 + HS_0001 (2body_dust_dust)
    kall(1578) = 1.00000000d+00*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.70250000d+03*invTd)*1.15537244d+12))*indns

    !D_0001 + SO2_0001 -> O2_0001 + DS_0001 (2body_dust_dust)
    kall(1579) = 1.00000000d+00*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.70250000d+03*invTd)*1.15537244d+12))*indns

    !HD_0001 + C_0001 -> CHD_0001 (2body_dust)
    kall(1580) = 4.53190249d+12*max(1.94534586d-14, exp(-2.50000000d+03 * invTd))&
        *(4.57259575d+12*max(1.94534586d-14, exp(-2.50000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !HD_0001 + C_0001 -> CHD_gas (2body_gas)
    kall(1581) = 4.06932548d+10*max(1.94534586d-14, exp(-2.50000000d+03 * invTd))&
        *(4.57259575d+12*max(1.94534586d-14, exp(-2.50000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !HD_0001 + C2_0001 -> D_0001 + CCH_0001 (2body_dust_dust)
    kall(1582) = 1.61665673d+12*max(1.85013492d-19, exp(-4.20000000d+03 * invTd))&
        *(3.23331346d+12*max(1.85013492d-19, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !HD_0001 + C2_0001 -> H_0001 + CCD_0001 (2body_dust_dust)
    kall(1583) = 1.61665673d+12*max(1.85013492d-19, exp(-4.20000000d+03 * invTd))&
        *(3.23331346d+12*max(1.85013492d-19, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !HD_0001 + CCH_0001 -> D_0001 + C2H2_0001 (2body_dust_dust)
    kall(1584) = 6.46016029d+11*max(1.68048905d-19, exp(-4.20000000d+03 * invTd))&
        *(1.93998808d+12*max(1.68048905d-19, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !HD_0001 + CCH_0001 -> H_0001 + C2HD_0001 (2body_dust_dust)
    kall(1585) = 1.29397205d+12*max(1.68048905d-19, exp(-4.20000000d+03 * invTd))&
        *(1.93998808d+12*max(1.68048905d-19, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !HD_0001 + CCD_0001 -> D_0001 + C2HD_0001 (2body_dust_dust)
    kall(1586) = 1.29397205d+12*max(1.53685134d-19, exp(-4.20000000d+03 * invTd))&
        *(1.93998808d+12*max(1.53685134d-19, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !HD_0001 + CCD_0001 -> H_0001 + C2D2_0001 (2body_dust_dust)
    kall(1587) = 6.46016029d+11*max(1.53685134d-19, exp(-4.20000000d+03 * invTd))&
        *(1.93998808d+12*max(1.53685134d-19, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !HD_0001 + C3_0001 -> D_0001 + l_C3H_0001 (2body_dust_dust)
    kall(1588) = 9.69994038d+11*max(8.13497389d-20, exp(-4.20000000d+03 * invTd))&
        *(1.93998808d+12*max(8.13497389d-20, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !HD_0001 + C3_0001 -> D_0001 + c_C3H_0001 (2body_dust_dust)
    kall(1589) = 9.69994038d+11*max(8.13497389d-20, exp(-4.20000000d+03 * invTd))&
        *(1.93998808d+12*max(8.13497389d-20, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !HD_0001 + C3_0001 -> H_0001 + l_C3D_0001 (2body_dust_dust)
    kall(1590) = 9.69994038d+11*max(8.13497389d-20, exp(-4.20000000d+03 * invTd))&
        *(1.93998808d+12*max(8.13497389d-20, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !HD_0001 + C3_0001 -> H_0001 + c_C3D_0001 (2body_dust_dust)
    kall(1591) = 9.69994038d+11*max(8.13497389d-20, exp(-4.20000000d+03 * invTd))&
        *(1.93998808d+12*max(8.13497389d-20, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !HD_0001 + l_C3H_0001 -> D_0001 + l_C3H2_0001 (2body_dust_dust)
    kall(1592) = 1.93998808d+12*max(7.77108150d-20, exp(-4.20000000d+03 * invTd))&
        *(1.93998808d+12*max(7.77108150d-20, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))*indns

    !HD_0001 + l_C3H_0001 -> H_0001 + l_C3HD_0001 (2body_dust_dust)
    kall(1593) = 1.93998808d+12*max(7.77108150d-20, exp(-4.20000000d+03 * invTd))&
        *(1.93998808d+12*max(7.77108150d-20, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))*indns

    !HD_0001 + l_C3D_0001 -> D_0001 + l_C3HD_0001 (2body_dust_dust)
    kall(1594) = 1.93998808d+12*max(7.44038528d-20, exp(-4.20000000d+03 * invTd))&
        *(1.93998808d+12*max(7.44038528d-20, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))*indns

    !HD_0001 + l_C3D_0001 -> H_0001 + l_C3D2_0001 (2body_dust_dust)
    kall(1595) = 1.93998808d+12*max(7.44038528d-20, exp(-4.20000000d+03 * invTd))&
        *(1.93998808d+12*max(7.44038528d-20, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))*indns

    !HD_0001 + CH2_0001 -> D_0001 + CH3_0001 (2body_dust_dust)
    kall(1596) = 4.84997019d+11*max(2.94906585d-17, exp(-3.53000000d+03 * invTd))&
        *(1.93998808d+12*max(2.94906585d-17, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !HD_0001 + CH2_0001 -> H_0001 + CH2D_0001 (2body_dust_dust)
    kall(1597) = 1.45499106d+12*max(2.94906585d-17, exp(-3.53000000d+03 * invTd))&
        *(1.93998808d+12*max(2.94906585d-17, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !HD_0001 + CHD_0001 -> D_0001 + CH2D_0001 (2body_dust_dust)
    kall(1598) = 9.69994038d+11*max(2.35277744d-17, exp(-3.53000000d+03 * invTd))&
        *(1.93998808d+12*max(2.35277744d-17, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !HD_0001 + CHD_0001 -> H_0001 + CHD2_0001 (2body_dust_dust)
    kall(1599) = 9.69994038d+11*max(2.35277744d-17, exp(-3.53000000d+03 * invTd))&
        *(1.93998808d+12*max(2.35277744d-17, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !HD_0001 + CD2_0001 -> H_0001 + CD3_0001 (2body_dust_dust)
    kall(1600) = 4.84997019d+11*max(1.92438497d-17, exp(-3.53000000d+03 * invTd))&
        *(1.93998808d+12*max(1.92438497d-17, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !HD_0001 + CD2_0001 -> D_0001 + CHD2_0001 (2body_dust_dust)
    kall(1601) = 1.45499106d+12*max(1.92438497d-17, exp(-3.53000000d+03 * invTd))&
        *(1.93998808d+12*max(1.92438497d-17, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !HD_0001 + CH3_0001 -> D_0001 + CH4_0001 (2body_dust_dust)
    kall(1602) = 3.87997615d+11*max(3.46869319d-23, exp(-6.44000000d+03 * invTd))&
        *(1.93998808d+12*max(3.46869319d-23, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !HD_0001 + CH3_0001 -> H_0001 + CH3D_0001 (2body_dust_dust)
    kall(1603) = 1.55199046d+12*max(3.46869319d-23, exp(-6.44000000d+03 * invTd))&
        *(1.93998808d+12*max(3.46869319d-23, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !HD_0001 + CD3_0001 -> D_0001 + CHD3_0001 (2body_dust_dust)
    kall(1604) = 1.55199046d+12*max(1.66560021d-23, exp(-6.44000000d+03 * invTd))&
        *(1.93998808d+12*max(1.66560021d-23, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !HD_0001 + CH2D_0001 -> D_0001 + CH3D_0001 (2body_dust_dust)
    kall(1605) = 7.75995231d+11*max(2.64402578d-23, exp(-6.44000000d+03 * invTd))&
        *(1.93998808d+12*max(2.64402578d-23, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !HD_0001 + CH2D_0001 -> H_0001 + CH2D2_0001 (2body_dust_dust)
    kall(1606) = 1.16399285d+12*max(2.64402578d-23, exp(-6.44000000d+03 * invTd))&
        *(1.93998808d+12*max(2.64402578d-23, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !HD_0001 + CHD2_0001 -> H_0001 + CHD3_0001 (2body_dust_dust)
    kall(1607) = 7.75995231d+11*max(2.07338259d-23, exp(-6.44000000d+03 * invTd))&
        *(1.93998808d+12*max(2.07338259d-23, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !HD_0001 + CHD2_0001 -> D_0001 + CH2D2_0001 (2body_dust_dust)
    kall(1608) = 1.16399285d+12*max(2.07338259d-23, exp(-6.44000000d+03 * invTd))&
        *(1.93998808d+12*max(2.07338259d-23, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !HD_0001 + CN_0001 -> D_0001 + HCN_0001 (2body_dust_dust)
    kall(1609) = 9.69994038d+11*max(6.19874085d-14, exp(-2.07000000d+03 * invTd))&
        *(1.93998808d+12*max(6.19874085d-14, exp(-2.07000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !HD_0001 + CN_0001 -> H_0001 + DCN_0001 (2body_dust_dust)
    kall(1610) = 9.69994038d+11*max(6.19874085d-14, exp(-2.07000000d+03 * invTd))&
        *(1.93998808d+12*max(6.19874085d-14, exp(-2.07000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !HD_0001 + NH2_0001 -> D_0001 + NH3_0001 (2body_dust_dust)
    kall(1611) = 2.01609475d+12*max(4.66685481d-23, exp(-6.30000000d+03 * invTd))&
        *(2.24010528d+12*max(4.66685481d-23, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !HD_0001 + NH2_0001 -> H_0001 + NH2D_0001 (2body_dust_dust)
    kall(1612) = 2.24010528d+11*max(4.66685481d-23, exp(-6.30000000d+03 * invTd))&
        *(2.24010528d+12*max(4.66685481d-23, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !HD_0001 + NHD_0001 -> D_0001 + NH2D_0001 (2body_dust_dust)
    kall(1613) = 1.95589920d+12*max(3.66937446d-23, exp(-6.30000000d+03 * invTd))&
        *(2.17322133d+12*max(3.66937446d-23, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !HD_0001 + NHD_0001 -> H_0001 + NHD2_0001 (2body_dust_dust)
    kall(1614) = 2.17322133d+11*max(3.66937446d-23, exp(-6.30000000d+03 * invTd))&
        *(2.17322133d+12*max(3.66937446d-23, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !HD_0001 + ND2_0001 -> H_0001 + ND3_0001 (2body_dust_dust)
    kall(1615) = 2.11199151d+11*max(2.95476417d-23, exp(-6.30000000d+03 * invTd))&
        *(2.11199151d+12*max(2.95476417d-23, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !HD_0001 + ND2_0001 -> D_0001 + NHD2_0001 (2body_dust_dust)
    kall(1616) = 1.90079236d+12*max(2.95476417d-23, exp(-6.30000000d+03 * invTd))&
        *(2.11199151d+12*max(2.95476417d-23, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !HD_0001 + OH_0001 -> D_0001 + H2O_0001 (2body_dust_dust)
    kall(1617) = 2.34504076d+12*max(1.11406664d-13, exp(-2.10000000d+03 * invTd))&
        *(2.60560084d+12*max(1.11406664d-13, exp(-2.10000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !HD_0001 + OH_0001 -> H_0001 + HDO_0001 (2body_dust_dust)
    kall(1618) = 2.60560084d+11*max(1.11406664d-13, exp(-2.10000000d+03 * invTd))&
        *(2.60560084d+12*max(1.11406664d-13, exp(-2.10000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !HD_0001 + OD_0001 -> D_0001 + HDO_0001 (2body_dust_dust)
    kall(1619) = 2.27896998d+12*max(9.83105867d-14, exp(-2.10000000d+03 * invTd))&
        *(2.53218886d+12*max(9.83105867d-14, exp(-2.10000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !HD_0001 + OD_0001 -> H_0001 + D2O_0001 (2body_dust_dust)
    kall(1620) = 2.53218886d+11*max(9.83105867d-14, exp(-2.10000000d+03 * invTd))&
        *(2.53218886d+12*max(9.83105867d-14, exp(-2.10000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.93998808d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !HCO_0001 + CH2OH_0001 -> HCOOCH3_0001 (2body_dust)
    kall(1621) = 9.99580227d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !HCO_0001 + CH2OH_0001 -> HCOOCH3_gas (2body_gas)
    kall(1622) = 4.19773404d-04*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !HCO_0001 + CHDOH_0001 -> DCOOCH3_0001 (2body_dust)
    kall(1623) = 2.49868528d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !HCO_0001 + CHDOH_0001 -> HCOOCH2D_0001 (2body_dust)
    kall(1624) = 7.49605584d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !HCO_0001 + CH2OD_0001 -> DCOOCH3_0001 (2body_dust)
    kall(1625) = 2.49895057d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !HCO_0001 + CH2OD_0001 -> HCOOCH2D_0001 (2body_dust)
    kall(1626) = 7.49685170d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !HCO_0001 + CD2OH_0001 -> HCOOCHD2_0001 (2body_dust)
    kall(1627) = 4.99749088d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !HCO_0001 + CD2OH_0001 -> DCOOCH2D_0001 (2body_dust)
    kall(1628) = 4.99749088d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !HCO_0001 + CHDOD_0001 -> HCOOCHD2_0001 (2body_dust)
    kall(1629) = 4.99749088d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !HCO_0001 + CHDOD_0001 -> DCOOCH2D_0001 (2body_dust)
    kall(1630) = 4.99749088d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !HCO_0001 + CD2OD_0001 -> HCOOCD3_0001 (2body_dust)
    kall(1631) = 2.50000000d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.80193994d+12))*indns

    !HCO_0001 + CD2OD_0001 -> DCOOCHD2_0001 (2body_dust)
    kall(1632) = 7.50000000d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.80193994d+12))*indns

    !DCO_0001 + CH2OH_0001 -> DCOOCH3_0001 (2body_dust)
    kall(1633) = 2.49898879d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !DCO_0001 + CH2OH_0001 -> HCOOCH2D_0001 (2body_dust)
    kall(1634) = 7.49696636d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !DCO_0001 + CHDOH_0001 -> HCOOCHD2_0001 (2body_dust)
    kall(1635) = 4.99745330d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !DCO_0001 + CHDOH_0001 -> DCOOCH2D_0001 (2body_dust)
    kall(1636) = 4.99745330d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !DCO_0001 + CH2OD_0001 -> HCOOCHD2_0001 (2body_dust)
    kall(1637) = 4.99797758d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !DCO_0001 + CH2OD_0001 -> DCOOCH2D_0001 (2body_dust)
    kall(1638) = 4.99797758d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !DCO_0001 + CD2OH_0001 -> HCOOCD3_0001 (2body_dust)
    kall(1639) = 2.50000000d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !DCO_0001 + CD2OH_0001 -> DCOOCHD2_0001 (2body_dust)
    kall(1640) = 7.50000000d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !DCO_0001 + CHDOD_0001 -> HCOOCD3_0001 (2body_dust)
    kall(1641) = 2.50000000d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !DCO_0001 + CHDOD_0001 -> DCOOCHD2_0001 (2body_dust)
    kall(1642) = 7.50000000d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !HCO_0001 + CH3O_0001 -> HCOOCH3_0001 (2body_dust)
    kall(1643) = 9.99392783d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !HCO_0001 + CH3O_0001 -> HCOOCH3_gas (2body_gas)
    kall(1644) = 6.07217125d-04*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !HCO_0001 + CH2DO_0001 -> DCOOCH3_0001 (2body_dust)
    kall(1645) = 2.49895057d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !HCO_0001 + CH2DO_0001 -> HCOOCH2D_0001 (2body_dust)
    kall(1646) = 7.49685170d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !HCO_0001 + CHD2O_0001 -> HCOOCHD2_0001 (2body_dust)
    kall(1647) = 4.99749088d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !HCO_0001 + CHD2O_0001 -> DCOOCH2D_0001 (2body_dust)
    kall(1648) = 4.99749088d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !HCO_0001 + CD3O_0001 -> HCOOCD3_0001 (2body_dust)
    kall(1649) = 2.50000000d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.80193994d+12))*indns

    !HCO_0001 + CD3O_0001 -> DCOOCHD2_0001 (2body_dust)
    kall(1650) = 7.50000000d-01*((exp(-1.20000000d+03*invTd)*1.44098697d+12) + (exp(-2.20000000d+03*invTd)*1.80193994d+12))*indns

    !DCO_0001 + CH3O_0001 -> DCOOCH3_0001 (2body_dust)
    kall(1651) = 2.49852520d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !DCO_0001 + CH3O_0001 -> HCOOCH2D_0001 (2body_dust)
    kall(1652) = 7.49557559d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.88711741d+12))*indns

    !DCO_0001 + CH2DO_0001 -> HCOOCHD2_0001 (2body_dust)
    kall(1653) = 4.99797758d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !DCO_0001 + CH2DO_0001 -> DCOOCH2D_0001 (2body_dust)
    kall(1654) = 4.99797758d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.85739717d+12))*indns

    !DCO_0001 + CHD2O_0001 -> HCOOCD3_0001 (2body_dust)
    kall(1655) = 2.50000000d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !DCO_0001 + CHD2O_0001 -> DCOOCHD2_0001 (2body_dust)
    kall(1656) = 7.50000000d-01*((exp(-1.20000000d+03*invTd)*1.41676697d+12) + (exp(-2.20000000d+03*invTd)*1.82903830d+12))*indns

    !N_0001 + C2_0001 -> CCN_0001 (2body_dust)
    kall(1657) = 9.91021537d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !N_0001 + C2_0001 -> CCN_gas (2body_gas)
    kall(1658) = 8.97846265d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !N_0001 + C3_0001 -> C3N_0001 (2body_dust)
    kall(1659) = 9.92374312d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !N_0001 + C3_0001 -> C3N_gas (2body_gas)
    kall(1660) = 7.62568783d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !N_0001 + l_C3H_0001 -> HC3N_0001 (2body_dust)
    kall(1661) = 9.93470909d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))*indns

    !N_0001 + l_C3H_0001 -> HC3N_gas (2body_gas)
    kall(1662) = 6.52909096d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))*indns

    !N_0001 + c_C3H_0001 -> HC3N_0001 (2body_dust)
    kall(1663) = 9.93470909d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))*indns

    !N_0001 + c_C3H_0001 -> HC3N_gas (2body_gas)
    kall(1664) = 6.52909096d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))*indns

    !N_0001 + l_C3D_0001 -> DC3N_0001 (2body_dust)
    kall(1665) = 9.93470909d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))*indns

    !N_0001 + l_C3D_0001 -> DC3N_gas (2body_gas)
    kall(1666) = 6.52909096d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))*indns

    !N_0001 + c_C3D_0001 -> DC3N_0001 (2body_dust)
    kall(1667) = 9.93470909d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))*indns

    !N_0001 + c_C3D_0001 -> DC3N_gas (2body_gas)
    kall(1668) = 6.52909096d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))*indns

    !CH_0001 + N_0001 -> HCN_0001 (2body_dust)
    kall(1669) = 9.91033632d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-3.60000000d+02*invTd)*1.13594070d+12))*indns

    !CH_0001 + N_0001 -> HCN_gas (2body_gas)
    kall(1670) = 8.96636802d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-3.60000000d+02*invTd)*1.13594070d+12))*indns

    !N_0001 + CD_0001 -> DCN_0001 (2body_dust)
    kall(1671) = 9.91037378d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-4.62500000d+02*invTd)*1.28753866d+12))*indns

    !N_0001 + CD_0001 -> DCN_gas (2body_gas)
    kall(1672) = 8.96262216d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-4.62500000d+02*invTd)*1.28753866d+12))*indns

    !N_0001 + CH2_0001 -> H2CN_0001 (2body_dust)
    kall(1673) = 9.91821907d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !N_0001 + CH2_0001 -> H2CN_gas (2body_gas)
    kall(1674) = 8.17809299d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !N_0001 + CHD_0001 -> HDCN_0001 (2body_dust)
    kall(1675) = 9.91846717d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !N_0001 + CHD_0001 -> HDCN_gas (2body_gas)
    kall(1676) = 8.15328327d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !N_0001 + CD2_0001 -> D2CN_0001 (2body_dust)
    kall(1677) = 9.91859785d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !N_0001 + CD2_0001 -> D2CN_gas (2body_gas)
    kall(1678) = 8.14021457d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !N_0001 + CH3_0001 -> CH2NH_0001 (2body_dust)
    kall(1679) = 9.95906799d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !N_0001 + CH3_0001 -> CH2NH_gas (2body_gas)
    kall(1680) = 4.09320083d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !N_0001 + CD3_0001 -> CD2ND_0001 (2body_dust)
    kall(1681) = 9.96066399d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !N_0001 + CD3_0001 -> CD2ND_gas (2body_gas)
    kall(1682) = 3.93360066d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !N_0001 + CH2D_0001 -> CHDNH_0001 (2body_dust)
    kall(1683) = 6.64322419d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !N_0001 + CH2D_0001 -> CHDNH_gas (2body_gas)
    kall(1684) = 2.67758135d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !N_0001 + CH2D_0001 -> CH2ND_0001 (2body_dust)
    kall(1685) = 3.31663217d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !N_0001 + CH2D_0001 -> CH2ND_gas (2body_gas)
    kall(1686) = 1.33678349d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !N_0001 + CHD2_0001 -> CD2NH_0001 (2body_dust)
    kall(1687) = 3.31676582d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !N_0001 + CHD2_0001 -> CD2NH_gas (2body_gas)
    kall(1688) = 1.32341765d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !N_0001 + CHD2_0001 -> CHDND_0001 (2body_dust)
    kall(1689) = 6.64349190d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !N_0001 + CHD2_0001 -> CHDND_gas (2body_gas)
    kall(1690) = 2.65080952d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !N_0001 + HS_0001 -> H_0001 + NS_0001 (2body_dust_dust)
    kall(1691) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.35000000d+03*invTd)*1.43277615d+12))*indns

    !N_0001 + DS_0001 -> D_0001 + NS_0001 (2body_dust_dust)
    kall(1692) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.35000000d+03*invTd)*1.41154866d+12))*indns

    !N_0001 + N_0001 -> N2_0001 (2body_dust)
    kall(1693) = 5d-1 * 9.90199825d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-3.60000000d+02*invTd)*1.13594070d+12))*indns

    !N_0001 + N_0001 -> N2_gas (2body_gas)
    kall(1694) = 5d-1 * 9.80017526d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-3.60000000d+02*invTd)*1.13594070d+12))*indns

    !N_0001 + NH_0001 -> H_0001 + N2_0001 (2body_dust_dust)
    kall(1695) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !N_0001 + ND_0001 -> D_0001 + N2_0001 (2body_dust_dust)
    kall(1696) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !N_0001 + NH2_0001 -> N2H2_0001 (2body_dust)
    kall(1697) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !N_0001 + NHD_0001 -> N2HD_0001 (2body_dust)
    kall(1698) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !N_0001 + ND2_0001 -> N2D2_0001 (2body_dust)
    kall(1699) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !N_0001 + NO_0001 -> O_0001 + N2_0001 (2body_dust_dust)
    kall(1700) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !N_0001 + NS_0001 -> N2_0001 + S_0001 (2body_dust_dust)
    kall(1701) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-9.50000000d+02*invTd)*1.01800829d+12))*indns

    !N_0001 + O_0001 -> NO_0001 (2body_dust)
    kall(1702) = 9.90305582d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !N_0001 + O_0001 -> NO_gas (2body_gas)
    kall(1703) = 9.69441806d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !N_0001 + O2H_0001 -> NH_0001 + O2_0001 (2body_dust_dust)
    kall(1704) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !N_0001 + O2D_0001 -> ND_0001 + O2_0001 (2body_dust_dust)
    kall(1705) = 1.00000000d+00*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !N_0001 + S_0001 -> NS_0001 (2body_dust)
    kall(1706) = 9.90421078d-01*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !N_0001 + S_0001 -> NS_gas (2body_gas)
    kall(1707) = 9.57892172d-03*((exp(-3.60000000d+02*invTd)*1.13594070d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !CH2_0001 + NH_0001 -> CH2NH_0001 (2body_dust)
    kall(1708) = 9.94953215d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !CH2_0001 + NH_0001 -> CH2NH_gas (2body_gas)
    kall(1709) = 5.04678540d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !NH_0001 + CHD_0001 -> CHDNH_0001 (2body_dust)
    kall(1710) = 6.63670343d-01*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !NH_0001 + CHD_0001 -> CHDNH_gas (2body_gas)
    kall(1711) = 3.32965713d-03*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !NH_0001 + CHD_0001 -> CH2ND_0001 (2body_dust)
    kall(1712) = 3.31337667d-01*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !NH_0001 + CHD_0001 -> CH2ND_gas (2body_gas)
    kall(1713) = 1.66233257d-03*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !NH_0001 + CD2_0001 -> CD2NH_0001 (2body_dust)
    kall(1714) = 3.31347226d-01*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !NH_0001 + CD2_0001 -> CD2NH_gas (2body_gas)
    kall(1715) = 1.65277416d-03*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !NH_0001 + CD2_0001 -> CHDND_0001 (2body_dust)
    kall(1716) = 6.63689488d-01*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !NH_0001 + CD2_0001 -> CHDND_gas (2body_gas)
    kall(1717) = 3.31051160d-03*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !CH2_0001 + ND_0001 -> CHDNH_0001 (2body_dust)
    kall(1718) = 6.63638461d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !CH2_0001 + ND_0001 -> CHDNH_gas (2body_gas)
    kall(1719) = 3.36153936d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !CH2_0001 + ND_0001 -> CH2ND_0001 (2body_dust)
    kall(1720) = 3.31321750d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !CH2_0001 + ND_0001 -> CH2ND_gas (2body_gas)
    kall(1721) = 1.67824978d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !CHD_0001 + ND_0001 -> CD2NH_0001 (2body_dust)
    kall(1722) = 3.31340047d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !CHD_0001 + ND_0001 -> CD2NH_gas (2body_gas)
    kall(1723) = 1.65995282d-03*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !CHD_0001 + ND_0001 -> CHDND_0001 (2body_dust)
    kall(1724) = 6.63675110d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !CHD_0001 + ND_0001 -> CHDND_gas (2body_gas)
    kall(1725) = 3.32489047d-03*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !ND_0001 + CD2_0001 -> CD2ND_0001 (2body_dust)
    kall(1726) = 9.95043940d-01*((exp(-1.30000000d+03*invTd)*2.01920361d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !ND_0001 + CD2_0001 -> CD2ND_gas (2body_gas)
    kall(1727) = 4.95606022d-03*((exp(-1.30000000d+03*invTd)*2.01920361d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !NH_0001 + CH3_0001 -> CH3NH_0001 (2body_dust)
    kall(1728) = 1.00000000d+00*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !NH_0001 + NH_0001 -> N2H2_0001 (2body_dust)
    kall(1729) = 5d-1 * 5.00000000d-01*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !NH_0001 + ND_0001 -> N2HD_0001 (2body_dust)
    kall(1730) = 5.00000000d-01*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !ND_0001 + ND_0001 -> N2D2_0001 (2body_dust)
    kall(1731) = 5d-1 * 5.00000000d-01*((exp(-1.30000000d+03*invTd)*2.01920361d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !NH_0001 + ND_0001 -> HD_0001 + N2_0001 (2body_dust_dust)
    kall(1732) = 5.00000000d-01*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !NH_0001 + NO_0001 -> H_0001 + O_0001 + N2_0001 (2body_dust_dust)
    kall(1733) = 1.00000000d+00*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !ND_0001 + NO_0001 -> D_0001 + O_0001 + N2_0001 (2body_dust_dust)
    kall(1734) = 1.00000000d+00*((exp(-1.30000000d+03*invTd)*2.01920361d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !CH3_0001 + NH2_0001 -> CH5N_0001 (2body_dust)
    kall(1735) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !NH2_0001 + HCO_0001 -> NH2CHO_0001 (2body_dust)
    kall(1736) = 9.98011631d-01*((exp(-1.60000000d+03*invTd)*2.24010528d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !NH2_0001 + HCO_0001 -> NH2CHO_gas (2body_gas)
    kall(1737) = 1.98836904d-03*((exp(-1.60000000d+03*invTd)*2.24010528d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !NH2_0001 + DCO_0001 -> NHDCHO_0001 (2body_dust)
    kall(1738) = 6.66690657d-01*((exp(-1.60000000d+03*invTd)*2.24010528d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !NH2_0001 + DCO_0001 -> NHDCHO_gas (2body_gas)
    kall(1739) = 1.30934285d-03*((exp(-1.60000000d+03*invTd)*2.24010528d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !NH2_0001 + DCO_0001 -> NH2CDO_0001 (2body_dust)
    kall(1740) = 3.32347289d-01*((exp(-1.60000000d+03*invTd)*2.24010528d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !NH2_0001 + DCO_0001 -> NH2CDO_gas (2body_gas)
    kall(1741) = 6.52711331d-04*((exp(-1.60000000d+03*invTd)*2.24010528d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !NHD_0001 + HCO_0001 -> NHDCHO_0001 (2body_dust)
    kall(1742) = 6.66685915d-01*((exp(-1.60000000d+03*invTd)*2.17322133d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !NHD_0001 + HCO_0001 -> NHDCHO_gas (2body_gas)
    kall(1743) = 1.31408463d-03*((exp(-1.60000000d+03*invTd)*2.17322133d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !NHD_0001 + HCO_0001 -> NH2CDO_0001 (2body_dust)
    kall(1744) = 3.32344925d-01*((exp(-1.60000000d+03*invTd)*2.17322133d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !NHD_0001 + HCO_0001 -> NH2CDO_gas (2body_gas)
    kall(1745) = 6.55075125d-04*((exp(-1.60000000d+03*invTd)*2.17322133d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !NHD_0001 + DCO_0001 -> ND2CHO_0001 (2body_dust)
    kall(1746) = 3.32354365d-01*((exp(-1.60000000d+03*invTd)*2.17322133d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !NHD_0001 + DCO_0001 -> ND2CHO_gas (2body_gas)
    kall(1747) = 6.45635487d-04*((exp(-1.60000000d+03*invTd)*2.17322133d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !NHD_0001 + DCO_0001 -> NHDCDO_0001 (2body_dust)
    kall(1748) = 6.66704851d-01*((exp(-1.60000000d+03*invTd)*2.17322133d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !NHD_0001 + DCO_0001 -> NHDCDO_gas (2body_gas)
    kall(1749) = 1.29514867d-03*((exp(-1.60000000d+03*invTd)*2.17322133d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !ND2_0001 + HCO_0001 -> ND2CHO_0001 (2body_dust)
    kall(1750) = 3.32351664d-01*((exp(-1.60000000d+03*invTd)*2.11199151d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !ND2_0001 + HCO_0001 -> ND2CHO_gas (2body_gas)
    kall(1751) = 6.48335768d-04*((exp(-1.60000000d+03*invTd)*2.11199151d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !ND2_0001 + HCO_0001 -> NHDCDO_0001 (2body_dust)
    kall(1752) = 6.66699435d-01*((exp(-1.60000000d+03*invTd)*2.11199151d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !ND2_0001 + HCO_0001 -> NHDCDO_gas (2body_gas)
    kall(1753) = 1.30056544d-03*((exp(-1.60000000d+03*invTd)*2.11199151d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !ND2_0001 + DCO_0001 -> ND2CDO_0001 (2body_dust)
    kall(1754) = 9.98081461d-01*((exp(-1.60000000d+03*invTd)*2.11199151d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !ND2_0001 + DCO_0001 -> ND2CDO_gas (2body_gas)
    kall(1755) = 1.91853938d-03*((exp(-1.60000000d+03*invTd)*2.11199151d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !NH2_0001 + NO_0001 -> H2O_0001 + N2_0001 (2body_dust_dust)
    kall(1756) = 1.00000000d+00*((exp(-1.60000000d+03*invTd)*2.24010528d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !NHD_0001 + NO_0001 -> HDO_0001 + N2_0001 (2body_dust_dust)
    kall(1757) = 1.00000000d+00*((exp(-1.60000000d+03*invTd)*2.17322133d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !ND2_0001 + NO_0001 -> D2O_0001 + N2_0001 (2body_dust_dust)
    kall(1758) = 1.00000000d+00*((exp(-1.60000000d+03*invTd)*2.11199151d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !O_0001 + C2_0001 -> CCO_0001 (2body_dust)
    kall(1759) = 9.90767971d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !O_0001 + C2_0001 -> CCO_gas (2body_gas)
    kall(1760) = 9.23202946d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !O_0001 + C3_0001 -> C3O_0001 (2body_dust)
    kall(1761) = 9.91896387d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !O_0001 + C3_0001 -> C3O_gas (2body_gas)
    kall(1762) = 8.10361300d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !CH_0001 + O_0001 -> HCO_0001 (2body_dust)
    kall(1763) = 9.90813727d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CH_0001 + O_0001 -> HCO_gas (2body_gas)
    kall(1764) = 9.18627322d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CD_0001 + O_0001 -> DCO_0001 (2body_dust)
    kall(1765) = 9.90815209d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CD_0001 + O_0001 -> DCO_gas (2body_gas)
    kall(1766) = 9.18479100d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CH2_0001 + O_0001 -> H2CO_0001 (2body_dust)
    kall(1767) = 9.92712201d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CH2_0001 + O_0001 -> H2CO_gas (2body_gas)
    kall(1768) = 7.28779859d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CHD_0001 + O_0001 -> HDCO_0001 (2body_dust)
    kall(1769) = 1.00000000d+00*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CHD_0001 + O_0001 -> HDCO_gas (2body_gas)
    kall(1770) = 0d0*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !O_0001 + CD2_0001 -> D2CO_0001 (2body_dust)
    kall(1771) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !O_0001 + CD2_0001 -> D2CO_gas (2body_gas)
    kall(1772) = 0d0*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !CH3_0001 + O_0001 -> CH3O_0001 (2body_dust)
    kall(1773) = 9.96089135d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !CH3_0001 + O_0001 -> CH3O_gas (2body_gas)
    kall(1774) = 3.91086518d-03*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !O_0001 + CD3_0001 -> CD3O_0001 (2body_dust)
    kall(1775) = 9.96023532d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !O_0001 + CD3_0001 -> CD3O_gas (2body_gas)
    kall(1776) = 3.97646814d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !O_0001 + CH2D_0001 -> CH2DO_0001 (2body_dust)
    kall(1777) = 9.95800862d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !O_0001 + CH2D_0001 -> CH2DO_gas (2body_gas)
    kall(1778) = 4.19913848d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !O_0001 + CHD2_0001 -> CHD2O_0001 (2body_dust)
    kall(1779) = 9.96023532d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !O_0001 + CHD2_0001 -> CHD2O_gas (2body_gas)
    kall(1780) = 3.97646814d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !O_0001 + CN_0001 -> OCN_0001 (2body_dust)
    kall(1781) = 9.91166678d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !O_0001 + CN_0001 -> OCN_gas (2body_gas)
    kall(1782) = 8.83332211d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !O_0001 + CO_0001 -> CO2_0001 (2body_dust)
    kall(1783) = 1.57014117d+12*max(4.71763203d-23, exp(-1.00000000d+03 * invTd))&
        *(1.58399363d+12*max(4.71763203d-23, exp(-1.00000000d+03 * invTd)) + (exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))**(-1d0) &
        *((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))*indns

    !O_0001 + CO_0001 -> CO2_gas (2body_gas)
    kall(1784) = 1.38524583d+10*max(4.71763203d-23, exp(-1.00000000d+03 * invTd))&
        *(1.58399363d+12*max(4.71763203d-23, exp(-1.00000000d+03 * invTd)) + (exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))**(-1d0) &
        *((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))*indns

    !O_0001 + CS_0001 -> OCS_0001 (2body_dust)
    kall(1785) = 9.91124604d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*1.35083431d+12))*indns

    !O_0001 + CS_0001 -> OCS_gas (2body_gas)
    kall(1786) = 8.87539584d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*1.35083431d+12))*indns

    !O_0001 + HCO_0001 -> H_0001 + CO2_0001 (2body_dust_dust)
    kall(1787) = 5.00000000d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !O_0001 + DCO_0001 -> D_0001 + CO2_0001 (2body_dust_dust)
    kall(1788) = 5.00000000d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !O_0001 + HCO_0001 -> OH_0001 + CO_0001 (2body_dust_dust)
    kall(1789) = 5.00000000d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !O_0001 + DCO_0001 -> OD_0001 + CO_0001 (2body_dust_dust)
    kall(1790) = 5.00000000d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !O_0001 + HNO_0001 -> OH_0001 + NO_0001 (2body_dust_dust)
    kall(1791) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))*indns

    !O_0001 + DNO_0001 -> OD_0001 + NO_0001 (2body_dust_dust)
    kall(1792) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))*indns

    !O_0001 + HS_0001 -> H_0001 + SO_0001 (2body_dust_dust)
    kall(1793) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.35000000d+03*invTd)*1.43277615d+12))*indns

    !O_0001 + DS_0001 -> D_0001 + SO_0001 (2body_dust_dust)
    kall(1794) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.35000000d+03*invTd)*1.41154866d+12))*indns

    !NH_0001 + O_0001 -> HNO_0001 (2body_dust)
    kall(1795) = 9.91455604d-01*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !NH_0001 + O_0001 -> HNO_gas (2body_gas)
    kall(1796) = 8.54439595d-03*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !O_0001 + ND_0001 -> DNO_0001 (2body_dust)
    kall(1797) = 9.91458742d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !O_0001 + ND_0001 -> DNO_gas (2body_gas)
    kall(1798) = 8.54125781d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !O_0001 + NH2_0001 -> H_0001 + HNO_0001 (2body_dust_dust)
    kall(1799) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !O_0001 + NHD_0001 -> D_0001 + HNO_0001 (2body_dust_dust)
    kall(1800) = 5.00000000d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !O_0001 + NHD_0001 -> H_0001 + DNO_0001 (2body_dust_dust)
    kall(1801) = 5.00000000d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !O_0001 + ND2_0001 -> D_0001 + DNO_0001 (2body_dust_dust)
    kall(1802) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !O_0001 + NO_0001 -> NO2_0001 (2body_dust)
    kall(1803) = 9.91913717d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !O_0001 + NO_0001 -> NO2_gas (2body_gas)
    kall(1804) = 8.08628286d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !O_0001 + NS_0001 -> NO_0001 + S_0001 (2body_dust_dust)
    kall(1805) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-9.50000000d+02*invTd)*1.01800829d+12))*indns

    !O_0001 + O_0001 -> O2_0001 (2body_dust)
    kall(1806) = 5d-1 * 9.90295305d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !O_0001 + O_0001 -> O2_gas (2body_gas)
    kall(1807) = 5d-1 * 9.70469491d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !O_0001 + O2_0001 -> O3_0001 (2body_dust)
    kall(1808) = 9.94189708d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))*indns

    !O_0001 + O2_0001 -> O3_gas (2body_gas)
    kall(1809) = 5.81029174d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-6.00000000d+02*invTd)*9.69994038d+11))*indns

    !O_0001 + O2H_0001 -> OH_0001 + O2_0001 (2body_dust_dust)
    kall(1810) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !O_0001 + O2D_0001 -> OD_0001 + O2_0001 (2body_dust_dust)
    kall(1811) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !O_0001 + OH_0001 -> O2H_0001 (2body_dust)
    kall(1812) = 9.93967933d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !O_0001 + OH_0001 -> O2H_gas (2body_gas)
    kall(1813) = 6.03206721d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !O_0001 + OD_0001 -> O2D_0001 (2body_dust)
    kall(1814) = 9.93973928d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !O_0001 + OD_0001 -> O2D_gas (2body_gas)
    kall(1815) = 6.02607230d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !O_0001 + S_0001 -> SO_0001 (2body_dust)
    kall(1816) = 9.90538854d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !O_0001 + S_0001 -> SO_gas (2body_gas)
    kall(1817) = 9.46114616d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !O_0001 + SO_0001 -> SO2_0001 (2body_dust)
    kall(1818) = 9.91542436d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.40000000d+03*invTd)*1.20979512d+12))*indns

    !O_0001 + SO_0001 -> SO2_gas (2body_gas)
    kall(1819) = 8.45756420d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.40000000d+03*invTd)*1.20979512d+12))*indns

    !CH2_0001 + OH_0001 -> CH2OH_0001 (2body_dust)
    kall(1820) = 9.95455131d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !CH2_0001 + OH_0001 -> CH2OH_gas (2body_gas)
    kall(1821) = 4.54486911d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !CHD_0001 + OH_0001 -> CHDOH_0001 (2body_dust)
    kall(1822) = 6.64160222d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !CHD_0001 + OH_0001 -> CHDOH_gas (2body_gas)
    kall(1823) = 2.83977765d-03*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !CHD_0001 + OH_0001 -> CH2OD_0001 (2body_dust)
    kall(1824) = 3.31514251d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !CHD_0001 + OH_0001 -> CH2OD_gas (2body_gas)
    kall(1825) = 1.48574948d-03*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !CD2_0001 + OH_0001 -> CD2OH_0001 (2body_dust)
    kall(1826) = 3.31582240d-01*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !CD2_0001 + OH_0001 -> CD2OH_gas (2body_gas)
    kall(1827) = 1.41776005d-03*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !CD2_0001 + OH_0001 -> CHDOD_0001 (2body_dust)
    kall(1828) = 6.64160222d-01*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !CD2_0001 + OH_0001 -> CHDOD_gas (2body_gas)
    kall(1829) = 2.83977765d-03*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !CH2_0001 + OD_0001 -> CHDOH_0001 (2body_dust)
    kall(1830) = 6.64103921d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CH2_0001 + OD_0001 -> CHDOH_gas (2body_gas)
    kall(1831) = 2.89607945d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CH2_0001 + OD_0001 -> CH2OD_0001 (2body_dust)
    kall(1832) = 3.31487967d-01*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CH2_0001 + OD_0001 -> CH2OD_gas (2body_gas)
    kall(1833) = 1.51203312d-03*((exp(-7.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CHD_0001 + OD_0001 -> CD2OH_0001 (2body_dust)
    kall(1834) = 3.31568395d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CHD_0001 + OD_0001 -> CD2OH_gas (2body_gas)
    kall(1835) = 1.43160528d-03*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CHD_0001 + OD_0001 -> CHDOD_0001 (2body_dust)
    kall(1836) = 6.64132490d-01*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CHD_0001 + OD_0001 -> CHDOD_gas (2body_gas)
    kall(1837) = 2.86750967d-03*((exp(-7.00000000d+02*invTd)*1.53028323d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CD2_0001 + OD_0001 -> CD2OD_0001 (2body_dust)
    kall(1838) = 9.95700885d-01*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CD2_0001 + OD_0001 -> CD2OD_gas (2body_gas)
    kall(1839) = 4.29911494d-03*((exp(-7.00000000d+02*invTd)*1.48169037d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CH3_0001 + OH_0001 -> CH3OH_0001 (2body_dust)
    kall(1840) = 9.97543081d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !CH3_0001 + OH_0001 -> CH3OH_gas (2body_gas)
    kall(1841) = 2.45691851d-03*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !OH_0001 + CD3_0001 -> CD3OH_0001 (2body_dust)
    kall(1842) = 2.49435143d-01*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !OH_0001 + CD3_0001 -> CD3OH_gas (2body_gas)
    kall(1843) = 5.64856999d-04*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !OH_0001 + CD3_0001 -> CHD2OD_0001 (2body_dust)
    kall(1844) = 7.48305429d-01*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !OH_0001 + CD3_0001 -> CHD2OD_gas (2body_gas)
    kall(1845) = 1.69457100d-03*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !CH2D_0001 + OH_0001 -> CH2DOH_0001 (2body_dust)
    kall(1846) = 7.48158108d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !CH2D_0001 + OH_0001 -> CH2DOH_gas (2body_gas)
    kall(1847) = 1.84189209d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !CH2D_0001 + OH_0001 -> CH3OD_0001 (2body_dust)
    kall(1848) = 2.49399404d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !CH2D_0001 + OH_0001 -> CH3OD_gas (2body_gas)
    kall(1849) = 6.00596456d-04*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !OH_0001 + CHD2_0001 -> CHD2OH_0001 (2body_dust)
    kall(1850) = 4.98820776d-01*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !OH_0001 + CHD2_0001 -> CHD2OH_gas (2body_gas)
    kall(1851) = 1.17922373d-03*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !OH_0001 + CHD2_0001 -> CH2DOD_0001 (2body_dust)
    kall(1852) = 4.98820776d-01*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !OH_0001 + CHD2_0001 -> CH2DOD_gas (2body_gas)
    kall(1853) = 1.17922373d-03*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !CH3_0001 + OD_0001 -> CH2DOH_0001 (2body_dust)
    kall(1854) = 7.48089850d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CH3_0001 + OD_0001 -> CH2DOH_gas (2body_gas)
    kall(1855) = 1.91015018d-03*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CH3_0001 + OD_0001 -> CH3OD_0001 (2body_dust)
    kall(1856) = 2.49376434d-01*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CH3_0001 + OD_0001 -> CH3OD_gas (2body_gas)
    kall(1857) = 6.23565883d-04*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CH2D_0001 + OD_0001 -> CHD2OH_0001 (2body_dust)
    kall(1858) = 4.98798758d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CH2D_0001 + OD_0001 -> CHD2OH_gas (2body_gas)
    kall(1859) = 1.20124166d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CH2D_0001 + OD_0001 -> CH2DOD_0001 (2body_dust)
    kall(1860) = 4.98798758d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CH2D_0001 + OD_0001 -> CH2DOD_gas (2body_gas)
    kall(1861) = 1.20124166d-03*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CHD2_0001 + OD_0001 -> CD3OH_0001 (2body_dust)
    kall(1862) = 2.49423954d-01*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CHD2_0001 + OD_0001 -> CD3OH_gas (2body_gas)
    kall(1863) = 5.76045765d-04*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CHD2_0001 + OD_0001 -> CHD2OD_0001 (2body_dust)
    kall(1864) = 7.48271863d-01*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CHD2_0001 + OD_0001 -> CHD2OD_gas (2body_gas)
    kall(1865) = 1.72813729d-03*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !OH_0001 + CO_0001 -> H_0001 + CO2_0001 (2body_dust_dust)
    kall(1866) = 2.60560084d+12*max(3.06297367d-12, exp(-4.00000000d+02 * invTd))&
        *(2.60560084d+12*max(3.06297367d-12, exp(-4.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))*indns

    !OD_0001 + CO_0001 -> D_0001 + CO2_0001 (2body_dust_dust)
    kall(1867) = 2.53218886d+12*max(1.91348773d-12, exp(-4.00000000d+02 * invTd))&
        *(2.53218886d+12*max(1.91348773d-12, exp(-4.00000000d+02 * invTd)) + (exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))**(-1d0) &
        *((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))*indns

    !OH_0001 + H2CO_0001 -> H2O_0001 + HCO_0001 (2body_dust_dust)
    kall(1868) = 1.00000000d+00*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !OH_0001 + HDCO_0001 -> HDO_0001 + HCO_0001 (2body_dust_dust)
    kall(1869) = 6.67000000d-01*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !OH_0001 + HDCO_0001 -> H2O_0001 + DCO_0001 (2body_dust_dust)
    kall(1870) = 3.33000000d-01*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !OH_0001 + D2CO_0001 -> D2O_0001 + HCO_0001 (2body_dust_dust)
    kall(1871) = 3.33000000d-01*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))*indns

    !OH_0001 + D2CO_0001 -> HDO_0001 + DCO_0001 (2body_dust_dust)
    kall(1872) = 6.67000000d-01*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))*indns

    !OD_0001 + H2CO_0001 -> HDO_0001 + HCO_0001 (2body_dust_dust)
    kall(1873) = 6.67000000d-01*((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !OD_0001 + H2CO_0001 -> H2O_0001 + DCO_0001 (2body_dust_dust)
    kall(1874) = 3.33000000d-01*((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !OD_0001 + HDCO_0001 -> D2O_0001 + HCO_0001 (2body_dust_dust)
    kall(1875) = 3.33000000d-01*((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !OD_0001 + HDCO_0001 -> HDO_0001 + DCO_0001 (2body_dust_dust)
    kall(1876) = 6.67000000d-01*((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !OD_0001 + D2CO_0001 -> D2O_0001 + DCO_0001 (2body_dust_dust)
    kall(1877) = 1.00000000d+00*((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))*indns

    !OH_0001 + HCO_0001 -> HCOOH_0001 (2body_dust)
    kall(1878) = 9.96181490d-01*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !OH_0001 + HCO_0001 -> HCOOH_gas (2body_gas)
    kall(1879) = 3.81851013d-03*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !OH_0001 + DCO_0001 -> DCOOH_0001 (2body_dust)
    kall(1880) = 4.98127371d-01*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !OH_0001 + DCO_0001 -> DCOOH_gas (2body_gas)
    kall(1881) = 1.87262930d-03*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !OH_0001 + DCO_0001 -> HCOOD_0001 (2body_dust)
    kall(1882) = 4.98127371d-01*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !OH_0001 + DCO_0001 -> HCOOD_gas (2body_gas)
    kall(1883) = 1.87262930d-03*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !OD_0001 + HCO_0001 -> DCOOH_0001 (2body_dust)
    kall(1884) = 4.98114850d-01*((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !OD_0001 + HCO_0001 -> DCOOH_gas (2body_gas)
    kall(1885) = 1.88515012d-03*((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !OD_0001 + HCO_0001 -> HCOOD_0001 (2body_dust)
    kall(1886) = 4.98114850d-01*((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !OD_0001 + HCO_0001 -> HCOOD_gas (2body_gas)
    kall(1887) = 1.88515012d-03*((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-1.20000000d+03*invTd)*1.44098697d+12))*indns

    !OD_0001 + DCO_0001 -> DCOOD_0001 (2body_dust)
    kall(1888) = 9.96302231d-01*((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !OD_0001 + DCO_0001 -> DCOOD_gas (2body_gas)
    kall(1889) = 3.69776938d-03*((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-1.20000000d+03*invTd)*1.41676697d+12))*indns

    !NH2_0001 + OH_0001 -> NH2OH_0001 (2body_dust)
    kall(1890) = 1.00000000d+00*((exp(-1.60000000d+03*invTd)*2.24010528d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !OH_0001 + NHD_0001 -> NH2OD_0001 (2body_dust)
    kall(1891) = 3.33000000d-01*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !NH2_0001 + OD_0001 -> NH2OD_0001 (2body_dust)
    kall(1892) = 3.33000000d-01*((exp(-1.60000000d+03*invTd)*2.24010528d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !OH_0001 + OH_0001 -> HOOH_0001 (2body_dust)
    kall(1893) = 5d-1 * 7.98514322d-01*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !OH_0001 + OH_0001 -> HOOH_gas (2body_gas)
    kall(1894) = 5d-1 * 1.48567773d-03*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !OH_0001 + OD_0001 -> HOOD_0001 (2body_dust)
    kall(1895) = 7.98505724d-01*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !OH_0001 + OD_0001 -> HOOD_gas (2body_gas)
    kall(1896) = 1.49427591d-03*((exp(-2.30000000d+03*invTd)*2.60560084d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !OD_0001 + OD_0001 -> DOOD_0001 (2body_dust)
    kall(1897) = 5d-1 * 7.98583600d-01*((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !OD_0001 + OD_0001 -> DOOD_gas (2body_gas)
    kall(1898) = 5d-1 * 1.41640005d-03*((exp(-2.30000000d+03*invTd)*2.53218886d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !CH_0001 + S_0001 -> HCS_0001 (2body_dust)
    kall(1899) = 9.91314167d-01*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !CH_0001 + S_0001 -> HCS_gas (2body_gas)
    kall(1900) = 8.68583320d-03*((exp(-4.62500000d+02*invTd)*1.33614202d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !CD_0001 + S_0001 -> DCS_0001 (2body_dust)
    kall(1901) = 9.91324930d-01*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !CD_0001 + S_0001 -> DCS_gas (2body_gas)
    kall(1902) = 8.67507004d-03*((exp(-4.62500000d+02*invTd)*1.28753866d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !CH3_0001 + S_0001 -> H_0001 + H2CS_0001 (2body_dust_dust)
    kall(1903) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.63594159d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !CD3_0001 + S_0001 -> D_0001 + D2CS_0001 (2body_dust_dust)
    kall(1904) = 1.00000000d+00*((exp(-8.00000000d+02*invTd)*1.49340352d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !CH2D_0001 + S_0001 -> D_0001 + H2CS_0001 (2body_dust_dust)
    kall(1905) = 3.33000000d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !CH2D_0001 + S_0001 -> H_0001 + HDCS_0001 (2body_dust_dust)
    kall(1906) = 6.67000000d-01*((exp(-8.00000000d+02*invTd)*1.58399363d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !CHD2_0001 + S_0001 -> D_0001 + HDCS_0001 (2body_dust_dust)
    kall(1907) = 6.67000000d-01*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !CHD2_0001 + S_0001 -> H_0001 + D2CS_0001 (2body_dust_dust)
    kall(1908) = 3.33000000d-01*((exp(-1.20000000d+03*invTd)*1.88206488d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !CO_0001 + S_0001 -> OCS_0001 (2body_dust)
    kall(1909) = 9.92252039d-01*((exp(-6.50000000d+02*invTd)*1.07930973d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !CO_0001 + S_0001 -> OCS_gas (2body_gas)
    kall(1910) = 7.74796148d-03*((exp(-6.50000000d+02*invTd)*1.07930973d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !NH_0001 + S_0001 -> H_0001 + NS_0001 (2body_dust_dust)
    kall(1911) = 1.00000000d+00*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !ND_0001 + S_0001 -> D_0001 + NS_0001 (2body_dust_dust)
    kall(1912) = 1.00000000d+00*((exp(-1.30000000d+03*invTd)*2.01920361d+12) + (exp(-1.30000000d+03*invTd)*1.42779256d+12))*indns

    !H_0001 + C2H6_0001 -> p_H2_0001 + C2H5_0001 (2body_dust_dust)
    kall(1913) = 2.50000000d-01 * 3.54191744d+12*max(6.65809812d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.65809812d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !H_0001 + C2H6_0001 -> o_H2_0001 + C2H5_0001 (2body_dust_dust)
    kall(1914) = 7.50000000d-01 * 3.54191744d+12*max(6.65809812d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.65809812d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !H_0001 + C2H5D_0001 -> p_H2_0001 + C2H4D_0001 (2body_dust_dust)
    kall(1915) = 2.50000000d-01 * 2.93979147d+12*max(6.56160128d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.56160128d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))*indns

    !H_0001 + C2H5D_0001 -> o_H2_0001 + C2H4D_0001 (2body_dust_dust)
    kall(1916) = 7.50000000d-01 * 2.93979147d+12*max(6.56160128d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.56160128d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))*indns

    !H_0001 + C2H4D2_0001 -> p_D2_0001 + C2H5_0001 (2body_dust_dust)
    kall(1917) = 3.33333333d-01 * 1.70012037d+11*max(6.47227185d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.47227185d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))*indns

    !H_0001 + C2H4D2_0001 -> o_D2_0001 + C2H5_0001 (2body_dust_dust)
    kall(1918) = 6.66666667d-01 * 1.70012037d+11*max(6.47227185d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.47227185d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))*indns

    !H_0001 + C2H4D2_0001 -> p_H2_0001 + C2H3D2_0001 (2body_dust_dust)
    kall(1919) = 2.50000000d-01 * 1.68595270d+12*max(6.47227185d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.47227185d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))*indns

    !H_0001 + C2H4D2_0001 -> o_H2_0001 + C2H3D2_0001 (2body_dust_dust)
    kall(1920) = 7.50000000d-01 * 1.68595270d+12*max(6.47227185d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.47227185d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))*indns

    !H_0001 + C2H3D3_0001 -> p_D2_0001 + C2H4D_0001 (2body_dust_dust)
    kall(1921) = 3.33333333d-01 * 5.06494193d+11*max(6.38934738d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.38934738d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))*indns

    !H_0001 + C2H3D3_0001 -> o_D2_0001 + C2H4D_0001 (2body_dust_dust)
    kall(1922) = 6.66666667d-01 * 5.06494193d+11*max(6.38934738d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.38934738d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))*indns

    !H_0001 + C2H3D3_0001 -> p_H2_0001 + C2H2D3_0001 (2body_dust_dust)
    kall(1923) = 2.50000000d-01 * 1.01298839d+12*max(6.38934738d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.38934738d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))*indns

    !H_0001 + C2H3D3_0001 -> o_H2_0001 + C2H2D3_0001 (2body_dust_dust)
    kall(1924) = 7.50000000d-01 * 1.01298839d+12*max(6.38934738d-13, exp(-4.89000000d+03 * invTd))&
        *(3.54191744d+12*max(6.38934738d-13, exp(-4.89000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))*indns

    !D_0001 + C2H6_0001 -> p_H2_0001 + C2H4D_0001 (2body_dust_dust)
    kall(1925) = 2.50000000d-01 * 1.82283276d+12*max(1.12414099d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.12414099d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !D_0001 + C2H6_0001 -> o_H2_0001 + C2H4D_0001 (2body_dust_dust)
    kall(1926) = 7.50000000d-01 * 1.82283276d+12*max(1.12414099d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.12414099d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.15678539d+12))*indns

    !D_0001 + C2H5D_0001 -> p_D2_0001 + C2H5_0001 (2body_dust_dust)
    kall(1927) = 3.33333333d-01 * 1.22715249d+11*max(1.08070967d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.08070967d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))*indns

    !D_0001 + C2H5D_0001 -> o_D2_0001 + C2H5_0001 (2body_dust_dust)
    kall(1928) = 6.66666667d-01 * 1.22715249d+11*max(1.08070967d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.08070967d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))*indns

    !D_0001 + C2H5D_0001 -> p_H2_0001 + C2H3D2_0001 (2body_dust_dust)
    kall(1929) = 2.50000000d-01 * 1.21692622d+12*max(1.08070967d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.08070967d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))*indns

    !D_0001 + C2H5D_0001 -> o_H2_0001 + C2H3D2_0001 (2body_dust_dust)
    kall(1930) = 7.50000000d-01 * 1.21692622d+12*max(1.08070967d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.08070967d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.13797462d+12))*indns

    !D_0001 + C2H4D2_0001 -> p_D2_0001 + C2H4D_0001 (2body_dust_dust)
    kall(1931) = 3.33333333d-01 * 3.65589179d+11*max(1.04140490d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.04140490d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))*indns

    !D_0001 + C2H4D2_0001 -> o_D2_0001 + C2H4D_0001 (2body_dust_dust)
    kall(1932) = 6.66666667d-01 * 3.65589179d+11*max(1.04140490d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.04140490d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))*indns

    !D_0001 + C2H4D2_0001 -> p_H2_0001 + C2H2D3_0001 (2body_dust_dust)
    kall(1933) = 2.50000000d-01 * 7.31178358d+11*max(1.04140490d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.04140490d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))*indns

    !D_0001 + C2H4D2_0001 -> o_H2_0001 + C2H2D3_0001 (2body_dust_dust)
    kall(1934) = 7.50000000d-01 * 7.31178358d+11*max(1.04140490d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.04140490d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.12005264d+12))*indns

    !D_0001 + C2H3D3_0001 -> p_D2_0001 + C2H3D2_0001 (2body_dust_dust)
    kall(1935) = 3.33333333d-01 * 7.31178358d+11*max(1.00568863d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.00568863d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))*indns

    !D_0001 + C2H3D3_0001 -> o_D2_0001 + C2H3D2_0001 (2body_dust_dust)
    kall(1936) = 6.66666667d-01 * 7.31178358d+11*max(1.00568863d-17, exp(-4.89000000d+03 * invTd))&
        *(2.55656768d+12*max(1.00568863d-17, exp(-4.89000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-8.00000000d+02*invTd)*1.10295159d+12))*indns

    !H_0001 + CH4_0001 -> p_H2_0001 + CH3_0001 (2body_dust_dust)
    kall(1937) = 2.50000000d-01 * 3.54191744d+12*max(1.90614046d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.90614046d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.22695619d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.22695619d+12))*indns

    !H_0001 + CH4_0001 -> o_H2_0001 + CH3_0001 (2body_dust_dust)
    kall(1938) = 7.50000000d-01 * 3.54191744d+12*max(1.90614046d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.90614046d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.22695619d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.22695619d+12))*indns

    !H_0001 + CH3D_0001 -> p_H2_0001 + CH2D_0001 (2body_dust_dust)
    kall(1939) = 2.50000000d-01 * 2.12515046d+12*max(1.69956299d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.69956299d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))*indns

    !H_0001 + CH3D_0001 -> o_H2_0001 + CH2D_0001 (2body_dust_dust)
    kall(1940) = 7.50000000d-01 * 2.12515046d+12*max(1.69956299d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.69956299d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))*indns

    !H_0001 + CHD3_0001 -> p_H2_0001 + CD3_0001 (2body_dust_dust)
    kall(1941) = 2.50000000d-01 * 3.54191744d+11*max(1.39908687d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.39908687d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))*indns

    !H_0001 + CHD3_0001 -> o_H2_0001 + CD3_0001 (2body_dust_dust)
    kall(1942) = 7.50000000d-01 * 3.54191744d+11*max(1.39908687d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.39908687d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))*indns

    !H_0001 + CHD3_0001 -> p_D2_0001 + CH2D_0001 (2body_dust_dust)
    kall(1943) = 3.33333333d-01 * 1.06257523d+12*max(1.39908687d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.39908687d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))*indns

    !H_0001 + CHD3_0001 -> o_D2_0001 + CH2D_0001 (2body_dust_dust)
    kall(1944) = 6.66666667d-01 * 1.06257523d+12*max(1.39908687d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.39908687d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))*indns

    !H_0001 + CH2D2_0001 -> p_D2_0001 + CH3_0001 (2body_dust_dust)
    kall(1945) = 3.33333333d-01 * 3.54191744d+11*max(1.53403970d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.53403970d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))*indns

    !H_0001 + CH2D2_0001 -> o_D2_0001 + CH3_0001 (2body_dust_dust)
    kall(1946) = 6.66666667d-01 * 3.54191744d+11*max(1.53403970d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.53403970d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))*indns

    !H_0001 + CH2D2_0001 -> p_H2_0001 + CHD2_0001 (2body_dust_dust)
    kall(1947) = 2.50000000d-01 * 1.06257523d+12*max(1.53403970d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.53403970d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))*indns

    !H_0001 + CH2D2_0001 -> o_H2_0001 + CHD2_0001 (2body_dust_dust)
    kall(1948) = 7.50000000d-01 * 1.06257523d+12*max(1.53403970d-29, exp(-5.94000000d+03 * invTd))&
        *(3.54191744d+12*max(1.53403970d-29, exp(-5.94000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))*indns

    !D_0001 + CH4_0001 -> p_H2_0001 + CH2D_0001 (2body_dust_dust)
    kall(1949) = 2.50000000d-01 * 1.53394061d+12*exp(-5.94000000d+03 * invTd)&
        *(2.55656768d+12*exp(-5.94000000d+03 * invTd) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.22695619d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.22695619d+12))*indns

    !D_0001 + CH4_0001 -> o_H2_0001 + CH2D_0001 (2body_dust_dust)
    kall(1950) = 7.50000000d-01 * 1.53394061d+12*exp(-5.94000000d+03 * invTd)&
        *(2.55656768d+12*exp(-5.94000000d+03 * invTd) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.22695619d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.22695619d+12))*indns

    !D_0001 + CH3D_0001 -> p_D2_0001 + CH3_0001 (2body_dust_dust)
    kall(1951) = 3.33333333d-01 * 2.55656768d+11*exp(-5.94000000d+03 * invTd)&
        *(2.55656768d+12*exp(-5.94000000d+03 * invTd) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))*indns

    !D_0001 + CH3D_0001 -> o_D2_0001 + CH3_0001 (2body_dust_dust)
    kall(1952) = 6.66666667d-01 * 2.55656768d+11*exp(-5.94000000d+03 * invTd)&
        *(2.55656768d+12*exp(-5.94000000d+03 * invTd) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))*indns

    !D_0001 + CH3D_0001 -> p_H2_0001 + CHD2_0001 (2body_dust_dust)
    kall(1953) = 2.50000000d-01 * 7.66970305d+11*exp(-5.94000000d+03 * invTd)&
        *(2.55656768d+12*exp(-5.94000000d+03 * invTd) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))*indns

    !D_0001 + CH3D_0001 -> o_H2_0001 + CHD2_0001 (2body_dust_dust)
    kall(1954) = 7.50000000d-01 * 7.66970305d+11*exp(-5.94000000d+03 * invTd)&
        *(2.55656768d+12*exp(-5.94000000d+03 * invTd) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.19032235d+12))*indns

    !D_0001 + CHD3_0001 -> p_D2_0001 + CHD2_0001 (2body_dust_dust)
    kall(1955) = 3.33333333d-01 * 1.53394061d+12*exp(-5.94000000d+03 * invTd)&
        *(2.55656768d+12*exp(-5.94000000d+03 * invTd) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))*indns

    !D_0001 + CHD3_0001 -> o_D2_0001 + CHD2_0001 (2body_dust_dust)
    kall(1956) = 6.66666667d-01 * 1.53394061d+12*exp(-5.94000000d+03 * invTd)&
        *(2.55656768d+12*exp(-5.94000000d+03 * invTd) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.12593222d+12))*indns

    !D_0001 + CH2D2_0001 -> p_H2_0001 + CD3_0001 (2body_dust_dust)
    kall(1957) = 2.50000000d-01 * 2.55656768d+11*exp(-5.94000000d+03 * invTd)&
        *(2.55656768d+12*exp(-5.94000000d+03 * invTd) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))*indns

    !D_0001 + CH2D2_0001 -> o_H2_0001 + CD3_0001 (2body_dust_dust)
    kall(1958) = 7.50000000d-01 * 2.55656768d+11*exp(-5.94000000d+03 * invTd)&
        *(2.55656768d+12*exp(-5.94000000d+03 * invTd) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))*indns

    !D_0001 + CH2D2_0001 -> p_D2_0001 + CH2D_0001 (2body_dust_dust)
    kall(1959) = 3.33333333d-01 * 7.66970305d+11*max(4.97698370d-19, exp(-5.94000000d+03 * invTd))&
        *(2.55656768d+12*max(4.97698370d-19, exp(-5.94000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))*indns

    !D_0001 + CH2D2_0001 -> o_D2_0001 + CH2D_0001 (2body_dust_dust)
    kall(1960) = 6.66666667d-01 * 7.66970305d+11*max(4.97698370d-19, exp(-5.94000000d+03 * invTd))&
        *(2.55656768d+12*max(4.97698370d-19, exp(-5.94000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-4.80000000d+02*invTd)*1.15678539d+12))*indns

    !H_0001 + H_0001 -> p_H2_0001 (2body_dust)
    kall(1961) = 2.50000000d-01 * 5d-1 * 9.90351539d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+02*invTd)*3.54191744d+12))*indns

    !H_0001 + H_0001 -> p_H2_gas (2body_gas)
    kall(1962) = 2.50000000d-01 * 5d-1 * 9.64846081d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+02*invTd)*3.54191744d+12))*indns

    !H_0001 + H_0001 -> o_H2_0001 (2body_dust)
    kall(1963) = 7.50000000d-01 * 5d-1 * 9.90351539d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+02*invTd)*3.54191744d+12))*indns

    !H_0001 + H_0001 -> o_H2_gas (2body_gas)
    kall(1964) = 7.50000000d-01 * 5d-1 * 9.64846081d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+02*invTd)*3.54191744d+12))*indns

    !D_0001 + D_0001 -> p_D2_0001 (2body_dust)
    kall(1965) = 3.33333333d-01 * 5d-1 * 9.90344808d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.60500000d+02*invTd)*2.55656768d+12))*indns

    !D_0001 + D_0001 -> p_D2_gas (2body_gas)
    kall(1966) = 3.33333333d-01 * 5d-1 * 9.65519157d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.60500000d+02*invTd)*2.55656768d+12))*indns

    !D_0001 + D_0001 -> o_D2_0001 (2body_dust)
    kall(1967) = 6.66666667d-01 * 5d-1 * 9.90344808d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.60500000d+02*invTd)*2.55656768d+12))*indns

    !D_0001 + D_0001 -> o_D2_gas (2body_gas)
    kall(1968) = 6.66666667d-01 * 5d-1 * 9.65519157d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.60500000d+02*invTd)*2.55656768d+12))*indns

    !H_0001 + HOOH_0001 -> p_H2_0001 + O2H_0001 (2body_dust_dust)
    kall(1969) = 2.50000000d-01 * 1.77095872d+12*max(2.48547454d-08, exp(-1.90000000d+03 * invTd))&
        *(3.54191744d+12*max(2.48547454d-08, exp(-1.90000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))*indns

    !H_0001 + HOOH_0001 -> o_H2_0001 + O2H_0001 (2body_dust_dust)
    kall(1970) = 7.50000000d-01 * 1.77095872d+12*max(2.48547454d-08, exp(-1.90000000d+03 * invTd))&
        *(3.54191744d+12*max(2.48547454d-08, exp(-1.90000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))*indns

    !H_0001 + HOOD_0001 -> p_H2_0001 + O2D_0001 (2body_dust_dust)
    kall(1971) = 2.50000000d-01 * 5.91500212d+11*max(2.46776330d-08, exp(-1.90000000d+03 * invTd))&
        *(3.54191744d+12*max(2.46776330d-08, exp(-1.90000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))*indns

    !H_0001 + HOOD_0001 -> o_H2_0001 + O2D_0001 (2body_dust_dust)
    kall(1972) = 7.50000000d-01 * 5.91500212d+11*max(2.46776330d-08, exp(-1.90000000d+03 * invTd))&
        *(3.54191744d+12*max(2.46776330d-08, exp(-1.90000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))*indns

    !H_0001 + DOOD_0001 -> p_D2_0001 + O2H_0001 (2body_dust_dust)
    kall(1973) = 3.33333333d-01 * 5.91500212d+11*max(2.45113219d-08, exp(-1.90000000d+03 * invTd))&
        *(3.54191744d+12*max(2.45113219d-08, exp(-1.90000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))*indns

    !H_0001 + DOOD_0001 -> o_D2_0001 + O2H_0001 (2body_dust_dust)
    kall(1974) = 6.66666667d-01 * 5.91500212d+11*max(2.45113219d-08, exp(-1.90000000d+03 * invTd))&
        *(3.54191744d+12*max(2.45113219d-08, exp(-1.90000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))*indns

    !D_0001 + HOOH_0001 -> p_H2_0001 + O2D_0001 (2body_dust_dust)
    kall(1975) = 2.50000000d-01 * 4.26946803d+11*max(2.48831255d-11, exp(-1.90000000d+03 * invTd))&
        *(2.55656768d+12*max(2.48831255d-11, exp(-1.90000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))*indns

    !D_0001 + HOOH_0001 -> o_H2_0001 + O2D_0001 (2body_dust_dust)
    kall(1976) = 7.50000000d-01 * 4.26946803d+11*max(2.48831255d-11, exp(-1.90000000d+03 * invTd))&
        *(2.55656768d+12*max(2.48831255d-11, exp(-1.90000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.10421251d+12))*indns

    !D_0001 + HOOD_0001 -> p_D2_0001 + O2H_0001 (2body_dust_dust)
    kall(1977) = 3.33333333d-01 * 4.26946803d+11*max(2.44050077d-11, exp(-1.90000000d+03 * invTd))&
        *(2.55656768d+12*max(2.44050077d-11, exp(-1.90000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))*indns

    !D_0001 + HOOD_0001 -> o_D2_0001 + O2H_0001 (2body_dust_dust)
    kall(1978) = 6.66666667d-01 * 4.26946803d+11*max(2.44050077d-11, exp(-1.90000000d+03 * invTd))&
        *(2.55656768d+12*max(2.44050077d-11, exp(-1.90000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.07393449d+12))*indns

    !D_0001 + DOOD_0001 -> p_D2_0001 + O2D_0001 (2body_dust_dust)
    kall(1979) = 3.33333333d-01 * 1.27828384d+12*max(2.39608714d-11, exp(-1.90000000d+03 * invTd))&
        *(2.55656768d+12*max(2.39608714d-11, exp(-1.90000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))*indns

    !D_0001 + DOOD_0001 -> o_D2_0001 + O2D_0001 (2body_dust_dust)
    kall(1980) = 6.66666667d-01 * 1.27828384d+12*max(2.39608714d-11, exp(-1.90000000d+03 * invTd))&
        *(2.55656768d+12*max(2.39608714d-11, exp(-1.90000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-3.00000000d+03*invTd)*2.04492698d+12))*indns

    !H_0001 + H2S_0001 -> p_H2_0001 + HS_0001 (2body_dust_dust)
    kall(1981) = 2.50000000d-01 * 3.54191744d+12*max(1.28626883d-07, exp(-1.56000000d+03 * invTd))&
        *(3.54191744d+12*max(1.28626883d-07, exp(-1.56000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.42274437d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.42274437d+12))*indns

    !H_0001 + H2S_0001 -> o_H2_0001 + HS_0001 (2body_dust_dust)
    kall(1982) = 7.50000000d-01 * 3.54191744d+12*max(1.28626883d-07, exp(-1.56000000d+03 * invTd))&
        *(3.54191744d+12*max(1.28626883d-07, exp(-1.56000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.42274437d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.42274437d+12))*indns

    !H_0001 + HDS_0001 -> p_H2_0001 + DS_0001 (2body_dust_dust)
    kall(1983) = 2.50000000d-01 * 1.17945851d+12*max(1.27796071d-07, exp(-1.56000000d+03 * invTd))&
        *(3.54191744d+12*max(1.27796071d-07, exp(-1.56000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.40227216d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.40227216d+12))*indns

    !H_0001 + HDS_0001 -> o_H2_0001 + DS_0001 (2body_dust_dust)
    kall(1984) = 7.50000000d-01 * 1.17945851d+12*max(1.27796071d-07, exp(-1.56000000d+03 * invTd))&
        *(3.54191744d+12*max(1.27796071d-07, exp(-1.56000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.40227216d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.40227216d+12))*indns

    !H_0001 + D2S_0001 -> p_D2_0001 + HS_0001 (2body_dust_dust)
    kall(1985) = 3.33333333d-01 * 1.17945851d+12*max(1.27015418d-07, exp(-1.56000000d+03 * invTd))&
        *(3.54191744d+12*max(1.27015418d-07, exp(-1.56000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.38265900d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.38265900d+12))*indns

    !H_0001 + D2S_0001 -> o_D2_0001 + HS_0001 (2body_dust_dust)
    kall(1986) = 6.66666667d-01 * 1.17945851d+12*max(1.27015418d-07, exp(-1.56000000d+03 * invTd))&
        *(3.54191744d+12*max(1.27015418d-07, exp(-1.56000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.38265900d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.37150000d+03*invTd)*1.38265900d+12))*indns

    !D_0001 + H2S_0001 -> p_H2_0001 + DS_0001 (2body_dust_dust)
    kall(1987) = 2.50000000d-01 * 8.51337039d+11*max(2.46276056d-10, exp(-1.56000000d+03 * invTd))&
        *(2.55656768d+12*max(2.46276056d-10, exp(-1.56000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.42274437d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.42274437d+12))*indns

    !D_0001 + H2S_0001 -> o_H2_0001 + DS_0001 (2body_dust_dust)
    kall(1988) = 7.50000000d-01 * 8.51337039d+11*max(2.46276056d-10, exp(-1.56000000d+03 * invTd))&
        *(2.55656768d+12*max(2.46276056d-10, exp(-1.56000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.42274437d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.42274437d+12))*indns

    !D_0001 + HDS_0001 -> p_D2_0001 + HS_0001 (2body_dust_dust)
    kall(1989) = 3.33333333d-01 * 8.51337039d+11*max(2.41984331d-10, exp(-1.56000000d+03 * invTd))&
        *(2.55656768d+12*max(2.41984331d-10, exp(-1.56000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.40227216d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.40227216d+12))*indns

    !D_0001 + HDS_0001 -> o_D2_0001 + HS_0001 (2body_dust_dust)
    kall(1990) = 6.66666667d-01 * 8.51337039d+11*max(2.41984331d-10, exp(-1.56000000d+03 * invTd))&
        *(2.55656768d+12*max(2.41984331d-10, exp(-1.56000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.40227216d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.40227216d+12))*indns

    !D_0001 + D2S_0001 -> p_D2_0001 + DS_0001 (2body_dust_dust)
    kall(1991) = 3.33333333d-01 * 2.55656768d+12*max(2.37990557d-10, exp(-1.56000000d+03 * invTd))&
        *(2.55656768d+12*max(2.37990557d-10, exp(-1.56000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.38265900d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.38265900d+12))*indns

    !D_0001 + D2S_0001 -> o_D2_0001 + DS_0001 (2body_dust_dust)
    kall(1992) = 6.66666667d-01 * 2.55656768d+12*max(2.37990557d-10, exp(-1.56000000d+03 * invTd))&
        *(2.55656768d+12*max(2.37990557d-10, exp(-1.56000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.38265900d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.37150000d+03*invTd)*1.38265900d+12))*indns

    !H_0001 + HNO_0001 -> p_H2_0001 + NO_0001 (2body_dust_dust)
    kall(1993) = 2.50000000d-01 * 3.54191744d+12*max(1.78840475d-07, exp(-1.50000000d+03 * invTd))&
        *(3.54191744d+12*max(1.78840475d-07, exp(-1.50000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))*indns

    !H_0001 + HNO_0001 -> o_H2_0001 + NO_0001 (2body_dust_dust)
    kall(1994) = 7.50000000d-01 * 3.54191744d+12*max(1.78840475d-07, exp(-1.50000000d+03 * invTd))&
        *(3.54191744d+12*max(1.78840475d-07, exp(-1.50000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))*indns

    !D_0001 + HNO_0001 -> HD_0001 + NO_0001 (2body_dust_dust)
    kall(1995) = 2.55656768d+12*max(4.01088790d-10, exp(-1.50000000d+03 * invTd))&
        *(2.55656768d+12*max(4.01088790d-10, exp(-1.50000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.50000000d+03*invTd)*1.55823592d+12))*indns

    !H_0001 + DNO_0001 -> HD_0001 + NO_0001 (2body_dust_dust)
    kall(1996) = 3.54191744d+12*max(1.77487881d-07, exp(-1.50000000d+03 * invTd))&
        *(3.54191744d+12*max(1.77487881d-07, exp(-1.50000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))*indns

    !D_0001 + DNO_0001 -> p_D2_0001 + NO_0001 (2body_dust_dust)
    kall(1997) = 3.33333333d-01 * 2.55656768d+12*max(3.92942875d-10, exp(-1.50000000d+03 * invTd))&
        *(2.55656768d+12*max(3.92942875d-10, exp(-1.50000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))*indns

    !D_0001 + DNO_0001 -> o_D2_0001 + NO_0001 (2body_dust_dust)
    kall(1998) = 6.66666667d-01 * 2.55656768d+12*max(3.92942875d-10, exp(-1.50000000d+03 * invTd))&
        *(2.55656768d+12*max(3.92942875d-10, exp(-1.50000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.50000000d+03*invTd)*1.53369524d+12))*indns

    !H_0001 + N2H2_0001 -> H_0001 + p_H2_0001 + N2_0001 (2body_dust_dust)
    kall(1999) = 3.54191744d+12*max(3.63530490d-05, exp(-6.50000000d+02 * invTd))&
        *(3.54191744d+12*max(3.63530490d-05, exp(-6.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.99440671d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.99440671d+12))*indns

    !H_0001 + N2H2_0001 -> H_0001 + o_H2_0001 + N2_0001 (2body_dust_dust)
    kall(2000) = 3.54191744d+12*max(3.63530490d-05, exp(-6.50000000d+02 * invTd))&
        *(3.54191744d+12*max(3.63530490d-05, exp(-6.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.99440671d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.99440671d+12))*indns

    !H_0001 + N2HD_0001 -> p_H2_0001 + D_0001 + N2_0001 (2body_dust_dust)
    kall(2001) = 2.50000000d-01 * 1.17945851d+12*max(3.61600669d-05, exp(-6.50000000d+02 * invTd))&
        *(3.54191744d+12*max(3.61600669d-05, exp(-6.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.96197517d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.96197517d+12))*indns

    !H_0001 + N2HD_0001 -> o_H2_0001 + D_0001 + N2_0001 (2body_dust_dust)
    kall(2002) = 7.50000000d-01 * 1.17945851d+12*max(3.61600669d-05, exp(-6.50000000d+02 * invTd))&
        *(3.54191744d+12*max(3.61600669d-05, exp(-6.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.96197517d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.96197517d+12))*indns

    !H_0001 + N2D2_0001 -> H_0001 + p_D2_0001 + N2_0001 (2body_dust_dust)
    kall(2003) = 1.17945851d+12*max(3.59798047d-05, exp(-6.50000000d+02 * invTd))&
        *(3.54191744d+12*max(3.59798047d-05, exp(-6.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.93107599d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.93107599d+12))*indns

    !H_0001 + N2D2_0001 -> H_0001 + o_D2_0001 + N2_0001 (2body_dust_dust)
    kall(2004) = 1.17945851d+12*max(3.59798047d-05, exp(-6.50000000d+02 * invTd))&
        *(3.54191744d+12*max(3.59798047d-05, exp(-6.50000000d+02 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.93107599d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.37800000d+03*invTd)*1.93107599d+12))*indns

    !D_0001 + N2H2_0001 -> p_H2_0001 + D_0001 + N2_0001 (2body_dust_dust)
    kall(2005) = 2.50000000d-01 * 8.51337039d+11*max(6.61506897d-07, exp(-6.50000000d+02 * invTd))&
        *(2.55656768d+12*max(6.61506897d-07, exp(-6.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.99440671d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.99440671d+12))*indns

    !D_0001 + N2H2_0001 -> o_H2_0001 + D_0001 + N2_0001 (2body_dust_dust)
    kall(2006) = 7.50000000d-01 * 8.51337039d+11*max(6.61506897d-07, exp(-6.50000000d+02 * invTd))&
        *(2.55656768d+12*max(6.61506897d-07, exp(-6.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.99440671d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.99440671d+12))*indns

    !D_0001 + N2HD_0001 -> H_0001 + p_D2_0001 + N2_0001 (2body_dust_dust)
    kall(2007) = 8.51337039d+11*max(6.52072134d-07, exp(-6.50000000d+02 * invTd))&
        *(2.55656768d+12*max(6.52072134d-07, exp(-6.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.96197517d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.96197517d+12))*indns

    !D_0001 + N2HD_0001 -> H_0001 + o_D2_0001 + N2_0001 (2body_dust_dust)
    kall(2008) = 8.51337039d+11*max(6.52072134d-07, exp(-6.50000000d+02 * invTd))&
        *(2.55656768d+12*max(6.52072134d-07, exp(-6.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.96197517d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.96197517d+12))*indns

    !D_0001 + N2D2_0001 -> D_0001 + p_D2_0001 + N2_0001 (2body_dust_dust)
    kall(2009) = 2.55656768d+12*max(6.43323821d-07, exp(-6.50000000d+02 * invTd))&
        *(2.55656768d+12*max(6.43323821d-07, exp(-6.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.93107599d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.93107599d+12))*indns

    !D_0001 + N2D2_0001 -> D_0001 + o_D2_0001 + N2_0001 (2body_dust_dust)
    kall(2010) = 2.55656768d+12*max(6.43323821d-07, exp(-6.50000000d+02 * invTd))&
        *(2.55656768d+12*max(6.43323821d-07, exp(-6.50000000d+02 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.93107599d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.37800000d+03*invTd)*1.93107599d+12))*indns

    !p_H2_0001 + C_0001 -> CH2_0001 (2body_dust)
    kall(2011) = 4.53203615d+12*max(2.58277949d-12, exp(-2.50000000d+03 * invTd))&
        *(4.57259575d+12*max(2.58277949d-12, exp(-2.50000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !p_H2_0001 + C_0001 -> CH2_gas (2body_gas)
    kall(2012) = 4.05595928d+10*max(2.58277949d-12, exp(-2.50000000d+03 * invTd))&
        *(4.57259575d+12*max(2.58277949d-12, exp(-2.50000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !o_H2_0001 + C_0001 -> CH2_0001 (2body_dust)
    kall(2013) = 4.53203615d+12*max(2.58277949d-12, exp(-2.50000000d+03 * invTd))&
        *(4.57259575d+12*max(2.58277949d-12, exp(-2.50000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !o_H2_0001 + C_0001 -> CH2_gas (2body_gas)
    kall(2014) = 4.05595928d+10*max(2.58277949d-12, exp(-2.50000000d+03 * invTd))&
        *(4.57259575d+12*max(2.58277949d-12, exp(-2.50000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !p_D2_0001 + C_0001 -> CD2_0001 (2body_dust)
    kall(2015) = 4.53183651d+12*max(4.68416764d-16, exp(-2.50000000d+03 * invTd))&
        *(4.57259575d+12*max(4.68416764d-16, exp(-2.50000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !p_D2_0001 + C_0001 -> CD2_gas (2body_gas)
    kall(2016) = 4.07592373d+10*max(4.68416764d-16, exp(-2.50000000d+03 * invTd))&
        *(4.57259575d+12*max(4.68416764d-16, exp(-2.50000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !o_D2_0001 + C_0001 -> CD2_0001 (2body_dust)
    kall(2017) = 4.53183651d+12*max(4.68416764d-16, exp(-2.50000000d+03 * invTd))&
        *(4.57259575d+12*max(4.68416764d-16, exp(-2.50000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !o_D2_0001 + C_0001 -> CD2_gas (2body_gas)
    kall(2018) = 4.07592373d+10*max(4.68416764d-16, exp(-2.50000000d+03 * invTd))&
        *(4.57259575d+12*max(4.68416764d-16, exp(-2.50000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-5.00000000d+03*invTd)*4.57259575d+12))*indns

    !p_H2_0001 + C2_0001 -> H_0001 + CCH_0001 (2body_dust_dust)
    kall(2019) = 3.23331346d+12*max(2.59042532d-16, exp(-4.20000000d+03 * invTd))&
        *(3.23331346d+12*max(2.59042532d-16, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !o_H2_0001 + C2_0001 -> H_0001 + CCH_0001 (2body_dust_dust)
    kall(2020) = 3.23331346d+12*max(2.59042532d-16, exp(-4.20000000d+03 * invTd))&
        *(3.23331346d+12*max(2.59042532d-16, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !p_D2_0001 + C2_0001 -> D_0001 + CCD_0001 (2body_dust_dust)
    kall(2021) = 3.23331346d+12*max(5.74123358d-22, exp(-4.20000000d+03 * invTd))&
        *(3.23331346d+12*max(5.74123358d-22, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !o_D2_0001 + C2_0001 -> D_0001 + CCD_0001 (2body_dust_dust)
    kall(2022) = 3.23331346d+12*max(5.74123358d-22, exp(-4.20000000d+03 * invTd))&
        *(3.23331346d+12*max(5.74123358d-22, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-5.00000000d+03*invTd)*3.23331346d+12))*indns

    !p_H2_0001 + CCH_0001 -> H_0001 + C2H2_0001 (2body_dust_dust)
    kall(2023) = 2.37599045d+12*max(2.45095968d-16, exp(-4.20000000d+03 * invTd))&
        *(2.37599045d+12*max(2.45095968d-16, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !o_H2_0001 + CCH_0001 -> H_0001 + C2H2_0001 (2body_dust_dust)
    kall(2024) = 2.37599045d+12*max(2.45095968d-16, exp(-4.20000000d+03 * invTd))&
        *(2.37599045d+12*max(2.45095968d-16, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !p_H2_0001 + CCD_0001 -> D_0001 + C2H2_0001 (2body_dust_dust)
    kall(2025) = 7.91204819d+11*max(2.32836534d-16, exp(-4.20000000d+03 * invTd))&
        *(2.37599045d+12*max(2.32836534d-16, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !o_H2_0001 + CCD_0001 -> D_0001 + C2H2_0001 (2body_dust_dust)
    kall(2026) = 7.91204819d+11*max(2.32836534d-16, exp(-4.20000000d+03 * invTd))&
        *(2.37599045d+12*max(2.32836534d-16, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !p_H2_0001 + CCD_0001 -> H_0001 + C2HD_0001 (2body_dust_dust)
    kall(2027) = 1.58478563d+12*max(2.32836534d-16, exp(-4.20000000d+03 * invTd))&
        *(2.37599045d+12*max(2.32836534d-16, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !o_H2_0001 + CCD_0001 -> H_0001 + C2HD_0001 (2body_dust_dust)
    kall(2028) = 1.58478563d+12*max(2.32836534d-16, exp(-4.20000000d+03 * invTd))&
        *(2.37599045d+12*max(2.32836534d-16, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !p_D2_0001 + CCH_0001 -> D_0001 + C2HD_0001 (2body_dust_dust)
    kall(2029) = 1.15736378d+12*max(4.98948041d-22, exp(-4.20000000d+03 * invTd))&
        *(1.73517809d+12*max(4.98948041d-22, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !o_D2_0001 + CCH_0001 -> D_0001 + C2HD_0001 (2body_dust_dust)
    kall(2030) = 1.15736378d+12*max(4.98948041d-22, exp(-4.20000000d+03 * invTd))&
        *(1.73517809d+12*max(4.98948041d-22, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !p_D2_0001 + CCH_0001 -> H_0001 + C2D2_0001 (2body_dust_dust)
    kall(2031) = 5.77814303d+11*max(4.98948041d-22, exp(-4.20000000d+03 * invTd))&
        *(1.73517809d+12*max(4.98948041d-22, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !o_D2_0001 + CCH_0001 -> H_0001 + C2D2_0001 (2body_dust_dust)
    kall(2032) = 5.77814303d+11*max(4.98948041d-22, exp(-4.20000000d+03 * invTd))&
        *(1.73517809d+12*max(4.98948041d-22, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.50000000d+03*invTd)*1.73517809d+12))*indns

    !p_D2_0001 + CCD_0001 -> D_0001 + C2D2_0001 (2body_dust_dust)
    kall(2033) = 1.70148210d+12*max(4.37850323d-22, exp(-4.20000000d+03 * invTd))&
        *(1.70148210d+12*max(4.37850323d-22, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !o_D2_0001 + CCD_0001 -> D_0001 + C2D2_0001 (2body_dust_dust)
    kall(2034) = 1.70148210d+12*max(4.37850323d-22, exp(-4.20000000d+03 * invTd))&
        *(1.70148210d+12*max(4.37850323d-22, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.50000000d+03*invTd)*1.70148210d+12))*indns

    !p_H2_0001 + C3_0001 -> H_0001 + l_C3H_0001 (2body_dust_dust)
    kall(2035) = 2.37599045d+12*max(1.62037541d-16, exp(-4.20000000d+03 * invTd))&
        *(2.37599045d+12*max(1.62037541d-16, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !p_D2_0001 + C3_0001 -> D_0001 + l_C3D_0001 (2body_dust_dust)
    kall(2036) = 1.68007896d+12*max(1.71576270d-22, exp(-4.20000000d+03 * invTd))&
        *(1.68007896d+12*max(1.71576270d-22, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.20000000d+03*invTd)*1.29332538d+12))*indns

    !p_H2_0001 + l_C3H_0001 -> H_0001 + l_C3H2_0001 (2body_dust_dust)
    kall(2037) = 2.37599045d+12*max(1.57896657d-16, exp(-4.20000000d+03 * invTd))&
        *(2.37599045d+12*max(1.57896657d-16, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))*indns

    !p_H2_0001 + l_C3D_0001 -> D_0001 + l_C3H2_0001 (2body_dust_dust)
    kall(2038) = 2.37599045d+12*max(1.54063502d-16, exp(-4.20000000d+03 * invTd))&
        *(2.37599045d+12*max(1.54063502d-16, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))*indns

    !p_H2_0001 + l_C3D_0001 -> H_0001 + l_C3HD_0001 (2body_dust_dust)
    kall(2039) = 2.37599045d+12*max(1.54063502d-16, exp(-4.20000000d+03 * invTd))&
        *(2.37599045d+12*max(1.54063502d-16, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))*indns

    !p_D2_0001 + l_C3H_0001 -> D_0001 + l_C3HD_0001 (2body_dust_dust)
    kall(2040) = 1.68007896d+12*max(1.60318844d-22, exp(-4.20000000d+03 * invTd))&
        *(1.68007896d+12*max(1.60318844d-22, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))*indns

    !p_D2_0001 + l_C3H_0001 -> H_0001 + l_C3D2_0001 (2body_dust_dust)
    kall(2041) = 1.68007896d+12*max(1.60318844d-22, exp(-4.20000000d+03 * invTd))&
        *(1.68007896d+12*max(1.60318844d-22, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-2.00000000d+03*invTd)*1.64695815d+12))*indns

    !p_D2_0001 + l_C3D_0001 -> D_0001 + l_C3D2_0001 (2body_dust_dust)
    kall(2042) = 1.68007896d+12*max(1.50297718d-22, exp(-4.20000000d+03 * invTd))&
        *(1.68007896d+12*max(1.50297718d-22, exp(-4.20000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-2.00000000d+03*invTd)*1.62514318d+12))*indns

    !p_H2_0001 + CH2_0001 -> H_0001 + CH3_0001 (2body_dust_dust)
    kall(2043) = 2.37599045d+12*max(1.22367585d-14, exp(-3.53000000d+03 * invTd))&
        *(2.37599045d+12*max(1.22367585d-14, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !o_H2_0001 + CH2_0001 -> H_0001 + CH3_0001 (2body_dust_dust)
    kall(2044) = 2.37599045d+12*max(1.22367585d-14, exp(-3.53000000d+03 * invTd))&
        *(2.37599045d+12*max(1.22367585d-14, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !p_H2_0001 + CHD_0001 -> D_0001 + CH3_0001 (2body_dust_dust)
    kall(2045) = 5.93997612d+11*max(1.06987613d-14, exp(-3.53000000d+03 * invTd))&
        *(2.37599045d+12*max(1.06987613d-14, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !o_H2_0001 + CHD_0001 -> D_0001 + CH3_0001 (2body_dust_dust)
    kall(2046) = 5.93997612d+11*max(1.06987613d-14, exp(-3.53000000d+03 * invTd))&
        *(2.37599045d+12*max(1.06987613d-14, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !p_H2_0001 + CHD_0001 -> H_0001 + CH2D_0001 (2body_dust_dust)
    kall(2047) = 1.78199284d+12*max(1.06987613d-14, exp(-3.53000000d+03 * invTd))&
        *(2.37599045d+12*max(1.06987613d-14, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !o_H2_0001 + CHD_0001 -> H_0001 + CH2D_0001 (2body_dust_dust)
    kall(2048) = 1.78199284d+12*max(1.06987613d-14, exp(-3.53000000d+03 * invTd))&
        *(2.37599045d+12*max(1.06987613d-14, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !p_H2_0001 + CD2_0001 -> D_0001 + CH2D_0001 (2body_dust_dust)
    kall(2049) = 1.18799522d+12*max(9.49917211d-15, exp(-3.53000000d+03 * invTd))&
        *(2.37599045d+12*max(9.49917211d-15, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !o_H2_0001 + CD2_0001 -> D_0001 + CH2D_0001 (2body_dust_dust)
    kall(2050) = 1.18799522d+12*max(9.49917211d-15, exp(-3.53000000d+03 * invTd))&
        *(2.37599045d+12*max(9.49917211d-15, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !p_H2_0001 + CD2_0001 -> H_0001 + CHD2_0001 (2body_dust_dust)
    kall(2051) = 1.18799522d+12*max(9.49917211d-15, exp(-3.53000000d+03 * invTd))&
        *(2.37599045d+12*max(9.49917211d-15, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !o_H2_0001 + CD2_0001 -> H_0001 + CHD2_0001 (2body_dust_dust)
    kall(2052) = 1.18799522d+12*max(9.49917211d-15, exp(-3.53000000d+03 * invTd))&
        *(2.37599045d+12*max(9.49917211d-15, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !p_D2_0001 + CH2_0001 -> D_0001 + CH2D_0001 (2body_dust_dust)
    kall(2053) = 8.40039479d+11*max(2.81982303d-19, exp(-3.53000000d+03 * invTd))&
        *(1.68007896d+12*max(2.81982303d-19, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !o_D2_0001 + CH2_0001 -> D_0001 + CH2D_0001 (2body_dust_dust)
    kall(2054) = 8.40039479d+11*max(2.81982303d-19, exp(-3.53000000d+03 * invTd))&
        *(1.68007896d+12*max(2.81982303d-19, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !p_D2_0001 + CH2_0001 -> H_0001 + CHD2_0001 (2body_dust_dust)
    kall(2055) = 8.40039479d+11*max(2.81982303d-19, exp(-3.53000000d+03 * invTd))&
        *(1.68007896d+12*max(2.81982303d-19, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !o_D2_0001 + CH2_0001 -> H_0001 + CHD2_0001 (2body_dust_dust)
    kall(2056) = 8.40039479d+11*max(2.81982303d-19, exp(-3.53000000d+03 * invTd))&
        *(1.68007896d+12*max(2.81982303d-19, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.58399363d+12))*indns

    !p_D2_0001 + CHD_0001 -> H_0001 + CD3_0001 (2body_dust_dust)
    kall(2057) = 4.20019739d+11*max(2.04771854d-19, exp(-3.53000000d+03 * invTd))&
        *(1.68007896d+12*max(2.04771854d-19, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !o_D2_0001 + CHD_0001 -> H_0001 + CD3_0001 (2body_dust_dust)
    kall(2058) = 4.20019739d+11*max(2.04771854d-19, exp(-3.53000000d+03 * invTd))&
        *(1.68007896d+12*max(2.04771854d-19, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !p_D2_0001 + CHD_0001 -> D_0001 + CHD2_0001 (2body_dust_dust)
    kall(2059) = 1.26005922d+12*max(2.04771854d-19, exp(-3.53000000d+03 * invTd))&
        *(1.68007896d+12*max(2.04771854d-19, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !o_D2_0001 + CHD_0001 -> D_0001 + CHD2_0001 (2body_dust_dust)
    kall(2060) = 1.26005922d+12*max(2.04771854d-19, exp(-3.53000000d+03 * invTd))&
        *(1.68007896d+12*max(2.04771854d-19, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.53028323d+12))*indns

    !p_D2_0001 + CD2_0001 -> D_0001 + CD3_0001 (2body_dust_dust)
    kall(2061) = 1.68007896d+12*max(1.53847800d-19, exp(-3.53000000d+03 * invTd))&
        *(1.68007896d+12*max(1.53847800d-19, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !o_D2_0001 + CD2_0001 -> D_0001 + CD3_0001 (2body_dust_dust)
    kall(2062) = 1.68007896d+12*max(1.53847800d-19, exp(-3.53000000d+03 * invTd))&
        *(1.68007896d+12*max(1.53847800d-19, exp(-3.53000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-7.00000000d+02*invTd)*1.48169037d+12))*indns

    !p_H2_0001 + CH3_0001 -> H_0001 + CH4_0001 (2body_dust_dust)
    kall(2063) = 2.37599045d+12*max(1.34885405d-19, exp(-6.44000000d+03 * invTd))&
        *(2.37599045d+12*max(1.34885405d-19, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !o_H2_0001 + CH3_0001 -> H_0001 + CH4_0001 (2body_dust_dust)
    kall(2064) = 2.37599045d+12*max(1.34885405d-19, exp(-6.44000000d+03 * invTd))&
        *(2.37599045d+12*max(1.34885405d-19, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !p_H2_0001 + CD3_0001 -> H_0001 + CHD3_0001 (2body_dust_dust)
    kall(2065) = 9.50396179d+11*max(8.75384382d-20, exp(-6.44000000d+03 * invTd))&
        *(2.37599045d+12*max(8.75384382d-20, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !o_H2_0001 + CD3_0001 -> H_0001 + CHD3_0001 (2body_dust_dust)
    kall(2066) = 9.50396179d+11*max(8.75384382d-20, exp(-6.44000000d+03 * invTd))&
        *(2.37599045d+12*max(8.75384382d-20, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !p_H2_0001 + CD3_0001 -> D_0001 + CH2D2_0001 (2body_dust_dust)
    kall(2067) = 1.42559427d+12*max(8.75384382d-20, exp(-6.44000000d+03 * invTd))&
        *(2.37599045d+12*max(8.75384382d-20, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !o_H2_0001 + CD3_0001 -> D_0001 + CH2D2_0001 (2body_dust_dust)
    kall(2068) = 1.42559427d+12*max(8.75384382d-20, exp(-6.44000000d+03 * invTd))&
        *(2.37599045d+12*max(8.75384382d-20, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.49340352d+12))*indns

    !p_H2_0001 + CH2D_0001 -> D_0001 + CH4_0001 (2body_dust_dust)
    kall(2069) = 4.75198089d+11*max(1.14869542d-19, exp(-6.44000000d+03 * invTd))&
        *(2.37599045d+12*max(1.14869542d-19, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !o_H2_0001 + CH2D_0001 -> D_0001 + CH4_0001 (2body_dust_dust)
    kall(2070) = 4.75198089d+11*max(1.14869542d-19, exp(-6.44000000d+03 * invTd))&
        *(2.37599045d+12*max(1.14869542d-19, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !p_H2_0001 + CH2D_0001 -> H_0001 + CH3D_0001 (2body_dust_dust)
    kall(2071) = 1.90079236d+12*max(1.14869542d-19, exp(-6.44000000d+03 * invTd))&
        *(2.37599045d+12*max(1.14869542d-19, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !o_H2_0001 + CH2D_0001 -> H_0001 + CH3D_0001 (2body_dust_dust)
    kall(2072) = 1.90079236d+12*max(1.14869542d-19, exp(-6.44000000d+03 * invTd))&
        *(2.37599045d+12*max(1.14869542d-19, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !p_H2_0001 + CHD2_0001 -> D_0001 + CH3D_0001 (2body_dust_dust)
    kall(2073) = 9.50396179d+11*max(9.95417030d-20, exp(-6.44000000d+03 * invTd))&
        *(2.37599045d+12*max(9.95417030d-20, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !o_H2_0001 + CHD2_0001 -> D_0001 + CH3D_0001 (2body_dust_dust)
    kall(2074) = 9.50396179d+11*max(9.95417030d-20, exp(-6.44000000d+03 * invTd))&
        *(2.37599045d+12*max(9.95417030d-20, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !p_H2_0001 + CHD2_0001 -> H_0001 + CH2D2_0001 (2body_dust_dust)
    kall(2075) = 1.42559427d+12*max(9.95417030d-20, exp(-6.44000000d+03 * invTd))&
        *(2.37599045d+12*max(9.95417030d-20, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !o_H2_0001 + CHD2_0001 -> H_0001 + CH2D2_0001 (2body_dust_dust)
    kall(2076) = 1.42559427d+12*max(9.95417030d-20, exp(-6.44000000d+03 * invTd))&
        *(2.37599045d+12*max(9.95417030d-20, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !p_D2_0001 + CH3_0001 -> D_0001 + CH3D_0001 (2body_dust_dust)
    kall(2077) = 6.72031583d+11*max(5.71907310d-26, exp(-6.44000000d+03 * invTd))&
        *(1.68007896d+12*max(5.71907310d-26, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !o_D2_0001 + CH3_0001 -> D_0001 + CH3D_0001 (2body_dust_dust)
    kall(2078) = 6.72031583d+11*max(5.71907310d-26, exp(-6.44000000d+03 * invTd))&
        *(1.68007896d+12*max(5.71907310d-26, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !p_D2_0001 + CH3_0001 -> H_0001 + CH2D2_0001 (2body_dust_dust)
    kall(2079) = 1.00804737d+12*max(5.71907310d-26, exp(-6.44000000d+03 * invTd))&
        *(1.68007896d+12*max(5.71907310d-26, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !o_D2_0001 + CH3_0001 -> H_0001 + CH2D2_0001 (2body_dust_dust)
    kall(2080) = 1.00804737d+12*max(5.71907310d-26, exp(-6.44000000d+03 * invTd))&
        *(1.68007896d+12*max(5.71907310d-26, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.63594159d+12))*indns

    !p_D2_0001 + CH2D_0001 -> H_0001 + CHD3_0001 (2body_dust_dust)
    kall(2081) = 6.72031583d+11*max(3.88685547d-26, exp(-6.44000000d+03 * invTd))&
        *(1.68007896d+12*max(3.88685547d-26, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !o_D2_0001 + CH2D_0001 -> H_0001 + CHD3_0001 (2body_dust_dust)
    kall(2082) = 6.72031583d+11*max(3.88685547d-26, exp(-6.44000000d+03 * invTd))&
        *(1.68007896d+12*max(3.88685547d-26, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !p_D2_0001 + CH2D_0001 -> D_0001 + CH2D2_0001 (2body_dust_dust)
    kall(2083) = 1.00804737d+12*max(3.88685547d-26, exp(-6.44000000d+03 * invTd))&
        *(1.68007896d+12*max(3.88685547d-26, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !o_D2_0001 + CH2D_0001 -> D_0001 + CH2D2_0001 (2body_dust_dust)
    kall(2084) = 1.00804737d+12*max(3.88685547d-26, exp(-6.44000000d+03 * invTd))&
        *(1.68007896d+12*max(3.88685547d-26, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-8.00000000d+02*invTd)*1.58399363d+12))*indns

    !p_D2_0001 + CHD2_0001 -> D_0001 + CHD3_0001 (2body_dust_dust)
    kall(2085) = 1.50565191d+12*max(2.74658787d-26, exp(-6.44000000d+03 * invTd))&
        *(1.88206488d+12*max(2.74658787d-26, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !o_D2_0001 + CHD2_0001 -> D_0001 + CHD3_0001 (2body_dust_dust)
    kall(2086) = 1.50565191d+12*max(2.74658787d-26, exp(-6.44000000d+03 * invTd))&
        *(1.88206488d+12*max(2.74658787d-26, exp(-6.44000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.20000000d+03*invTd)*1.88206488d+12))*indns

    !p_H2_0001 + CN_0001 -> H_0001 + HCN_0001 (2body_dust_dust)
    kall(2087) = 2.37599045d+12*max(1.05944443d-11, exp(-2.07000000d+03 * invTd))&
        *(2.37599045d+12*max(1.05944443d-11, exp(-2.07000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !o_H2_0001 + CN_0001 -> H_0001 + HCN_0001 (2body_dust_dust)
    kall(2088) = 2.37599045d+12*max(1.05944443d-11, exp(-2.07000000d+03 * invTd))&
        *(2.37599045d+12*max(1.05944443d-11, exp(-2.07000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !p_D2_0001 + CN_0001 -> D_0001 + DCN_0001 (2body_dust_dust)
    kall(2089) = 1.68007896d+12*max(1.01251671d-15, exp(-2.07000000d+03 * invTd))&
        *(1.68007896d+12*max(1.01251671d-15, exp(-2.07000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !o_D2_0001 + CN_0001 -> D_0001 + DCN_0001 (2body_dust_dust)
    kall(2090) = 1.68007896d+12*max(1.01251671d-15, exp(-2.07000000d+03 * invTd))&
        *(1.68007896d+12*max(1.01251671d-15, exp(-2.07000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.40000000d+03*invTd)*1.64378788d+12))*indns

    !p_H2_0001 + NH2_0001 -> H_0001 + NH3_0001 (2body_dust_dust)
    kall(2091) = 2.37599045d+12*max(1.85013492d-19, exp(-6.30000000d+03 * invTd))&
        *(2.37599045d+12*max(1.85013492d-19, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !o_H2_0001 + NH2_0001 -> H_0001 + NH3_0001 (2body_dust_dust)
    kall(2092) = 2.37599045d+12*max(1.85013492d-19, exp(-6.30000000d+03 * invTd))&
        *(2.37599045d+12*max(1.85013492d-19, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !p_H2_0001 + NHD_0001 -> H_0001 + NH2D_0001 (2body_dust_dust)
    kall(2093) = 2.37599045d+12*max(1.60577014d-19, exp(-6.30000000d+03 * invTd))&
        *(2.37599045d+12*max(1.60577014d-19, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !o_H2_0001 + NHD_0001 -> H_0001 + NH2D_0001 (2body_dust_dust)
    kall(2094) = 2.37599045d+12*max(1.60577014d-19, exp(-6.30000000d+03 * invTd))&
        *(2.37599045d+12*max(1.60577014d-19, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !p_H2_0001 + ND2_0001 -> H_0001 + NHD2_0001 (2body_dust_dust)
    kall(2095) = 2.37599045d+12*max(1.41412249d-19, exp(-6.30000000d+03 * invTd))&
        *(2.37599045d+12*max(1.41412249d-19, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !o_H2_0001 + ND2_0001 -> H_0001 + NHD2_0001 (2body_dust_dust)
    kall(2096) = 2.37599045d+12*max(1.41412249d-19, exp(-6.30000000d+03 * invTd))&
        *(2.37599045d+12*max(1.41412249d-19, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !p_D2_0001 + NH2_0001 -> D_0001 + NH2D_0001 (2body_dust_dust)
    kall(2097) = 2.24010528d+12*max(7.36743033d-26, exp(-6.30000000d+03 * invTd))&
        *(2.24010528d+12*max(7.36743033d-26, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !o_D2_0001 + NH2_0001 -> D_0001 + NH2D_0001 (2body_dust_dust)
    kall(2098) = 2.24010528d+12*max(7.36743033d-26, exp(-6.30000000d+03 * invTd))&
        *(2.24010528d+12*max(7.36743033d-26, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.60000000d+03*invTd)*2.24010528d+12))*indns

    !p_D2_0001 + NHD_0001 -> D_0001 + NHD2_0001 (2body_dust_dust)
    kall(2099) = 2.17322133d+12*max(5.22587881d-26, exp(-6.30000000d+03 * invTd))&
        *(2.17322133d+12*max(5.22587881d-26, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !o_D2_0001 + NHD_0001 -> D_0001 + NHD2_0001 (2body_dust_dust)
    kall(2100) = 2.17322133d+12*max(5.22587881d-26, exp(-6.30000000d+03 * invTd))&
        *(2.17322133d+12*max(5.22587881d-26, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.60000000d+03*invTd)*2.17322133d+12))*indns

    !p_D2_0001 + ND2_0001 -> D_0001 + ND3_0001 (2body_dust_dust)
    kall(2101) = 2.11199151d+12*max(3.83108594d-26, exp(-6.30000000d+03 * invTd))&
        *(2.11199151d+12*max(3.83108594d-26, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !o_D2_0001 + ND2_0001 -> D_0001 + ND3_0001 (2body_dust_dust)
    kall(2102) = 2.11199151d+12*max(3.83108594d-26, exp(-6.30000000d+03 * invTd))&
        *(2.11199151d+12*max(3.83108594d-26, exp(-6.30000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-1.60000000d+03*invTd)*2.11199151d+12))*indns

    !p_H2_0001 + OH_0001 -> H_0001 + H2O_0001 (2body_dust_dust)
    kall(2103) = 2.60560084d+12*max(1.40960113d-11, exp(-2.10000000d+03 * invTd))&
        *(2.60560084d+12*max(1.40960113d-11, exp(-2.10000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !o_H2_0001 + OH_0001 -> H_0001 + H2O_0001 (2body_dust_dust)
    kall(2104) = 2.60560084d+12*max(1.40960113d-11, exp(-2.10000000d+03 * invTd))&
        *(2.60560084d+12*max(1.40960113d-11, exp(-2.10000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !p_H2_0001 + OD_0001 -> H_0001 + HDO_0001 (2body_dust_dust)
    kall(2105) = 2.53218886d+12*max(1.30987125d-11, exp(-2.10000000d+03 * invTd))&
        *(2.53218886d+12*max(1.30987125d-11, exp(-2.10000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !o_H2_0001 + OD_0001 -> H_0001 + HDO_0001 (2body_dust_dust)
    kall(2106) = 2.53218886d+12*max(1.30987125d-11, exp(-2.10000000d+03 * invTd))&
        *(2.53218886d+12*max(1.30987125d-11, exp(-2.10000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*2.37599045d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !p_D2_0001 + OH_0001 -> D_0001 + HDO_0001 (2body_dust_dust)
    kall(2107) = 2.60560084d+12*max(2.53234274d-15, exp(-2.10000000d+03 * invTd))&
        *(2.60560084d+12*max(2.53234274d-15, exp(-2.10000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !o_D2_0001 + OH_0001 -> D_0001 + HDO_0001 (2body_dust_dust)
    kall(2108) = 2.60560084d+12*max(2.53234274d-15, exp(-2.10000000d+03 * invTd))&
        *(2.60560084d+12*max(2.53234274d-15, exp(-2.10000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-2.30000000d+03*invTd)*2.60560084d+12))*indns

    !p_D2_0001 + OD_0001 -> D_0001 + D2O_0001 (2body_dust_dust)
    kall(2109) = 2.53218886d+12*max(2.11677186d-15, exp(-2.10000000d+03 * invTd))&
        *(2.53218886d+12*max(2.11677186d-15, exp(-2.10000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !o_D2_0001 + OD_0001 -> D_0001 + D2O_0001 (2body_dust_dust)
    kall(2110) = 2.53218886d+12*max(2.11677186d-15, exp(-2.10000000d+03 * invTd))&
        *(2.53218886d+12*max(2.11677186d-15, exp(-2.10000000d+03 * invTd)) + (exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))**(-1d0) &
        *((exp(-2.25000000d+02*invTd)*1.68007896d+12) + (exp(-2.30000000d+03*invTd)*2.53218886d+12))*indns

    !NH_0001 + NH_0001 -> p_H2_0001 + N2_0001 (2body_dust_dust)
    kall(2111) = 2.50000000d-01 * 5d-1 * 5.00000000d-01*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !NH_0001 + NH_0001 -> o_H2_0001 + N2_0001 (2body_dust_dust)
    kall(2112) = 7.50000000d-01 * 5d-1 * 5.00000000d-01*((exp(-1.30000000d+03*invTd)*2.08542452d+12) + (exp(-1.30000000d+03*invTd)*2.08542452d+12))*indns

    !ND_0001 + ND_0001 -> p_D2_0001 + N2_0001 (2body_dust_dust)
    kall(2113) = 3.33333333d-01 * 5d-1 * 5.00000000d-01*((exp(-1.30000000d+03*invTd)*2.01920361d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !ND_0001 + ND_0001 -> o_D2_0001 + N2_0001 (2body_dust_dust)
    kall(2114) = 6.66666667d-01 * 5d-1 * 5.00000000d-01*((exp(-1.30000000d+03*invTd)*2.01920361d+12) + (exp(-1.30000000d+03*invTd)*2.01920361d+12))*indns

    !H_0001 + CO_0001 -> HCO_0001 (2body_dust)
    kall(2115) = 3.53134247d+12*max(4.04608501d-18, exp(-2.50000000d+03 * invTd))&
        *(3.54191744d+12*max(4.04608501d-18, exp(-2.50000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))*indns

    !H_0001 + CO_0001 -> HCO_gas (2body_gas)
    kall(2116) = 1.05749673d+10*max(4.04608501d-18, exp(-2.50000000d+03 * invTd))&
        *(3.54191744d+12*max(4.04608501d-18, exp(-2.50000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))*indns

    !D_0001 + CO_0001 -> DCO_0001 (2body_dust)
    kall(2117) = 2.54730901d+12*max(5.93729953d-19, exp(-1.42000000d+03 * invTd))&
        *(2.55656768d+12*max(5.93729953d-19, exp(-1.42000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))*indns

    !D_0001 + CO_0001 -> DCO_gas (2body_gas)
    kall(2118) = 9.25867501d+09*max(5.93729953d-19, exp(-1.42000000d+03 * invTd))&
        *(2.55656768d+12*max(5.93729953d-19, exp(-1.42000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-6.50000000d+02*invTd)*1.07930973d+12))*indns

    !H_0001 + H2CO_0001 -> CH2OH_0001 (2body_dust)
    kall(2119) = 3.54051267d+12*max(2.04369635d-19, exp(-2.88000000d+03 * invTd))&
        *(3.54191744d+12*max(2.04369635d-19, exp(-2.88000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !H_0001 + H2CO_0001 -> CH2OH_gas (2body_gas)
    kall(2120) = 1.40476776d+09*max(2.04369635d-19, exp(-2.88000000d+03 * invTd))&
        *(3.54191744d+12*max(2.04369635d-19, exp(-2.88000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !H_0001 + HDCO_0001 -> CHDOH_0001 (2body_dust)
    kall(2121) = 3.54191744d+12*max(1.99841069d-19, exp(-2.88000000d+03 * invTd))&
        *(3.54191744d+12*max(1.99841069d-19, exp(-2.88000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !H_0001 + HDCO_0001 -> CHDOH_gas (2body_gas)
    kall(2122) = 0d0*max(1.99841069d-19, exp(-2.88000000d+03 * invTd))&
        *(3.54191744d+12*max(1.99841069d-19, exp(-2.88000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !H_0001 + D2CO_0001 -> CD2OH_0001 (2body_dust)
    kall(2123) = 3.54191744d+12*max(1.95680490d-19, exp(-2.88000000d+03 * invTd))&
        *(3.54191744d+12*max(1.95680490d-19, exp(-2.88000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))*indns

    !H_0001 + D2CO_0001 -> CD2OH_gas (2body_gas)
    kall(2124) = 0d0*max(1.95680490d-19, exp(-2.88000000d+03 * invTd))&
        *(3.54191744d+12*max(1.95680490d-19, exp(-2.88000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))*indns

    !D_0001 + H2CO_0001 -> CH2OD_0001 (2body_dust)
    kall(2125) = 2.55535387d+12*max(2.68394514d-20, exp(-1.63000000d+03 * invTd))&
        *(2.55656768d+12*max(2.68394514d-20, exp(-1.63000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !D_0001 + H2CO_0001 -> CH2OD_gas (2body_gas)
    kall(2126) = 1.21381583d+09*max(2.68394514d-20, exp(-1.63000000d+03 * invTd))&
        *(2.55656768d+12*max(2.68394514d-20, exp(-1.63000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !D_0001 + HDCO_0001 -> CHDOD_0001 (2body_dust)
    kall(2127) = 2.55656768d+12*max(2.94579278d-20, exp(-1.62000000d+03 * invTd))&
        *(2.55656768d+12*max(2.94579278d-20, exp(-1.62000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !D_0001 + HDCO_0001 -> CHDOD_gas (2body_gas)
    kall(2128) = 0d0*max(2.94579278d-20, exp(-1.62000000d+03 * invTd))&
        *(2.55656768d+12*max(2.94579278d-20, exp(-1.62000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !D_0001 + D2CO_0001 -> CD2OD_0001 (2body_dust)
    kall(2129) = 2.55656768d+12*max(2.82280454d-20, exp(-1.62000000d+03 * invTd))&
        *(2.55656768d+12*max(2.82280454d-20, exp(-1.62000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))*indns

    !D_0001 + D2CO_0001 -> CD2OD_gas (2body_gas)
    kall(2130) = 0d0*max(2.82280454d-20, exp(-1.62000000d+03 * invTd))&
        *(2.55656768d+12*max(2.82280454d-20, exp(-1.62000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))*indns

    !H_0001 + H2CO_0001 -> CH3O_0001 (2body_dust)
    kall(2131) = 3.54173888d+12*max(2.04453519d-18, exp(-2.58000000d+03 * invTd))&
        *(3.54191744d+12*max(2.04453519d-18, exp(-2.58000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !H_0001 + H2CO_0001 -> CH3O_gas (2body_gas)
    kall(2132) = 1.78550856d+08*max(2.04453519d-18, exp(-2.58000000d+03 * invTd))&
        *(3.54191744d+12*max(2.04453519d-18, exp(-2.58000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !H_0001 + HDCO_0001 -> CH2DO_0001 (2body_dust)
    kall(2133) = 3.54191744d+12*max(2.34485777d-18, exp(-2.56000000d+03 * invTd))&
        *(3.54191744d+12*max(2.34485777d-18, exp(-2.56000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !H_0001 + HDCO_0001 -> CH2DO_gas (2body_gas)
    kall(2134) = 0d0*max(2.34485777d-18, exp(-2.56000000d+03 * invTd))&
        *(3.54191744d+12*max(2.34485777d-18, exp(-2.56000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !H_0001 + D2CO_0001 -> CHD2O_0001 (2body_dust)
    kall(2135) = 3.54191744d+12*max(2.69486700d-18, exp(-2.54000000d+03 * invTd))&
        *(3.54191744d+12*max(2.69486700d-18, exp(-2.54000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))*indns

    !H_0001 + D2CO_0001 -> CHD2O_gas (2body_gas)
    kall(2136) = 0d0*max(2.69486700d-18, exp(-2.54000000d+03 * invTd))&
        *(3.54191744d+12*max(2.69486700d-18, exp(-2.54000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))*indns

    !D_0001 + H2CO_0001 -> CH2DO_0001 (2body_dust)
    kall(2137) = 2.55535387d+12*max(2.59495341d-19, exp(-1.47000000d+03 * invTd))&
        *(2.55656768d+12*max(2.59495341d-19, exp(-1.47000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !D_0001 + H2CO_0001 -> CH2DO_gas (2body_gas)
    kall(2138) = 1.21381583d+09*max(2.59495341d-19, exp(-1.47000000d+03 * invTd))&
        *(2.55656768d+12*max(2.59495341d-19, exp(-1.47000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !D_0001 + HDCO_0001 -> CHD2O_0001 (2body_dust)
    kall(2139) = 2.55656768d+12*max(3.32936191d-19, exp(-1.45000000d+03 * invTd))&
        *(2.55656768d+12*max(3.32936191d-19, exp(-1.45000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !D_0001 + HDCO_0001 -> CHD2O_gas (2body_gas)
    kall(2140) = 0d0*max(3.32936191d-19, exp(-1.45000000d+03 * invTd))&
        *(2.55656768d+12*max(3.32936191d-19, exp(-1.45000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !D_0001 + D2CO_0001 -> CD3O_0001 (2body_dust)
    kall(2141) = 2.55656768d+12*max(3.70446169d-19, exp(-1.44000000d+03 * invTd))&
        *(2.55656768d+12*max(3.70446169d-19, exp(-1.44000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))*indns

    !D_0001 + D2CO_0001 -> CD3O_gas (2body_gas)
    kall(2142) = 0d0*max(3.70446169d-19, exp(-1.44000000d+03 * invTd))&
        *(2.55656768d+12*max(3.70446169d-19, exp(-1.44000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))*indns

    !H_0001 + H2CO_0001 -> p_H2_0001 + HCO_0001 (2body_dust_dust)
    kall(2143) = 2.50000000d-01 * 3.54191744d+12*max(3.86374587d-18, exp(-2.50000000d+03 * invTd))&
        *(3.54191744d+12*max(3.86374587d-18, exp(-2.50000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !H_0001 + H2CO_0001 -> o_H2_0001 + HCO_0001 (2body_dust_dust)
    kall(2144) = 7.50000000d-01 * 3.54191744d+12*max(3.86374587d-18, exp(-2.50000000d+03 * invTd))&
        *(3.54191744d+12*max(3.86374587d-18, exp(-2.50000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !D_0001 + H2CO_0001 -> HD_0001 + HCO_0001 (2body_dust_dust)
    kall(2145) = 2.55656768d+12*max(2.26241778d-18, exp(-1.32500000d+03 * invTd))&
        *(2.55656768d+12*max(2.26241778d-18, exp(-1.32500000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !H_0001 + HDCO_0001 -> HD_0001 + HCO_0001 (2body_dust_dust)
    kall(2146) = 3.54191744d+12*max(1.51963541d-18, exp(-2.61500000d+03 * invTd))&
        *(3.54191744d+12*max(1.51963541d-18, exp(-2.61500000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !D_0001 + HDCO_0001 -> p_D2_0001 + HCO_0001 (2body_dust_dust)
    kall(2147) = 3.33333333d-01 * 2.55656768d+12*max(6.98140596d-18, exp(-1.25000000d+03 * invTd))&
        *(2.55656768d+12*max(6.98140596d-18, exp(-1.25000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !D_0001 + HDCO_0001 -> o_D2_0001 + HCO_0001 (2body_dust_dust)
    kall(2148) = 6.66666667d-01 * 2.55656768d+12*max(6.98140596d-18, exp(-1.25000000d+03 * invTd))&
        *(2.55656768d+12*max(6.98140596d-18, exp(-1.25000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !H_0001 + HDCO_0001 -> p_H2_0001 + DCO_0001 (2body_dust_dust)
    kall(2149) = 2.50000000d-01 * 3.54191744d+12*max(3.78391741d-18, exp(-2.50000000d+03 * invTd))&
        *(3.54191744d+12*max(3.78391741d-18, exp(-2.50000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !H_0001 + HDCO_0001 -> o_H2_0001 + DCO_0001 (2body_dust_dust)
    kall(2150) = 7.50000000d-01 * 3.54191744d+12*max(3.78391741d-18, exp(-2.50000000d+03 * invTd))&
        *(3.54191744d+12*max(3.78391741d-18, exp(-2.50000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !D_0001 + HDCO_0001 -> HD_0001 + DCO_0001 (2body_dust_dust)
    kall(2151) = 2.55656768d+12*max(8.56155516d-25, exp(-2.46000000d+03 * invTd))&
        *(2.55656768d+12*max(8.56155516d-25, exp(-2.46000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !H_0001 + D2CO_0001 -> HD_0001 + DCO_0001 (2body_dust_dust)
    kall(2152) = 3.54191744d+12*max(1.48947324d-18, exp(-2.61500000d+03 * invTd))&
        *(3.54191744d+12*max(1.48947324d-18, exp(-2.61500000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))*indns

    !D_0001 + D2CO_0001 -> p_D2_0001 + DCO_0001 (2body_dust_dust)
    kall(2153) = 3.33333333d-01 * 2.55656768d+12*max(6.72470926d-18, exp(-1.25000000d+03 * invTd))&
        *(2.55656768d+12*max(6.72470926d-18, exp(-1.25000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))*indns

    !D_0001 + D2CO_0001 -> o_D2_0001 + DCO_0001 (2body_dust_dust)
    kall(2154) = 6.66666667d-01 * 2.55656768d+12*max(6.72470926d-18, exp(-1.25000000d+03 * invTd))&
        *(2.55656768d+12*max(6.72470926d-18, exp(-1.25000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.87838538d+12))*indns

    !D_0001 + H2CO_0001 -> H_0001 + HDCO_0001 (2body_dust_dust)
    kall(2155) = 2.55656768d+12*max(2.84929143d-18, exp(-1.31000000d+03 * invTd))&
        *(2.55656768d+12*max(2.84929143d-18, exp(-1.31000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.93998808d+12))*indns

    !D_0001 + HDCO_0001 -> H_0001 + D2CO_0001 (2body_dust_dust)
    kall(2156) = 2.55656768d+12*max(2.73541564d-18, exp(-1.31000000d+03 * invTd))&
        *(2.55656768d+12*max(2.73541564d-18, exp(-1.31000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.25000000d+03*invTd)*1.90844145d+12))*indns

    !H_0001 + CH3OH_0001 -> p_H2_0001 + CH2OH_0001 (2body_dust_dust)
    kall(2157) = 2.50000000d-01 * 1.77095872d+12*max(3.71046678d-18, exp(-2.50000000d+03 * invTd))&
        *(3.54191744d+12*max(3.71046678d-18, exp(-2.50000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))*indns

    !H_0001 + CH3OH_0001 -> o_H2_0001 + CH2OH_0001 (2body_dust_dust)
    kall(2158) = 7.50000000d-01 * 1.77095872d+12*max(3.71046678d-18, exp(-2.50000000d+03 * invTd))&
        *(3.54191744d+12*max(3.71046678d-18, exp(-2.50000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))*indns

    !D_0001 + CH3OH_0001 -> HD_0001 + CH2OH_0001 (2body_dust_dust)
    kall(2159) = 2.55656768d+12*max(3.58839236d-24, exp(-2.33000000d+03 * invTd))&
        *(2.55656768d+12*max(3.58839236d-24, exp(-2.33000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))*indns

    !D_0001 + CH3OH_0001 -> HD_0001 + CH3O_0001 (2body_dust_dust)
    kall(2160) = 2.55656768d+12*max(1.14043774d-24, exp(-2.43000000d+03 * invTd))&
        *(2.55656768d+12*max(1.14043774d-24, exp(-2.43000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))*indns

    !H_0001 + CH3OH_0001 -> p_H2_0001 + CH3O_0001 (2body_dust_dust)
    kall(2161) = 2.50000000d-01 * 1.77095872d+12*max(2.48877568d-18, exp(-2.55000000d+03 * invTd))&
        *(3.54191744d+12*max(2.48877568d-18, exp(-2.55000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))*indns

    !H_0001 + CH3OH_0001 -> o_H2_0001 + CH3O_0001 (2body_dust_dust)
    kall(2162) = 7.50000000d-01 * 1.77095872d+12*max(2.48877568d-18, exp(-2.55000000d+03 * invTd))&
        *(3.54191744d+12*max(2.48877568d-18, exp(-2.55000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.97999204d+12))*indns

    !H_0001 + CH2DOH_0001 -> p_H2_0001 + CHDOH_0001 (2body_dust_dust)
    kall(2163) = 2.50000000d-01 * 3.54191744d+12*max(2.44285412d-18, exp(-2.55000000d+03 * invTd))&
        *(3.54191744d+12*max(2.44285412d-18, exp(-2.55000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !H_0001 + CH2DOH_0001 -> o_H2_0001 + CHDOH_0001 (2body_dust_dust)
    kall(2164) = 7.50000000d-01 * 3.54191744d+12*max(2.44285412d-18, exp(-2.55000000d+03 * invTd))&
        *(3.54191744d+12*max(2.44285412d-18, exp(-2.55000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !D_0001 + CH2DOH_0001 -> HD_0001 + CHDOH_0001 (2body_dust_dust)
    kall(2165) = 2.55656768d+12*max(2.15480603d-24, exp(-2.37000000d+03 * invTd))&
        *(2.55656768d+12*max(2.15480603d-24, exp(-2.37000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !H_0001 + CHD2OH_0001 -> p_H2_0001 + CD2OH_0001 (2body_dust_dust)
    kall(2166) = 2.50000000d-01 * 3.54191744d+12*max(1.89204958d-18, exp(-2.58000000d+03 * invTd))&
        *(3.54191744d+12*max(1.89204958d-18, exp(-2.58000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !H_0001 + CHD2OH_0001 -> o_H2_0001 + CD2OH_0001 (2body_dust_dust)
    kall(2167) = 7.50000000d-01 * 3.54191744d+12*max(1.89204958d-18, exp(-2.58000000d+03 * invTd))&
        *(3.54191744d+12*max(1.89204958d-18, exp(-2.58000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !D_0001 + CHD2OH_0001 -> HD_0001 + CD2OH_0001 (2body_dust_dust)
    kall(2168) = 2.55656768d+12*max(1.45898905d-24, exp(-2.40000000d+03 * invTd))&
        *(2.55656768d+12*max(1.45898905d-24, exp(-2.40000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !H_0001 + CH3OD_0001 -> p_H2_0001 + CH2OD_0001 (2body_dust_dust)
    kall(2169) = 2.50000000d-01 * 3.54191744d+12*max(3.64267156d-18, exp(-2.50000000d+03 * invTd))&
        *(3.54191744d+12*max(3.64267156d-18, exp(-2.50000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !H_0001 + CH3OD_0001 -> o_H2_0001 + CH2OD_0001 (2body_dust_dust)
    kall(2170) = 7.50000000d-01 * 3.54191744d+12*max(3.64267156d-18, exp(-2.50000000d+03 * invTd))&
        *(3.54191744d+12*max(3.64267156d-18, exp(-2.50000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !D_0001 + CH3OD_0001 -> HD_0001 + CH2OD_0001 (2body_dust_dust)
    kall(2171) = 2.55656768d+12*max(3.41960651d-24, exp(-2.33000000d+03 * invTd))&
        *(2.55656768d+12*max(3.41960651d-24, exp(-2.33000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.94976138d+12))*indns

    !H_0001 + CH2DOD_0001 -> p_H2_0001 + CHDOD_0001 (2body_dust_dust)
    kall(2172) = 2.50000000d-01 * 3.54191744d+12*max(2.40035179d-18, exp(-2.55000000d+03 * invTd))&
        *(3.54191744d+12*max(2.40035179d-18, exp(-2.55000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !H_0001 + CH2DOD_0001 -> o_H2_0001 + CHDOD_0001 (2body_dust_dust)
    kall(2173) = 7.50000000d-01 * 3.54191744d+12*max(2.40035179d-18, exp(-2.55000000d+03 * invTd))&
        *(3.54191744d+12*max(2.40035179d-18, exp(-2.55000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !D_0001 + CH2DOD_0001 -> HD_0001 + CHDOD_0001 (2body_dust_dust)
    kall(2174) = 2.55656768d+12*max(2.05823607d-24, exp(-2.37000000d+03 * invTd))&
        *(2.55656768d+12*max(2.05823607d-24, exp(-2.37000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.92087443d+12))*indns

    !H_0001 + CHD2OD_0001 -> p_H2_0001 + CD2OD_0001 (2body_dust_dust)
    kall(2175) = 2.50000000d-01 * 3.54191744d+12*max(1.86077634d-18, exp(-2.58000000d+03 * invTd))&
        *(3.54191744d+12*max(1.86077634d-18, exp(-2.58000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))*indns

    !H_0001 + CHD2OD_0001 -> o_H2_0001 + CD2OD_0001 (2body_dust_dust)
    kall(2176) = 7.50000000d-01 * 3.54191744d+12*max(1.86077634d-18, exp(-2.58000000d+03 * invTd))&
        *(3.54191744d+12*max(1.86077634d-18, exp(-2.58000000d+03 * invTd)) + (exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))**(-1d0) &
        *((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))*indns

    !D_0001 + CHD2OD_0001 -> HD_0001 + CD2OD_0001 (2body_dust_dust)
    kall(2177) = 2.55656768d+12*max(1.39672872d-24, exp(-2.40000000d+03 * invTd))&
        *(2.55656768d+12*max(1.39672872d-24, exp(-2.40000000d+03 * invTd)) + (exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))**(-1d0) &
        *((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-2.50000000d+03*invTd)*1.89323451d+12))*indns

    !H_0001 + Cl_0001 -> HCl_0001 (2body_dust)
    kall(2178) = 9.92776162d-01*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.50000000d+03*invTd)*1.46649314d+12))*indns

    !H_0001 + Cl_0001 -> HCl_gas (2body_gas)
    kall(2179) = 7.22383809d-03*((exp(-2.50000000d+02*invTd)*3.54191744d+12) + (exp(-1.50000000d+03*invTd)*1.46649314d+12))*indns

    !D_0001 + Cl_0001 -> DCl_0001 (2body_dust)
    kall(2180) = 9.92743027d-01*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.50000000d+03*invTd)*1.46649314d+12))*indns

    !D_0001 + Cl_0001 -> DCl_gas (2body_gas)
    kall(2181) = 7.25697297d-03*((exp(-2.60500000d+02*invTd)*2.55656768d+12) + (exp(-1.50000000d+03*invTd)*1.46649314d+12))*indns

    !H2O_0001 -> OH_0001 + H_gas (photodesorption)
    kall(2182) = Ffuva * 1.19000000d-02

    !H2O_0001 -> OH_gas + H_0001 (photodesorption)
    kall(2183) = Ffuva * 4.81000000d-04

    !H2O_0001 -> OH_gas + H_gas (photodesorption)
    kall(2184) = Ffuva * 1.66000000d-04

    !H2O_0001 -> H2O_gas (photodesorption)
    kall(2185) = Ffuva * 1.17000000d-04

    !H2O_0001 -> OH_0001 + H_0001 (photodesorption)
    kall(2186) = Ffuva * 3.96000000d-03

    !HDO_0001 -> OH_0001 + D_gas (photodesorption)
    kall(2187) = Ffuva * 3.66000000d-03

    !HDO_0001 -> OD_0001 + H_gas (photodesorption)
    kall(2188) = Ffuva * 7.93000000d-03

    !HDO_0001 -> OH_gas + D_0001 (photodesorption)
    kall(2189) = Ffuva * 4.05000000d-05

    !HDO_0001 -> OD_gas + H_0001 (photodesorption)
    kall(2190) = Ffuva * 3.14000000d-05

    !HDO_0001 -> OD_gas + H_gas (photodesorption)
    kall(2191) = Ffuva * 1.40000000d-04

    !HDO_0001 -> OH_gas + D_gas (photodesorption)
    kall(2192) = Ffuva * 1.63000000d-04

    !HDO_0001 -> OD_0001 + H_0001 (photodesorption)
    kall(2193) = Ffuva * 3.11000000d-03

    !HDO_0001 -> OH_0001 + D_0001 (photodesorption)
    kall(2194) = Ffuva * 1.70000000d-03

    !HDO_0001 -> HDO_gas (photodesorption)
    kall(2195) = Ffuva * 6.81000000d-05

    !D2O_0001 -> OD_0001 + D_gas (photodesorption)
    kall(2196) = Ffuva * 1.10000000d-02

    !D2O_0001 -> OD_gas + D_0001 (photodesorption)
    kall(2197) = Ffuva * 1.02000000d-04

    !D2O_0001 -> OD_gas + D_gas (photodesorption)
    kall(2198) = Ffuva * 4.13000000d-04

    !D2O_0001 -> D2O_gas (photodesorption)
    kall(2199) = Ffuva * 8.46000000d-05

    !D2O_0001 -> OD_0001 + D_0001 (photodesorption)
    kall(2200) = Ffuva * 3.68000000d-03

    !CO_0001 -> CO_gas (photodesorption)
    kall(2201) = Ffuva_CO * Y_CO

    !CO2_0001 -> CO2_gas (photodesorption)
    kall(2202) = Ffuva * 5.00000000d-04

    !CO2_0001 -> CO_gas + O_gas (photodesorption)
    kall(2203) = Ffuva * 5.00000000d-04

    !N2_0001 -> N2_gas (photodesorption)
    kall(2204) = Ffuva_N2 * 5.00000000d-04

    !H_0001 -> H_gas (photodesorption)
    kall(2205) = Ffuva * 1.00000000d-03

    !D_0001 -> D_gas (photodesorption)
    kall(2206) = Ffuva * 1.00000000d-03

    !p_H2_0001 -> p_H2_gas (photodesorption)
    kall(2207) = Ffuva_H2 * 1.00000000d-03

    !p_D2_0001 -> p_D2_gas (photodesorption)
    kall(2208) = Ffuva * 1.00000000d-03

    !o_H2_0001 -> o_H2_gas (photodesorption)
    kall(2209) = Ffuva_H2 * 1.00000000d-03

    !o_D2_0001 -> o_D2_gas (photodesorption)
    kall(2210) = Ffuva * 1.00000000d-03

    !HD_0001 -> HD_gas (photodesorption)
    kall(2211) = Ffuva_HD * 1.00000000d-03

    !He_0001 -> He_gas (photodesorption)
    kall(2212) = Ffuva * 1.00000000d-03

    !O_0001 -> O_gas (photodesorption)
    kall(2213) = Ffuva * 1.00000000d-03

    !O2_0001 -> O2_gas (photodesorption)
    kall(2214) = Ffuva * 1.00000000d-03

    !O3_0001 -> O3_gas (photodesorption)
    kall(2215) = Ffuva * 1.00000000d-03

    !OH_0001 -> OH_gas (photodesorption)
    kall(2216) = Ffuva * 1.00000000d-03

    !OD_0001 -> OD_gas (photodesorption)
    kall(2217) = Ffuva * 1.00000000d-03

    !O2H_0001 -> O2H_gas (photodesorption)
    kall(2218) = Ffuva * 1.00000000d-03

    !O2D_0001 -> O2D_gas (photodesorption)
    kall(2219) = Ffuva * 1.00000000d-03

    !HOOH_0001 -> HOOH_gas (photodesorption)
    kall(2220) = Ffuva * 1.00000000d-03

    !HOOD_0001 -> HOOD_gas (photodesorption)
    kall(2221) = Ffuva * 1.00000000d-03

    !DOOD_0001 -> DOOD_gas (photodesorption)
    kall(2222) = Ffuva * 1.00000000d-03

    !Fe_0001 -> Fe_gas (photodesorption)
    kall(2223) = Ffuva * 1.00000000d-03

    !FeH_0001 -> FeH_gas (photodesorption)
    kall(2224) = Ffuva * 1.00000000d-03

    !N_0001 -> N_gas (photodesorption)
    kall(2225) = Ffuva * 1.00000000d-03

    !S_0001 -> S_gas (photodesorption)
    kall(2226) = Ffuva * 1.00000000d-03

    !C_0001 -> C_gas (photodesorption)
    kall(2227) = Ffuva * 1.00000000d-03

    !C2_0001 -> C2_gas (photodesorption)
    kall(2228) = Ffuva * 1.00000000d-03

    !HCO_0001 -> HCO_gas (photodesorption)
    kall(2229) = Ffuva * 1.00000000d-03

    !DCO_0001 -> DCO_gas (photodesorption)
    kall(2230) = Ffuva * 1.00000000d-03

    !H2CO_0001 -> H2CO_gas (photodesorption)
    kall(2231) = Ffuva * 1.00000000d-03

    !HDCO_0001 -> HDCO_gas (photodesorption)
    kall(2232) = Ffuva * 1.00000000d-03

    !D2CO_0001 -> D2CO_gas (photodesorption)
    kall(2233) = Ffuva * 1.00000000d-03

    !CH2OH_0001 -> CH2OH_gas (photodesorption)
    kall(2234) = Ffuva * 1.00000000d-03

    !CD2OD_0001 -> CD2OD_gas (photodesorption)
    kall(2235) = Ffuva * 1.00000000d-03

    !CH2OD_0001 -> CH2OD_gas (photodesorption)
    kall(2236) = Ffuva * 1.00000000d-03

    !CHDOH_0001 -> CHDOH_gas (photodesorption)
    kall(2237) = Ffuva * 1.00000000d-03

    !CHDOD_0001 -> CHDOD_gas (photodesorption)
    kall(2238) = Ffuva * 1.00000000d-03

    !CD2OH_0001 -> CD2OH_gas (photodesorption)
    kall(2239) = Ffuva * 1.00000000d-03

    !CH3O_0001 -> CH3O_gas (photodesorption)
    kall(2240) = Ffuva * 1.00000000d-03

    !CHD2O_0001 -> CHD2O_gas (photodesorption)
    kall(2241) = Ffuva * 1.00000000d-03

    !CH2DO_0001 -> CH2DO_gas (photodesorption)
    kall(2242) = Ffuva * 1.00000000d-03

    !CD3O_0001 -> CD3O_gas (photodesorption)
    kall(2243) = Ffuva * 1.00000000d-03

    !CH3OH_0001 -> CH3OH_gas (photodesorption)
    kall(2244) = Ffuva * 1.00000000d-03

    !CH3OD_0001 -> CH3OD_gas (photodesorption)
    kall(2245) = Ffuva * 1.00000000d-03

    !CHD2OH_0001 -> CHD2OH_gas (photodesorption)
    kall(2246) = Ffuva * 1.00000000d-03

    !CHD2OD_0001 -> CHD2OD_gas (photodesorption)
    kall(2247) = Ffuva * 1.00000000d-03

    !CH2DOH_0001 -> CH2DOH_gas (photodesorption)
    kall(2248) = Ffuva * 1.00000000d-03

    !CH2DOD_0001 -> CH2DOD_gas (photodesorption)
    kall(2249) = Ffuva * 1.00000000d-03

    !CD3OD_0001 -> CD3OD_gas (photodesorption)
    kall(2250) = Ffuva * 1.00000000d-03

    !CD3OH_0001 -> CD3OH_gas (photodesorption)
    kall(2251) = Ffuva * 1.00000000d-03

    !CH_0001 -> CH_gas (photodesorption)
    kall(2252) = Ffuva * 1.00000000d-03

    !CD_0001 -> CD_gas (photodesorption)
    kall(2253) = Ffuva * 1.00000000d-03

    !CH2_0001 -> CH2_gas (photodesorption)
    kall(2254) = Ffuva * 1.00000000d-03

    !CHD_0001 -> CHD_gas (photodesorption)
    kall(2255) = Ffuva * 1.00000000d-03

    !CD2_0001 -> CD2_gas (photodesorption)
    kall(2256) = Ffuva * 1.00000000d-03

    !CH3_0001 -> CH3_gas (photodesorption)
    kall(2257) = Ffuva * 1.00000000d-03

    !CH2D_0001 -> CH2D_gas (photodesorption)
    kall(2258) = Ffuva * 1.00000000d-03

    !CHD2_0001 -> CHD2_gas (photodesorption)
    kall(2259) = Ffuva * 1.00000000d-03

    !CD3_0001 -> CD3_gas (photodesorption)
    kall(2260) = Ffuva * 1.00000000d-03

    !CH4_0001 -> CH4_gas (photodesorption)
    kall(2261) = Ffuva * 1.00000000d-03

    !CH3D_0001 -> CH3D_gas (photodesorption)
    kall(2262) = Ffuva * 1.00000000d-03

    !CH2D2_0001 -> CH2D2_gas (photodesorption)
    kall(2263) = Ffuva * 1.00000000d-03

    !CHD3_0001 -> CHD3_gas (photodesorption)
    kall(2264) = Ffuva * 1.00000000d-03

    !CD4_0001 -> CD4_gas (photodesorption)
    kall(2265) = Ffuva * 1.00000000d-03

    !HCOOH_0001 -> HCOOH_gas (photodesorption)
    kall(2266) = Ffuva * 1.00000000d-03

    !HCOOD_0001 -> HCOOD_gas (photodesorption)
    kall(2267) = Ffuva * 1.00000000d-03

    !DCOOH_0001 -> DCOOH_gas (photodesorption)
    kall(2268) = Ffuva * 1.00000000d-03

    !DCOOD_0001 -> DCOOD_gas (photodesorption)
    kall(2269) = Ffuva * 1.00000000d-03

    !HOCO_0001 -> HOCO_gas (photodesorption)
    kall(2270) = Ffuva * 1.00000000d-03

    !DOCO_0001 -> DOCO_gas (photodesorption)
    kall(2271) = Ffuva * 1.00000000d-03

    !NH_0001 -> NH_gas (photodesorption)
    kall(2272) = Ffuva * 1.00000000d-03

    !ND_0001 -> ND_gas (photodesorption)
    kall(2273) = Ffuva * 1.00000000d-03

    !NH2_0001 -> NH2_gas (photodesorption)
    kall(2274) = Ffuva * 1.00000000d-03

    !NHD_0001 -> NHD_gas (photodesorption)
    kall(2275) = Ffuva * 1.00000000d-03

    !ND2_0001 -> ND2_gas (photodesorption)
    kall(2276) = Ffuva * 1.00000000d-03

    !NH3_0001 -> NH3_gas (photodesorption)
    kall(2277) = Ffuva * 1.00000000d-03

    !NH2D_0001 -> NH2D_gas (photodesorption)
    kall(2278) = Ffuva * 1.00000000d-03

    !NHD2_0001 -> NHD2_gas (photodesorption)
    kall(2279) = Ffuva * 1.00000000d-03

    !ND3_0001 -> ND3_gas (photodesorption)
    kall(2280) = Ffuva * 1.00000000d-03

    !C10_0001 -> C10_gas (photodesorption)
    kall(2281) = Ffuva * 1.00000000d-03

    !C10H_0001 -> C10H_gas (photodesorption)
    kall(2282) = Ffuva * 1.00000000d-03

    !C10H2_0001 -> C10H2_gas (photodesorption)
    kall(2283) = Ffuva * 1.00000000d-03

    !C10N_0001 -> C10N_gas (photodesorption)
    kall(2284) = Ffuva * 1.00000000d-03

    !C11_0001 -> C11_gas (photodesorption)
    kall(2285) = Ffuva * 1.00000000d-03

    !C2H2_0001 -> C2H2_gas (photodesorption)
    kall(2286) = Ffuva * 1.00000000d-03

    !C2HD_0001 -> C2HD_gas (photodesorption)
    kall(2287) = Ffuva * 1.00000000d-03

    !C2D2_0001 -> C2D2_gas (photodesorption)
    kall(2288) = Ffuva * 1.00000000d-03

    !C2H3_0001 -> C2H3_gas (photodesorption)
    kall(2289) = Ffuva * 1.00000000d-03

    !C2H2D_0001 -> C2H2D_gas (photodesorption)
    kall(2290) = Ffuva * 1.00000000d-03

    !C2HD2_0001 -> C2HD2_gas (photodesorption)
    kall(2291) = Ffuva * 1.00000000d-03

    !C2D3_0001 -> C2D3_gas (photodesorption)
    kall(2292) = Ffuva * 1.00000000d-03

    !C2H4_0001 -> C2H4_gas (photodesorption)
    kall(2293) = Ffuva * 1.00000000d-03

    !C2H3D_0001 -> C2H3D_gas (photodesorption)
    kall(2294) = Ffuva * 1.00000000d-03

    !C2H2D2_0001 -> C2H2D2_gas (photodesorption)
    kall(2295) = Ffuva * 1.00000000d-03

    !C2HD3_0001 -> C2HD3_gas (photodesorption)
    kall(2296) = Ffuva * 1.00000000d-03

    !C2D4_0001 -> C2D4_gas (photodesorption)
    kall(2297) = Ffuva * 1.00000000d-03

    !C2H5_0001 -> C2H5_gas (photodesorption)
    kall(2298) = Ffuva * 1.00000000d-03

    !C2H6_0001 -> C2H6_gas (photodesorption)
    kall(2299) = Ffuva * 1.00000000d-03

    !C3_0001 -> C3_gas (photodesorption)
    kall(2300) = Ffuva * 1.00000000d-03

    !C3N_0001 -> C3N_gas (photodesorption)
    kall(2301) = Ffuva * 1.00000000d-03

    !C3O_0001 -> C3O_gas (photodesorption)
    kall(2302) = Ffuva * 1.00000000d-03

    !C3S_0001 -> C3S_gas (photodesorption)
    kall(2303) = Ffuva * 1.00000000d-03

    !C4_0001 -> C4_gas (photodesorption)
    kall(2304) = Ffuva * 1.00000000d-03

    !C4H_0001 -> C4H_gas (photodesorption)
    kall(2305) = Ffuva * 1.00000000d-03

    !C4D_0001 -> C4D_gas (photodesorption)
    kall(2306) = Ffuva * 1.00000000d-03

    !C4H2_0001 -> C4H2_gas (photodesorption)
    kall(2307) = Ffuva * 1.00000000d-03

    !C4HD_0001 -> C4HD_gas (photodesorption)
    kall(2308) = Ffuva * 1.00000000d-03

    !C4D2_0001 -> C4D2_gas (photodesorption)
    kall(2309) = Ffuva * 1.00000000d-03

    !C4H3_0001 -> C4H3_gas (photodesorption)
    kall(2310) = Ffuva * 1.00000000d-03

    !C4N_0001 -> C4N_gas (photodesorption)
    kall(2311) = Ffuva * 1.00000000d-03

    !C4S_0001 -> C4S_gas (photodesorption)
    kall(2312) = Ffuva * 1.00000000d-03

    !C5_0001 -> C5_gas (photodesorption)
    kall(2313) = Ffuva * 1.00000000d-03

    !C5H_0001 -> C5H_gas (photodesorption)
    kall(2314) = Ffuva * 1.00000000d-03

    !C5D_0001 -> C5D_gas (photodesorption)
    kall(2315) = Ffuva * 1.00000000d-03

    !C5H2_0001 -> C5H2_gas (photodesorption)
    kall(2316) = Ffuva * 1.00000000d-03

    !C5H3_0001 -> C5H3_gas (photodesorption)
    kall(2317) = Ffuva * 1.00000000d-03

    !C5H4_0001 -> C5H4_gas (photodesorption)
    kall(2318) = Ffuva * 1.00000000d-03

    !C5N_0001 -> C5N_gas (photodesorption)
    kall(2319) = Ffuva * 1.00000000d-03

    !C5O_0001 -> C5O_gas (photodesorption)
    kall(2320) = Ffuva * 1.00000000d-03

    !C6_0001 -> C6_gas (photodesorption)
    kall(2321) = Ffuva * 1.00000000d-03

    !C6H_0001 -> C6H_gas (photodesorption)
    kall(2322) = Ffuva * 1.00000000d-03

    !C6H2_0001 -> C6H2_gas (photodesorption)
    kall(2323) = Ffuva * 1.00000000d-03

    !C6H3_0001 -> C6H3_gas (photodesorption)
    kall(2324) = Ffuva * 1.00000000d-03

    !C6H4_0001 -> C6H4_gas (photodesorption)
    kall(2325) = Ffuva * 1.00000000d-03

    !C6H6_0001 -> C6H6_gas (photodesorption)
    kall(2326) = Ffuva * 1.00000000d-03

    !C6N_0001 -> C6N_gas (photodesorption)
    kall(2327) = Ffuva * 1.00000000d-03

    !C7_0001 -> C7_gas (photodesorption)
    kall(2328) = Ffuva * 1.00000000d-03

    !C7H_0001 -> C7H_gas (photodesorption)
    kall(2329) = Ffuva * 1.00000000d-03

    !C7H2_0001 -> C7H2_gas (photodesorption)
    kall(2330) = Ffuva * 1.00000000d-03

    !C7H3_0001 -> C7H3_gas (photodesorption)
    kall(2331) = Ffuva * 1.00000000d-03

    !C7H4_0001 -> C7H4_gas (photodesorption)
    kall(2332) = Ffuva * 1.00000000d-03

    !C7N_0001 -> C7N_gas (photodesorption)
    kall(2333) = Ffuva * 1.00000000d-03

    !C7O_0001 -> C7O_gas (photodesorption)
    kall(2334) = Ffuva * 1.00000000d-03

    !C8_0001 -> C8_gas (photodesorption)
    kall(2335) = Ffuva * 1.00000000d-03

    !C8H_0001 -> C8H_gas (photodesorption)
    kall(2336) = Ffuva * 1.00000000d-03

    !C8H2_0001 -> C8H2_gas (photodesorption)
    kall(2337) = Ffuva * 1.00000000d-03

    !C8H3_0001 -> C8H3_gas (photodesorption)
    kall(2338) = Ffuva * 1.00000000d-03

    !C8H4_0001 -> C8H4_gas (photodesorption)
    kall(2339) = Ffuva * 1.00000000d-03

    !C8N_0001 -> C8N_gas (photodesorption)
    kall(2340) = Ffuva * 1.00000000d-03

    !C9_0001 -> C9_gas (photodesorption)
    kall(2341) = Ffuva * 1.00000000d-03

    !C9H_0001 -> C9H_gas (photodesorption)
    kall(2342) = Ffuva * 1.00000000d-03

    !C9H2_0001 -> C9H2_gas (photodesorption)
    kall(2343) = Ffuva * 1.00000000d-03

    !C9H3_0001 -> C9H3_gas (photodesorption)
    kall(2344) = Ffuva * 1.00000000d-03

    !C9H4_0001 -> C9H4_gas (photodesorption)
    kall(2345) = Ffuva * 1.00000000d-03

    !C9N_0001 -> C9N_gas (photodesorption)
    kall(2346) = Ffuva * 1.00000000d-03

    !C9O_0001 -> C9O_gas (photodesorption)
    kall(2347) = Ffuva * 1.00000000d-03

    !CCH_0001 -> CCH_gas (photodesorption)
    kall(2348) = Ffuva * 1.00000000d-03

    !CCD_0001 -> CCD_gas (photodesorption)
    kall(2349) = Ffuva * 1.00000000d-03

    !CCN_0001 -> CCN_gas (photodesorption)
    kall(2350) = Ffuva * 1.00000000d-03

    !CCO_0001 -> CCO_gas (photodesorption)
    kall(2351) = Ffuva * 1.00000000d-03

    !CCS_0001 -> CCS_gas (photodesorption)
    kall(2352) = Ffuva * 1.00000000d-03

    !CD3CN_0001 -> CD3CN_gas (photodesorption)
    kall(2353) = Ffuva * 1.00000000d-03

    !CH2CCH_0001 -> CH2CCH_gas (photodesorption)
    kall(2354) = Ffuva * 1.00000000d-03

    !CH2CCD_0001 -> CH2CCD_gas (photodesorption)
    kall(2355) = Ffuva * 1.00000000d-03

    !CHDCCH_0001 -> CHDCCH_gas (photodesorption)
    kall(2356) = Ffuva * 1.00000000d-03

    !CHDCCD_0001 -> CHDCCD_gas (photodesorption)
    kall(2357) = Ffuva * 1.00000000d-03

    !CD2CCH_0001 -> CD2CCH_gas (photodesorption)
    kall(2358) = Ffuva * 1.00000000d-03

    !CD2CCD_0001 -> CD2CCD_gas (photodesorption)
    kall(2359) = Ffuva * 1.00000000d-03

    !CH2CHC2H_0001 -> CH2CHC2H_gas (photodesorption)
    kall(2360) = Ffuva * 1.00000000d-03

    !CH2CHCHCH2_0001 -> CH2CHCHCH2_gas (photodesorption)
    kall(2361) = Ffuva * 1.00000000d-03

    !CH2CHCN_0001 -> CH2CHCN_gas (photodesorption)
    kall(2362) = Ffuva * 1.00000000d-03

    !CH2NH_0001 -> CH2NH_gas (photodesorption)
    kall(2363) = Ffuva * 1.00000000d-03

    !CHDNH_0001 -> CHDNH_gas (photodesorption)
    kall(2364) = Ffuva * 1.00000000d-03

    !CHDND_0001 -> CHDND_gas (photodesorption)
    kall(2365) = Ffuva * 1.00000000d-03

    !CH2ND_0001 -> CH2ND_gas (photodesorption)
    kall(2366) = Ffuva * 1.00000000d-03

    !CD2NH_0001 -> CD2NH_gas (photodesorption)
    kall(2367) = Ffuva * 1.00000000d-03

    !CD2ND_0001 -> CD2ND_gas (photodesorption)
    kall(2368) = Ffuva * 1.00000000d-03

    !CH2NH2_0001 -> CH2NH2_gas (photodesorption)
    kall(2369) = Ffuva * 1.00000000d-03

    !CH2NHD_0001 -> CH2NHD_gas (photodesorption)
    kall(2370) = Ffuva * 1.00000000d-03

    !CHDNH2_0001 -> CHDNH2_gas (photodesorption)
    kall(2371) = Ffuva * 1.00000000d-03

    !CHDNHD_0001 -> CHDNHD_gas (photodesorption)
    kall(2372) = Ffuva * 1.00000000d-03

    !CHDND2_0001 -> CHDND2_gas (photodesorption)
    kall(2373) = Ffuva * 1.00000000d-03

    !CH2ND2_0001 -> CH2ND2_gas (photodesorption)
    kall(2374) = Ffuva * 1.00000000d-03

    !CD2NH2_0001 -> CD2NH2_gas (photodesorption)
    kall(2375) = Ffuva * 1.00000000d-03

    !CD2NHD_0001 -> CD2NHD_gas (photodesorption)
    kall(2376) = Ffuva * 1.00000000d-03

    !CD2ND2_0001 -> CD2ND2_gas (photodesorption)
    kall(2377) = Ffuva * 1.00000000d-03

    !CH3C3N_0001 -> CH3C3N_gas (photodesorption)
    kall(2378) = Ffuva * 1.00000000d-03

    !CH3C4H_0001 -> CH3C4H_gas (photodesorption)
    kall(2379) = Ffuva * 1.00000000d-03

    !CH3C5N_0001 -> CH3C5N_gas (photodesorption)
    kall(2380) = Ffuva * 1.00000000d-03

    !CH3C6H_0001 -> CH3C6H_gas (photodesorption)
    kall(2381) = Ffuva * 1.00000000d-03

    !CH3C7N_0001 -> CH3C7N_gas (photodesorption)
    kall(2382) = Ffuva * 1.00000000d-03

    !CH3CCH_0001 -> CH3CCH_gas (photodesorption)
    kall(2383) = Ffuva * 1.00000000d-03

    !CH3CH2OH_0001 -> CH3CH2OH_gas (photodesorption)
    kall(2384) = Ffuva * 1.00000000d-03

    !CH3CHCH2_0001 -> CH3CHCH2_gas (photodesorption)
    kall(2385) = Ffuva * 1.00000000d-03

    !CH3CHO_0001 -> CH3CHO_gas (photodesorption)
    kall(2386) = Ffuva * 1.00000000d-03

    !CH3CN_0001 -> CH3CN_gas (photodesorption)
    kall(2387) = Ffuva * 1.00000000d-03

    !CH2DCN_0001 -> CH2DCN_gas (photodesorption)
    kall(2388) = Ffuva * 1.00000000d-03

    !CHD2CN_0001 -> CHD2CN_gas (photodesorption)
    kall(2389) = Ffuva * 1.00000000d-03

    !CH3COCH3_0001 -> CH3COCH3_gas (photodesorption)
    kall(2390) = Ffuva * 1.00000000d-03

    !CH3NH2_0001 -> CH3NH2_gas (photodesorption)
    kall(2391) = Ffuva * 1.00000000d-03

    !CH3OCH2_0001 -> CH3OCH2_gas (photodesorption)
    kall(2392) = Ffuva * 1.00000000d-03

    !CH3OCH3_0001 -> CH3OCH3_gas (photodesorption)
    kall(2393) = Ffuva * 1.00000000d-03

    !CN_0001 -> CN_gas (photodesorption)
    kall(2394) = Ffuva * 1.00000000d-03

    !CS_0001 -> CS_gas (photodesorption)
    kall(2395) = Ffuva * 1.00000000d-03

    !H2CCN_0001 -> H2CCN_gas (photodesorption)
    kall(2396) = Ffuva * 1.00000000d-03

    !HDCCN_0001 -> HDCCN_gas (photodesorption)
    kall(2397) = Ffuva * 1.00000000d-03

    !D2CCN_0001 -> D2CCN_gas (photodesorption)
    kall(2398) = Ffuva * 1.00000000d-03

    !H2CCO_0001 -> H2CCO_gas (photodesorption)
    kall(2399) = Ffuva * 1.00000000d-03

    !HDCCO_0001 -> HDCCO_gas (photodesorption)
    kall(2400) = Ffuva * 1.00000000d-03

    !D2CCO_0001 -> D2CCO_gas (photodesorption)
    kall(2401) = Ffuva * 1.00000000d-03

    !H2CN_0001 -> H2CN_gas (photodesorption)
    kall(2402) = Ffuva * 1.00000000d-03

    !HDCN_0001 -> HDCN_gas (photodesorption)
    kall(2403) = Ffuva * 1.00000000d-03

    !D2CN_0001 -> D2CN_gas (photodesorption)
    kall(2404) = Ffuva * 1.00000000d-03

    !H2CS_0001 -> H2CS_gas (photodesorption)
    kall(2405) = Ffuva * 1.00000000d-03

    !HDCS_0001 -> HDCS_gas (photodesorption)
    kall(2406) = Ffuva * 1.00000000d-03

    !D2CS_0001 -> D2CS_gas (photodesorption)
    kall(2407) = Ffuva * 1.00000000d-03

    !H2S_0001 -> H2S_gas (photodesorption)
    kall(2408) = Ffuva * 1.00000000d-03

    !HDS_0001 -> HDS_gas (photodesorption)
    kall(2409) = Ffuva * 1.00000000d-03

    !D2S_0001 -> D2S_gas (photodesorption)
    kall(2410) = Ffuva * 1.00000000d-03

    !HC2O_0001 -> HC2O_gas (photodesorption)
    kall(2411) = Ffuva * 1.00000000d-03

    !DC2O_0001 -> DC2O_gas (photodesorption)
    kall(2412) = Ffuva * 1.00000000d-03

    !HC3N_0001 -> HC3N_gas (photodesorption)
    kall(2413) = Ffuva * 1.00000000d-03

    !DC3N_0001 -> DC3N_gas (photodesorption)
    kall(2414) = Ffuva * 1.00000000d-03

    !HC4N_0001 -> HC4N_gas (photodesorption)
    kall(2415) = Ffuva * 1.00000000d-03

    !DC4N_0001 -> DC4N_gas (photodesorption)
    kall(2416) = Ffuva * 1.00000000d-03

    !HC5N_0001 -> HC5N_gas (photodesorption)
    kall(2417) = Ffuva * 1.00000000d-03

    !HC6N_0001 -> HC6N_gas (photodesorption)
    kall(2418) = Ffuva * 1.00000000d-03

    !HC7N_0001 -> HC7N_gas (photodesorption)
    kall(2419) = Ffuva * 1.00000000d-03

    !HC8N_0001 -> HC8N_gas (photodesorption)
    kall(2420) = Ffuva * 1.00000000d-03

    !HC9N_0001 -> HC9N_gas (photodesorption)
    kall(2421) = Ffuva * 1.00000000d-03

    !HCCNC_0001 -> HCCNC_gas (photodesorption)
    kall(2422) = Ffuva * 1.00000000d-03

    !DCCNC_0001 -> DCCNC_gas (photodesorption)
    kall(2423) = Ffuva * 1.00000000d-03

    !HCN_0001 -> HCN_gas (photodesorption)
    kall(2424) = Ffuva * 1.00000000d-03

    !DCN_0001 -> DCN_gas (photodesorption)
    kall(2425) = Ffuva * 1.00000000d-03

    !HCNCC_0001 -> HCNCC_gas (photodesorption)
    kall(2426) = Ffuva * 1.00000000d-03

    !DCNCC_0001 -> DCNCC_gas (photodesorption)
    kall(2427) = Ffuva * 1.00000000d-03

    !HCOOCH3_0001 -> HCOOCH3_gas (photodesorption)
    kall(2428) = Ffuva * 1.00000000d-03

    !HCS_0001 -> HCS_gas (photodesorption)
    kall(2429) = Ffuva * 1.00000000d-03

    !DCS_0001 -> DCS_gas (photodesorption)
    kall(2430) = Ffuva * 1.00000000d-03

    !HNC_0001 -> HNC_gas (photodesorption)
    kall(2431) = Ffuva * 1.00000000d-03

    !DNC_0001 -> DNC_gas (photodesorption)
    kall(2432) = Ffuva * 1.00000000d-03

    !HNCCC_0001 -> HNCCC_gas (photodesorption)
    kall(2433) = Ffuva * 1.00000000d-03

    !DNCCC_0001 -> DNCCC_gas (photodesorption)
    kall(2434) = Ffuva * 1.00000000d-03

    !HNCO_0001 -> HNCO_gas (photodesorption)
    kall(2435) = Ffuva * 1.00000000d-03

    !DNCO_0001 -> DNCO_gas (photodesorption)
    kall(2436) = Ffuva * 1.00000000d-03

    !HNO_0001 -> HNO_gas (photodesorption)
    kall(2437) = Ffuva * 1.00000000d-03

    !DNO_0001 -> DNO_gas (photodesorption)
    kall(2438) = Ffuva * 1.00000000d-03

    !HS_0001 -> HS_gas (photodesorption)
    kall(2439) = Ffuva * 1.00000000d-03

    !DS_0001 -> DS_gas (photodesorption)
    kall(2440) = Ffuva * 1.00000000d-03

    !N2O_0001 -> N2O_gas (photodesorption)
    kall(2441) = Ffuva * 1.00000000d-03

    !NC4N_0001 -> NC4N_gas (photodesorption)
    kall(2442) = Ffuva * 1.00000000d-03

    !NC6N_0001 -> NC6N_gas (photodesorption)
    kall(2443) = Ffuva * 1.00000000d-03

    !NC8N_0001 -> NC8N_gas (photodesorption)
    kall(2444) = Ffuva * 1.00000000d-03

    !NH2CHO_0001 -> NH2CHO_gas (photodesorption)
    kall(2445) = Ffuva * 1.00000000d-03

    !NH2CDO_0001 -> NH2CDO_gas (photodesorption)
    kall(2446) = Ffuva * 1.00000000d-03

    !NHDCHO_0001 -> NHDCHO_gas (photodesorption)
    kall(2447) = Ffuva * 1.00000000d-03

    !NHDCDO_0001 -> NHDCDO_gas (photodesorption)
    kall(2448) = Ffuva * 1.00000000d-03

    !ND2CHO_0001 -> ND2CHO_gas (photodesorption)
    kall(2449) = Ffuva * 1.00000000d-03

    !ND2CDO_0001 -> ND2CDO_gas (photodesorption)
    kall(2450) = Ffuva * 1.00000000d-03

    !NH2CN_0001 -> NH2CN_gas (photodesorption)
    kall(2451) = Ffuva * 1.00000000d-03

    !NHDCN_0001 -> NHDCN_gas (photodesorption)
    kall(2452) = Ffuva * 1.00000000d-03

    !ND2CN_0001 -> ND2CN_gas (photodesorption)
    kall(2453) = Ffuva * 1.00000000d-03

    !HSS_0001 -> HSS_gas (photodesorption)
    kall(2454) = Ffuva * 1.00000000d-03

    !DSS_0001 -> DSS_gas (photodesorption)
    kall(2455) = Ffuva * 1.00000000d-03

    !HSSH_0001 -> HSSH_gas (photodesorption)
    kall(2456) = Ffuva * 1.00000000d-03

    !HSSD_0001 -> HSSD_gas (photodesorption)
    kall(2457) = Ffuva * 1.00000000d-03

    !DSSH_0001 -> DSSH_gas (photodesorption)
    kall(2458) = Ffuva * 1.00000000d-03

    !DSSD_0001 -> DSSD_gas (photodesorption)
    kall(2459) = Ffuva * 1.00000000d-03

    !NO_0001 -> NO_gas (photodesorption)
    kall(2460) = Ffuva * 1.00000000d-03

    !NO2_0001 -> NO2_gas (photodesorption)
    kall(2461) = Ffuva * 1.00000000d-03

    !NS_0001 -> NS_gas (photodesorption)
    kall(2462) = Ffuva * 1.00000000d-03

    !OCN_0001 -> OCN_gas (photodesorption)
    kall(2463) = Ffuva * 1.00000000d-03

    !OCS_0001 -> OCS_gas (photodesorption)
    kall(2464) = Ffuva * 1.00000000d-03

    !S2_0001 -> S2_gas (photodesorption)
    kall(2465) = Ffuva * 1.00000000d-03

    !SO_0001 -> SO_gas (photodesorption)
    kall(2466) = Ffuva * 1.00000000d-03

    !SO2_0001 -> SO2_gas (photodesorption)
    kall(2467) = Ffuva * 1.00000000d-03

    !CH3OCHO_0001 -> CH3OCHO_gas (photodesorption)
    kall(2468) = Ffuva * 1.00000000d-03

    !Si_0001 -> Si_gas (photodesorption)
    kall(2469) = Ffuva * 1.00000000d-03

    !SiS_0001 -> SiS_gas (photodesorption)
    kall(2470) = Ffuva * 1.00000000d-03

    !SiN_0001 -> SiN_gas (photodesorption)
    kall(2471) = Ffuva * 1.00000000d-03

    !SiC_0001 -> SiC_gas (photodesorption)
    kall(2472) = Ffuva * 1.00000000d-03

    !SiH_0001 -> SiH_gas (photodesorption)
    kall(2473) = Ffuva * 1.00000000d-03

    !SiH2_0001 -> SiH2_gas (photodesorption)
    kall(2474) = Ffuva * 1.00000000d-03

    !SiH3_0001 -> SiH3_gas (photodesorption)
    kall(2475) = Ffuva * 1.00000000d-03

    !SiH4_0001 -> SiH4_gas (photodesorption)
    kall(2476) = Ffuva * 1.00000000d-03

    !SiC2CH3_0001 -> SiC2CH3_gas (photodesorption)
    kall(2477) = Ffuva * 1.00000000d-03

    !SiC3H_0001 -> SiC3H_gas (photodesorption)
    kall(2478) = Ffuva * 1.00000000d-03

    !SiC3H5_0001 -> SiC3H5_gas (photodesorption)
    kall(2479) = Ffuva * 1.00000000d-03

    !SiC4_0001 -> SiC4_gas (photodesorption)
    kall(2480) = Ffuva * 1.00000000d-03

    !SiC4H_0001 -> SiC4H_gas (photodesorption)
    kall(2481) = Ffuva * 1.00000000d-03

    !SiC6H_0001 -> SiC6H_gas (photodesorption)
    kall(2482) = Ffuva * 1.00000000d-03

    !SiC8H_0001 -> SiC8H_gas (photodesorption)
    kall(2483) = Ffuva * 1.00000000d-03

    !c_HCCHSi_0001 -> c_HCCHSi_gas (photodesorption)
    kall(2484) = Ffuva * 1.00000000d-03

    !c_SiC2_0001 -> c_SiC2_gas (photodesorption)
    kall(2485) = Ffuva * 1.00000000d-03

    !l_C3H_0001 -> l_C3H_gas (photodesorption)
    kall(2486) = Ffuva * 1.00000000d-03

    !l_C3D_0001 -> l_C3D_gas (photodesorption)
    kall(2487) = Ffuva * 1.00000000d-03

    !c_C3H_0001 -> c_C3H_gas (photodesorption)
    kall(2488) = Ffuva * 1.00000000d-03

    !c_C3D_0001 -> c_C3D_gas (photodesorption)
    kall(2489) = Ffuva * 1.00000000d-03

    !l_SiC3_0001 -> l_SiC3_gas (photodesorption)
    kall(2490) = Ffuva * 1.00000000d-03

    !l_C3H2_0001 -> l_C3H2_gas (photodesorption)
    kall(2491) = Ffuva * 1.00000000d-03

    !l_C3HD_0001 -> l_C3HD_gas (photodesorption)
    kall(2492) = Ffuva * 1.00000000d-03

    !l_C3D2_0001 -> l_C3D2_gas (photodesorption)
    kall(2493) = Ffuva * 1.00000000d-03

    !c_C3H2_0001 -> c_C3H2_gas (photodesorption)
    kall(2494) = Ffuva * 1.00000000d-03

    !c_C3HD_0001 -> c_C3HD_gas (photodesorption)
    kall(2495) = Ffuva * 1.00000000d-03

    !c_C3D2_0001 -> c_C3D2_gas (photodesorption)
    kall(2496) = Ffuva * 1.00000000d-03

    !Mg_0001 -> Mg_gas (photodesorption)
    kall(2497) = Ffuva * 1.00000000d-03

    !MgH_0001 -> MgH_gas (photodesorption)
    kall(2498) = Ffuva * 1.00000000d-03

    !MgH2_0001 -> MgH2_gas (photodesorption)
    kall(2499) = Ffuva * 1.00000000d-03

    !Na_0001 -> Na_gas (photodesorption)
    kall(2500) = Ffuva * 1.00000000d-03

    !NaH_0001 -> NaH_gas (photodesorption)
    kall(2501) = Ffuva * 1.00000000d-03

    !F_0001 -> F_gas (photodesorption)
    kall(2502) = Ffuva * 1.00000000d-03

    !HF_0001 -> HF_gas (photodesorption)
    kall(2503) = Ffuva * 1.00000000d-03

    !DF_0001 -> DF_gas (photodesorption)
    kall(2504) = Ffuva * 1.00000000d-03

    !MgD_0001 -> MgD_gas (photodesorption)
    kall(2505) = Ffuva * 1.00000000d-03

    !MgHD_0001 -> MgHD_gas (photodesorption)
    kall(2506) = Ffuva * 1.00000000d-03

    !MgD2_0001 -> MgD2_gas (photodesorption)
    kall(2507) = Ffuva * 1.00000000d-03

    !NaD_0001 -> NaD_gas (photodesorption)
    kall(2508) = Ffuva * 1.00000000d-03

    !SiD_0001 -> SiD_gas (photodesorption)
    kall(2509) = Ffuva * 1.00000000d-03

    !Cl_0001 -> Cl_gas (photodesorption)
    kall(2510) = Ffuva * 1.00000000d-03

    !HCl_0001 -> HCl_gas (photodesorption)
    kall(2511) = Ffuva * 1.00000000d-03

    !DCl_0001 -> DCl_gas (photodesorption)
    kall(2512) = Ffuva * 1.00000000d-03

    !H2Cl_0001 -> H2Cl_gas (photodesorption)
    kall(2513) = Ffuva * 1.00000000d-03

    !HDCl_0001 -> HDCl_gas (photodesorption)
    kall(2514) = Ffuva * 1.00000000d-03

    !D2Cl_0001 -> D2Cl_gas (photodesorption)
    kall(2515) = Ffuva * 1.00000000d-03

    !CCl_0001 -> CCl_gas (photodesorption)
    kall(2516) = Ffuva * 1.00000000d-03

    !ClO_0001 -> ClO_gas (photodesorption)
    kall(2517) = Ffuva * 1.00000000d-03

    !C3H4_0001 -> C3H4_gas (photodesorption)
    kall(2518) = Ffuva * 1.00000000d-03

    !FeD_0001 -> FeD_gas (photodesorption)
    kall(2519) = Ffuva * 1.00000000d-03

    !P_0001 -> P_gas (photodesorption)
    kall(2520) = Ffuva * 1.00000000d-03

    !PO_0001 -> PO_gas (photodesorption)
    kall(2521) = Ffuva * 1.00000000d-03

    !PH_0001 -> PH_gas (photodesorption)
    kall(2522) = Ffuva * 1.00000000d-03

    !PD_0001 -> PD_gas (photodesorption)
    kall(2523) = Ffuva * 1.00000000d-03

    !PH2_0001 -> PH2_gas (photodesorption)
    kall(2524) = Ffuva * 1.00000000d-03

    !PHD_0001 -> PHD_gas (photodesorption)
    kall(2525) = Ffuva * 1.00000000d-03

    !PD2_0001 -> PD2_gas (photodesorption)
    kall(2526) = Ffuva * 1.00000000d-03

    !PN_0001 -> PN_gas (photodesorption)
    kall(2527) = Ffuva * 1.00000000d-03

    !CP_0001 -> CP_gas (photodesorption)
    kall(2528) = Ffuva * 1.00000000d-03

    !CCP_0001 -> CCP_gas (photodesorption)
    kall(2529) = Ffuva * 1.00000000d-03

    !C3P_0001 -> C3P_gas (photodesorption)
    kall(2530) = Ffuva * 1.00000000d-03

    !C4P_0001 -> C4P_gas (photodesorption)
    kall(2531) = Ffuva * 1.00000000d-03

    !CH2PH_0001 -> CH2PH_gas (photodesorption)
    kall(2532) = Ffuva * 1.00000000d-03

    !CH2PD_0001 -> CH2PD_gas (photodesorption)
    kall(2533) = Ffuva * 1.00000000d-03

    !CHDPD_0001 -> CHDPD_gas (photodesorption)
    kall(2534) = Ffuva * 1.00000000d-03

    !HCP_0001 -> HCP_gas (photodesorption)
    kall(2535) = Ffuva * 1.00000000d-03

    !DCP_0001 -> DCP_gas (photodesorption)
    kall(2536) = Ffuva * 1.00000000d-03

    !HCCP_0001 -> HCCP_gas (photodesorption)
    kall(2537) = Ffuva * 1.00000000d-03

    !DCCP_0001 -> DCCP_gas (photodesorption)
    kall(2538) = Ffuva * 1.00000000d-03

    !HPO_0001 -> HPO_gas (photodesorption)
    kall(2539) = Ffuva * 1.00000000d-03

    !DPO_0001 -> DPO_gas (photodesorption)
    kall(2540) = Ffuva * 1.00000000d-03

    !H_0001 -> H_gas (CRdesorption)
    kall(2541) = variable_crflux*6.80579278d+07

    !D_0001 -> D_gas (CRdesorption)
    kall(2542) = variable_crflux*3.63922778d+07

    !p_D2_0001 -> p_D2_gas (CRdesorption)
    kall(2543) = variable_crflux*6.59447650d+07

    !o_D2_0001 -> o_D2_gas (CRdesorption)
    kall(2544) = variable_crflux*6.59447650d+07

    !HD_0001 -> HD_gas (CRdesorption)
    kall(2545) = variable_crflux*7.61464556d+07

    !He_0001 -> He_gas (CRdesorption)
    kall(2546) = variable_crflux*4.61366946d+09

    !O_0001 -> O_gas (CRdesorption)
    kall(2547) = variable_crflux*4.55791421d+00

    !O2_0001 -> O2_gas (CRdesorption)
    kall(2548) = variable_crflux*8.46183464d+02

    !O3_0001 -> O3_gas (CRdesorption)
    kall(2549) = variable_crflux*2.38315195d-03

    !OH_0001 -> OH_gas (CRdesorption)
    kall(2550) = variable_crflux*1.82936324d-18

    !OD_0001 -> OD_gas (CRdesorption)
    kall(2551) = variable_crflux*1.77782151d-18

    !H2O_0001 -> H2O_gas (CRdesorption)
    kall(2552) = variable_crflux*1.22573439d-24

    !HDO_0001 -> HDO_gas (CRdesorption)
    kall(2553) = variable_crflux*1.19304225d-24

    !D2O_0001 -> D2O_gas (CRdesorption)
    kall(2554) = variable_crflux*1.16283374d-24

    !O2H_0001 -> O2H_gas (CRdesorption)
    kall(2555) = variable_crflux*4.51534327d-21

    !O2D_0001 -> O2D_gas (CRdesorption)
    kall(2556) = variable_crflux*4.44844559d-21

    !HOOH_0001 -> HOOH_gas (CRdesorption)
    kall(2557) = variable_crflux*3.04503313d-27

    !HOOD_0001 -> HOOD_gas (CRdesorption)
    kall(2558) = variable_crflux*3.00121742d-27

    !DOOD_0001 -> DOOD_gas (CRdesorption)
    kall(2559) = variable_crflux*2.95924028d-27

    !Fe_0001 -> Fe_gas (CRdesorption)
    kall(2560) = variable_crflux*2.91983882d-16

    !FeH_0001 -> FeH_gas (CRdesorption)
    kall(2561) = variable_crflux*4.91727161d-19

    !N_0001 -> N_gas (CRdesorption)
    kall(2562) = variable_crflux*9.42042367d+05

    !S_0001 -> S_gas (CRdesorption)
    kall(2563) = variable_crflux*2.56726646d-06

    !C_0001 -> C_gas (CRdesorption)
    kall(2564) = variable_crflux*1.00887578d-51

    !C2_0001 -> C2_gas (CRdesorption)
    kall(2565) = variable_crflux*7.13382904d-52

    !CO_0001 -> CO_gas (CRdesorption)
    kall(2566) = variable_crflux*2.25642484d+02

    !HCO_0001 -> HCO_gas (CRdesorption)
    kall(2567) = variable_crflux*4.51135768d-05

    !DCO_0001 -> DCO_gas (CRdesorption)
    kall(2568) = variable_crflux*4.43553114d-05

    !H2CO_0001 -> H2CO_gas (CRdesorption)
    kall(2569) = variable_crflux*5.68344696d-18

    !HDCO_0001 -> HDCO_gas (CRdesorption)
    kall(2570) = variable_crflux*5.59102703d-18

    !D2CO_0001 -> D2CO_gas (CRdesorption)
    kall(2571) = variable_crflux*5.50297386d-18

    !CH2OH_0001 -> CH2OH_gas (CRdesorption)
    kall(2572) = variable_crflux*2.30691908d-17

    !CD2OD_0001 -> CD2OD_gas (CRdesorption)
    kall(2573) = variable_crflux*2.20279332d-17

    !CH2OD_0001 -> CH2OD_gas (CRdesorption)
    kall(2574) = variable_crflux*2.27058738d-17

    !CHDOH_0001 -> CHDOH_gas (CRdesorption)
    kall(2575) = variable_crflux*2.27058738d-17

    !CHDOD_0001 -> CHDOD_gas (CRdesorption)
    kall(2576) = variable_crflux*2.23591989d-17

    !CD2OH_0001 -> CD2OH_gas (CRdesorption)
    kall(2577) = variable_crflux*2.23591989d-17

    !CH3O_0001 -> CH3O_gas (CRdesorption)
    kall(2578) = variable_crflux*2.30691908d-17

    !CHD2O_0001 -> CHD2O_gas (CRdesorption)
    kall(2579) = variable_crflux*2.23591989d-17

    !CH2DO_0001 -> CH2DO_gas (CRdesorption)
    kall(2580) = variable_crflux*2.27058738d-17

    !CD3O_0001 -> CD3O_gas (CRdesorption)
    kall(2581) = variable_crflux*2.20279332d-17

    !CH3OH_0001 -> CH3OH_gas (CRdesorption)
    kall(2582) = variable_crflux*4.58535276d-21

    !CH3OD_0001 -> CH3OD_gas (CRdesorption)
    kall(2583) = variable_crflux*4.51534327d-21

    !CHD2OH_0001 -> CHD2OH_gas (CRdesorption)
    kall(2584) = variable_crflux*4.44844559d-21

    !CHD2OD_0001 -> CHD2OD_gas (CRdesorption)
    kall(2585) = variable_crflux*4.38443585d-21

    !CH2DOH_0001 -> CH2DOH_gas (CRdesorption)
    kall(2586) = variable_crflux*4.51534327d-21

    !CH2DOD_0001 -> CH2DOD_gas (CRdesorption)
    kall(2587) = variable_crflux*4.44844559d-21

    !CD3OD_0001 -> CD3OD_gas (CRdesorption)
    kall(2588) = variable_crflux*4.32311204d-21

    !CD3OH_0001 -> CD3OH_gas (CRdesorption)
    kall(2589) = variable_crflux*4.38443585d-21

    !CH_0001 -> CH_gas (CRdesorption)
    kall(2590) = variable_crflux*5.92522687d+04

    !CD_0001 -> CD_gas (CRdesorption)
    kall(2591) = variable_crflux*5.70969148d+04

    !CH2_0001 -> CH2_gas (CRdesorption)
    kall(2592) = variable_crflux*7.93610716d+01

    !CHD_0001 -> CHD_gas (CRdesorption)
    kall(2593) = variable_crflux*7.66700792d+01

    !CD2_0001 -> CD2_gas (CRdesorption)
    kall(2594) = variable_crflux*7.42354849d+01

    !CH3_0001 -> CH3_gas (CRdesorption)
    kall(2595) = variable_crflux*4.70739355d+00

    !CH2D_0001 -> CH2D_gas (CRdesorption)
    kall(2596) = variable_crflux*4.55791421d+00

    !CHD2_0001 -> CHD2_gas (CRdesorption)
    kall(2597) = variable_crflux*5.89225860d-05

    !CD3_0001 -> CD3_gas (CRdesorption)
    kall(2598) = variable_crflux*4.29724273d+00

    !CH4_0001 -> CH4_gas (CRdesorption)
    kall(2599) = variable_crflux*3.30015378d+04

    !CH3D_0001 -> CH3D_gas (CRdesorption)
    kall(2600) = variable_crflux*3.20161944d+04

    !CH2D2_0001 -> CH2D2_gas (CRdesorption)
    kall(2601) = variable_crflux*3.11141482d+04

    !CHD3_0001 -> CHD3_gas (CRdesorption)
    kall(2602) = variable_crflux*3.02842881d+04

    !CD4_0001 -> CD4_gas (CRdesorption)
    kall(2603) = variable_crflux*2.95174728d+04

    !CO2_0001 -> CO2_gas (CRdesorption)
    kall(2604) = variable_crflux*2.18937219d-06

    !HCOOH_0001 -> HCOOH_gas (CRdesorption)
    kall(2605) = variable_crflux*1.17385241d-24

    !HCOOD_0001 -> HCOOD_gas (CRdesorption)
    kall(2606) = variable_crflux*1.16129748d-24

    !DCOOH_0001 -> DCOOH_gas (CRdesorption)
    kall(2607) = variable_crflux*1.16129748d-24

    !DCOOD_0001 -> DCOOD_gas (CRdesorption)
    kall(2608) = variable_crflux*1.14913696d-24

    !HOCO_0001 -> HOCO_gas (CRdesorption)
    kall(2609) = variable_crflux*1.00228708d-02

    !DOCO_0001 -> DOCO_gas (CRdesorption)
    kall(2610) = variable_crflux*9.91332798d-03

    !NH_0001 -> NH_gas (CRdesorption)
    kall(2611) = variable_crflux*3.74973267d-06

    !ND_0001 -> ND_gas (CRdesorption)
    kall(2612) = variable_crflux*3.63066305d-06

    !NH2_0001 -> NH2_gas (CRdesorption)
    kall(2613) = variable_crflux*7.63044962d-10

    !NHD_0001 -> NHD_gas (CRdesorption)
    kall(2614) = variable_crflux*7.40262347d-10

    !ND2_0001 -> ND2_gas (CRdesorption)
    kall(2615) = variable_crflux*7.19405690d-10

    !NH3_0001 -> NH3_gas (CRdesorption)
    kall(2616) = variable_crflux*5.21574297d-24

    !NH2D_0001 -> NH2D_gas (CRdesorption)
    kall(2617) = variable_crflux*5.06879106d-24

    !NHD2_0001 -> NHD2_gas (CRdesorption)
    kall(2618) = variable_crflux*4.93359894d-24

    !ND3_0001 -> ND3_gas (CRdesorption)
    kall(2619) = variable_crflux*4.80867742d-24

    !C10_0001 -> C10_gas (CRdesorption)
    kall(2620) = variable_crflux*7.30796498d-40

    !C10H_0001 -> C10H_gas (CRdesorption)
    kall(2621) = variable_crflux*1.20777117d-42

    !C10H2_0001 -> C10H2_gas (CRdesorption)
    kall(2622) = variable_crflux*2.89955561d-43

    !C10N_0001 -> C10N_gas (CRdesorption)
    kall(2623) = variable_crflux*7.89160821d-45

    !C11_0001 -> C11_gas (CRdesorption)
    kall(2624) = variable_crflux*9.03565672d-50

    !C2H2_0001 -> C2H2_gas (CRdesorption)
    kall(2625) = variable_crflux*3.42078274d-06

    !C2HD_0001 -> C2HD_gas (CRdesorption)
    kall(2626) = variable_crflux*3.35683724d-06

    !C2D2_0001 -> C2D2_gas (CRdesorption)
    kall(2627) = variable_crflux*3.29634873d-06

    !C2H3_0001 -> C2H3_gas (CRdesorption)
    kall(2628) = variable_crflux*1.66576993d-07

    !C2H2D_0001 -> C2H2D_gas (CRdesorption)
    kall(2629) = variable_crflux*1.63575360d-07

    !C2HD2_0001 -> C2HD2_gas (CRdesorption)
    kall(2630) = variable_crflux*1.60730354d-07

    !C2D3_0001 -> C2D3_gas (CRdesorption)
    kall(2631) = variable_crflux*1.58028811d-07

    !C2H4_0001 -> C2H4_gas (CRdesorption)
    kall(2632) = variable_crflux*1.12297720d-05

    !C2H3D_0001 -> C2H3D_gas (CRdesorption)
    kall(2633) = variable_crflux*1.10344567d-05

    !C2H2D2_0001 -> C2H2D2_gas (CRdesorption)
    kall(2634) = variable_crflux*2.28863239d-09

    !C2HD3_0001 -> C2HD3_gas (CRdesorption)
    kall(2635) = variable_crflux*1.06725724d-05

    !C2D4_0001 -> C2D4_gas (CRdesorption)
    kall(2636) = variable_crflux*1.05044899d-05

    !C2H5_0001 -> C2H5_gas (CRdesorption)
    kall(2637) = variable_crflux*2.32775714d-09

    !C2H6_0001 -> C2H6_gas (CRdesorption)
    kall(2638) = variable_crflux*3.32862990d+00

    !C3_0001 -> C3_gas (CRdesorption)
    kall(2639) = variable_crflux*4.04906743d-05

    !C3N_0001 -> C3N_gas (CRdesorption)
    kall(2640) = variable_crflux*4.31643414d-10

    !C3O_0001 -> C3O_gas (CRdesorption)
    kall(2641) = variable_crflux*2.42992520d-07

    !C3S_0001 -> C3S_gas (CRdesorption)
    kall(2642) = variable_crflux*5.32785789d-12

    !C4_0001 -> C4_gas (CRdesorption)
    kall(2643) = variable_crflux*4.40544214d-10

    !C4H_0001 -> C4H_gas (CRdesorption)
    kall(2644) = variable_crflux*2.19553085d-13

    !C4D_0001 -> C4D_gas (CRdesorption)
    kall(2645) = variable_crflux*2.17346465d-13

    !C4H2_0001 -> C4H2_gas (CRdesorption)
    kall(2646) = variable_crflux*3.71491957d-16

    !C4HD_0001 -> C4HD_gas (CRdesorption)
    kall(2647) = variable_crflux*3.67831848d-16

    !C4D2_0001 -> C4D2_gas (CRdesorption)
    kall(2648) = variable_crflux*3.64277834d-16

    !C4H3_0001 -> C4H3_gas (CRdesorption)
    kall(2649) = variable_crflux*6.25062290d-19

    !C4N_0001 -> C4N_gas (CRdesorption)
    kall(2650) = variable_crflux*4.71523970d-15

    !C4S_0001 -> C4S_gas (CRdesorption)
    kall(2651) = variable_crflux*5.92375097d-17

    !C5_0001 -> C5_gas (CRdesorption)
    kall(2652) = variable_crflux*4.79318283d-15

    !C5H_0001 -> C5H_gas (CRdesorption)
    kall(2653) = variable_crflux*2.35901118d-18

    !C5D_0001 -> C5D_gas (CRdesorption)
    kall(2654) = variable_crflux*2.33990956d-18

    !C5H2_0001 -> C5H2_gas (CRdesorption)
    kall(2655) = variable_crflux*3.96133224d-21

    !C5H3_0001 -> C5H3_gas (CRdesorption)
    kall(2656) = variable_crflux*6.62572983d-24

    !C5H4_0001 -> C5H4_gas (CRdesorption)
    kall(2657) = variable_crflux*1.10455734d-26

    !C5N_0001 -> C5N_gas (CRdesorption)
    kall(2658) = variable_crflux*5.14409376d-20

    !C5O_0001 -> C5O_gas (CRdesorption)
    kall(2659) = variable_crflux*2.99250612d-17

    !C6_0001 -> C6_gas (CRdesorption)
    kall(2660) = variable_crflux*5.21505013d-20

    !C6H_0001 -> C6H_gas (CRdesorption)
    kall(2661) = variable_crflux*2.54467423d-23

    !C6H2_0001 -> C6H2_gas (CRdesorption)
    kall(2662) = variable_crflux*4.24974433d-26

    !C6H3_0001 -> C6H3_gas (CRdesorption)
    kall(2663) = variable_crflux*7.07646013d-29

    !C6H4_0001 -> C6H4_gas (CRdesorption)
    kall(2664) = variable_crflux*1.17536997d-31

    !C6H6_0001 -> C6H6_gas (CRdesorption)
    kall(2665) = variable_crflux*3.22230964d-37

    !C6N_0001 -> C6N_gas (CRdesorption)
    kall(2666) = variable_crflux*5.60768225d-25

    !C7_0001 -> C7_gas (CRdesorption)
    kall(2667) = variable_crflux*5.67404766d-25

    !C7H_0001 -> C7H_gas (CRdesorption)
    kall(2668) = variable_crflux*2.75136438d-28

    !C7H2_0001 -> C7H2_gas (CRdesorption)
    kall(2669) = variable_crflux*4.57594737d-31

    !C7H3_0001 -> C7H3_gas (CRdesorption)
    kall(2670) = variable_crflux*7.59324592d-34

    !C7H4_0001 -> C7H4_gas (CRdesorption)
    kall(2671) = variable_crflux*1.25751391d-36

    !C7N_0001 -> C7N_gas (CRdesorption)
    kall(2672) = variable_crflux*6.11012438d-30

    !C7O_0001 -> C7O_gas (CRdesorption)
    kall(2673) = variable_crflux*3.61180822d-27

    !C8_0001 -> C8_gas (CRdesorption)
    kall(2674) = variable_crflux*6.17344342d-30

    !C8H_0001 -> C8H_gas (CRdesorption)
    kall(2675) = variable_crflux*2.97929827d-33

    !C8H2_0001 -> C8H2_gas (CRdesorption)
    kall(2676) = variable_crflux*4.93903289d-36

    !C8H3_0001 -> C8H3_gas (CRdesorption)
    kall(2677) = variable_crflux*8.17306843d-39

    !C8H4_0001 -> C8H4_gas (CRdesorption)
    kall(2678) = variable_crflux*1.35030968d-41

    !C8N_0001 -> C8N_gas (CRdesorption)
    kall(2679) = variable_crflux*6.65545115d-35

    !C9_0001 -> C9_gas (CRdesorption)
    kall(2680) = variable_crflux*6.71679301d-35

    !C9H_0001 -> C9H_gas (CRdesorption)
    kall(2681) = variable_crflux*3.22939781d-38

    !C9H2_0001 -> C9H2_gas (CRdesorption)
    kall(2682) = variable_crflux*5.33975649d-41

    !C9H3_0001 -> C9H3_gas (CRdesorption)
    kall(2683) = variable_crflux*8.81621701d-44

    !C9H4_0001 -> C9H4_gas (CRdesorption)
    kall(2684) = variable_crflux*1.45368549d-46

    !C9N_0001 -> C9N_gas (CRdesorption)
    kall(2685) = variable_crflux*7.24781609d-40

    !C9O_0001 -> C9O_gas (CRdesorption)
    kall(2686) = variable_crflux*4.32511658d-37

    !CCH_0001 -> CCH_gas (CRdesorption)
    kall(2687) = variable_crflux*1.02912264d-08

    !CCD_0001 -> CCD_gas (CRdesorption)
    kall(2688) = variable_crflux*1.00913777d-08

    !CCN_0001 -> CCN_gas (CRdesorption)
    kall(2689) = variable_crflux*3.94107284d-05

    !CCO_0001 -> CCO_gas (CRdesorption)
    kall(2690) = variable_crflux*2.14427779d-02

    !CCS_0001 -> CCS_gas (CRdesorption)
    kall(2691) = variable_crflux*4.73943359d-07

    !CD3CN_0001 -> CD3CN_gas (CRdesorption)
    kall(2692) = variable_crflux*3.65767851d-19

    !CH2CCH_0001 -> CH2CCH_gas (CRdesorption)
    kall(2693) = variable_crflux*1.18943061d-10

    !CH2CCD_0001 -> CH2CCD_gas (CRdesorption)
    kall(2694) = variable_crflux*1.17446862d-10

    !CHDCCH_0001 -> CHDCCH_gas (CRdesorption)
    kall(2695) = variable_crflux*1.17446862d-10

    !CHDCCD_0001 -> CHDCCD_gas (CRdesorption)
    kall(2696) = variable_crflux*1.16005742d-10

    !CD2CCH_0001 -> CD2CCH_gas (CRdesorption)
    kall(2697) = variable_crflux*1.16005742d-10

    !CD2CCD_0001 -> CD2CCD_gas (CRdesorption)
    kall(2698) = variable_crflux*1.14616401d-10

    !CH2CHC2H_0001 -> CH2CHC2H_gas (CRdesorption)
    kall(2699) = variable_crflux*1.04694999d-21

    !CH2CHCHCH2_0001 -> CH2CHCHCH2_gas (CRdesorption)
    kall(2700) = variable_crflux*2.90615179d-27

    !CH2CHCN_0001 -> CH2CHCN_gas (CRdesorption)
    kall(2701) = variable_crflux*3.92369987d-24

    !CH2NH_0001 -> CH2NH_gas (CRdesorption)
    kall(2702) = variable_crflux*2.46454551d-24

    !CHDNH_0001 -> CHDNH_gas (CRdesorption)
    kall(2703) = variable_crflux*2.42312163d-24

    !CHDND_0001 -> CHDND_gas (CRdesorption)
    kall(2704) = variable_crflux*2.38371865d-24

    !CH2ND_0001 -> CH2ND_gas (CRdesorption)
    kall(2705) = variable_crflux*2.42312163d-24

    !CD2NH_0001 -> CD2NH_gas (CRdesorption)
    kall(2706) = variable_crflux*2.38371865d-24

    !CD2ND_0001 -> CD2ND_gas (CRdesorption)
    kall(2707) = variable_crflux*2.34617743d-24

    !CH2NH2_0001 -> CH2NH2_gas (CRdesorption)
    kall(2708) = variable_crflux*2.42312163d-24

    !CH2NHD_0001 -> CH2NHD_gas (CRdesorption)
    kall(2709) = variable_crflux*2.38371865d-24

    !CHDNH2_0001 -> CHDNH2_gas (CRdesorption)
    kall(2710) = variable_crflux*2.38371865d-24

    !CHDNHD_0001 -> CHDNHD_gas (CRdesorption)
    kall(2711) = variable_crflux*2.34617743d-24

    !CHDND2_0001 -> CHDND2_gas (CRdesorption)
    kall(2712) = variable_crflux*2.31035582d-24

    !CH2ND2_0001 -> CH2ND2_gas (CRdesorption)
    kall(2713) = variable_crflux*2.34617743d-24

    !CD2NH2_0001 -> CD2NH2_gas (CRdesorption)
    kall(2714) = variable_crflux*2.34617743d-24

    !CD2NHD_0001 -> CD2NHD_gas (CRdesorption)
    kall(2715) = variable_crflux*2.31035582d-24

    !CD2ND2_0001 -> CD2ND2_gas (CRdesorption)
    kall(2716) = variable_crflux*2.27612643d-24

    !CH3C3N_0001 -> CH3C3N_gas (CRdesorption)
    kall(2717) = variable_crflux*2.40750599d-30

    !CH3C4H_0001 -> CH3C4H_gas (CRdesorption)
    kall(2718) = variable_crflux*1.10455734d-26

    !CH3C5N_0001 -> CH3C5N_gas (CRdesorption)
    kall(2719) = variable_crflux*4.67643625d-39

    !CH3C6H_0001 -> CH3C6H_gas (CRdesorption)
    kall(2720) = variable_crflux*1.25751391d-36

    !CH3C7N_0001 -> CH3C7N_gas (CRdesorption)
    kall(2721) = variable_crflux*5.38866077d-49

    !CH3CCH_0001 -> CH3CCH_gas (CRdesorption)
    kall(2722) = variable_crflux*9.96260252d-14

    !CH3CH2OH_0001 -> CH3CH2OH_gas (CRdesorption)
    kall(2723) = variable_crflux*1.31098569d-23

    !CH3CHCH2_0001 -> CH3CHCH2_gas (CRdesorption)
    kall(2724) = variable_crflux*1.93424740d-09

    !CH3CHO_0001 -> CH3CHO_gas (CRdesorption)
    kall(2725) = variable_crflux*1.34044972d-23

    !CH3CN_0001 -> CH3CN_gas (CRdesorption)
    kall(2726) = variable_crflux*3.78913380d-19

    !CH2DCN_0001 -> CH2DCN_gas (CRdesorption)
    kall(2727) = variable_crflux*3.74375331d-19

    !CHD2CN_0001 -> CHD2CN_gas (CRdesorption)
    kall(2728) = variable_crflux*3.69996522d-19

    !CH3COCH3_0001 -> CH3COCH3_gas (CRdesorption)
    kall(2729) = variable_crflux*5.76890107d-12

    !CH3NH2_0001 -> CH3NH2_gas (CRdesorption)
    kall(2730) = variable_crflux*7.95358816d-31

    !CH3OCH2_0001 -> CH3OCH2_gas (CRdesorption)
    kall(2731) = variable_crflux*6.54938969d-12

    !CH3OCH3_0001 -> CH3OCH3_gas (CRdesorption)
    kall(2732) = variable_crflux*9.12056697d-10

    !CN_0001 -> CN_gas (CRdesorption)
    kall(2733) = variable_crflux*1.69750174d-07

    !CS_0001 -> CS_gas (CRdesorption)
    kall(2734) = variable_crflux*4.60133425d-10

    !H2CCN_0001 -> H2CCN_gas (CRdesorption)
    kall(2735) = variable_crflux*2.25861511d-16

    !HDCCN_0001 -> HDCCN_gas (CRdesorption)
    kall(2736) = variable_crflux*2.23090099d-16

    !D2CCN_0001 -> D2CCN_gas (CRdesorption)
    kall(2737) = variable_crflux*2.20418265d-16

    !H2CCO_0001 -> H2CCO_gas (CRdesorption)
    kall(2738) = variable_crflux*1.33558722d-07

    !HDCCO_0001 -> HDCCO_gas (CRdesorption)
    kall(2739) = variable_crflux*1.31996578d-07

    !D2CCO_0001 -> D2CCO_gas (CRdesorption)
    kall(2740) = variable_crflux*1.30487996d-07

    !H2CN_0001 -> H2CN_gas (CRdesorption)
    kall(2741) = variable_crflux*4.59121091d-05

    !HDCN_0001 -> HDCN_gas (CRdesorption)
    kall(2742) = variable_crflux*4.51135768d-05

    !D2CN_0001 -> D2CN_gas (CRdesorption)
    kall(2743) = variable_crflux*4.43553114d-05

    !H2CS_0001 -> H2CS_gas (CRdesorption)
    kall(2744) = variable_crflux*1.89380079d-17

    !HDCS_0001 -> HDCS_gas (CRdesorption)
    kall(2745) = variable_crflux*1.87354565d-17

    !D2CS_0001 -> D2CS_gas (CRdesorption)
    kall(2746) = variable_crflux*1.85392683d-17

    !H2S_0001 -> H2S_gas (CRdesorption)
    kall(2747) = variable_crflux*3.31689041d-07

    !HDS_0001 -> HDS_gas (CRdesorption)
    kall(2748) = variable_crflux*3.26916288d-07

    !D2S_0001 -> D2S_gas (CRdesorption)
    kall(2749) = variable_crflux*3.22343807d-07

    !HC2O_0001 -> HC2O_gas (CRdesorption)
    kall(2750) = variable_crflux*3.79414856d-05

    !DC2O_0001 -> DC2O_gas (CRdesorption)
    kall(2751) = variable_crflux*3.74870801d-05

    !HC3N_0001 -> HC3N_gas (CRdesorption)
    kall(2752) = variable_crflux*1.40241737d-18

    !DC3N_0001 -> DC3N_gas (CRdesorption)
    kall(2753) = variable_crflux*1.38886713d-18

    !HC4N_0001 -> HC4N_gas (CRdesorption)
    kall(2754) = variable_crflux*1.48793852d-23

    !DC4N_0001 -> DC4N_gas (CRdesorption)
    kall(2755) = variable_crflux*1.47626823d-23

    !HC5N_0001 -> HC5N_gas (CRdesorption)
    kall(2756) = variable_crflux*1.59023835d-28

    !HC6N_0001 -> HC6N_gas (CRdesorption)
    kall(2757) = variable_crflux*1.96109112d-38

    !HC7N_0001 -> HC7N_gas (CRdesorption)
    kall(2758) = variable_crflux*1.83839912d-38

    !HC8N_0001 -> HC8N_gas (CRdesorption)
    kall(2759) = variable_crflux*2.25671396d-48

    !HC9N_0001 -> HC9N_gas (CRdesorption)
    kall(2760) = variable_crflux*2.14380584d-48

    !HCCNC_0001 -> HCCNC_gas (CRdesorption)
    kall(2761) = variable_crflux*1.40241737d-18

    !DCCNC_0001 -> DCCNC_gas (CRdesorption)
    kall(2762) = variable_crflux*1.38886713d-18

    !HCN_0001 -> HCN_gas (CRdesorption)
    kall(2763) = variable_crflux*4.99287472d-13

    !DCN_0001 -> DCN_gas (CRdesorption)
    kall(2764) = variable_crflux*4.90290565d-13

    !HCNCC_0001 -> HCNCC_gas (CRdesorption)
    kall(2765) = variable_crflux*1.40241737d-18

    !DCNCC_0001 -> DCNCC_gas (CRdesorption)
    kall(2766) = variable_crflux*1.38886713d-18

    !HCOOCH3_0001 -> HCOOCH3_gas (CRdesorption)
    kall(2767) = variable_crflux*3.47085976d-29

    !HCS_0001 -> HCS_gas (CRdesorption)
    kall(2768) = variable_crflux*3.14695070d-08

    !DCS_0001 -> DCS_gas (CRdesorption)
    kall(2769) = variable_crflux*3.11255676d-08

    !HNC_0001 -> HNC_gas (CRdesorption)
    kall(2770) = variable_crflux*1.21260936d-13

    !DNC_0001 -> DNC_gas (CRdesorption)
    kall(2771) = variable_crflux*1.19075876d-13

    !HNCCC_0001 -> HNCCC_gas (CRdesorption)
    kall(2772) = variable_crflux*1.40241737d-18

    !DNCCC_0001 -> DNCCC_gas (CRdesorption)
    kall(2773) = variable_crflux*1.38886713d-18

    !HNCO_0001 -> HNCO_gas (CRdesorption)
    kall(2774) = variable_crflux*1.95874987d-17

    !DNCO_0001 -> DNCO_gas (CRdesorption)
    kall(2775) = variable_crflux*1.93636342d-17

    !HNO_0001 -> HNO_gas (CRdesorption)
    kall(2776) = variable_crflux*9.24179411d-09

    !DNO_0001 -> DNO_gas (CRdesorption)
    kall(2777) = variable_crflux*9.09624495d-09

    !HS_0001 -> HS_gas (CRdesorption)
    kall(2778) = variable_crflux*6.17395537d-07

    !DS_0001 -> DS_gas (CRdesorption)
    kall(2779) = variable_crflux*6.08248431d-07

    !N2_0001 -> N2_gas (CRdesorption)
    kall(2780) = variable_crflux*1.37116280d+03

    !N2O_0001 -> N2O_gas (CRdesorption)
    kall(2781) = variable_crflux*3.66251930d-05

    !NC4N_0001 -> NC4N_gas (CRdesorption)
    kall(2782) = variable_crflux*2.50637527d-16

    !NC6N_0001 -> NC6N_gas (CRdesorption)
    kall(2783) = variable_crflux*3.03956356d-26

    !NC8N_0001 -> NC8N_gas (CRdesorption)
    kall(2784) = variable_crflux*3.64981822d-36

    !NH2CHO_0001 -> NH2CHO_gas (CRdesorption)
    kall(2785) = variable_crflux*3.73299844d-29

    !NH2CDO_0001 -> NH2CDO_gas (CRdesorption)
    kall(2786) = variable_crflux*3.69219941d-29

    !NHDCHO_0001 -> NHDCHO_gas (CRdesorption)
    kall(2787) = variable_crflux*3.69219941d-29

    !NHDCDO_0001 -> NHDCDO_gas (CRdesorption)
    kall(2788) = variable_crflux*3.65270951d-29

    !ND2CHO_0001 -> ND2CHO_gas (CRdesorption)
    kall(2789) = variable_crflux*3.65270951d-29

    !ND2CDO_0001 -> ND2CDO_gas (CRdesorption)
    kall(2790) = variable_crflux*3.61446019d-29

    !NH2CN_0001 -> NH2CN_gas (CRdesorption)
    kall(2791) = variable_crflux*1.49858087d-24

    !NHDCN_0001 -> NHDCN_gas (CRdesorption)
    kall(2792) = variable_crflux*1.48105301d-24

    !ND2CN_0001 -> ND2CN_gas (CRdesorption)
    kall(2793) = variable_crflux*1.46412613d-24

    !HSS_0001 -> HSS_gas (CRdesorption)
    kall(2794) = variable_crflux*8.90256548d-07

    !DSS_0001 -> DSS_gas (CRdesorption)
    kall(2795) = variable_crflux*8.83486438d-07

    !HSSH_0001 -> HSSH_gas (CRdesorption)
    kall(2796) = variable_crflux*1.54299564d-09

    !HSSD_0001 -> HSSD_gas (CRdesorption)
    kall(2797) = variable_crflux*1.53143746d-09

    !DSSH_0001 -> DSSH_gas (CRdesorption)
    kall(2798) = variable_crflux*1.53143746d-09

    !DSSD_0001 -> DSSD_gas (CRdesorption)
    kall(2799) = variable_crflux*1.52013518d-09

    !NO_0001 -> NO_gas (CRdesorption)
    kall(2800) = variable_crflux*3.32862990d+00

    !NO2_0001 -> NO2_gas (CRdesorption)
    kall(2801) = variable_crflux*3.58201454d-05

    !NS_0001 -> NS_gas (CRdesorption)
    kall(2802) = variable_crflux*4.03182771d-02

    !OCN_0001 -> OCN_gas (CRdesorption)
    kall(2803) = variable_crflux*3.74870801d-05

    !OCS_0001 -> OCS_gas (CRdesorption)
    kall(2804) = variable_crflux*3.22827594d-08

    !S2_0001 -> S2_gas (CRdesorption)
    kall(2805) = variable_crflux*5.06248057d-04

    !SO_0001 -> SO_gas (CRdesorption)
    kall(2806) = variable_crflux*1.24932745d-07

    !SO2_0001 -> SO2_gas (CRdesorption)
    kall(2807) = variable_crflux*2.10446298d-11

    !CH3OCHO_0001 -> CH3OCHO_gas (CRdesorption)
    kall(2808) = variable_crflux*3.47085976d-29

    !Si_0001 -> Si_gas (CRdesorption)
    kall(2809) = variable_crflux*8.42068230d-62

    !SiS_0001 -> SiS_gas (CRdesorption)
    kall(2810) = variable_crflux*8.13443090d-14

    !SiN_0001 -> SiN_gas (CRdesorption)
    kall(2811) = variable_crflux*6.77926239d-12

    !SiC_0001 -> SiC_gas (CRdesorption)
    kall(2812) = variable_crflux*6.94667679d-12

    !SiH_0001 -> SiH_gas (CRdesorption)
    kall(2813) = variable_crflux*1.80542870d-70

    !SiH2_0001 -> SiH2_gas (CRdesorption)
    kall(2814) = variable_crflux*1.94958871d-12

    !SiH3_0001 -> SiH3_gas (CRdesorption)
    kall(2815) = variable_crflux*3.28477740d-15

    !SiH4_0001 -> SiH4_gas (CRdesorption)
    kall(2816) = variable_crflux*5.50297386d-18

    !SiC2CH3_0001 -> SiC2CH3_gas (CRdesorption)
    kall(2817) = variable_crflux*3.74644301d-24

    !SiC3H_0001 -> SiC3H_gas (CRdesorption)
    kall(2818) = variable_crflux*1.31171283d-24

    !SiC3H5_0001 -> SiC3H5_gas (CRdesorption)
    kall(2819) = variable_crflux*6.26048641d-25

    !SiC4_0001 -> SiC4_gas (CRdesorption)
    kall(2820) = variable_crflux*8.42744164d-27

    !SiC4H_0001 -> SiC4H_gas (CRdesorption)
    kall(2821) = variable_crflux*4.07487076d-30

    !SiC6H_0001 -> SiC6H_gas (CRdesorption)
    kall(2822) = variable_crflux*4.70622546d-40

    !SiC8H_0001 -> SiC8H_gas (CRdesorption)
    kall(2823) = variable_crflux*4.65133496d-25

    !c_HCCHSi_0001 -> c_HCCHSi_gas (CRdesorption)
    kall(2824) = variable_crflux*2.06740782d-22

    !c_SiC2_0001 -> c_SiC2_gas (CRdesorption)
    kall(2825) = variable_crflux*7.34750880d-17

    !l_C3H_0001 -> l_C3H_gas (CRdesorption)
    kall(2826) = variable_crflux*6.10377840d-15

    !l_C3D_0001 -> l_C3D_gas (CRdesorption)
    kall(2827) = variable_crflux*6.02293008d-15

    !c_C3H_0001 -> c_C3H_gas (CRdesorption)
    kall(2828) = variable_crflux*6.10377840d-15

    !c_C3D_0001 -> c_C3D_gas (CRdesorption)
    kall(2829) = variable_crflux*6.02293008d-15

    !l_SiC3_0001 -> l_SiC3_gas (CRdesorption)
    kall(2830) = variable_crflux*7.84760531d-22

    !l_C3H2_0001 -> l_C3H2_gas (CRdesorption)
    kall(2831) = variable_crflux*3.52260741d-11

    !l_C3HD_0001 -> l_C3HD_gas (CRdesorption)
    kall(2832) = variable_crflux*3.47715251d-11

    !l_C3D2_0001 -> l_C3D2_gas (CRdesorption)
    kall(2833) = variable_crflux*3.43341300d-11

    !c_C3H2_0001 -> c_C3H2_gas (CRdesorption)
    kall(2834) = variable_crflux*3.52260741d-11

    !c_C3HD_0001 -> c_C3HD_gas (CRdesorption)
    kall(2835) = variable_crflux*3.47715251d-11

    !c_C3D2_0001 -> c_C3D2_gas (CRdesorption)
    kall(2836) = variable_crflux*3.43341300d-11

    !Mg_0001 -> Mg_gas (CRdesorption)
    kall(2837) = variable_crflux*7.50296617d-23

    !MgH_0001 -> MgH_gas (CRdesorption)
    kall(2838) = variable_crflux*1.23643558d-25

    !MgH2_0001 -> MgH2_gas (CRdesorption)
    kall(2839) = variable_crflux*2.03293573d-28

    !Na_0001 -> Na_gas (CRdesorption)
    kall(2840) = variable_crflux*5.38186763d-63

    !NaH_0001 -> NaH_gas (CRdesorption)
    kall(2841) = variable_crflux*8.66812672d-66

    !F_0001 -> F_gas (CRdesorption)
    kall(2842) = variable_crflux*2.71831541d+05

    !HF_0001 -> HF_gas (CRdesorption)
    kall(2843) = variable_crflux*2.19260885d-36

    !DF_0001 -> DF_gas (CRdesorption)
    kall(2844) = variable_crflux*2.13976714d-36

    !MgD_0001 -> MgD_gas (CRdesorption)
    kall(2845) = variable_crflux*1.21242484d-25

    !MgHD_0001 -> MgHD_gas (CRdesorption)
    kall(2846) = variable_crflux*1.99493358d-28

    !MgD2_0001 -> MgD2_gas (CRdesorption)
    kall(2847) = variable_crflux*1.95898589d-28

    !NaD_0001 -> NaD_gas (CRdesorption)
    kall(2848) = variable_crflux*8.49299500d-66

    !SiD_0001 -> SiD_gas (CRdesorption)
    kall(2849) = variable_crflux*1.77508320d-70

    !Cl_0001 -> Cl_gas (CRdesorption)
    kall(2850) = variable_crflux*8.69767376d-09

    !HCl_0001 -> HCl_gas (CRdesorption)
    kall(2851) = variable_crflux*3.66178384d-22

    !DCl_0001 -> DCl_gas (CRdesorption)
    kall(2852) = variable_crflux*3.61196133d-22

    !H2Cl_0001 -> H2Cl_gas (CRdesorption)
    kall(2853) = variable_crflux*8.54932369d-25

    !HDCl_0001 -> HDCl_gas (CRdesorption)
    kall(2854) = variable_crflux*8.43608262d-25

    !D2Cl_0001 -> D2Cl_gas (CRdesorption)
    kall(2855) = variable_crflux*8.32722538d-25

    !CCl_0001 -> CCl_gas (CRdesorption)
    kall(2856) = variable_crflux*3.98870532d-02

    !ClO_0001 -> ClO_gas (CRdesorption)
    kall(2857) = variable_crflux*3.82909193d-02

    !C3H4_0001 -> C3H4_gas (CRdesorption)
    kall(2858) = variable_crflux*9.96260252d-14

    !FeD_0001 -> FeD_gas (CRdesorption)
    kall(2859) = variable_crflux*4.87469703d-19

    !P_0001 -> P_gas (CRdesorption)
    kall(2860) = variable_crflux*3.43467053d+03

    !PO_0001 -> PO_gas (CRdesorption)
    kall(2861) = variable_crflux*3.98870532d-02

    !PH_0001 -> PH_gas (CRdesorption)
    kall(2862) = variable_crflux*6.47988578d+00

    !PD_0001 -> PD_gas (CRdesorption)
    kall(2863) = variable_crflux*6.38095042d+00

    !PH2_0001 -> PH2_gas (CRdesorption)
    kall(2864) = variable_crflux*1.17041915d-02

    !PHD_0001 -> PHD_gas (CRdesorption)
    kall(2865) = variable_crflux*1.15307865d-02

    !PD2_0001 -> PD2_gas (CRdesorption)
    kall(2866) = variable_crflux*1.13648673d-02

    !PN_0001 -> PN_gas (CRdesorption)
    kall(2867) = variable_crflux*4.07637965d-02

    !CP_0001 -> CP_gas (CRdesorption)
    kall(2868) = variable_crflux*4.17010177d-02

    !CCP_0001 -> CCP_gas (CRdesorption)
    kall(2869) = variable_crflux*7.14431248d-17

    !C3P_0001 -> C3P_gas (CRdesorption)
    kall(2870) = variable_crflux*8.97563428d-27

    !C4P_0001 -> C4P_gas (CRdesorption)
    kall(2871) = variable_crflux*1.10322124d-36

    !CH2PH_0001 -> CH2PH_gas (CRdesorption)
    kall(2872) = variable_crflux*2.14124825d-06

    !CH2PD_0001 -> CH2PD_gas (CRdesorption)
    kall(2873) = variable_crflux*2.11834654d-06

    !CHDPD_0001 -> CHDPD_gas (CRdesorption)
    kall(2874) = variable_crflux*2.09616429d-06

    !HCP_0001 -> HCP_gas (CRdesorption)
    kall(2875) = variable_crflux*7.40318456d-05

    !DCP_0001 -> DCP_gas (CRdesorption)
    kall(2876) = variable_crflux*7.32046481d-05

    !HCCP_0001 -> HCCP_gas (CRdesorption)
    kall(2877) = variable_crflux*1.20162029d-19

    !DCCP_0001 -> DCCP_gas (CRdesorption)
    kall(2878) = variable_crflux*1.19103312d-19

    !HPO_0001 -> HPO_gas (CRdesorption)
    kall(2879) = variable_crflux*7.08800958d-05

    !DPO_0001 -> DPO_gas (CRdesorption)
    kall(2880) = variable_crflux*7.01531012d-05

    !C4H2_0001 -> H_0001 + C4H_0001 (photodissociation)
    kall(2881) = 5.40000000d-03* (1.90000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.640000 * variable_Av) + 4.70000000d-03* (1.90000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4H2_0001 -> CCH_0001 + CCH_0001 (photodissociation)
    kall(2882) = 5.40000000d-03* (1.90000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.640000 * variable_Av) + 4.70000000d-03* (1.90000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C5H_0001 -> CCH_0001 + C3_0001 (photodissociation)
    kall(2883) = 5.40000000d-03* (4.29000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.760000 * variable_Av) + 4.70000000d-03* (4.29000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C5H_0001 -> C2_0001 + c_C3H_0001 (photodissociation)
    kall(2884) = 5.40000000d-03* (2.21000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.760000 * variable_Av) + 4.70000000d-03* (2.21000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C5H_0001 -> C2_0001 + l_C3H_0001 (photodissociation)
    kall(2885) = 5.40000000d-03* (2.21000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.760000 * variable_Av) + 4.70000000d-03* (2.21000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C5H_0001 -> H_0001 + C5_0001 (photodissociation)
    kall(2886) = 5.40000000d-03* (4.29000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.760000 * variable_Av) + 4.70000000d-03* (4.29000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C5N_0001 -> CN_0001 + C4_0001 (photodissociation)
    kall(2887) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3CHO_0001 -> HCO_0001 + CH3_0001 (photodissociation)
    kall(2888) = 5.40000000d-03* (9.95000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.460000 * variable_Av) + 4.70000000d-03* (9.95000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3CHO_0001 -> CO_0001 + CH4_0001 (photodissociation)
    kall(2889) = 5.40000000d-03* (9.95000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.460000 * variable_Av) + 4.70000000d-03* (9.95000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3CCH_0001 -> p_H2_0001 + c_C3H2_0001 (photodissociation)
    kall(2890) = 5.40000000d-03* (5.03000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.690000 * variable_Av) + 4.70000000d-03* (5.03000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3CCH_0001 -> o_H2_0001 + c_C3H2_0001 (photodissociation)
    kall(2891) = 5.40000000d-03* (1.31000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.690000 * variable_Av) + 4.70000000d-03* (1.31000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3CCH_0001 -> p_H2_0001 + l_C3H2_0001 (photodissociation)
    kall(2892) = 5.40000000d-03* (5.03000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.690000 * variable_Av) + 4.70000000d-03* (5.03000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3CCH_0001 -> o_H2_0001 + l_C3H2_0001 (photodissociation)
    kall(2893) = 5.40000000d-03* (1.31000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.690000 * variable_Av) + 4.70000000d-03* (1.31000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3CCH_0001 -> H_0001 + CH2CCH_0001 (photodissociation)
    kall(2894) = 5.40000000d-03* (2.94000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.720000 * variable_Av) + 4.70000000d-03* (2.94000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CHC2H_0001 -> C_0001 + CH3CCH_0001 (photodissociation)
    kall(2895) = 5.40000000d-03* (1.13000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.640000 * variable_Av) + 4.70000000d-03* (1.13000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CHCHCH2_0001 -> C_0001 + CH3CHCH2_0001 (photodissociation)
    kall(2896) = 5.40000000d-03* (1.13000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.640000 * variable_Av) + 4.70000000d-03* (1.13000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C5H2_0001 -> CCH_0001 + c_C3H_0001 (photodissociation)
    kall(2897) = 5.40000000d-03* (5.00000000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C5H2_0001 -> CCH_0001 + l_C3H_0001 (photodissociation)
    kall(2898) = 5.40000000d-03* (5.00000000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C5H2_0001 -> H_0001 + C5H_0001 (photodissociation)
    kall(2899) = 5.40000000d-03* (1.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C6H_0001 -> C3_0001 + l_C3H_0001 (photodissociation)
    kall(2900) = 5.40000000d-03* (2.50000000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.50000000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C6H_0001 -> CCH_0001 + C4_0001 (photodissociation)
    kall(2901) = 5.40000000d-03* (5.00000000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C6H_0001 -> C3_0001 + c_C3H_0001 (photodissociation)
    kall(2902) = 5.40000000d-03* (2.50000000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.50000000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3NH2_0001 -> H_0001 + CN_0001 + p_H2_0001 + o_H2_0001 (photodissociation)
    kall(2903) = 5.40000000d-03* (1.53000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.370000 * variable_Av) + 4.70000000d-03* (1.53000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3NH2_0001 -> H_0001 + CN_0001 + o_H2_0001 + p_H2_0001 (photodissociation)
    kall(2904) = 5.40000000d-03* (1.53000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.370000 * variable_Av) + 4.70000000d-03* (1.53000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3NH2_0001 -> H_0001 + CN_0001 + o_H2_0001 + o_H2_0001 (photodissociation)
    kall(2905) = 5.40000000d-03* (6.44500000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.370000 * variable_Av) + 4.70000000d-03* (6.44500000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3NH2_0001 -> H_0001 + CN_0001 + p_H2_0001 + p_H2_0001 (photodissociation)
    kall(2906) = 5.40000000d-03* (6.81500000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.370000 * variable_Av) + 4.70000000d-03* (6.81500000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3NH2_0001 -> H_0001 + H_0001 + p_H2_0001 + HCN_0001 (photodissociation)
    kall(2907) = 5.40000000d-03* (8.63000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.370000 * variable_Av) + 4.70000000d-03* (8.63000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3NH2_0001 -> H_0001 + H_0001 + o_H2_0001 + HCN_0001 (photodissociation)
    kall(2908) = 5.40000000d-03* (3.12000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.370000 * variable_Av) + 4.70000000d-03* (3.12000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3NH2_0001 -> NH2_0001 + CH3_0001 (photodissociation)
    kall(2909) = 5.40000000d-03* (1.86000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.370000 * variable_Av) + 4.70000000d-03* (1.86000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3NH2_0001 -> H_0001 + H_0001 + CH2NH_0001 (photodissociation)
    kall(2910) = 5.40000000d-03* (4.04000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.370000 * variable_Av) + 4.70000000d-03* (4.04000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HC5N_0001 -> H_0001 + C5N_0001 (photodissociation)
    kall(2911) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.830000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HC5N_0001 -> CN_0001 + C4H_0001 (photodissociation)
    kall(2912) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.830000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3C3N_0001 -> C3N_0001 + CH3_0001 (photodissociation)
    kall(2913) = 5.40000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NC4N_0001 -> C2_0001 + CN_0001 + CN_0001 (photodissociation)
    kall(2914) = 5.40000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NC6N_0001 -> CN_0001 + CN_0001 + C4_0001 (photodissociation)
    kall(2915) = 5.40000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NC8N_0001 -> CN_0001 + CN_0001 + C6_0001 (photodissociation)
    kall(2916) = 5.40000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HC4N_0001 -> C2_0001 + CH_0001 + CN_0001 (photodissociation)
    kall(2917) = 5.40000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HC6N_0001 -> CH_0001 + CN_0001 + C4_0001 (photodissociation)
    kall(2918) = 5.40000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HC8N_0001 -> CH_0001 + CN_0001 + C6_0001 (photodissociation)
    kall(2919) = 5.40000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiC2CH3_0001 -> p_H2_0001 + SiC3H_0001 (photodissociation)
    kall(2920) = 5.40000000d-03* (2.50000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (2.50000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiC2CH3_0001 -> o_H2_0001 + SiC3H_0001 (photodissociation)
    kall(2921) = 5.40000000d-03* (7.50000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (7.50000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiC3H5_0001 -> o_H2_0001 + p_H2_0001 + SiC3H_0001 (photodissociation)
    kall(2922) = 5.40000000d-03* (1.50000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (1.50000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiC3H5_0001 -> p_H2_0001 + p_H2_0001 + SiC3H_0001 (photodissociation)
    kall(2923) = 5.40000000d-03* (6.67000000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (6.67000000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiC3H5_0001 -> p_H2_0001 + o_H2_0001 + SiC3H_0001 (photodissociation)
    kall(2924) = 5.40000000d-03* (1.50000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (1.50000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiC3H5_0001 -> o_H2_0001 + o_H2_0001 + SiC3H_0001 (photodissociation)
    kall(2925) = 5.40000000d-03* (6.33000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (6.33000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiC4H_0001 -> Si_0001 + C4H_0001 (photodissociation)
    kall(2926) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiC6H_0001 -> Si_0001 + C6H_0001 (photodissociation)
    kall(2927) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiC8H_0001 -> Si_0001 + C8H_0001 (photodissociation)
    kall(2928) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3CH2OH_0001 -> H2O_0001 + C2H4_0001 (photodissociation)
    kall(2929) = 5.40000000d-03* (1.38000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.730000 * variable_Av) + 4.70000000d-03* (1.38000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3CH2OH_0001 -> OH_0001 + C2H5_0001 (photodissociation)
    kall(2930) = 5.40000000d-03* (2.90000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.730000 * variable_Av) + 4.70000000d-03* (2.90000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3C4H_0001 -> CH3_0001 + C4H_0001 (photodissociation)
    kall(2931) = 5.40000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3OCH3_0001 -> H2CO_0001 + CH4_0001 (photodissociation)
    kall(2932) = 5.40000000d-03* (1.50000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.600000 * variable_Av) + 4.70000000d-03* (1.50000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2_0001 -> C_0001 + C_0001 (photodissociation)
    kall(2933) = 5.40000000d-03* (2.35000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.040000 * variable_Av) + 4.70000000d-03* (2.35000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CCl_0001 -> C_0001 + Cl_0001 (photodissociation)
    kall(2934) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.000000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH_0001 -> C_0001 + H_0001 (photodissociation)
    kall(2935) = 5.40000000d-03* (9.14000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.120000 * variable_Av) + 4.70000000d-03* (9.14000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !ClO_0001 -> Cl_0001 + O_0001 (photodissociation)
    kall(2936) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.000000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CN_0001 -> C_0001 + N_0001 (photodissociation)
    kall(2937) = 5.40000000d-03* (5.19000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.500000 * variable_Av) + 4.70000000d-03* (5.19000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CO_0001 -> C_0001 + O_0001 (photodissociation)
    kall(2938) = 5.40000000d-03* ss_CO * (2.43000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.880000 * variable_Av) + 4.70000000d-03* (2.43000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CS_0001 -> C_0001 + S_0001 (photodissociation)
    kall(2939) = 5.40000000d-03* (9.49000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.770000 * variable_Av) + 4.70000000d-03* (9.49000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !p_H2_0001 -> H_0001 + H_0001 (photodissociation)
    kall(2940) = 5.40000000d-03* ss_H2 * (5.68000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-4.180000 * variable_Av) + 4.70000000d-03* (5.68000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !o_H2_0001 -> H_0001 + H_0001 (photodissociation)
    kall(2941) = 5.40000000d-03* ss_H2 * (5.68000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-4.180000 * variable_Av) + 4.70000000d-03* (5.68000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HCl_0001 -> Cl_0001 + H_0001 (photodissociation)
    kall(2942) = 5.40000000d-03* (1.73000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.880000 * variable_Av) + 4.70000000d-03* (1.73000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HF_0001 -> H_0001 + F_0001 (photodissociation)
    kall(2943) = 5.40000000d-03* (1.17000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.210000 * variable_Av) + 4.70000000d-03* (1.17000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HS_0001 -> H_0001 + S_0001 (photodissociation)
    kall(2944) = 5.40000000d-03* (1.23000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.400000 * variable_Av) + 4.70000000d-03* (1.23000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !N2_0001 -> N_0001 + N_0001 (photodissociation)
    kall(2945) = 5.40000000d-03* ss_N2 * (1.70000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-4.250000 * variable_Av) + 4.70000000d-03* (1.70000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NH_0001 -> H_0001 + N_0001 (photodissociation)
    kall(2946) = 5.40000000d-03* (5.75000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.630000 * variable_Av) + 4.70000000d-03* (5.75000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NO_0001 -> N_0001 + O_0001 (photodissociation)
    kall(2947) = 5.40000000d-03* (3.84000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.560000 * variable_Av) + 4.70000000d-03* (3.84000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NS_0001 -> N_0001 + S_0001 (photodissociation)
    kall(2948) = 5.40000000d-03* (1.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.000000 * variable_Av) + 4.70000000d-03* (1.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !O2_0001 -> O_0001 + O_0001 (photodissociation)
    kall(2949) = 5.40000000d-03* (7.71000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.450000 * variable_Av) + 4.70000000d-03* (7.71000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !OH_0001 -> H_0001 + O_0001 (photodissociation)
    kall(2950) = 5.40000000d-03* (3.79000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.660000 * variable_Av) + 4.70000000d-03* (3.79000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiC_0001 -> C_0001 + Si_0001 (photodissociation)
    kall(2951) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiH_0001 -> H_0001 + Si_0001 (photodissociation)
    kall(2952) = 5.40000000d-03* (2.71000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.950000 * variable_Av) + 4.70000000d-03* (2.71000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiS_0001 -> S_0001 + Si_0001 (photodissociation)
    kall(2953) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SO_0001 -> O_0001 + S_0001 (photodissociation)
    kall(2954) = 5.40000000d-03* (4.25000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (4.25000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CCH_0001 -> H_0001 + C2_0001 (photodissociation)
    kall(2955) = 5.40000000d-03* (1.58000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.670000 * variable_Av) + 4.70000000d-03* (1.58000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CCN_0001 -> N_0001 + C2_0001 (photodissociation)
    kall(2956) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CCN_0001 -> C_0001 + CN_0001 (photodissociation)
    kall(2957) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CCS_0001 -> S_0001 + C2_0001 (photodissociation)
    kall(2958) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.000000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C3_0001 -> C_0001 + C2_0001 (photodissociation)
    kall(2959) = 5.40000000d-03* (4.97000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.390000 * variable_Av) + 4.70000000d-03* (4.97000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2_0001 -> H_0001 + CH_0001 (photodissociation)
    kall(2960) = 5.40000000d-03* (5.80000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.350000 * variable_Av) + 4.70000000d-03* (5.80000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CO2_0001 -> O_0001 + CO_0001 (photodissociation)
    kall(2961) = (Fnot * Gnot* kph_factor * exp(-gamma_CO2 * variable_Av) + kph_factor * F_cr ) * 2.40000000d-03

    !H2S_0001 -> H_0001 + HS_0001 (photodissociation)
    kall(2962) = 5.40000000d-03* (3.12000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.640000 * variable_Av) + 4.70000000d-03* (3.12000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HCN_0001 -> H_0001 + CN_0001 (photodissociation)
    kall(2963) = 5.40000000d-03* (1.64000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.120000 * variable_Av) + 4.70000000d-03* (1.64000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HCO_0001 -> H_0001 + CO_0001 (photodissociation)
    kall(2964) = 5.40000000d-03* (1.10000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.430000 * variable_Av) + 4.70000000d-03* (1.10000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HNC_0001 -> H_0001 + CN_0001 (photodissociation)
    kall(2965) = 5.40000000d-03* (1.64000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.120000 * variable_Av) + 4.70000000d-03* (1.64000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HNO_0001 -> H_0001 + NO_0001 (photodissociation)
    kall(2966) = 5.40000000d-03* (1.70000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-0.530000 * variable_Av) + 4.70000000d-03* (1.70000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NH2_0001 -> H_0001 + NH_0001 (photodissociation)
    kall(2967) = 5.40000000d-03* (9.51000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.310000 * variable_Av) + 4.70000000d-03* (9.51000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NO2_0001 -> O_0001 + NO_0001 (photodissociation)
    kall(2968) = 5.40000000d-03* (1.39000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.500000 * variable_Av) + 4.70000000d-03* (1.39000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !OCN_0001 -> O_0001 + CN_0001 (photodissociation)
    kall(2969) = 5.40000000d-03* (1.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.000000 * variable_Av) + 4.70000000d-03* (1.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !OCS_0001 -> S_0001 + CO_0001 (photodissociation)
    kall(2970) = 5.40000000d-03* (4.67000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.460000 * variable_Av) + 4.70000000d-03* (4.67000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SO2_0001 -> O_0001 + SO_0001 (photodissociation)
    kall(2971) = 5.40000000d-03* (2.43000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.780000 * variable_Av) + 4.70000000d-03* (2.43000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H2_0001 -> H_0001 + CCH_0001 (photodissociation)
    kall(2972) = 5.40000000d-03* (2.41000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.640000 * variable_Av) + 4.70000000d-03* (2.41000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !c_C3H_0001 -> H_0001 + C3_0001 (photodissociation)
    kall(2973) = 5.40000000d-03* (1.06000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.150000 * variable_Av) + 4.70000000d-03* (1.06000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !l_C3H_0001 -> H_0001 + C3_0001 (photodissociation)
    kall(2974) = 5.40000000d-03* (1.77000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.080000 * variable_Av) + 4.70000000d-03* (1.77000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C3N_0001 -> C2_0001 + CN_0001 (photodissociation)
    kall(2975) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C3O_0001 -> C2_0001 + CO_0001 (photodissociation)
    kall(2976) = 5.40000000d-03* (7.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.580000 * variable_Av) + 4.70000000d-03* (7.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C3S_0001 -> C2_0001 + CS_0001 (photodissociation)
    kall(2977) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.000000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4_0001 -> C2_0001 + C2_0001 (photodissociation)
    kall(2978) = 5.40000000d-03* (1.27200000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (1.27200000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4_0001 -> C_0001 + C3_0001 (photodissociation)
    kall(2979) = 5.40000000d-03* (7.20800000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (7.20800000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3_0001 -> CH_0001 + p_H2_0001 (photodissociation)
    kall(2980) = 5.40000000d-03* (5.11500000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.500000 * variable_Av) + 4.70000000d-03* (5.11500000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3_0001 -> CH_0001 + o_H2_0001 (photodissociation)
    kall(2981) = 5.40000000d-03* (1.53500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.500000 * variable_Av) + 4.70000000d-03* (1.53500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3_0001 -> H_0001 + CH2_0001 (photodissociation)
    kall(2982) = 5.40000000d-03* (2.05000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.500000 * variable_Av) + 4.70000000d-03* (2.05000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !H2CO_0001 -> H_0001 + H_0001 + CO_0001 (photodissociation)
    kall(2983) = 5.40000000d-03* (7.10000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.540000 * variable_Av) + 4.70000000d-03* (7.10000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !H2CO_0001 -> CO_0001 + p_H2_0001 (photodissociation)
    kall(2984) = 5.40000000d-03* (3.55000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.540000 * variable_Av) + 4.70000000d-03* (3.55000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !H2CO_0001 -> CO_0001 + o_H2_0001 (photodissociation)
    kall(2985) = 5.40000000d-03* (3.55000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.540000 * variable_Av) + 4.70000000d-03* (3.55000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !H2CS_0001 -> CS_0001 + p_H2_0001 (photodissociation)
    kall(2986) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !H2CS_0001 -> CS_0001 + o_H2_0001 (photodissociation)
    kall(2987) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NH3_0001 -> H_0001 + NH2_0001 (photodissociation)
    kall(2988) = 5.40000000d-03* (8.30000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.420000 * variable_Av) + 4.70000000d-03* (8.30000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NH3_0001 -> p_H2_0001 + NH_0001 (photodissociation)
    kall(2989) = 5.40000000d-03* (1.52500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.200000 * variable_Av) + 4.70000000d-03* (1.52500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NH3_0001 -> o_H2_0001 + NH_0001 (photodissociation)
    kall(2990) = 5.40000000d-03* (4.57500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.200000 * variable_Av) + 4.70000000d-03* (4.57500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !H2CCO_0001 -> CO_0001 + CH2_0001 (photodissociation)
    kall(2991) = 5.40000000d-03* (1.40000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.580000 * variable_Av) + 4.70000000d-03* (1.40000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H3_0001 -> H_0001 + C2H2_0001 (photodissociation)
    kall(2992) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !c_C3H2_0001 -> p_H2_0001 + C3_0001 (photodissociation)
    kall(2993) = 5.40000000d-03* (3.47500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.260000 * variable_Av) + 4.70000000d-03* (3.47500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !c_C3H2_0001 -> o_H2_0001 + C3_0001 (photodissociation)
    kall(2994) = 5.40000000d-03* (3.47500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.260000 * variable_Av) + 4.70000000d-03* (3.47500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !c_C3H2_0001 -> H_0001 + c_C3H_0001 (photodissociation)
    kall(2995) = 5.40000000d-03* (6.95000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.260000 * variable_Av) + 4.70000000d-03* (6.95000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !l_C3H2_0001 -> p_H2_0001 + C3_0001 (photodissociation)
    kall(2996) = 5.40000000d-03* (1.03000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.510000 * variable_Av) + 4.70000000d-03* (1.03000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !l_C3H2_0001 -> o_H2_0001 + C3_0001 (photodissociation)
    kall(2997) = 5.40000000d-03* (1.03000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.510000 * variable_Av) + 4.70000000d-03* (1.03000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !l_C3H2_0001 -> H_0001 + l_C3H_0001 (photodissociation)
    kall(2998) = 5.40000000d-03* (2.06000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.510000 * variable_Av) + 4.70000000d-03* (2.06000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4H_0001 -> C2_0001 + CCH_0001 (photodissociation)
    kall(2999) = 5.40000000d-03* (1.83000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.360000 * variable_Av) + 4.70000000d-03* (1.83000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4H_0001 -> H_0001 + C4_0001 (photodissociation)
    kall(3000) = 5.40000000d-03* (1.83000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.360000 * variable_Av) + 4.70000000d-03* (1.83000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4N_0001 -> CN_0001 + C3_0001 (photodissociation)
    kall(3001) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4S_0001 -> CS_0001 + C3_0001 (photodissociation)
    kall(3002) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.000000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C5_0001 -> C_0001 + C4_0001 (photodissociation)
    kall(3003) = 5.40000000d-03* (1.50000000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (1.50000000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C5_0001 -> C2_0001 + C3_0001 (photodissociation)
    kall(3004) = 5.40000000d-03* (8.50000000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (8.50000000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HCOOH_0001 -> OH_0001 + HCO_0001 (photodissociation)
    kall(3005) = 5.40000000d-03* (4.10000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (4.10000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2NH_0001 -> p_H2_0001 + HCN_0001 (photodissociation)
    kall(3006) = 5.40000000d-03* (8.75000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (8.75000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2NH_0001 -> o_H2_0001 + HCN_0001 (photodissociation)
    kall(3007) = 5.40000000d-03* (2.62000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (2.62000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH4_0001 -> H_0001 + CH_0001 + o_H2_0001 (photodissociation)
    kall(3008) = 5.40000000d-03* (2.19000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.080000 * variable_Av) + 4.70000000d-03* (2.19000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH4_0001 -> H_0001 + CH_0001 + p_H2_0001 (photodissociation)
    kall(3009) = 5.40000000d-03* (8.44900000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.080000 * variable_Av) + 4.70000000d-03* (8.44900000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH4_0001 -> p_H2_0001 + CH2_0001 (photodissociation)
    kall(3010) = 5.40000000d-03* (6.58500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.080000 * variable_Av) + 4.70000000d-03* (6.58500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH4_0001 -> o_H2_0001 + CH2_0001 (photodissociation)
    kall(3011) = 5.40000000d-03* (2.53300000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.080000 * variable_Av) + 4.70000000d-03* (2.53300000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH4_0001 -> H_0001 + CH3_0001 (photodissociation)
    kall(3012) = 5.40000000d-03* (3.04000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.080000 * variable_Av) + 4.70000000d-03* (3.04000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HC3N_0001 -> CN_0001 + CCH_0001 (photodissociation)
    kall(3013) = 5.40000000d-03* (7.13000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (7.13000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3CN_0001 -> CN_0001 + CH3_0001 (photodissociation)
    kall(3014) = 5.40000000d-03* (2.95000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.070000 * variable_Av) + 4.70000000d-03* (2.95000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H4_0001 -> p_H2_0001 + C2H2_0001 (photodissociation)
    kall(3015) = 5.40000000d-03* (8.49700000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (8.49700000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H4_0001 -> o_H2_0001 + C2H2_0001 (photodissociation)
    kall(3016) = 5.40000000d-03* (2.21300000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (2.21300000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCH_0001 -> p_H2_0001 + c_C3H_0001 (photodissociation)
    kall(3017) = 5.40000000d-03* (1.25000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.25000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCH_0001 -> o_H2_0001 + c_C3H_0001 (photodissociation)
    kall(3018) = 5.40000000d-03* (3.75000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.75000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCH_0001 -> p_H2_0001 + l_C3H_0001 (photodissociation)
    kall(3019) = 5.40000000d-03* (1.25000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.25000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCH_0001 -> o_H2_0001 + l_C3H_0001 (photodissociation)
    kall(3020) = 5.40000000d-03* (3.75000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.75000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCH_0001 -> H_0001 + c_C3H2_0001 (photodissociation)
    kall(3021) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCH_0001 -> H_0001 + l_C3H2_0001 (photodissociation)
    kall(3022) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C6_0001 -> C_0001 + C5_0001 (photodissociation)
    kall(3023) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C6_0001 -> C2_0001 + C4_0001 (photodissociation)
    kall(3024) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C6_0001 -> C3_0001 + C3_0001 (photodissociation)
    kall(3025) = 5.40000000d-03* (8.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C7_0001 -> C2_0001 + C5_0001 (photodissociation)
    kall(3026) = 5.40000000d-03* (2.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C7_0001 -> C3_0001 + C4_0001 (photodissociation)
    kall(3027) = 5.40000000d-03* (8.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C8_0001 -> C_0001 + C7_0001 (photodissociation)
    kall(3028) = 5.40000000d-03* (5.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C8_0001 -> C3_0001 + C5_0001 (photodissociation)
    kall(3029) = 5.40000000d-03* (9.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (9.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C8_0001 -> C4_0001 + C4_0001 (photodissociation)
    kall(3030) = 5.40000000d-03* (5.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C9_0001 -> C4_0001 + C5_0001 (photodissociation)
    kall(3031) = 5.40000000d-03* (3.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C9_0001 -> C2_0001 + C7_0001 (photodissociation)
    kall(3032) = 5.40000000d-03* (5.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C9_0001 -> C3_0001 + C6_0001 (photodissociation)
    kall(3033) = 5.40000000d-03* (6.50000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (6.50000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C10_0001 -> C3_0001 + C7_0001 (photodissociation)
    kall(3034) = 5.40000000d-03* (7.50000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (7.50000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C10_0001 -> C5_0001 + C5_0001 (photodissociation)
    kall(3035) = 5.40000000d-03* (2.50000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.50000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C10H2_0001 -> H_0001 + C10H_0001 (photodissociation)
    kall(3036) = 5.40000000d-03* (5.27000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.100000 * variable_Av) + 4.70000000d-03* (5.27000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C10H2_0001 -> H_0001 + H_0001 + C10_0001 (photodissociation)
    kall(3037) = 5.40000000d-03* (3.10000000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-5.200000 * variable_Av) + 4.70000000d-03* (3.10000000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C10H2_0001 -> H_0001 + C3_0001 + C7H_0001 (photodissociation)
    kall(3038) = 5.40000000d-03* (6.14000000d-13/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-5.200000 * variable_Av) + 4.70000000d-03* (6.14000000d-13/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C11_0001 -> C3_0001 + C8_0001 (photodissociation)
    kall(3039) = 5.40000000d-03* (4.41000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.700000 * variable_Av) + 4.70000000d-03* (4.41000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C11_0001 -> C5_0001 + C6_0001 (photodissociation)
    kall(3040) = 5.40000000d-03* (5.28000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.100000 * variable_Av) + 4.70000000d-03* (5.28000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C11_0001 -> C3_0001 + C3_0001 + C5_0001 (photodissociation)
    kall(3041) = 5.40000000d-03* (3.83000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-4.500000 * variable_Av) + 4.70000000d-03* (3.83000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C11_0001 -> C4_0001 + C7_0001 (photodissociation)
    kall(3042) = 5.40000000d-03* (2.31000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.300000 * variable_Av) + 4.70000000d-03* (2.31000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C11_0001 -> C2_0001 + C9_0001 (photodissociation)
    kall(3043) = 5.40000000d-03* (2.86000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.800000 * variable_Av) + 4.70000000d-03* (2.86000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C10H_0001 -> C5_0001 + C5H_0001 (photodissociation)
    kall(3044) = 5.40000000d-03* (3.19000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.000000 * variable_Av) + 4.70000000d-03* (3.19000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C10H_0001 -> H_0001 + C3_0001 + C7_0001 (photodissociation)
    kall(3045) = 5.40000000d-03* (5.95000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.900000 * variable_Av) + 4.70000000d-03* (5.95000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C10H_0001 -> H_0001 + C5_0001 + C5_0001 (photodissociation)
    kall(3046) = 5.40000000d-03* (2.30000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.800000 * variable_Av) + 4.70000000d-03* (2.30000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C10H_0001 -> C_0001 + H_0001 + C9_0001 (photodissociation)
    kall(3047) = 5.40000000d-03* (2.07000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-4.600000 * variable_Av) + 4.70000000d-03* (2.07000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C10H_0001 -> H_0001 + C2_0001 + C8_0001 (photodissociation)
    kall(3048) = 5.40000000d-03* (1.42000000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-5.000000 * variable_Av) + 4.70000000d-03* (1.42000000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C10H_0001 -> H_0001 + C10_0001 (photodissociation)
    kall(3049) = 5.40000000d-03* (7.36000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.600000 * variable_Av) + 4.70000000d-03* (7.36000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C10H_0001 -> C3_0001 + C7H_0001 (photodissociation)
    kall(3050) = 5.40000000d-03* (1.86000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.600000 * variable_Av) + 4.70000000d-03* (1.86000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3C6H_0001 -> CH3_0001 + C6H_0001 (photodissociation)
    kall(3051) = 5.40000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-5.000000 * variable_Av) + 4.70000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3C5N_0001 -> CH3_0001 + C5N_0001 (photodissociation)
    kall(3052) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3C7N_0001 -> CH3_0001 + C7N_0001 (photodissociation)
    kall(3053) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CCO_0001 -> O_0001 + C2_0001 (photodissociation)
    kall(3054) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CCO_0001 -> C_0001 + CO_0001 (photodissociation)
    kall(3055) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HOOH_0001 -> OH_0001 + OH_0001 (photodissociation)
    kall(3056) = 5.40000000d-03* (8.08000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.610000 * variable_Av) + 4.70000000d-03* (8.08000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H5_0001 -> o_H2_0001 + C2H3_0001 (photodissociation)
    kall(3057) = 5.40000000d-03* (7.83000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (7.83000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H5_0001 -> p_H2_0001 + C2H3_0001 (photodissociation)
    kall(3058) = 5.40000000d-03* (2.17000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.17000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CHCN_0001 -> CN_0001 + C2H3_0001 (photodissociation)
    kall(3059) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3COCH3_0001 -> CO_0001 + CH3_0001 + CH3_0001 (photodissociation)
    kall(3060) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C6H2_0001 -> H_0001 + C6H_0001 (photodissociation)
    kall(3061) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C7H_0001 -> H_0001 + C7_0001 (photodissociation)
    kall(3062) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C7H2_0001 -> H_0001 + C7H_0001 (photodissociation)
    kall(3063) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C7N_0001 -> CN_0001 + C6_0001 (photodissociation)
    kall(3064) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C8H_0001 -> H_0001 + C8_0001 (photodissociation)
    kall(3065) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C8H2_0001 -> H_0001 + C8H_0001 (photodissociation)
    kall(3066) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C9H_0001 -> H_0001 + C9_0001 (photodissociation)
    kall(3067) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C9H2_0001 -> H_0001 + C9H_0001 (photodissociation)
    kall(3068) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C9N_0001 -> CN_0001 + C8_0001 (photodissociation)
    kall(3069) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HC7N_0001 -> CN_0001 + C6H_0001 (photodissociation)
    kall(3070) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HC9N_0001 -> CN_0001 + C8H_0001 (photodissociation)
    kall(3071) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CP_0001 -> C_0001 + P_0001 (photodissociation)
    kall(3072) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.800000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !MgH_0001 -> H_0001 + Mg_0001 (photodissociation)
    kall(3073) = 5.40000000d-03* (5.07000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (5.07000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NaH_0001 -> H_0001 + Na_0001 (photodissociation)
    kall(3074) = 5.40000000d-03* (7.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.770000 * variable_Av) + 4.70000000d-03* (7.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !PH_0001 -> H_0001 + P_0001 (photodissociation)
    kall(3075) = 5.40000000d-03* (5.78000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.480000 * variable_Av) + 4.70000000d-03* (5.78000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !PN_0001 -> N_0001 + P_0001 (photodissociation)
    kall(3076) = 5.40000000d-03* (5.00000000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.000000 * variable_Av) + 4.70000000d-03* (5.00000000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !PO_0001 -> O_0001 + P_0001 (photodissociation)
    kall(3077) = 5.40000000d-03* (3.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.000000 * variable_Av) + 4.70000000d-03* (3.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !S2_0001 -> S_0001 + S_0001 (photodissociation)
    kall(3078) = 5.40000000d-03* (6.56000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.900000 * variable_Av) + 4.70000000d-03* (6.56000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiN_0001 -> N_0001 + Si_0001 (photodissociation)
    kall(3079) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.800000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CCP_0001 -> P_0001 + C2_0001 (photodissociation)
    kall(3080) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CCP_0001 -> C_0001 + CP_0001 (photodissociation)
    kall(3081) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HCP_0001 -> H_0001 + CP_0001 (photodissociation)
    kall(3082) = 5.40000000d-03* (5.48000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.000000 * variable_Av) + 4.70000000d-03* (5.48000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HPO_0001 -> H_0001 + PO_0001 (photodissociation)
    kall(3083) = 5.40000000d-03* (1.70000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-0.530000 * variable_Av) + 4.70000000d-03* (1.70000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HSS_0001 -> S_0001 + HS_0001 (photodissociation)
    kall(3084) = 5.40000000d-03* (3.30000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.400000 * variable_Av) + 4.70000000d-03* (3.30000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !N2O_0001 -> O_0001 + N2_0001 (photodissociation)
    kall(3085) = 5.40000000d-03* (1.87000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.810000 * variable_Av) + 4.70000000d-03* (1.87000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !PH2_0001 -> H_0001 + PH_0001 (photodissociation)
    kall(3086) = 5.40000000d-03* (2.11000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.520000 * variable_Av) + 4.70000000d-03* (2.11000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !c_SiC2_0001 -> Si_0001 + C2_0001 (photodissociation)
    kall(3087) = 5.40000000d-03* (2.60000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.280000 * variable_Av) + 4.70000000d-03* (2.60000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiH2_0001 -> H_0001 + SiH_0001 (photodissociation)
    kall(3088) = 5.40000000d-03* (5.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C3P_0001 -> C2_0001 + CP_0001 (photodissociation)
    kall(3089) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HSSH_0001 -> HS_0001 + HS_0001 (photodissociation)
    kall(3090) = 5.40000000d-03* (8.30000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (8.30000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !l_SiC3_0001 -> C_0001 + c_SiC2_0001 (photodissociation)
    kall(3091) = 5.40000000d-03* (2.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (2.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !l_SiC3_0001 -> C2_0001 + SiC_0001 (photodissociation)
    kall(3092) = 5.40000000d-03* (2.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (2.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiH3_0001 -> p_H2_0001 + SiH_0001 (photodissociation)
    kall(3093) = 5.40000000d-03* (7.50000000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (7.50000000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiH3_0001 -> o_H2_0001 + SiH_0001 (photodissociation)
    kall(3094) = 5.40000000d-03* (2.25000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.25000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiH3_0001 -> H_0001 + SiH2_0001 (photodissociation)
    kall(3095) = 5.40000000d-03* (3.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !H2CCN_0001 -> CN_0001 + CH2_0001 (photodissociation)
    kall(3096) = 5.40000000d-03* (1.56000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.950000 * variable_Av) + 4.70000000d-03* (1.56000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3CHCH2_0001 -> CH2_0001 + C2H4_0001 (photodissociation)
    kall(3097) = 5.40000000d-03* (1.13000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.640000 * variable_Av) + 4.70000000d-03* (1.13000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4P_0001 -> CP_0001 + C3_0001 (photodissociation)
    kall(3098) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2PH_0001 -> PH_0001 + CH2_0001 (photodissociation)
    kall(3099) = 5.40000000d-03* (9.54000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.830000 * variable_Av) + 4.70000000d-03* (9.54000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HCNCC_0001 -> C2_0001 + HCN_0001 (photodissociation)
    kall(3100) = 5.40000000d-03* (9.54000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.830000 * variable_Av) + 4.70000000d-03* (9.54000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HCCNC_0001 -> CN_0001 + CCH_0001 (photodissociation)
    kall(3101) = 5.40000000d-03* (9.54000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.830000 * variable_Av) + 4.70000000d-03* (9.54000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HNCCC_0001 -> H_0001 + C3N_0001 (photodissociation)
    kall(3102) = 5.40000000d-03* (9.54000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.830000 * variable_Av) + 4.70000000d-03* (9.54000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiC3H_0001 -> H_0001 + l_SiC3_0001 (photodissociation)
    kall(3103) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiC4_0001 -> Si_0001 + C4_0001 (photodissociation)
    kall(3104) = 5.40000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.300000 * variable_Av) + 4.70000000d-03* (1.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiH4_0001 -> H_0001 + o_H2_0001 + SiH_0001 (photodissociation)
    kall(3105) = 5.40000000d-03* (1.16000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.200000 * variable_Av) + 4.70000000d-03* (1.16000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiH4_0001 -> H_0001 + p_H2_0001 + SiH_0001 (photodissociation)
    kall(3106) = 5.40000000d-03* (4.44000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.200000 * variable_Av) + 4.70000000d-03* (4.44000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiH4_0001 -> p_H2_0001 + SiH2_0001 (photodissociation)
    kall(3107) = 5.40000000d-03* (1.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.200000 * variable_Av) + 4.70000000d-03* (1.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiH4_0001 -> o_H2_0001 + SiH2_0001 (photodissociation)
    kall(3108) = 5.40000000d-03* (3.47000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.200000 * variable_Av) + 4.70000000d-03* (3.47000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiH4_0001 -> H_0001 + SiH3_0001 (photodissociation)
    kall(3109) = 5.40000000d-03* (1.60000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.200000 * variable_Av) + 4.70000000d-03* (1.60000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4H3_0001 -> p_H2_0001 + C4H_0001 (photodissociation)
    kall(3110) = 5.40000000d-03* (2.50000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.50000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4H3_0001 -> o_H2_0001 + C4H_0001 (photodissociation)
    kall(3111) = 5.40000000d-03* (7.50000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (7.50000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4H3_0001 -> H_0001 + C4H2_0001 (photodissociation)
    kall(3112) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HCOOCH3_0001 -> H2CO_0001 + H2CO_0001 (photodissociation)
    kall(3113) = 5.40000000d-03* (1.38000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.730000 * variable_Av) + 4.70000000d-03* (1.38000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C6H6_0001 -> C2H2_0001 + CH2CHC2H_0001 (photodissociation)
    kall(3114) = 5.40000000d-03* (1.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C6H6_0001 -> CH3CCH_0001 + l_C3H2_0001 (photodissociation)
    kall(3115) = 5.40000000d-03* (5.00000000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C6H6_0001 -> c_C3H2_0001 + CH3CCH_0001 (photodissociation)
    kall(3116) = 5.40000000d-03* (5.00000000d-12/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-12/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4HD_0001 -> H_0001 + C4D_0001 (photodissociation)
    kall(3117) = 5.40000000d-03* (9.50000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.640000 * variable_Av) + 4.70000000d-03* (9.50000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4HD_0001 -> D_0001 + C4H_0001 (photodissociation)
    kall(3118) = 5.40000000d-03* (9.50000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.640000 * variable_Av) + 4.70000000d-03* (9.50000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4D2_0001 -> D_0001 + C4D_0001 (photodissociation)
    kall(3119) = 5.40000000d-03* (1.90000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.640000 * variable_Av) + 4.70000000d-03* (1.90000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4HD_0001 -> CCH_0001 + CCD_0001 (photodissociation)
    kall(3120) = 5.40000000d-03* (9.50000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.640000 * variable_Av) + 4.70000000d-03* (9.50000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4D2_0001 -> CCD_0001 + CCD_0001 (photodissociation)
    kall(3121) = 5.40000000d-03* (1.90000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.640000 * variable_Av) + 4.70000000d-03* (1.90000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C5D_0001 -> CCD_0001 + C3_0001 (photodissociation)
    kall(3122) = 5.40000000d-03* (4.29000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.760000 * variable_Av) + 4.70000000d-03* (4.29000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C5D_0001 -> C2_0001 + c_C3D_0001 (photodissociation)
    kall(3123) = 5.40000000d-03* (2.21000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.760000 * variable_Av) + 4.70000000d-03* (2.21000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C5D_0001 -> C2_0001 + l_C3D_0001 (photodissociation)
    kall(3124) = 5.40000000d-03* (2.21000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.760000 * variable_Av) + 4.70000000d-03* (2.21000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C5D_0001 -> D_0001 + C5_0001 (photodissociation)
    kall(3125) = 5.40000000d-03* (4.29000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.760000 * variable_Av) + 4.70000000d-03* (4.29000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2DOD_0001 -> p_H2_0001 + D2CO_0001 (photodissociation)
    kall(3126) = 5.40000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2DOD_0001 -> o_H2_0001 + D2CO_0001 (photodissociation)
    kall(3127) = 5.40000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2DOD_0001 -> HD_0001 + HDCO_0001 (photodissociation)
    kall(3128) = 5.40000000d-03* (4.95200000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (4.95200000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2DOD_0001 -> o_D2_0001 + H2CO_0001 (photodissociation)
    kall(3129) = 5.40000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2DOD_0001 -> p_D2_0001 + H2CO_0001 (photodissociation)
    kall(3130) = 5.40000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2OH_0001 -> p_H2_0001 + D2CO_0001 (photodissociation)
    kall(3131) = 5.40000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2OH_0001 -> o_H2_0001 + D2CO_0001 (photodissociation)
    kall(3132) = 5.40000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2OH_0001 -> HD_0001 + HDCO_0001 (photodissociation)
    kall(3133) = 5.40000000d-03* (4.95200000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (4.95200000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2OH_0001 -> o_D2_0001 + H2CO_0001 (photodissociation)
    kall(3134) = 5.40000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2OH_0001 -> p_D2_0001 + H2CO_0001 (photodissociation)
    kall(3135) = 5.40000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (6.19600000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2OD_0001 -> HD_0001 + D2CO_0001 (photodissociation)
    kall(3136) = 5.40000000d-03* (3.71500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (3.71500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2OD_0001 -> p_D2_0001 + HDCO_0001 (photodissociation)
    kall(3137) = 5.40000000d-03* (1.85800000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (1.85800000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2OD_0001 -> o_D2_0001 + HDCO_0001 (photodissociation)
    kall(3138) = 5.40000000d-03* (1.85800000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (1.85800000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD3OH_0001 -> HD_0001 + D2CO_0001 (photodissociation)
    kall(3139) = 5.40000000d-03* (3.71500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (3.71500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD3OH_0001 -> p_D2_0001 + HDCO_0001 (photodissociation)
    kall(3140) = 5.40000000d-03* (1.85800000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (1.85800000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD3OH_0001 -> o_D2_0001 + HDCO_0001 (photodissociation)
    kall(3141) = 5.40000000d-03* (1.85800000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (1.85800000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD3OD_0001 -> o_D2_0001 + D2CO_0001 (photodissociation)
    kall(3142) = 5.40000000d-03* (4.64400000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (4.64400000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD3OD_0001 -> p_D2_0001 + D2CO_0001 (photodissociation)
    kall(3143) = 5.40000000d-03* (2.78600000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (2.78600000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2DOD_0001 -> OH_0001 + CHD2_0001 (photodissociation)
    kall(3144) = 5.40000000d-03* (3.18500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (3.18500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2DOD_0001 -> OD_0001 + CH2D_0001 (photodissociation)
    kall(3145) = 5.40000000d-03* (3.18500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (3.18500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2OH_0001 -> OH_0001 + CHD2_0001 (photodissociation)
    kall(3146) = 5.40000000d-03* (3.18500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (3.18500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2OH_0001 -> OD_0001 + CH2D_0001 (photodissociation)
    kall(3147) = 5.40000000d-03* (3.18500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (3.18500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2OD_0001 -> OH_0001 + CD3_0001 (photodissociation)
    kall(3148) = 5.40000000d-03* (1.59300000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (1.59300000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2OD_0001 -> OD_0001 + CHD2_0001 (photodissociation)
    kall(3149) = 5.40000000d-03* (4.77800000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (4.77800000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD3OH_0001 -> OH_0001 + CD3_0001 (photodissociation)
    kall(3150) = 5.40000000d-03* (1.59300000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (1.59300000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD3OH_0001 -> OD_0001 + CHD2_0001 (photodissociation)
    kall(3151) = 5.40000000d-03* (4.77800000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (4.77800000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD3OD_0001 -> OD_0001 + CD3_0001 (photodissociation)
    kall(3152) = 5.40000000d-03* (6.37000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (6.37000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DC4N_0001 -> C2_0001 + CD_0001 + CN_0001 (photodissociation)
    kall(3153) = 5.40000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD_0001 -> C_0001 + D_0001 (photodissociation)
    kall(3154) = 5.40000000d-03* (9.14000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.120000 * variable_Av) + 4.70000000d-03* (9.14000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HD_0001 -> H_0001 + D_0001 (photodissociation)
    kall(3155) = 5.40000000d-03* ss_HD * (5.68000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-4.180000 * variable_Av) + 4.70000000d-03* (5.68000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !p_D2_0001 -> D_0001 + D_0001 (photodissociation)
    kall(3156) = 5.40000000d-03* (5.68000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-4.180000 * variable_Av) + 4.70000000d-03* (5.68000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !o_D2_0001 -> D_0001 + D_0001 (photodissociation)
    kall(3157) = 5.40000000d-03* (5.68000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-4.180000 * variable_Av) + 4.70000000d-03* (5.68000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DCl_0001 -> Cl_0001 + D_0001 (photodissociation)
    kall(3158) = 5.40000000d-03* (1.10000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (1.10000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DF_0001 -> D_0001 + F_0001 (photodissociation)
    kall(3159) = 5.40000000d-03* (1.17000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.210000 * variable_Av) + 4.70000000d-03* (1.17000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DS_0001 -> D_0001 + S_0001 (photodissociation)
    kall(3160) = 5.40000000d-03* (9.80000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.040000 * variable_Av) + 4.70000000d-03* (9.80000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !ND_0001 -> D_0001 + N_0001 (photodissociation)
    kall(3161) = 5.40000000d-03* (5.75000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.630000 * variable_Av) + 4.70000000d-03* (5.75000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !OD_0001 -> D_0001 + O_0001 (photodissociation)
    kall(3162) = 5.40000000d-03* (3.79000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.660000 * variable_Av) + 4.70000000d-03* (3.79000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !SiD_0001 -> D_0001 + Si_0001 (photodissociation)
    kall(3163) = 5.40000000d-03* (2.80000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.590000 * variable_Av) + 4.70000000d-03* (2.80000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CCD_0001 -> D_0001 + C2_0001 (photodissociation)
    kall(3164) = 5.40000000d-03* (1.58000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.670000 * variable_Av) + 4.70000000d-03* (1.58000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD_0001 -> H_0001 + CD_0001 (photodissociation)
    kall(3165) = 5.40000000d-03* (2.90000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.350000 * variable_Av) + 4.70000000d-03* (2.90000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD_0001 -> D_0001 + CH_0001 (photodissociation)
    kall(3166) = 5.40000000d-03* (2.90000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.350000 * variable_Av) + 4.70000000d-03* (2.90000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2_0001 -> D_0001 + CD_0001 (photodissociation)
    kall(3167) = 5.40000000d-03* (5.80000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.350000 * variable_Av) + 4.70000000d-03* (5.80000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HDS_0001 -> H_0001 + DS_0001 (photodissociation)
    kall(3168) = 5.40000000d-03* (1.55000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (1.55000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HDS_0001 -> D_0001 + HS_0001 (photodissociation)
    kall(3169) = 5.40000000d-03* (1.55000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (1.55000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !D2S_0001 -> D_0001 + DS_0001 (photodissociation)
    kall(3170) = 5.40000000d-03* (3.10000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (3.10000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DCN_0001 -> D_0001 + CN_0001 (photodissociation)
    kall(3171) = 5.40000000d-03* (1.60000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.690000 * variable_Av) + 4.70000000d-03* (1.60000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DCO_0001 -> D_0001 + CO_0001 (photodissociation)
    kall(3172) = 5.40000000d-03* (1.10000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.090000 * variable_Av) + 4.70000000d-03* (1.10000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DNC_0001 -> D_0001 + CN_0001 (photodissociation)
    kall(3173) = 5.40000000d-03* (1.60000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.690000 * variable_Av) + 4.70000000d-03* (1.60000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DNO_0001 -> D_0001 + NO_0001 (photodissociation)
    kall(3174) = 5.40000000d-03* (1.70000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-0.530000 * variable_Av) + 4.70000000d-03* (1.70000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NHD_0001 -> H_0001 + ND_0001 (photodissociation)
    kall(3175) = 5.40000000d-03* (4.75500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.310000 * variable_Av) + 4.70000000d-03* (4.75500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NHD_0001 -> D_0001 + NH_0001 (photodissociation)
    kall(3176) = 5.40000000d-03* (4.75500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.310000 * variable_Av) + 4.70000000d-03* (4.75500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !ND2_0001 -> D_0001 + ND_0001 (photodissociation)
    kall(3177) = 5.40000000d-03* (9.51000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.310000 * variable_Av) + 4.70000000d-03* (9.51000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H2_0001 -> H_0001 + CCH_0001 (photodissociation)
    kall(3178) = 5.40000000d-03* (2.41000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.640000 * variable_Av) + 4.70000000d-03* (2.41000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2HD_0001 -> H_0001 + CCD_0001 (photodissociation)
    kall(3179) = 5.40000000d-03* (1.20500000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.640000 * variable_Av) + 4.70000000d-03* (1.20500000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2HD_0001 -> D_0001 + CCH_0001 (photodissociation)
    kall(3180) = 5.40000000d-03* (1.20500000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.640000 * variable_Av) + 4.70000000d-03* (1.20500000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2D2_0001 -> D_0001 + CCD_0001 (photodissociation)
    kall(3181) = 5.40000000d-03* (2.41000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.640000 * variable_Av) + 4.70000000d-03* (2.41000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !c_C3D_0001 -> D_0001 + C3_0001 (photodissociation)
    kall(3182) = 5.40000000d-03* (1.06000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.150000 * variable_Av) + 4.70000000d-03* (1.06000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !l_C3D_0001 -> D_0001 + C3_0001 (photodissociation)
    kall(3183) = 5.40000000d-03* (1.77000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.080000 * variable_Av) + 4.70000000d-03* (1.77000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D_0001 -> CH_0001 + HD_0001 (photodissociation)
    kall(3184) = 5.40000000d-03* (9.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (9.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D_0001 -> CD_0001 + p_H2_0001 (photodissociation)
    kall(3185) = 5.40000000d-03* (2.25000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (2.25000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D_0001 -> CD_0001 + o_H2_0001 (photodissociation)
    kall(3186) = 5.40000000d-03* (2.25000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (2.25000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2_0001 -> CH_0001 + p_D2_0001 (photodissociation)
    kall(3187) = 5.40000000d-03* (2.25000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (2.25000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2_0001 -> CH_0001 + o_D2_0001 (photodissociation)
    kall(3188) = 5.40000000d-03* (2.25000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (2.25000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2_0001 -> CD_0001 + HD_0001 (photodissociation)
    kall(3189) = 5.40000000d-03* (9.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (9.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD3_0001 -> CD_0001 + p_D2_0001 (photodissociation)
    kall(3190) = 5.40000000d-03* (6.75000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (6.75000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD3_0001 -> CD_0001 + o_D2_0001 (photodissociation)
    kall(3191) = 5.40000000d-03* (6.75000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (6.75000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D_0001 -> H_0001 + CHD_0001 (photodissociation)
    kall(3192) = 5.40000000d-03* (9.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (9.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D_0001 -> D_0001 + CH2_0001 (photodissociation)
    kall(3193) = 5.40000000d-03* (4.50000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (4.50000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2_0001 -> H_0001 + CD2_0001 (photodissociation)
    kall(3194) = 5.40000000d-03* (4.50000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (4.50000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2_0001 -> D_0001 + CHD_0001 (photodissociation)
    kall(3195) = 5.40000000d-03* (9.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (9.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD3_0001 -> D_0001 + CD2_0001 (photodissociation)
    kall(3196) = 5.40000000d-03* (1.35000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.270000 * variable_Av) + 4.70000000d-03* (1.35000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HDCO_0001 -> H_0001 + D_0001 + CO_0001 (photodissociation)
    kall(3197) = 5.40000000d-03* (7.10000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.540000 * variable_Av) + 4.70000000d-03* (7.10000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !D2CO_0001 -> D_0001 + D_0001 + CO_0001 (photodissociation)
    kall(3198) = 5.40000000d-03* (7.10000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.540000 * variable_Av) + 4.70000000d-03* (7.10000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HDCO_0001 -> CO_0001 + HD_0001 (photodissociation)
    kall(3199) = 5.40000000d-03* (7.10000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.540000 * variable_Av) + 4.70000000d-03* (7.10000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !D2CO_0001 -> CO_0001 + p_D2_0001 (photodissociation)
    kall(3200) = 5.40000000d-03* (3.55000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.540000 * variable_Av) + 4.70000000d-03* (3.55000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !D2CO_0001 -> CO_0001 + o_D2_0001 (photodissociation)
    kall(3201) = 5.40000000d-03* (3.55000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.540000 * variable_Av) + 4.70000000d-03* (3.55000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HDCS_0001 -> CS_0001 + HD_0001 (photodissociation)
    kall(3202) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !D2CS_0001 -> CS_0001 + p_D2_0001 (photodissociation)
    kall(3203) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !D2CS_0001 -> CS_0001 + o_D2_0001 (photodissociation)
    kall(3204) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NH2D_0001 -> H_0001 + NHD_0001 (photodissociation)
    kall(3205) = 5.40000000d-03* (5.53300000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.420000 * variable_Av) + 4.70000000d-03* (5.53300000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NH2D_0001 -> D_0001 + NH2_0001 (photodissociation)
    kall(3206) = 5.40000000d-03* (2.76700000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.420000 * variable_Av) + 4.70000000d-03* (2.76700000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NHD2_0001 -> H_0001 + ND2_0001 (photodissociation)
    kall(3207) = 5.40000000d-03* (2.76700000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.420000 * variable_Av) + 4.70000000d-03* (2.76700000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NHD2_0001 -> D_0001 + NHD_0001 (photodissociation)
    kall(3208) = 5.40000000d-03* (5.53300000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.420000 * variable_Av) + 4.70000000d-03* (5.53300000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !ND3_0001 -> D_0001 + ND2_0001 (photodissociation)
    kall(3209) = 5.40000000d-03* (8.30000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.420000 * variable_Av) + 4.70000000d-03* (8.30000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NH2D_0001 -> p_H2_0001 + ND_0001 (photodissociation)
    kall(3210) = 5.40000000d-03* (1.01700000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.200000 * variable_Av) + 4.70000000d-03* (1.01700000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NH2D_0001 -> o_H2_0001 + ND_0001 (photodissociation)
    kall(3211) = 5.40000000d-03* (1.01700000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.200000 * variable_Av) + 4.70000000d-03* (1.01700000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NH2D_0001 -> HD_0001 + NH_0001 (photodissociation)
    kall(3212) = 5.40000000d-03* (4.06700000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.200000 * variable_Av) + 4.70000000d-03* (4.06700000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NHD2_0001 -> HD_0001 + ND_0001 (photodissociation)
    kall(3213) = 5.40000000d-03* (4.06700000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.200000 * variable_Av) + 4.70000000d-03* (4.06700000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NHD2_0001 -> p_D2_0001 + NH_0001 (photodissociation)
    kall(3214) = 5.40000000d-03* (1.01700000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.200000 * variable_Av) + 4.70000000d-03* (1.01700000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NHD2_0001 -> o_D2_0001 + NH_0001 (photodissociation)
    kall(3215) = 5.40000000d-03* (1.01700000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.200000 * variable_Av) + 4.70000000d-03* (1.01700000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !ND3_0001 -> p_D2_0001 + ND_0001 (photodissociation)
    kall(3216) = 5.40000000d-03* (3.05000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.200000 * variable_Av) + 4.70000000d-03* (3.05000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !ND3_0001 -> o_D2_0001 + ND_0001 (photodissociation)
    kall(3217) = 5.40000000d-03* (3.05000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.200000 * variable_Av) + 4.70000000d-03* (3.05000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HDCCO_0001 -> CO_0001 + CHD_0001 (photodissociation)
    kall(3218) = 5.40000000d-03* (1.40000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.580000 * variable_Av) + 4.70000000d-03* (1.40000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !D2CCO_0001 -> CO_0001 + CD2_0001 (photodissociation)
    kall(3219) = 5.40000000d-03* (1.40000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.580000 * variable_Av) + 4.70000000d-03* (1.40000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H2D_0001 -> H_0001 + C2HD_0001 (photodissociation)
    kall(3220) = 5.40000000d-03* (6.67000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (6.67000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H2D_0001 -> D_0001 + C2H2_0001 (photodissociation)
    kall(3221) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2HD2_0001 -> H_0001 + C2D2_0001 (photodissociation)
    kall(3222) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2HD2_0001 -> D_0001 + C2HD_0001 (photodissociation)
    kall(3223) = 5.40000000d-03* (6.67000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (6.67000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2D3_0001 -> D_0001 + C2D2_0001 (photodissociation)
    kall(3224) = 5.40000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.00000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !c_C3HD_0001 -> HD_0001 + C3_0001 (photodissociation)
    kall(3225) = 5.40000000d-03* (6.95000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.260000 * variable_Av) + 4.70000000d-03* (6.95000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !c_C3D2_0001 -> p_D2_0001 + C3_0001 (photodissociation)
    kall(3226) = 5.40000000d-03* (3.47500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.260000 * variable_Av) + 4.70000000d-03* (3.47500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !c_C3D2_0001 -> o_D2_0001 + C3_0001 (photodissociation)
    kall(3227) = 5.40000000d-03* (3.47500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.260000 * variable_Av) + 4.70000000d-03* (3.47500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !c_C3HD_0001 -> H_0001 + c_C3D_0001 (photodissociation)
    kall(3228) = 5.40000000d-03* (3.47500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.260000 * variable_Av) + 4.70000000d-03* (3.47500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !c_C3HD_0001 -> D_0001 + c_C3H_0001 (photodissociation)
    kall(3229) = 5.40000000d-03* (3.47500000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.260000 * variable_Av) + 4.70000000d-03* (3.47500000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !c_C3D2_0001 -> D_0001 + c_C3D_0001 (photodissociation)
    kall(3230) = 5.40000000d-03* (6.95000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.260000 * variable_Av) + 4.70000000d-03* (6.95000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !l_C3HD_0001 -> HD_0001 + C3_0001 (photodissociation)
    kall(3231) = 5.40000000d-03* (2.06000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.510000 * variable_Av) + 4.70000000d-03* (2.06000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !l_C3D2_0001 -> p_D2_0001 + C3_0001 (photodissociation)
    kall(3232) = 5.40000000d-03* (1.03000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.510000 * variable_Av) + 4.70000000d-03* (1.03000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !l_C3D2_0001 -> o_D2_0001 + C3_0001 (photodissociation)
    kall(3233) = 5.40000000d-03* (1.03000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.510000 * variable_Av) + 4.70000000d-03* (1.03000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !l_C3HD_0001 -> H_0001 + l_C3D_0001 (photodissociation)
    kall(3234) = 5.40000000d-03* (1.03000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.510000 * variable_Av) + 4.70000000d-03* (1.03000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !l_C3HD_0001 -> D_0001 + l_C3H_0001 (photodissociation)
    kall(3235) = 5.40000000d-03* (1.03000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.510000 * variable_Av) + 4.70000000d-03* (1.03000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !l_C3D2_0001 -> D_0001 + l_C3D_0001 (photodissociation)
    kall(3236) = 5.40000000d-03* (2.06000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.510000 * variable_Av) + 4.70000000d-03* (2.06000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4D_0001 -> C2_0001 + CCD_0001 (photodissociation)
    kall(3237) = 5.40000000d-03* (1.83000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.360000 * variable_Av) + 4.70000000d-03* (1.83000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C4D_0001 -> D_0001 + C4_0001 (photodissociation)
    kall(3238) = 5.40000000d-03* (1.83000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.360000 * variable_Av) + 4.70000000d-03* (1.83000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HCOOD_0001 -> OH_0001 + DCO_0001 (photodissociation)
    kall(3239) = 5.40000000d-03* (2.05000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (2.05000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HCOOD_0001 -> OD_0001 + HCO_0001 (photodissociation)
    kall(3240) = 5.40000000d-03* (2.05000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (2.05000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DCOOH_0001 -> OH_0001 + DCO_0001 (photodissociation)
    kall(3241) = 5.40000000d-03* (2.05000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (2.05000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DCOOH_0001 -> OD_0001 + HCO_0001 (photodissociation)
    kall(3242) = 5.40000000d-03* (2.05000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (2.05000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DCOOD_0001 -> OD_0001 + DCO_0001 (photodissociation)
    kall(3243) = 5.40000000d-03* (4.10000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (4.10000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2ND_0001 -> p_H2_0001 + DCN_0001 (photodissociation)
    kall(3244) = 5.40000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2ND_0001 -> o_H2_0001 + DCN_0001 (photodissociation)
    kall(3245) = 5.40000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2ND_0001 -> HD_0001 + HCN_0001 (photodissociation)
    kall(3246) = 5.40000000d-03* (2.33000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (2.33000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDNH_0001 -> p_H2_0001 + DCN_0001 (photodissociation)
    kall(3247) = 5.40000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDNH_0001 -> o_H2_0001 + DCN_0001 (photodissociation)
    kall(3248) = 5.40000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDNH_0001 -> HD_0001 + HCN_0001 (photodissociation)
    kall(3249) = 5.40000000d-03* (2.33000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (2.33000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDND_0001 -> HD_0001 + DCN_0001 (photodissociation)
    kall(3250) = 5.40000000d-03* (2.33000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (2.33000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDND_0001 -> p_D2_0001 + HCN_0001 (photodissociation)
    kall(3251) = 5.40000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDND_0001 -> o_D2_0001 + HCN_0001 (photodissociation)
    kall(3252) = 5.40000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2NH_0001 -> HD_0001 + DCN_0001 (photodissociation)
    kall(3253) = 5.40000000d-03* (2.33000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (2.33000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2NH_0001 -> p_D2_0001 + HCN_0001 (photodissociation)
    kall(3254) = 5.40000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2NH_0001 -> o_D2_0001 + HCN_0001 (photodissociation)
    kall(3255) = 5.40000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (5.85000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2ND_0001 -> p_D2_0001 + DCN_0001 (photodissociation)
    kall(3256) = 5.40000000d-03* (1.75000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (1.75000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2ND_0001 -> o_D2_0001 + DCN_0001 (photodissociation)
    kall(3257) = 5.40000000d-03* (1.75000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.630000 * variable_Av) + 4.70000000d-03* (1.75000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3D_0001 -> H_0001 + CH_0001 + HD_0001 (photodissociation)
    kall(3258) = 5.40000000d-03* (1.20000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (1.20000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3D_0001 -> H_0001 + CD_0001 + o_H2_0001 (photodissociation)
    kall(3259) = 5.40000000d-03* (4.50000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (4.50000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3D_0001 -> H_0001 + CD_0001 + p_H2_0001 (photodissociation)
    kall(3260) = 5.40000000d-03* (1.50000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (1.50000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3D_0001 -> D_0001 + CH_0001 + o_H2_0001 (photodissociation)
    kall(3261) = 5.40000000d-03* (4.50000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (4.50000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3D_0001 -> D_0001 + CH_0001 + p_H2_0001 (photodissociation)
    kall(3262) = 5.40000000d-03* (1.50000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (1.50000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D2_0001 -> H_0001 + CH_0001 + o_D2_0001 (photodissociation)
    kall(3263) = 5.40000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D2_0001 -> H_0001 + CH_0001 + p_D2_0001 (photodissociation)
    kall(3264) = 5.40000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D2_0001 -> H_0001 + CD_0001 + HD_0001 (photodissociation)
    kall(3265) = 5.40000000d-03* (8.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (8.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D2_0001 -> D_0001 + CH_0001 + HD_0001 (photodissociation)
    kall(3266) = 5.40000000d-03* (8.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (8.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D2_0001 -> D_0001 + CD_0001 + p_H2_0001 (photodissociation)
    kall(3267) = 5.40000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D2_0001 -> D_0001 + CD_0001 + o_H2_0001 (photodissociation)
    kall(3268) = 5.40000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (2.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD3_0001 -> H_0001 + CD_0001 + p_D2_0001 (photodissociation)
    kall(3269) = 5.40000000d-03* (3.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (3.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD3_0001 -> H_0001 + CD_0001 + o_D2_0001 (photodissociation)
    kall(3270) = 5.40000000d-03* (3.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (3.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD3_0001 -> D_0001 + CH_0001 + p_D2_0001 (photodissociation)
    kall(3271) = 5.40000000d-03* (3.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (3.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD3_0001 -> D_0001 + CH_0001 + o_D2_0001 (photodissociation)
    kall(3272) = 5.40000000d-03* (3.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (3.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD3_0001 -> D_0001 + CD_0001 + HD_0001 (photodissociation)
    kall(3273) = 5.40000000d-03* (1.20000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (1.20000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD4_0001 -> D_0001 + CD_0001 + p_D2_0001 (photodissociation)
    kall(3274) = 5.40000000d-03* (9.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (9.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD4_0001 -> D_0001 + CD_0001 + o_D2_0001 (photodissociation)
    kall(3275) = 5.40000000d-03* (1.50000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (1.50000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3D_0001 -> p_H2_0001 + CHD_0001 (photodissociation)
    kall(3276) = 5.40000000d-03* (9.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (9.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3D_0001 -> o_H2_0001 + CHD_0001 (photodissociation)
    kall(3277) = 5.40000000d-03* (2.70000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (2.70000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3D_0001 -> HD_0001 + CH2_0001 (photodissociation)
    kall(3278) = 5.40000000d-03* (3.60000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (3.60000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D2_0001 -> p_H2_0001 + CD2_0001 (photodissociation)
    kall(3279) = 5.40000000d-03* (6.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (6.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D2_0001 -> o_H2_0001 + CD2_0001 (photodissociation)
    kall(3280) = 5.40000000d-03* (6.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (6.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D2_0001 -> HD_0001 + CHD_0001 (photodissociation)
    kall(3281) = 5.40000000d-03* (4.80000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (4.80000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D2_0001 -> o_D2_0001 + CH2_0001 (photodissociation)
    kall(3282) = 5.40000000d-03* (6.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (6.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D2_0001 -> p_D2_0001 + CH2_0001 (photodissociation)
    kall(3283) = 5.40000000d-03* (6.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (6.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD3_0001 -> HD_0001 + CD2_0001 (photodissociation)
    kall(3284) = 5.40000000d-03* (3.60000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (3.60000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD3_0001 -> p_D2_0001 + CHD_0001 (photodissociation)
    kall(3285) = 5.40000000d-03* (1.80000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (1.80000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD3_0001 -> o_D2_0001 + CHD_0001 (photodissociation)
    kall(3286) = 5.40000000d-03* (1.80000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (1.80000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD4_0001 -> o_D2_0001 + CD2_0001 (photodissociation)
    kall(3287) = 5.40000000d-03* (4.50000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (4.50000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD4_0001 -> p_D2_0001 + CD2_0001 (photodissociation)
    kall(3288) = 5.40000000d-03* (2.70000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (2.70000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3D_0001 -> H_0001 + CH2D_0001 (photodissociation)
    kall(3289) = 5.40000000d-03* (1.80000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (1.80000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3D_0001 -> D_0001 + CH3_0001 (photodissociation)
    kall(3290) = 5.40000000d-03* (6.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (6.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D2_0001 -> H_0001 + CHD2_0001 (photodissociation)
    kall(3291) = 5.40000000d-03* (1.20000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (1.20000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2D2_0001 -> D_0001 + CH2D_0001 (photodissociation)
    kall(3292) = 5.40000000d-03* (1.20000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (1.20000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD3_0001 -> H_0001 + CD3_0001 (photodissociation)
    kall(3293) = 5.40000000d-03* (6.00000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (6.00000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD3_0001 -> D_0001 + CHD2_0001 (photodissociation)
    kall(3294) = 5.40000000d-03* (1.80000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (1.80000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD4_0001 -> D_0001 + CD3_0001 (photodissociation)
    kall(3295) = 5.40000000d-03* (2.40000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.590000 * variable_Av) + 4.70000000d-03* (2.40000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DC3N_0001 -> CN_0001 + CCD_0001 (photodissociation)
    kall(3296) = 5.40000000d-03* (5.60000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.160000 * variable_Av) + 4.70000000d-03* (5.60000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2DCN_0001 -> CN_0001 + CH2D_0001 (photodissociation)
    kall(3297) = 5.40000000d-03* (2.95000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.070000 * variable_Av) + 4.70000000d-03* (2.95000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHD2CN_0001 -> CN_0001 + CHD2_0001 (photodissociation)
    kall(3298) = 5.40000000d-03* (2.95000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.070000 * variable_Av) + 4.70000000d-03* (2.95000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD3CN_0001 -> CN_0001 + CD3_0001 (photodissociation)
    kall(3299) = 5.40000000d-03* (2.95000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-3.070000 * variable_Av) + 4.70000000d-03* (2.95000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H4_0001 -> p_H2_0001 + C2H2_0001 (photodissociation)
    kall(3300) = 5.40000000d-03* (8.49700000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (8.49700000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H4_0001 -> o_H2_0001 + C2H2_0001 (photodissociation)
    kall(3301) = 5.40000000d-03* (2.21300000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (2.21300000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H3D_0001 -> p_H2_0001 + C2HD_0001 (photodissociation)
    kall(3302) = 5.40000000d-03* (3.86300000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (3.86300000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H3D_0001 -> o_H2_0001 + C2HD_0001 (photodissociation)
    kall(3303) = 5.40000000d-03* (1.15900000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (1.15900000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H3D_0001 -> HD_0001 + C2H2_0001 (photodissociation)
    kall(3304) = 5.40000000d-03* (1.54500000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (1.54500000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H2D2_0001 -> p_H2_0001 + C2D2_0001 (photodissociation)
    kall(3305) = 5.40000000d-03* (2.55000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (2.55000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H2D2_0001 -> o_H2_0001 + C2D2_0001 (photodissociation)
    kall(3306) = 5.40000000d-03* (2.55000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (2.55000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H2D2_0001 -> HD_0001 + C2HD_0001 (photodissociation)
    kall(3307) = 5.40000000d-03* (2.04000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (2.04000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H2D2_0001 -> o_D2_0001 + C2H2_0001 (photodissociation)
    kall(3308) = 5.40000000d-03* (2.55000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (2.55000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H2D2_0001 -> p_D2_0001 + C2H2_0001 (photodissociation)
    kall(3309) = 5.40000000d-03* (2.55000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (2.55000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2HD3_0001 -> HD_0001 + C2D2_0001 (photodissociation)
    kall(3310) = 5.40000000d-03* (1.54500000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (1.54500000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2HD3_0001 -> p_D2_0001 + C2HD_0001 (photodissociation)
    kall(3311) = 5.40000000d-03* (7.65000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (7.65000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2HD3_0001 -> o_D2_0001 + C2HD_0001 (photodissociation)
    kall(3312) = 5.40000000d-03* (7.65000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (7.65000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2D4_0001 -> o_D2_0001 + C2D2_0001 (photodissociation)
    kall(3313) = 5.40000000d-03* (1.90700000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (1.90700000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2D4_0001 -> p_D2_0001 + C2D2_0001 (photodissociation)
    kall(3314) = 5.40000000d-03* (1.14200000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.490000 * variable_Av) + 4.70000000d-03* (1.14200000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCD_0001 -> p_H2_0001 + c_C3D_0001 (photodissociation)
    kall(3315) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCD_0001 -> o_H2_0001 + c_C3D_0001 (photodissociation)
    kall(3316) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCD_0001 -> HD_0001 + c_C3H_0001 (photodissociation)
    kall(3317) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCH_0001 -> p_H2_0001 + c_C3D_0001 (photodissociation)
    kall(3318) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCH_0001 -> o_H2_0001 + c_C3D_0001 (photodissociation)
    kall(3319) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCH_0001 -> HD_0001 + c_C3H_0001 (photodissociation)
    kall(3320) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCD_0001 -> HD_0001 + c_C3D_0001 (photodissociation)
    kall(3321) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCD_0001 -> p_D2_0001 + c_C3H_0001 (photodissociation)
    kall(3322) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCD_0001 -> o_D2_0001 + c_C3H_0001 (photodissociation)
    kall(3323) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCH_0001 -> HD_0001 + c_C3D_0001 (photodissociation)
    kall(3324) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCH_0001 -> p_D2_0001 + c_C3H_0001 (photodissociation)
    kall(3325) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCH_0001 -> o_D2_0001 + c_C3H_0001 (photodissociation)
    kall(3326) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCD_0001 -> p_D2_0001 + c_C3D_0001 (photodissociation)
    kall(3327) = 5.40000000d-03* (2.50000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.50000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCD_0001 -> o_D2_0001 + c_C3D_0001 (photodissociation)
    kall(3328) = 5.40000000d-03* (2.50000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.50000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCD_0001 -> p_H2_0001 + l_C3D_0001 (photodissociation)
    kall(3329) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCD_0001 -> o_H2_0001 + l_C3D_0001 (photodissociation)
    kall(3330) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCD_0001 -> HD_0001 + l_C3H_0001 (photodissociation)
    kall(3331) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCH_0001 -> p_H2_0001 + l_C3D_0001 (photodissociation)
    kall(3332) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCH_0001 -> o_H2_0001 + l_C3D_0001 (photodissociation)
    kall(3333) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCH_0001 -> HD_0001 + l_C3H_0001 (photodissociation)
    kall(3334) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCD_0001 -> HD_0001 + l_C3D_0001 (photodissociation)
    kall(3335) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCD_0001 -> p_D2_0001 + l_C3H_0001 (photodissociation)
    kall(3336) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCD_0001 -> o_D2_0001 + l_C3H_0001 (photodissociation)
    kall(3337) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCH_0001 -> HD_0001 + l_C3D_0001 (photodissociation)
    kall(3338) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCH_0001 -> p_D2_0001 + l_C3H_0001 (photodissociation)
    kall(3339) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCH_0001 -> o_D2_0001 + l_C3H_0001 (photodissociation)
    kall(3340) = 5.40000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (8.35000000d-11/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCD_0001 -> p_D2_0001 + l_C3D_0001 (photodissociation)
    kall(3341) = 5.40000000d-03* (2.50000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.50000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCD_0001 -> o_D2_0001 + l_C3D_0001 (photodissociation)
    kall(3342) = 5.40000000d-03* (2.50000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (2.50000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCD_0001 -> H_0001 + c_C3HD_0001 (photodissociation)
    kall(3343) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCD_0001 -> D_0001 + c_C3H2_0001 (photodissociation)
    kall(3344) = 5.40000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCH_0001 -> H_0001 + c_C3HD_0001 (photodissociation)
    kall(3345) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCH_0001 -> D_0001 + c_C3H2_0001 (photodissociation)
    kall(3346) = 5.40000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCD_0001 -> H_0001 + c_C3D2_0001 (photodissociation)
    kall(3347) = 5.40000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCD_0001 -> D_0001 + c_C3HD_0001 (photodissociation)
    kall(3348) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCH_0001 -> H_0001 + c_C3D2_0001 (photodissociation)
    kall(3349) = 5.40000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCH_0001 -> D_0001 + c_C3HD_0001 (photodissociation)
    kall(3350) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCD_0001 -> D_0001 + c_C3D2_0001 (photodissociation)
    kall(3351) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCD_0001 -> H_0001 + l_C3HD_0001 (photodissociation)
    kall(3352) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2CCD_0001 -> D_0001 + l_C3H2_0001 (photodissociation)
    kall(3353) = 5.40000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCH_0001 -> H_0001 + l_C3HD_0001 (photodissociation)
    kall(3354) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCH_0001 -> D_0001 + l_C3H2_0001 (photodissociation)
    kall(3355) = 5.40000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCD_0001 -> H_0001 + l_C3D2_0001 (photodissociation)
    kall(3356) = 5.40000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDCCD_0001 -> D_0001 + l_C3HD_0001 (photodissociation)
    kall(3357) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCH_0001 -> H_0001 + l_C3D2_0001 (photodissociation)
    kall(3358) = 5.40000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (1.67000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCH_0001 -> D_0001 + l_C3HD_0001 (photodissociation)
    kall(3359) = 5.40000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (3.33000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CD2CCD_0001 -> D_0001 + l_C3D2_0001 (photodissociation)
    kall(3360) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.700000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HOOD_0001 -> OH_0001 + OD_0001 (photodissociation)
    kall(3361) = 5.40000000d-03* (4.15000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (4.15000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DOOD_0001 -> OD_0001 + OD_0001 (photodissociation)
    kall(3362) = 5.40000000d-03* (8.30000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (8.30000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !MgD_0001 -> D_0001 + Mg_0001 (photodissociation)
    kall(3363) = 5.40000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.450000 * variable_Av) + 4.70000000d-03* (5.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NaD_0001 -> D_0001 + Na_0001 (photodissociation)
    kall(3364) = 5.40000000d-03* (7.30000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.130000 * variable_Av) + 4.70000000d-03* (7.30000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !PD_0001 -> D_0001 + P_0001 (photodissociation)
    kall(3365) = 5.40000000d-03* (4.00000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.500000 * variable_Av) + 4.70000000d-03* (4.00000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DCP_0001 -> D_0001 + CP_0001 (photodissociation)
    kall(3366) = 5.40000000d-03* (5.48000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.000000 * variable_Av) + 4.70000000d-03* (5.48000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DPO_0001 -> D_0001 + PO_0001 (photodissociation)
    kall(3367) = 5.40000000d-03* (1.70000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-0.530000 * variable_Av) + 4.70000000d-03* (1.70000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DSS_0001 -> S_0001 + DS_0001 (photodissociation)
    kall(3368) = 5.40000000d-03* (3.30000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.400000 * variable_Av) + 4.70000000d-03* (3.30000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !PHD_0001 -> H_0001 + PD_0001 (photodissociation)
    kall(3369) = 5.40000000d-03* (1.05000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.520000 * variable_Av) + 4.70000000d-03* (1.05000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !PHD_0001 -> D_0001 + PH_0001 (photodissociation)
    kall(3370) = 5.40000000d-03* (1.05000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.520000 * variable_Av) + 4.70000000d-03* (1.05000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !PD2_0001 -> D_0001 + PD_0001 (photodissociation)
    kall(3371) = 5.40000000d-03* (2.11000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.520000 * variable_Av) + 4.70000000d-03* (2.11000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HSSD_0001 -> HS_0001 + DS_0001 (photodissociation)
    kall(3372) = 5.40000000d-03* (4.15000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (4.15000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DSSH_0001 -> HS_0001 + DS_0001 (photodissociation)
    kall(3373) = 5.40000000d-03* (4.15000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (4.15000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DSSD_0001 -> DS_0001 + DS_0001 (photodissociation)
    kall(3374) = 5.40000000d-03* (8.30000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.800000 * variable_Av) + 4.70000000d-03* (8.30000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HDCCN_0001 -> CN_0001 + CHD_0001 (photodissociation)
    kall(3375) = 5.40000000d-03* (1.56000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.950000 * variable_Av) + 4.70000000d-03* (1.56000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !D2CCN_0001 -> CN_0001 + CD2_0001 (photodissociation)
    kall(3376) = 5.40000000d-03* (1.56000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.950000 * variable_Av) + 4.70000000d-03* (1.56000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2PD_0001 -> PH_0001 + CHD_0001 (photodissociation)
    kall(3377) = 5.40000000d-03* (6.36000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.830000 * variable_Av) + 4.70000000d-03* (6.36000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH2PD_0001 -> PD_0001 + CH2_0001 (photodissociation)
    kall(3378) = 5.40000000d-03* (3.18000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.830000 * variable_Av) + 4.70000000d-03* (3.18000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDPD_0001 -> PH_0001 + CD2_0001 (photodissociation)
    kall(3379) = 5.40000000d-03* (3.18000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.830000 * variable_Av) + 4.70000000d-03* (3.18000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CHDPD_0001 -> PD_0001 + CHD_0001 (photodissociation)
    kall(3380) = 5.40000000d-03* (6.36000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.830000 * variable_Av) + 4.70000000d-03* (6.36000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DCNCC_0001 -> C2_0001 + DCN_0001 (photodissociation)
    kall(3381) = 5.40000000d-03* (9.54000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.830000 * variable_Av) + 4.70000000d-03* (9.54000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DCCNC_0001 -> CN_0001 + CCD_0001 (photodissociation)
    kall(3382) = 5.40000000d-03* (9.54000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.830000 * variable_Av) + 4.70000000d-03* (9.54000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DNCCC_0001 -> D_0001 + C3N_0001 (photodissociation)
    kall(3383) = 5.40000000d-03* (9.54000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-1.830000 * variable_Av) + 4.70000000d-03* (9.54000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !MgH2_0001 -> MgH_0001 + H_0001 (photodissociation)
    kall(3384) = 5.40000000d-03* (3.00000000d+03/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-0.000000 * variable_Av) + 4.70000000d-03* (3.00000000d+03/k_H2O_ph_0) &
        * F_cr * kph_factor

    !O3_0001 -> O2_0001 + O_0001 (photodissociation)
    kall(3385) = 5.40000000d-03* (1.83000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.250000 * variable_Av) + 4.70000000d-03* (1.83000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !C2H6_0001 -> C2H5_0001 + H_0001 (photodissociation)
    kall(3386) = 5.40000000d-03* (2.13000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.940000 * variable_Av) + 4.70000000d-03* (2.13000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !O2H_0001 -> H_0001 + O2_0001 (photodissociation)
    kall(3387) = 5.40000000d-03* (6.69000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.460000 * variable_Av) + 4.70000000d-03* (6.69000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HOOH_0001 -> OH_0001 + OH_0001 (photodissociation)
    kall(3388) = 5.40000000d-03* (8.08000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.610000 * variable_Av) + 4.70000000d-03* (8.08000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !HOOD_0001 -> OD_0001 + OH_0001 (photodissociation)
    kall(3389) = 5.40000000d-03* (8.08000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.610000 * variable_Av) + 4.70000000d-03* (8.08000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !DOOD_0001 -> OD_0001 + OD_0001 (photodissociation)
    kall(3390) = 5.40000000d-03* (8.08000000d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.610000 * variable_Av) + 4.70000000d-03* (8.08000000d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !NH2CHO_0001 -> CO_0001 + NH3_0001 (photodissociation)
    kall(3391) = 5.40000000d-03* (2.70000000d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.400000 * variable_Av) + 4.70000000d-03* (2.70000000d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3OH_0001 -> CH2OH_0001 + H_0001 (photodissociation)
    kall(3392) = 5.40000000d-03* (1.20714286d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (1.20714286d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3OH_0001 -> CH3O_0001 + H_0001 (photodissociation)
    kall(3393) = 5.40000000d-03* (2.41428571d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (2.41428571d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3OH_0001 -> CH3_0001 + OH_0001 (photodissociation)
    kall(3394) = 5.40000000d-03* (2.41428571d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (2.41428571d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3OD_0001 -> CH2OD_0001 + H_0001 (photodissociation)
    kall(3395) = 5.40000000d-03* (1.20714286d-09/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (1.20714286d-09/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3OD_0001 -> CH3O_0001 + D_0001 (photodissociation)
    kall(3396) = 5.40000000d-03* (2.41428571d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (2.41428571d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3OD_0001 -> CH3_0001 + OD_0001 (photodissociation)
    kall(3397) = 5.40000000d-03* (2.41428571d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (2.41428571d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3OH_0001 -> CHDOD_0001 + H_0001 (photodissociation)
    kall(3398) = 5.40000000d-03* (8.04761905d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (8.04761905d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3OH_0001 -> CH2DO_0001 + H_0001 (photodissociation)
    kall(3399) = 5.40000000d-03* (2.41428571d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (2.41428571d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    !CH3OH_0001 -> CH2D_0001 + OH_0001 (photodissociation)
    kall(3400) = 5.40000000d-03* (2.41428571d-10/k_H2O_ph_0) &
        * Fnot * Gnot * kph_factor * exp(-2.760000 * variable_Av) + 4.70000000d-03* (2.41428571d-10/k_H2O_ph_0) &
        * F_cr * kph_factor

    ! monolayer reaction  (dummy)
    kall(52129) = 0.0000005305164770 / (ngas*xdust)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_RATES


  end subroutine computeDustRates


end module kemimo_dust_rates
