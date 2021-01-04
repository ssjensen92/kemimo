# -*- coding: utf-8 -*-
import numpy as np
import ctypes
import numpy.ctypeslib as npctypes

# define aliases for complicated variable types
int_byref = ctypes.POINTER(ctypes.c_int)
dble_byref = ctypes.POINTER(ctypes.c_double)
array_1d_int = npctypes.ndpointer(dtype=np.int,ndim=1,flags='CONTIGUOUS')
array_1d_double = npctypes.ndpointer(dtype=np.float,ndim=1,flags='FORTRAN')
array_2d_double = npctypes.ndpointer(dtype=np.float,ndim=2,flags='CONTIGUOUS')
array_pointer = ctypes.POINTER(ctypes.c_int)

class Pykemimo(object):

    def __init__(self):
        fortran = npctypes.load_library('libkemimo.so', loader_path = '.')
        self.lib = fortran
        # constants
        self.spy = 365.0 * 24.0 * 3600.0
        self.pmass = 1.6726219e-24 # g
        self.gravity = 6.67259e-8  # cgs
        self.au2cm = 1.49597871e13 # AU -> cm
        self.kb = 1.38064852e-16   # erg/K

        ##BEGIN_ARRAYSIZE
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # NOTE: This block is auto-generated
        # WHEN: 2020-12-22 21:21:49
        # CHANGESET: xxxxxxx
        # BY: unknown@unknown

        #number of species
        self.nmols=1865
        #number of unique species
        self.nmolsu=1412
        #number of dust species
        self.nmols_dust=906
        #number of reactions (nlayer*dust+gas)
        self.nrea=52129
        #number of unique reactions (dust+gas)
        #number of dust-phase reactions
        self.nreadust=3401
        #number of gas-phase reactions
        self.nreagas=self.nrea-self.nreadust

        #monolayer thickness of each layer in model:
        self.layerThickness=4
        #idx for mantle and surface species:
        self.surface_start=957
        self.surface_end=1409
        self.mantle_start=1413
        self.mantle_end=1865
        self. CO_desorption_idx=2201

        # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        ##END_ARRAYSIZE

        ##BEGIN_SPECIESNAMES
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # NOTE: This block is auto-generated
        # WHEN: 2020-12-22 21:21:49
        # CHANGESET: xxxxxxx
        # BY: unknown@unknown

        self.speciesNames = ("GRAIN0_gas     ",
        "GRAIN-_gas     ",
        "E_gas          ",
        "H_gas          ",
        "H+_gas         ",
        "H-_gas         ",
        "D_gas          ",
        "D+_gas         ",
        "D-_gas         ",
        "p_H2_gas       ",
        "o_H2_gas       ",
        "p_H2+_gas      ",
        "o_H2+_gas      ",
        "p_D2_gas       ",
        "o_D2_gas       ",
        "p_D2+_gas      ",
        "o_D2+_gas      ",
        "p_D3+_gas      ",
        "o_D3+_gas      ",
        "m_D3+_gas      ",
        "p_D2H+_gas     ",
        "o_D2H+_gas     ",
        "p_H2D+_gas     ",
        "o_H2D+_gas     ",
        "HD_gas         ",
        "HD+_gas        ",
        "He_gas         ",
        "He+_gas        ",
        "O_gas          ",
        "O+_gas         ",
        "O-_gas         ",
        "O2_gas         ",
        "O2+_gas        ",
        "O3_gas         ",
        "OH_gas         ",
        "OH+_gas        ",
        "OH-_gas        ",
        "OD_gas         ",
        "OD+_gas        ",
        "OD-_gas        ",
        "H2O_gas        ",
        "H2O+_gas       ",
        "HDO_gas        ",
        "HDO+_gas       ",
        "D2O_gas        ",
        "D2O+_gas       ",
        "O2H_gas        ",
        "O2D_gas        ",
        "HO2+_gas       ",
        "DO2+_gas       ",
        "HOOH_gas       ",
        "HOOD_gas       ",
        "DOOH_gas       ",
        "DOOD_gas       ",
        "p_H3+_gas      ",
        "o_H3+_gas      ",
        "H3O+_gas       ",
        "H2DO+_gas      ",
        "HD2O+_gas      ",
        "D3O+_gas       ",
        "F_gas          ",
        "F+_gas         ",
        "Cl_gas         ",
        "Cl+_gas        ",
        "CF+_gas        ",
        "C_gas          ",
        "C+_gas         ",
        "C-_gas         ",
        "C2_gas         ",
        "C2+_gas        ",
        "C2N+_gas       ",
        "CO_gas         ",
        "CO+_gas        ",
        "CO2_gas        ",
        "CO2+_gas       ",
        "N_gas          ",
        "N+_gas         ",
        "N2_gas         ",
        "N2+_gas        ",
        "N2H+_gas       ",
        "N2D+_gas       ",
        "HCO_gas        ",
        "DCO_gas        ",
        "HCO+_gas       ",
        "DCO+_gas       ",
        "H2CO_gas       ",
        "HDCO_gas       ",
        "D2CO_gas       ",
        "CH2OH_gas      ",
        "CD2OD_gas      ",
        "CH2OD_gas      ",
        "CHDOH_gas      ",
        "CHDOD_gas      ",
        "CD2OH_gas      ",
        "CH3O_gas       ",
        "CHD2O_gas      ",
        "CH2DO_gas      ",
        "CD3O_gas       ",
        "CH3OH_gas      ",
        "CH3OD_gas      ",
        "CHD2OH_gas     ",
        "CHD2OD_gas     ",
        "CH2DOH_gas     ",
        "CH2DOD_gas     ",
        "CD3OD_gas      ",
        "CD3OH_gas      ",
        "CH_gas         ",
        "CD_gas         ",
        "CH2_gas        ",
        "CHD_gas        ",
        "CD2_gas        ",
        "CH3_gas        ",
        "CH2D_gas       ",
        "CHD2_gas       ",
        "CD3_gas        ",
        "CH4_gas        ",
        "CH3D_gas       ",
        "CH2D2_gas      ",
        "CHD3_gas       ",
        "CD4_gas        ",
        "CH+_gas        ",
        "CD+_gas        ",
        "CH2+_gas       ",
        "CHD+_gas       ",
        "CD2+_gas       ",
        "CH3+_gas       ",
        "CH2D+_gas      ",
        "CHD2+_gas      ",
        "CD3+_gas       ",
        "CH4+_gas       ",
        "CH3D+_gas      ",
        "CH2D2+_gas     ",
        "CHD3+_gas      ",
        "CD4+_gas       ",
        "CH5+_gas       ",
        "CH4D+_gas      ",
        "CH3D2+_gas     ",
        "CH2D3+_gas     ",
        "CHD4+_gas      ",
        "CD5+_gas       ",
        "HCOOH_gas      ",
        "HCOOD_gas      ",
        "DCOOH_gas      ",
        "DCOOD_gas      ",
        "HOCO_gas       ",
        "DOCO_gas       ",
        "NH_gas         ",
        "NH+_gas        ",
        "ND_gas         ",
        "ND+_gas        ",
        "NH2_gas        ",
        "NH2+_gas       ",
        "NHD_gas        ",
        "NHD+_gas       ",
        "ND2_gas        ",
        "ND2+_gas       ",
        "NH3_gas        ",
        "NH3+_gas       ",
        "NH2D_gas       ",
        "NH2D+_gas      ",
        "NHD2_gas       ",
        "NHD2+_gas      ",
        "ND3_gas        ",
        "ND3+_gas       ",
        "NH4+_gas       ",
        "NH3D+_gas      ",
        "NH2D2+_gas     ",
        "NHD3+_gas      ",
        "ND4+_gas       ",
        "Fe_gas         ",
        "Fe+_gas        ",
        "S_gas          ",
        "S+_gas         ",
        "S-_gas         ",
        "C10_gas        ",
        "C10+_gas       ",
        "C10-_gas       ",
        "C10H_gas       ",
        "C10H+_gas      ",
        "C10H-_gas      ",
        "C10H2_gas      ",
        "C10H2+_gas     ",
        "C10H3+_gas     ",
        "C10N_gas       ",
        "C10N+_gas      ",
        "C11_gas        ",
        "C11+_gas       ",
        "C2H+_gas       ",
        "C2D+_gas       ",
        "C2H2_gas       ",
        "C2H2+_gas      ",
        "C2HD_gas       ",
        "C2HD+_gas      ",
        "C2D2_gas       ",
        "C2D2+_gas      ",
        "C2DO+_gas      ",
        "C2H3_gas       ",
        "C2H3+_gas      ",
        "C2H2D_gas      ",
        "C2H2D+_gas     ",
        "C2HD2_gas      ",
        "C2HD2+_gas     ",
        "C2D3_gas       ",
        "C2D3+_gas      ",
        "C2H4_gas       ",
        "C2H4+_gas      ",
        "C2H3D_gas      ",
        "C2H3D+_gas     ",
        "C2H2D2_gas     ",
        "C2H2D2+_gas    ",
        "C2HD3_gas      ",
        "C2HD3+_gas     ",
        "C2D4_gas       ",
        "C2D4+_gas      ",
        "C2H4O+_gas     ",
        "C2H5_gas       ",
        "C2H5+_gas      ",
        "C2H5OH+_gas    ",
        "C2H5OH2+_gas   ",
        "C2H6_gas       ",
        "C2H6+_gas      ",
        "C2H6CO+_gas    ",
        "C2H7+_gas      ",
        "C2HO+_gas      ",
        "C2N2+_gas      ",
        "C2O+_gas       ",
        "C2S+_gas       ",
        "C3_gas         ",
        "C3+_gas        ",
        "C3-_gas        ",
        "C3H+_gas       ",
        "C3D+_gas       ",
        "C3H3N+_gas     ",
        "C3H3NH+_gas    ",
        "C3H4_gas       ",
        "C3H4+_gas      ",
        "C3H5+_gas      ",
        "C3H6OH+_gas    ",
        "C3N_gas        ",
        "C3N+_gas       ",
        "C3N-_gas       ",
        "C3O_gas        ",
        "C3O+_gas       ",
        "C3S_gas        ",
        "C3S+_gas       ",
        "C4_gas         ",
        "C4+_gas        ",
        "C4-_gas        ",
        "C4H_gas        ",
        "C4H+_gas       ",
        "C4H-_gas       ",
        "C4D_gas        ",
        "C4D+_gas       ",
        "C4D-_gas       ",
        "C4H2_gas       ",
        "C4H2+_gas      ",
        "C4HD_gas       ",
        "C4HD+_gas      ",
        "C4D2_gas       ",
        "C4D2+_gas      ",
        "C4H3_gas       ",
        "C4H3+_gas      ",
        "C4H4+_gas      ",
        "C4H5+_gas      ",
        "C4H7+_gas      ",
        "C4N_gas        ",
        "C4N+_gas       ",
        "C4S_gas        ",
        "C4S+_gas       ",
        "C5_gas         ",
        "C5+_gas        ",
        "C5-_gas        ",
        "C5H_gas        ",
        "C5H+_gas       ",
        "C5H-_gas       ",
        "C5D_gas        ",
        "C5D+_gas       ",
        "C5D-_gas       ",
        "C5H2_gas       ",
        "C5H2+_gas      ",
        "C5H3_gas       ",
        "C5H3+_gas      ",
        "C5H3N+_gas     ",
        "C5H4_gas       ",
        "C5H4+_gas      ",
        "C5H4N+_gas     ",
        "C5H5+_gas      ",
        "C5N_gas        ",
        "C5N+_gas       ",
        "C5O_gas        ",
        "C6_gas         ",
        "C6+_gas        ",
        "C6-_gas        ",
        "C6H_gas        ",
        "C6H+_gas       ",
        "C6H-_gas       ",
        "C6H2_gas       ",
        "C6H2+_gas      ",
        "C6H3_gas       ",
        "C6H3+_gas      ",
        "C6H4_gas       ",
        "C6H4+_gas      ",
        "C6H5+_gas      ",
        "C6H6_gas       ",
        "C6H7+_gas      ",
        "C6N_gas        ",
        "C6N+_gas       ",
        "C7_gas         ",
        "C7+_gas        ",
        "C7-_gas        ",
        "C7H_gas        ",
        "C7H+_gas       ",
        "C7H-_gas       ",
        "C7H2_gas       ",
        "C7H2+_gas      ",
        "C7H2N+_gas     ",
        "C7H3_gas       ",
        "C7H3+_gas      ",
        "C7H4_gas       ",
        "C7H4+_gas      ",
        "C7H5+_gas      ",
        "C7N_gas        ",
        "C7N+_gas       ",
        "C7O_gas        ",
        "C8_gas         ",
        "C8+_gas        ",
        "C8-_gas        ",
        "C8H_gas        ",
        "C8H+_gas       ",
        "C8H-_gas       ",
        "C8H2_gas       ",
        "C8H2+_gas      ",
        "C8H3_gas       ",
        "C8H3+_gas      ",
        "C8H4_gas       ",
        "C8H4+_gas      ",
        "C8H4N+_gas     ",
        "C8H5+_gas      ",
        "C8N_gas        ",
        "C8N+_gas       ",
        "C9_gas         ",
        "C9+_gas        ",
        "C9-_gas        ",
        "C9H_gas        ",
        "C9H+_gas       ",
        "C9H-_gas       ",
        "C9H2_gas       ",
        "C9H2+_gas      ",
        "C9H2N+_gas     ",
        "C9H3_gas       ",
        "C9H3+_gas      ",
        "C9H3N+_gas     ",
        "C9H4_gas       ",
        "C9H4+_gas      ",
        "C9H5+_gas      ",
        "C9HN+_gas      ",
        "C9N_gas        ",
        "C9N+_gas       ",
        "C9O_gas        ",
        "CCH_gas        ",
        "CCD_gas        ",
        "CCN_gas        ",
        "CCO_gas        ",
        "CCP_gas        ",
        "CCP+_gas       ",
        "CCS_gas        ",
        "CD2ND2_gas     ",
        "CD2ND2+_gas    ",
        "CD2NH2_gas     ",
        "CD2NH2+_gas    ",
        "CD2NHD_gas     ",
        "CD2NHD+_gas    ",
        "CD3CN_gas      ",
        "CD3CN+_gas     ",
        "CD3CO+_gas     ",
        "CD3O2+_gas     ",
        "CD3OH+_gas     ",
        "CD3OD+_gas     ",
        "CH2CCH_gas     ",
        "CH2CCD_gas     ",
        "CD2CCH_gas     ",
        "CD2CCD_gas     ",
        "CHDCCH_gas     ",
        "CHDCCD_gas     ",
        "CH2CHC2H_gas   ",
        "CH2CHCHCH2_gas ",
        "CH2CHCN_gas    ",
        "CH2CN+_gas     ",
        "CHDCN+_gas     ",
        "CD2CN+_gas     ",
        "CH2DCN_gas     ",
        "CH2DCN+_gas    ",
        "CH2DCO+_gas    ",
        "CH2DO2+_gas    ",
        "CH2DOD+_gas    ",
        "CH2ND2_gas     ",
        "CH2ND2+_gas    ",
        "CH2NH_gas      ",
        "CH2ND_gas      ",
        "CHDNH_gas      ",
        "CHDND_gas      ",
        "CD2NH_gas      ",
        "CD2ND_gas      ",
        "CH2NH2_gas     ",
        "CH2NH2+_gas    ",
        "CH2NHD_gas     ",
        "CH2NHD+_gas    ",
        "CH3C3N_gas     ",
        "CH3C4H_gas     ",
        "CH3C5N_gas     ",
        "CH3C6H_gas     ",
        "CH3C7N_gas     ",
        "CH3CCH_gas     ",
        "CH3CH2OH_gas   ",
        "CH3CHCH2_gas   ",
        "CH3CHO_gas     ",
        "CH3CHOH+_gas   ",
        "CH3CN_gas      ",
        "CH3CN+_gas     ",
        "CH3CNH+_gas    ",
        "CH3CO+_gas     ",
        "CH3COCH3_gas   ",
        "CH3NH2_gas     ",
        "CH3NH2+_gas    ",
        "CH3NH3+_gas    ",
        "CH3O2+_gas     ",
        "CH3OCH2_gas    ",
        "CH3OCH3_gas    ",
        "CH3OCH3+_gas   ",
        "CH3OCH4+_gas   ",
        "CH3OH+_gas     ",
        "CH2DOH+_gas    ",
        "CH3OD+_gas     ",
        "CH3OH2+_gas    ",
        "CHD2CN_gas     ",
        "CHD2CN+_gas    ",
        "CHD2CO+_gas    ",
        "CHD2O2+_gas    ",
        "CHD2OH+_gas    ",
        "CHD2OD+_gas    ",
        "CHDNH2_gas     ",
        "CHDNH2+_gas    ",
        "CHDNHD_gas     ",
        "CHDNHD+_gas    ",
        "CHDND2_gas     ",
        "CHDND2+_gas    ",
        "CN_gas         ",
        "CN+_gas        ",
        "CN-_gas        ",
        "CNC+_gas       ",
        "COOCH4+_gas    ",
        "CS_gas         ",
        "CS+_gas        ",
        "D2C3O+_gas     ",
        "DC2NCH+_gas    ",
        "DC2NCD+_gas    ",
        "DCND+_gas      ",
        "DCNH+_gas      ",
        "DCS_gas        ",
        "DCS+_gas       ",
        "DN2O+_gas      ",
        "DNC_gas        ",
        "DNC+_gas       ",
        "DNCCC_gas      ",
        "DNCO_gas       ",
        "DNCO+_gas      ",
        "DNO_gas        ",
        "DNO+_gas       ",
        "DNS+_gas       ",
        "DOC+_gas       ",
        "DOCO+_gas      ",
        "DOCS+_gas      ",
        "HS_gas         ",
        "HS+_gas        ",
        "DS_gas         ",
        "DS+_gas        ",
        "FeH_gas        ",
        "FeD_gas        ",
        "H2C10N+_gas    ",
        "H2C3O+_gas     ",
        "H2C4N+_gas     ",
        "H2C5N+_gas     ",
        "H2C6N+_gas     ",
        "H2C8N+_gas     ",
        "H2CCN_gas      ",
        "HDCCN_gas      ",
        "D2CCN_gas      ",
        "H2CCO_gas      ",
        "H2CCO+_gas     ",
        "HDCCO_gas      ",
        "HDCCO+_gas     ",
        "D2CCO_gas      ",
        "D2CCO+_gas     ",
        "H2CN_gas       ",
        "H2CN+_gas      ",
        "HDCN_gas       ",
        "D2CN_gas       ",
        "H2CO+_gas      ",
        "HDCO+_gas      ",
        "D2CO+_gas      ",
        "H2COH+_gas     ",
        "H2COD+_gas     ",
        "HDCOH+_gas     ",
        "HDCOD+_gas     ",
        "D2COH+_gas     ",
        "D2COD+_gas     ",
        "H2CS_gas       ",
        "H2CS+_gas      ",
        "HDCS_gas       ",
        "HDCS+_gas      ",
        "D2CS_gas       ",
        "D2CS+_gas      ",
        "H2NC+_gas      ",
        "HDNC+_gas      ",
        "D2NC+_gas      ",
        "H2NO+_gas      ",
        "HDNO+_gas      ",
        "D2NO+_gas      ",
        "H2S_gas        ",
        "H2S+_gas       ",
        "HDS_gas        ",
        "HDS+_gas       ",
        "D2S_gas        ",
        "D2S+_gas       ",
        "H2S2+_gas      ",
        "HDS2+_gas      ",
        "D2S2+_gas      ",
        "H3C3O+_gas     ",
        "H3C4N+_gas     ",
        "H3C4NH+_gas    ",
        "H3C6NH+_gas    ",
        "H3C7N+_gas     ",
        "H3CS+_gas      ",
        "H2DCS+_gas     ",
        "HD2CS+_gas     ",
        "D3CS+_gas      ",
        "H3S+_gas       ",
        "H2DS+_gas      ",
        "HD2S+_gas      ",
        "D3S+_gas       ",
        "H3S2+_gas      ",
        "H2DS2+_gas     ",
        "HD2S2+_gas     ",
        "D3S2+_gas      ",
        "H5C2O2+_gas    ",
        "HC10N+_gas     ",
        "HC2N_gas       ",
        "HC2N+_gas      ",
        "DC2N_gas       ",
        "DC2N+_gas      ",
        "HC2NCD+_gas    ",
        "HC2NCH+_gas    ",
        "HC2O_gas       ",
        "DC2O_gas       ",
        "HC2S+_gas      ",
        "DC2S+_gas      ",
        "HC3N_gas       ",
        "HC3N+_gas      ",
        "DC3N_gas       ",
        "DC3N+_gas      ",
        "HC3NH+_gas     ",
        "HC3ND+_gas     ",
        "DC3NH+_gas     ",
        "DC3ND+_gas     ",
        "HC3O+_gas      ",
        "DC3O+_gas      ",
        "HC3S+_gas      ",
        "DC3S+_gas      ",
        "HC4N_gas       ",
        "HC4N+_gas      ",
        "DC4N_gas       ",
        "DC4N+_gas      ",
        "HC4O+_gas      ",
        "DC4O+_gas      ",
        "HC4S+_gas      ",
        "DC4S+_gas      ",
        "HC5N_gas       ",
        "HC5N+_gas      ",
        "HC5O+_gas      ",
        "HC6N_gas       ",
        "HC6N+_gas      ",
        "HC7N_gas       ",
        "HC7N+_gas      ",
        "HC7O+_gas      ",
        "HC8N_gas       ",
        "HC8N+_gas      ",
        "HC9N_gas       ",
        "HC9O+_gas      ",
        "HCCNC_gas      ",
        "DCCNC_gas      ",
        "HCN_gas        ",
        "HCN+_gas       ",
        "DCN_gas        ",
        "DCN+_gas       ",
        "HCNCC_gas      ",
        "DCNCC_gas      ",
        "HCNH+_gas      ",
        "HCND+_gas      ",
        "HCOOCH3_gas    ",
        "HCOOH+_gas     ",
        "HCOOD+_gas     ",
        "DCOOH+_gas     ",
        "DCOOD+_gas     ",
        "HCS_gas        ",
        "HCS+_gas       ",
        "HDC3O+_gas     ",
        "HN2O+_gas      ",
        "HNC_gas        ",
        "HNC+_gas       ",
        "HNCCC_gas      ",
        "HNCO_gas       ",
        "HNCO+_gas      ",
        "HNO_gas        ",
        "HNO+_gas       ",
        "HNS+_gas       ",
        "HOC+_gas       ",
        "HOCO+_gas      ",
        "HOCS+_gas      ",
        "HS2+_gas       ",
        "DS2+_gas       ",
        "HSO+_gas       ",
        "DSO+_gas       ",
        "HSO2+_gas      ",
        "DSO2+_gas      ",
        "HSS_gas        ",
        "DSS_gas        ",
        "HSSH_gas       ",
        "HSSD_gas       ",
        "DSSH_gas       ",
        "DSSD_gas       ",
        "N2O_gas        ",
        "NC4N_gas       ",
        "NC6N_gas       ",
        "NC8N_gas       ",
        "NCO+_gas       ",
        "NH2CH2O+_gas   ",
        "NH2CHO_gas     ",
        "NH2CDO_gas     ",
        "NHDCHO_gas     ",
        "NHDCDO_gas     ",
        "ND2CHO_gas     ",
        "ND2CDO_gas     ",
        "NH2CN_gas      ",
        "NHDCN_gas      ",
        "ND2CN_gas      ",
        "NH2CND+_gas    ",
        "NHDCND+_gas    ",
        "ND2CND+_gas    ",
        "NH2CNH+_gas    ",
        "NHDCNH+_gas    ",
        "ND2CNH+_gas    ",
        "NO_gas         ",
        "NO+_gas        ",
        "NO2_gas        ",
        "NO2+_gas       ",
        "NS_gas         ",
        "NS+_gas        ",
        "OCN_gas        ",
        "OCS_gas        ",
        "OCS+_gas       ",
        "S2_gas         ",
        "S2+_gas        ",
        "SO_gas         ",
        "SO+_gas        ",
        "SO2_gas        ",
        "SO2+_gas       ",
        "CH3OCHO_gas    ",
        "Si_gas         ",
        "Si+_gas        ",
        "SiO_gas        ",
        "SiO+_gas       ",
        "SiO2_gas       ",
        "SiS_gas        ",
        "SiS+_gas       ",
        "SiN_gas        ",
        "SiN+_gas       ",
        "SiC_gas        ",
        "SiC+_gas       ",
        "SiC2_gas       ",
        "SiC2+_gas      ",
        "SiH_gas        ",
        "SiD_gas        ",
        "SiH2_gas       ",
        "SiH2+_gas      ",
        "SiHD_gas       ",
        "SiHD+_gas      ",
        "SiD2_gas       ",
        "SiD2+_gas      ",
        "SiH3_gas       ",
        "SiH4_gas       ",
        "P_gas          ",
        "P+_gas         ",
        "PH_gas         ",
        "PH+_gas        ",
        "PD_gas         ",
        "PD+_gas        ",
        "PO_gas         ",
        "PO+_gas        ",
        "PH2_gas        ",
        "PH2+_gas       ",
        "PHD_gas        ",
        "PHD+_gas       ",
        "PD2_gas        ",
        "PD2+_gas       ",
        "C3P_gas        ",
        "C4P_gas        ",
        "C4P+_gas       ",
        "CH2PH_gas      ",
        "CH2PD_gas      ",
        "CHDPH_gas      ",
        "CHDPD_gas      ",
        "CD2PH_gas      ",
        "CD2PD_gas      ",
        "CH2Si+_gas     ",
        "CHDSi+_gas     ",
        "CD2Si+_gas     ",
        "CHSi+_gas      ",
        "CDSi+_gas      ",
        "CP_gas         ",
        "CP+_gas        ",
        "D2PO+_gas      ",
        "D3SiO+_gas     ",
        "DCP_gas        ",
        "DCP+_gas       ",
        "DF_gas         ",
        "DF+_gas        ",
        "DPO_gas        ",
        "DPO+_gas       ",
        "DSiNH+_gas     ",
        "DSiND+_gas     ",
        "H2CSiCH_gas    ",
        "H2CSiCD_gas    ",
        "HDCSiCD_gas    ",
        "HDCSiCH_gas    ",
        "D2CSiCH_gas    ",
        "D2CSiCD_gas    ",
        "H2DSiO+_gas    ",
        "H2F+_gas       ",
        "HDF+_gas       ",
        "D2F+_gas       ",
        "H2PO+_gas      ",
        "H2SiO_gas      ",
        "H2SiO+_gas     ",
        "HDSiO_gas      ",
        "HDSiO+_gas     ",
        "D2SiO_gas      ",
        "D2SiO+_gas     ",
        "H3SiO+_gas     ",
        "HCCP_gas       ",
        "DCCP_gas       ",
        "HCCSi_gas      ",
        "DCCSi_gas      ",
        "HCP_gas        ",
        "HCP+_gas       ",
        "HCSi_gas       ",
        "DCSi_gas       ",
        "HD2SiO+_gas    ",
        "HDPO+_gas      ",
        "HF_gas         ",
        "HF+_gas        ",
        "HNSi_gas       ",
        "HNSi+_gas      ",
        "DNSi_gas       ",
        "DNSi+_gas      ",
        "HPN+_gas       ",
        "DPN+_gas       ",
        "HPO_gas        ",
        "HPO+_gas       ",
        "HSiNH+_gas     ",
        "HSiND+_gas     ",
        "HSiO+_gas      ",
        "DSiO+_gas      ",
        "HSiO2+_gas     ",
        "DSiO2+_gas     ",
        "HSiS+_gas      ",
        "DSiS+_gas      ",
        "HeH+_gas       ",
        "HeD+_gas       ",
        "PC2H+_gas      ",
        "PC2D+_gas      ",
        "PC2H2+_gas     ",
        "PC2HD+_gas     ",
        "PC2D2+_gas     ",
        "PC2H3+_gas     ",
        "PC2H2D+_gas    ",
        "PC2HD2+_gas    ",
        "PC2D3+_gas     ",
        "PC2H4+_gas     ",
        "PC3H+_gas      ",
        "PC3D+_gas      ",
        "PC4H+_gas      ",
        "PC4D+_gas      ",
        "PC4H2+_gas     ",
        "PCH2+_gas      ",
        "PCHD+_gas      ",
        "PCD2+_gas      ",
        "PCH3+_gas      ",
        "PCH2D+_gas     ",
        "PCHD2+_gas     ",
        "PCD3+_gas      ",
        "PCH4+_gas      ",
        "PCH3D+_gas     ",
        "PCH2D2+_gas    ",
        "PCHD3+_gas     ",
        "PCD4+_gas      ",
        "PH3+_gas       ",
        "PH2D+_gas      ",
        "PHD2+_gas      ",
        "PD3+_gas       ",
        "PN_gas         ",
        "PN+_gas        ",
        "PNH2+_gas      ",
        "PNHD+_gas      ",
        "PND2+_gas      ",
        "PNH3+_gas      ",
        "PNHD2+_gas     ",
        "PNH2D+_gas     ",
        "PND3+_gas      ",
        "SiC2CH3_gas    ",
        "SiC2H+_gas     ",
        "SiC2D+_gas     ",
        "SiC2H2+_gas    ",
        "SiC2HD+_gas    ",
        "SiC2D2+_gas    ",
        "SiC2H3+_gas    ",
        "SiC2H2D+_gas   ",
        "SiC2HD2+_gas   ",
        "SiC2D3+_gas    ",
        "SiC3D_gas      ",
        "SiC3D+_gas     ",
        "SiC3H_gas      ",
        "SiC3H+_gas     ",
        "SiC3H2+_gas    ",
        "SiC3HD+_gas    ",
        "SiC3D2+_gas    ",
        "SiC3H5_gas     ",
        "SiC4_gas       ",
        "SiC4+_gas      ",
        "SiC4D_gas      ",
        "SiC4D+_gas     ",
        "SiC4H_gas      ",
        "SiC4H+_gas     ",
        "SiC6H_gas      ",
        "SiC8H_gas      ",
        "SiCD3_gas      ",
        "SiCD3+_gas     ",
        "SiCH2_gas      ",
        "SiCHD_gas      ",
        "SiCD2_gas      ",
        "SiCH2D_gas     ",
        "SiCH2D+_gas    ",
        "SiCH3_gas      ",
        "SiCH3+_gas     ",
        "SiCH4+_gas     ",
        "SiCH3D+_gas    ",
        "SiCH2D2+_gas   ",
        "SiCHD3+_gas    ",
        "SiCD4+_gas     ",
        "SiCHD2_gas     ",
        "SiCHD2+_gas    ",
        "SiD3_gas       ",
        "SiD3+_gas      ",
        "SiF+_gas       ",
        "SiH+_gas       ",
        "SiD+_gas       ",
        "SiH2D_gas      ",
        "SiH2D+_gas     ",
        "SiH2D2_gas     ",
        "SiH2D2+_gas    ",
        "SiH3+_gas      ",
        "SiH3D_gas      ",
        "SiH3D+_gas     ",
        "SiH4+_gas      ",
        "SiD4_gas       ",
        "SiD4+_gas      ",
        "SiH5+_gas      ",
        "SiH4D+_gas     ",
        "SiH3D2+_gas    ",
        "SiH2D3+_gas    ",
        "SiHD4+_gas     ",
        "SiD5+_gas      ",
        "SiHD2_gas      ",
        "SiHD2+_gas     ",
        "SiHD3_gas      ",
        "SiHD3+_gas     ",
        "SiNC_gas       ",
        "SiNC+_gas      ",
        "SiNCH+_gas     ",
        "SiNCD+_gas     ",
        "c_DCCHSi_gas   ",
        "c_DCCDSi_gas   ",
        "c_HCCHSi_gas   ",
        "c_HCCDSi_gas   ",
        "c_SiC2_gas     ",
        "l_SiC3_gas     ",
        "l_SiC3+_gas    ",
        "l_C3H_gas      ",
        "l_C3D_gas      ",
        "c_C3H_gas      ",
        "c_C3D_gas      ",
        "l_C3H2_gas     ",
        "l_C3HD_gas     ",
        "l_C3D2_gas     ",
        "c_C3H2_gas     ",
        "c_C3HD_gas     ",
        "c_C3D2_gas     ",
        "l_C3H2+_gas    ",
        "l_C3HD+_gas    ",
        "l_C3D2+_gas    ",
        "c_C3H2+_gas    ",
        "c_C3HD+_gas    ",
        "c_C3D2+_gas    ",
        "l_C3H3+_gas    ",
        "l_C3H2D+_gas   ",
        "l_C3HD2+_gas   ",
        "l_C3D3+_gas    ",
        "c_C3H3+_gas    ",
        "c_C3H2D+_gas   ",
        "c_C3HD2+_gas   ",
        "c_C3D3+_gas    ",
        "Mg_gas         ",
        "Mg+_gas        ",
        "MgH_gas        ",
        "MgD_gas        ",
        "MgH2_gas       ",
        "MgHD_gas       ",
        "MgD2_gas       ",
        "HCl_gas        ",
        "HCl+_gas       ",
        "DCl_gas        ",
        "DCl+_gas       ",
        "H2Cl_gas       ",
        "H2Cl+_gas      ",
        "HDCl_gas       ",
        "HDCl+_gas      ",
        "D2Cl_gas       ",
        "D2Cl+_gas      ",
        "CCl_gas        ",
        "CCl+_gas       ",
        "ClO_gas        ",
        "ClO+_gas       ",
        "H2CCl+_gas     ",
        "HDCCl+_gas     ",
        "D2CCl+_gas     ",
        "Na_gas         ",
        "Na+_gas        ",
        "NaH_gas        ",
        "NaD_gas        ",
        "NaOH_gas       ",
        "NaOD_gas       ",
        "NaH2+_gas      ",
        "NaHD+_gas      ",
        "NaD2+_gas      ",
        "NaH2O+_gas     ",
        "NaHDO+_gas     ",
        "NaD2O+_gas     ",
        "H_0001         ",
        "D_0001         ",
        "p_H2_0001      ",
        "p_D2_0001      ",
        "o_H2_0001      ",
        "o_D2_0001      ",
        "HD_0001        ",
        "He_0001        ",
        "O_0001         ",
        "O2_0001        ",
        "O3_0001        ",
        "OH_0001        ",
        "OD_0001        ",
        "H2O_0001       ",
        "HDO_0001       ",
        "D2O_0001       ",
        "O2H_0001       ",
        "O2D_0001       ",
        "HOOH_0001      ",
        "HOOD_0001      ",
        "DOOD_0001      ",
        "Fe_0001        ",
        "FeH_0001       ",
        "N_0001         ",
        "S_0001         ",
        "C_0001         ",
        "C2_0001        ",
        "CO_0001        ",
        "HCO_0001       ",
        "DCO_0001       ",
        "H2CO_0001      ",
        "HDCO_0001      ",
        "D2CO_0001      ",
        "CH2OH_0001     ",
        "CD2OD_0001     ",
        "CH2OD_0001     ",
        "CHDOH_0001     ",
        "CHDOD_0001     ",
        "CD2OH_0001     ",
        "CH3O_0001      ",
        "CHD2O_0001     ",
        "CH2DO_0001     ",
        "CD3O_0001      ",
        "CH3OH_0001     ",
        "CH3OD_0001     ",
        "CHD2OH_0001    ",
        "CHD2OD_0001    ",
        "CH2DOH_0001    ",
        "CH2DOD_0001    ",
        "CD3OD_0001     ",
        "CD3OH_0001     ",
        "CH_0001        ",
        "CD_0001        ",
        "CH2_0001       ",
        "CHD_0001       ",
        "CD2_0001       ",
        "CH3_0001       ",
        "CH2D_0001      ",
        "CHD2_0001      ",
        "CD3_0001       ",
        "CH4_0001       ",
        "CH3D_0001      ",
        "CH2D2_0001     ",
        "CHD3_0001      ",
        "CD4_0001       ",
        "CO2_0001       ",
        "HCOOH_0001     ",
        "HCOOD_0001     ",
        "DCOOH_0001     ",
        "DCOOD_0001     ",
        "HOCO_0001      ",
        "DOCO_0001      ",
        "NH_0001        ",
        "ND_0001        ",
        "NH2_0001       ",
        "NHD_0001       ",
        "ND2_0001       ",
        "NH3_0001       ",
        "NH2D_0001      ",
        "NHD2_0001      ",
        "ND3_0001       ",
        "C10_0001       ",
        "C10H_0001      ",
        "C10H2_0001     ",
        "C10N_0001      ",
        "C11_0001       ",
        "C2H2_0001      ",
        "C2HD_0001      ",
        "C2D2_0001      ",
        "C2H3_0001      ",
        "C2H2D_0001     ",
        "C2HD2_0001     ",
        "C2D3_0001      ",
        "C2H4_0001      ",
        "C2H3D_0001     ",
        "C2H2D2_0001    ",
        "C2HD3_0001     ",
        "C2D4_0001      ",
        "C2H5_0001      ",
        "C2H6_0001      ",
        "C3_0001        ",
        "C3N_0001       ",
        "C3O_0001       ",
        "C3S_0001       ",
        "C4_0001        ",
        "C4H_0001       ",
        "C4D_0001       ",
        "C4H2_0001      ",
        "C4HD_0001      ",
        "C4D2_0001      ",
        "C4H3_0001      ",
        "C4N_0001       ",
        "C4S_0001       ",
        "C5_0001        ",
        "C5H_0001       ",
        "C5D_0001       ",
        "C5H2_0001      ",
        "C5H3_0001      ",
        "C5H4_0001      ",
        "C5N_0001       ",
        "C5O_0001       ",
        "C6_0001        ",
        "C6H_0001       ",
        "C6H2_0001      ",
        "C6H3_0001      ",
        "C6H4_0001      ",
        "C6H6_0001      ",
        "C6N_0001       ",
        "C7_0001        ",
        "C7H_0001       ",
        "C7H2_0001      ",
        "C7H3_0001      ",
        "C7H4_0001      ",
        "C7N_0001       ",
        "C7O_0001       ",
        "C8_0001        ",
        "C8H_0001       ",
        "C8H2_0001      ",
        "C8H3_0001      ",
        "C8H4_0001      ",
        "C8N_0001       ",
        "C9_0001        ",
        "C9H_0001       ",
        "C9H2_0001      ",
        "C9H3_0001      ",
        "C9H4_0001      ",
        "C9N_0001       ",
        "C9O_0001       ",
        "CCH_0001       ",
        "CCD_0001       ",
        "CCN_0001       ",
        "CCO_0001       ",
        "CCS_0001       ",
        "CD3CN_0001     ",
        "CH2CCH_0001    ",
        "CH2CCD_0001    ",
        "CHDCCH_0001    ",
        "CHDCCD_0001    ",
        "CD2CCH_0001    ",
        "CD2CCD_0001    ",
        "CH2CHC2H_0001  ",
        "CH2CHCHCH2_0001",
        "CH2CHCN_0001   ",
        "CH2NH_0001     ",
        "CHDNH_0001     ",
        "CHDND_0001     ",
        "CH2ND_0001     ",
        "CD2NH_0001     ",
        "CD2ND_0001     ",
        "CH2NH2_0001    ",
        "CH2NHD_0001    ",
        "CHDNH2_0001    ",
        "CHDNHD_0001    ",
        "CHDND2_0001    ",
        "CH2ND2_0001    ",
        "CD2NH2_0001    ",
        "CD2NHD_0001    ",
        "CD2ND2_0001    ",
        "CH3C3N_0001    ",
        "CH3C4H_0001    ",
        "CH3C5N_0001    ",
        "CH3C6H_0001    ",
        "CH3C7N_0001    ",
        "CH3CCH_0001    ",
        "CH3CH2OH_0001  ",
        "CH3CHCH2_0001  ",
        "CH3CHO_0001    ",
        "CH3CN_0001     ",
        "CH2DCN_0001    ",
        "CHD2CN_0001    ",
        "CH3COCH3_0001  ",
        "CH3NH2_0001    ",
        "CH3OCH2_0001   ",
        "CH3OCH3_0001   ",
        "CN_0001        ",
        "CS_0001        ",
        "H2CCN_0001     ",
        "HDCCN_0001     ",
        "D2CCN_0001     ",
        "H2CCO_0001     ",
        "HDCCO_0001     ",
        "D2CCO_0001     ",
        "H2CN_0001      ",
        "HDCN_0001      ",
        "D2CN_0001      ",
        "H2CS_0001      ",
        "HDCS_0001      ",
        "D2CS_0001      ",
        "H2S_0001       ",
        "HDS_0001       ",
        "D2S_0001       ",
        "HC2O_0001      ",
        "DC2O_0001      ",
        "HC3N_0001      ",
        "DC3N_0001      ",
        "HC4N_0001      ",
        "DC4N_0001      ",
        "HC5N_0001      ",
        "HC6N_0001      ",
        "HC7N_0001      ",
        "HC8N_0001      ",
        "HC9N_0001      ",
        "HCCNC_0001     ",
        "DCCNC_0001     ",
        "HCN_0001       ",
        "DCN_0001       ",
        "HCNCC_0001     ",
        "DCNCC_0001     ",
        "HCOOCH3_0001   ",
        "HCS_0001       ",
        "DCS_0001       ",
        "HNC_0001       ",
        "DNC_0001       ",
        "HNCCC_0001     ",
        "DNCCC_0001     ",
        "HNCO_0001      ",
        "DNCO_0001      ",
        "HNO_0001       ",
        "DNO_0001       ",
        "HS_0001        ",
        "DS_0001        ",
        "N2_0001        ",
        "N2O_0001       ",
        "NC4N_0001      ",
        "NC6N_0001      ",
        "NC8N_0001      ",
        "NH2CHO_0001    ",
        "NH2CDO_0001    ",
        "NHDCHO_0001    ",
        "NHDCDO_0001    ",
        "ND2CHO_0001    ",
        "ND2CDO_0001    ",
        "NH2CN_0001     ",
        "NHDCN_0001     ",
        "ND2CN_0001     ",
        "HSO_0001       ",
        "DSO_0001       ",
        "HSS_0001       ",
        "DSS_0001       ",
        "HSSH_0001      ",
        "HSSD_0001      ",
        "DSSH_0001      ",
        "DSSD_0001      ",
        "NO_0001        ",
        "NO2_0001       ",
        "NS_0001        ",
        "OCN_0001       ",
        "OCS_0001       ",
        "S2_0001        ",
        "SO_0001        ",
        "SO2_0001       ",
        "CH3OCHO_0001   ",
        "Si_0001        ",
        "SiS_0001       ",
        "SiN_0001       ",
        "SiC_0001       ",
        "SiH_0001       ",
        "SiH2_0001      ",
        "SiH3_0001      ",
        "SiH4_0001      ",
        "SiC2CH3_0001   ",
        "SiC3H_0001     ",
        "SiC3H5_0001    ",
        "SiC4_0001      ",
        "SiC4H_0001     ",
        "SiC6H_0001     ",
        "SiC8H_0001     ",
        "c_HCCHSi_0001  ",
        "c_SiC2_0001    ",
        "l_C3H_0001     ",
        "l_C3D_0001     ",
        "c_C3H_0001     ",
        "c_C3D_0001     ",
        "l_SiC3_0001    ",
        "l_C3H2_0001    ",
        "l_C3HD_0001    ",
        "l_C3D2_0001    ",
        "c_C3H2_0001    ",
        "c_C3HD_0001    ",
        "c_C3D2_0001    ",
        "Mg_0001        ",
        "MgH_0001       ",
        "MgH2_0001      ",
        "Na_0001        ",
        "NaH_0001       ",
        "F_0001         ",
        "HF_0001        ",
        "DF_0001        ",
        "MgD_0001       ",
        "MgHD_0001      ",
        "MgD2_0001      ",
        "NaD_0001       ",
        "SiD_0001       ",
        "Cl_0001        ",
        "HCl_0001       ",
        "DCl_0001       ",
        "H2Cl_0001      ",
        "HDCl_0001      ",
        "D2Cl_0001      ",
        "CCl_0001       ",
        "ClO_0001       ",
        "C3H3_0001      ",
        "C3H2D_0001     ",
        "C3HD2_0001     ",
        "C3D3_0001      ",
        "C3H4_0001      ",
        "C3H3D_0001     ",
        "C3H2D2_0001    ",
        "HCCN_0001      ",
        "DCCN_0001      ",
        "C2H2N_0001     ",
        "C2H5OH_0001    ",
        "C2H3D3_0001    ",
        "C2H5D_0001     ",
        "C2H4D2_0001    ",
        "C2H3N_0001     ",
        "C2H4O_0001     ",
        "CH2DCHO_0001   ",
        "CH3CDO_0001    ",
        "CD3CHO_0001    ",
        "CHD2CDO_0001   ",
        "CHD2CHO_0001   ",
        "CH2DCDO_0001   ",
        "C2H4D_0001     ",
        "C2H2D3_0001    ",
        "C2H3D2_0001    ",
        "C3H3N_0001     ",
        "H4C3N_0001     ",
        "C3D3N_0001     ",
        "HD3C3N_0001    ",
        "C3H2DN_0001    ",
        "H3DC3N_0001    ",
        "C3HD2N_0001    ",
        "H2D2C3N_0001   ",
        "HC3O_0001      ",
        "DC3O_0001      ",
        "C4H4_0001      ",
        "CH5N_0001      ",
        "CH3NH_0001     ",
        "FeD_0001       ",
        "H2C3N_0001     ",
        "H2C5N_0001     ",
        "H3C5N_0001     ",
        "H2C7N_0001     ",
        "H3C7N_0001     ",
        "H2C9N_0001     ",
        "H3C9N_0001     ",
        "H5C3N_0001     ",
        "C2H2O_0001     ",
        "C2HDO_0001     ",
        "C2D2O_0001     ",
        "HDC3N_0001     ",
        "D2C3N_0001     ",
        "H2C3O_0001     ",
        "HDC3O_0001     ",
        "D2C3O_0001     ",
        "C2HDN_0001     ",
        "C2D2N_0001     ",
        "CH3N_0001      ",
        "N2H2_0001      ",
        "N2HD_0001      ",
        "N2D2_0001      ",
        "NH2OH_0001     ",
        "NH2OD_0001     ",
        "N2H_0001       ",
        "N2D_0001       ",
        "CNH2_0001      ",
        "CHNH2_0001     ",
        "HON_0001       ",
        "DON_0001       ",
        "NHNO_0001      ",
        "CH2DN_0001     ",
        "CHD2N_0001     ",
        "CD3N_0001      ",
        "NH2NO_0001     ",
        "CH3CO_0001     ",
        "CH2DOCH3_0001  ",
        "CHD2OCH3_0001  ",
        "CH2DOCH2D_0001 ",
        "CD3OCH3_0001   ",
        "CHD2OCH2D_0001 ",
        "DCOOCH3_0001   ",
        "HCOOCH2D_0001  ",
        "DCOOCH2D_0001  ",
        "HCOOCHD2_0001  ",
        "DCOOCHD2_0001  ",
        "HCOOCD3_0001   ",
        "DCOOCD3_0001   ",
        "NH2CO_0001     ",
        "NHDCO_0001     ",
        "ND2CO_0001     ",
        "CH2DOCHO_0001  ",
        "CH3OCDO_0001   ",
        "CHD2OCHO_0001  ",
        "CH2DOCDO_0001  ",
        "CD3OCHO_0001   ",
        "CHD2OCDO_0001  ",
        "CHDOCDO_0001   ",
        "CHDOCHDO_0001  ",
        "CH3OCHD_0001   ",
        "CH2DOCH2_0001  ",
        "CH2DOCHD_0001  ",
        "CHD2OCH2_0001  ",
        "CH3OCD2_0001   ",
        "CH2DOCD2_0001  ",
        "CH3OCHD2_0001  ",
        "CH3OCH2D_0001  ",
        "CH2DCO_0001    ",
        "CH3NHD_0001    ",
        "CH2DNH2_0001   ",
        "P_0001         ",
        "PO_0001        ",
        "PH_0001        ",
        "PD_0001        ",
        "PH2_0001       ",
        "PHD_0001       ",
        "PD2_0001       ",
        "PN_0001        ",
        "CP_0001        ",
        "CCP_0001       ",
        "C3P_0001       ",
        "C4P_0001       ",
        "CH2PH_0001     ",
        "CHD2PH_0001    ",
        "CH2PD_0001     ",
        "CHDPD_0001     ",
        "CD2P2_0001     ",
        "HCP_0001       ",
        "DCP_0001       ",
        "HCCP_0001      ",
        "DCCP_0001      ",
        "HPO_0001       ",
        "DPO_0001       ",
        "surface_mask   ",
        "mantle_mask    ",
        "dummy          ")
        # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        ##END_SPECIESNAMES

        ##BEGIN_IDXLIST
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # NOTE: This block is auto-generated
        # WHEN: 2020-12-22 21:21:49
        # CHANGESET: xxxxxxx
        # BY: unknown@unknown

        self.idx_GRAIN0_gas=1
        self.idx_GRAINk_gas=2
        self.idx_E_gas=3
        self.idx_H_gas=4
        self.idx_Hj_gas=5
        self.idx_Hk_gas=6
        self.idx_D_gas=7
        self.idx_Dj_gas=8
        self.idx_Dk_gas=9
        self.idx_p_H2_gas=10
        self.idx_o_H2_gas=11
        self.idx_p_H2j_gas=12
        self.idx_o_H2j_gas=13
        self.idx_p_D2_gas=14
        self.idx_o_D2_gas=15
        self.idx_p_D2j_gas=16
        self.idx_o_D2j_gas=17
        self.idx_p_D3j_gas=18
        self.idx_o_D3j_gas=19
        self.idx_m_D3j_gas=20
        self.idx_p_D2Hj_gas=21
        self.idx_o_D2Hj_gas=22
        self.idx_p_H2Dj_gas=23
        self.idx_o_H2Dj_gas=24
        self.idx_HD_gas=25
        self.idx_HDj_gas=26
        self.idx_He_gas=27
        self.idx_Hej_gas=28
        self.idx_O_gas=29
        self.idx_Oj_gas=30
        self.idx_Ok_gas=31
        self.idx_O2_gas=32
        self.idx_O2j_gas=33
        self.idx_O3_gas=34
        self.idx_OH_gas=35
        self.idx_OHj_gas=36
        self.idx_OHk_gas=37
        self.idx_OD_gas=38
        self.idx_ODj_gas=39
        self.idx_ODk_gas=40
        self.idx_H2O_gas=41
        self.idx_H2Oj_gas=42
        self.idx_HDO_gas=43
        self.idx_HDOj_gas=44
        self.idx_D2O_gas=45
        self.idx_D2Oj_gas=46
        self.idx_O2H_gas=47
        self.idx_O2D_gas=48
        self.idx_HO2j_gas=49
        self.idx_DO2j_gas=50
        self.idx_HOOH_gas=51
        self.idx_HOOD_gas=52
        self.idx_DOOH_gas=53
        self.idx_DOOD_gas=54
        self.idx_p_H3j_gas=55
        self.idx_o_H3j_gas=56
        self.idx_H3Oj_gas=57
        self.idx_H2DOj_gas=58
        self.idx_HD2Oj_gas=59
        self.idx_D3Oj_gas=60
        self.idx_F_gas=61
        self.idx_Fj_gas=62
        self.idx_Cl_gas=63
        self.idx_Clj_gas=64
        self.idx_CFj_gas=65
        self.idx_C_gas=66
        self.idx_Cj_gas=67
        self.idx_Ck_gas=68
        self.idx_C2_gas=69
        self.idx_C2j_gas=70
        self.idx_C2Nj_gas=71
        self.idx_CO_gas=72
        self.idx_COj_gas=73
        self.idx_CO2_gas=74
        self.idx_CO2j_gas=75
        self.idx_N_gas=76
        self.idx_Nj_gas=77
        self.idx_N2_gas=78
        self.idx_N2j_gas=79
        self.idx_N2Hj_gas=80
        self.idx_N2Dj_gas=81
        self.idx_HCO_gas=82
        self.idx_DCO_gas=83
        self.idx_HCOj_gas=84
        self.idx_DCOj_gas=85
        self.idx_H2CO_gas=86
        self.idx_HDCO_gas=87
        self.idx_D2CO_gas=88
        self.idx_CH2OH_gas=89
        self.idx_CD2OD_gas=90
        self.idx_CH2OD_gas=91
        self.idx_CHDOH_gas=92
        self.idx_CHDOD_gas=93
        self.idx_CD2OH_gas=94
        self.idx_CH3O_gas=95
        self.idx_CHD2O_gas=96
        self.idx_CH2DO_gas=97
        self.idx_CD3O_gas=98
        self.idx_CH3OH_gas=99
        self.idx_CH3OD_gas=100
        self.idx_CHD2OH_gas=101
        self.idx_CHD2OD_gas=102
        self.idx_CH2DOH_gas=103
        self.idx_CH2DOD_gas=104
        self.idx_CD3OD_gas=105
        self.idx_CD3OH_gas=106
        self.idx_CH_gas=107
        self.idx_CD_gas=108
        self.idx_CH2_gas=109
        self.idx_CHD_gas=110
        self.idx_CD2_gas=111
        self.idx_CH3_gas=112
        self.idx_CH2D_gas=113
        self.idx_CHD2_gas=114
        self.idx_CD3_gas=115
        self.idx_CH4_gas=116
        self.idx_CH3D_gas=117
        self.idx_CH2D2_gas=118
        self.idx_CHD3_gas=119
        self.idx_CD4_gas=120
        self.idx_CHj_gas=121
        self.idx_CDj_gas=122
        self.idx_CH2j_gas=123
        self.idx_CHDj_gas=124
        self.idx_CD2j_gas=125
        self.idx_CH3j_gas=126
        self.idx_CH2Dj_gas=127
        self.idx_CHD2j_gas=128
        self.idx_CD3j_gas=129
        self.idx_CH4j_gas=130
        self.idx_CH3Dj_gas=131
        self.idx_CH2D2j_gas=132
        self.idx_CHD3j_gas=133
        self.idx_CD4j_gas=134
        self.idx_CH5j_gas=135
        self.idx_CH4Dj_gas=136
        self.idx_CH3D2j_gas=137
        self.idx_CH2D3j_gas=138
        self.idx_CHD4j_gas=139
        self.idx_CD5j_gas=140
        self.idx_HCOOH_gas=141
        self.idx_HCOOD_gas=142
        self.idx_DCOOH_gas=143
        self.idx_DCOOD_gas=144
        self.idx_HOCO_gas=145
        self.idx_DOCO_gas=146
        self.idx_NH_gas=147
        self.idx_NHj_gas=148
        self.idx_ND_gas=149
        self.idx_NDj_gas=150
        self.idx_NH2_gas=151
        self.idx_NH2j_gas=152
        self.idx_NHD_gas=153
        self.idx_NHDj_gas=154
        self.idx_ND2_gas=155
        self.idx_ND2j_gas=156
        self.idx_NH3_gas=157
        self.idx_NH3j_gas=158
        self.idx_NH2D_gas=159
        self.idx_NH2Dj_gas=160
        self.idx_NHD2_gas=161
        self.idx_NHD2j_gas=162
        self.idx_ND3_gas=163
        self.idx_ND3j_gas=164
        self.idx_NH4j_gas=165
        self.idx_NH3Dj_gas=166
        self.idx_NH2D2j_gas=167
        self.idx_NHD3j_gas=168
        self.idx_ND4j_gas=169
        self.idx_Fe_gas=170
        self.idx_Fej_gas=171
        self.idx_S_gas=172
        self.idx_Sj_gas=173
        self.idx_Sk_gas=174
        self.idx_C10_gas=175
        self.idx_C10j_gas=176
        self.idx_C10k_gas=177
        self.idx_C10H_gas=178
        self.idx_C10Hj_gas=179
        self.idx_C10Hk_gas=180
        self.idx_C10H2_gas=181
        self.idx_C10H2j_gas=182
        self.idx_C10H3j_gas=183
        self.idx_C10N_gas=184
        self.idx_C10Nj_gas=185
        self.idx_C11_gas=186
        self.idx_C11j_gas=187
        self.idx_C2Hj_gas=188
        self.idx_C2Dj_gas=189
        self.idx_C2H2_gas=190
        self.idx_C2H2j_gas=191
        self.idx_C2HD_gas=192
        self.idx_C2HDj_gas=193
        self.idx_C2D2_gas=194
        self.idx_C2D2j_gas=195
        self.idx_C2DOj_gas=196
        self.idx_C2H3_gas=197
        self.idx_C2H3j_gas=198
        self.idx_C2H2D_gas=199
        self.idx_C2H2Dj_gas=200
        self.idx_C2HD2_gas=201
        self.idx_C2HD2j_gas=202
        self.idx_C2D3_gas=203
        self.idx_C2D3j_gas=204
        self.idx_C2H4_gas=205
        self.idx_C2H4j_gas=206
        self.idx_C2H3D_gas=207
        self.idx_C2H3Dj_gas=208
        self.idx_C2H2D2_gas=209
        self.idx_C2H2D2j_gas=210
        self.idx_C2HD3_gas=211
        self.idx_C2HD3j_gas=212
        self.idx_C2D4_gas=213
        self.idx_C2D4j_gas=214
        self.idx_C2H4Oj_gas=215
        self.idx_C2H5_gas=216
        self.idx_C2H5j_gas=217
        self.idx_C2H5OHj_gas=218
        self.idx_C2H5OH2j_gas=219
        self.idx_C2H6_gas=220
        self.idx_C2H6j_gas=221
        self.idx_C2H6COj_gas=222
        self.idx_C2H7j_gas=223
        self.idx_C2HOj_gas=224
        self.idx_C2N2j_gas=225
        self.idx_C2Oj_gas=226
        self.idx_C2Sj_gas=227
        self.idx_C3_gas=228
        self.idx_C3j_gas=229
        self.idx_C3k_gas=230
        self.idx_C3Hj_gas=231
        self.idx_C3Dj_gas=232
        self.idx_C3H3Nj_gas=233
        self.idx_C3H3NHj_gas=234
        self.idx_C3H4_gas=235
        self.idx_C3H4j_gas=236
        self.idx_C3H5j_gas=237
        self.idx_C3H6OHj_gas=238
        self.idx_C3N_gas=239
        self.idx_C3Nj_gas=240
        self.idx_C3Nk_gas=241
        self.idx_C3O_gas=242
        self.idx_C3Oj_gas=243
        self.idx_C3S_gas=244
        self.idx_C3Sj_gas=245
        self.idx_C4_gas=246
        self.idx_C4j_gas=247
        self.idx_C4k_gas=248
        self.idx_C4H_gas=249
        self.idx_C4Hj_gas=250
        self.idx_C4Hk_gas=251
        self.idx_C4D_gas=252
        self.idx_C4Dj_gas=253
        self.idx_C4Dk_gas=254
        self.idx_C4H2_gas=255
        self.idx_C4H2j_gas=256
        self.idx_C4HD_gas=257
        self.idx_C4HDj_gas=258
        self.idx_C4D2_gas=259
        self.idx_C4D2j_gas=260
        self.idx_C4H3_gas=261
        self.idx_C4H3j_gas=262
        self.idx_C4H4j_gas=263
        self.idx_C4H5j_gas=264
        self.idx_C4H7j_gas=265
        self.idx_C4N_gas=266
        self.idx_C4Nj_gas=267
        self.idx_C4S_gas=268
        self.idx_C4Sj_gas=269
        self.idx_C5_gas=270
        self.idx_C5j_gas=271
        self.idx_C5k_gas=272
        self.idx_C5H_gas=273
        self.idx_C5Hj_gas=274
        self.idx_C5Hk_gas=275
        self.idx_C5D_gas=276
        self.idx_C5Dj_gas=277
        self.idx_C5Dk_gas=278
        self.idx_C5H2_gas=279
        self.idx_C5H2j_gas=280
        self.idx_C5H3_gas=281
        self.idx_C5H3j_gas=282
        self.idx_C5H3Nj_gas=283
        self.idx_C5H4_gas=284
        self.idx_C5H4j_gas=285
        self.idx_C5H4Nj_gas=286
        self.idx_C5H5j_gas=287
        self.idx_C5N_gas=288
        self.idx_C5Nj_gas=289
        self.idx_C5O_gas=290
        self.idx_C6_gas=291
        self.idx_C6j_gas=292
        self.idx_C6k_gas=293
        self.idx_C6H_gas=294
        self.idx_C6Hj_gas=295
        self.idx_C6Hk_gas=296
        self.idx_C6H2_gas=297
        self.idx_C6H2j_gas=298
        self.idx_C6H3_gas=299
        self.idx_C6H3j_gas=300
        self.idx_C6H4_gas=301
        self.idx_C6H4j_gas=302
        self.idx_C6H5j_gas=303
        self.idx_C6H6_gas=304
        self.idx_C6H7j_gas=305
        self.idx_C6N_gas=306
        self.idx_C6Nj_gas=307
        self.idx_C7_gas=308
        self.idx_C7j_gas=309
        self.idx_C7k_gas=310
        self.idx_C7H_gas=311
        self.idx_C7Hj_gas=312
        self.idx_C7Hk_gas=313
        self.idx_C7H2_gas=314
        self.idx_C7H2j_gas=315
        self.idx_C7H2Nj_gas=316
        self.idx_C7H3_gas=317
        self.idx_C7H3j_gas=318
        self.idx_C7H4_gas=319
        self.idx_C7H4j_gas=320
        self.idx_C7H5j_gas=321
        self.idx_C7N_gas=322
        self.idx_C7Nj_gas=323
        self.idx_C7O_gas=324
        self.idx_C8_gas=325
        self.idx_C8j_gas=326
        self.idx_C8k_gas=327
        self.idx_C8H_gas=328
        self.idx_C8Hj_gas=329
        self.idx_C8Hk_gas=330
        self.idx_C8H2_gas=331
        self.idx_C8H2j_gas=332
        self.idx_C8H3_gas=333
        self.idx_C8H3j_gas=334
        self.idx_C8H4_gas=335
        self.idx_C8H4j_gas=336
        self.idx_C8H4Nj_gas=337
        self.idx_C8H5j_gas=338
        self.idx_C8N_gas=339
        self.idx_C8Nj_gas=340
        self.idx_C9_gas=341
        self.idx_C9j_gas=342
        self.idx_C9k_gas=343
        self.idx_C9H_gas=344
        self.idx_C9Hj_gas=345
        self.idx_C9Hk_gas=346
        self.idx_C9H2_gas=347
        self.idx_C9H2j_gas=348
        self.idx_C9H2Nj_gas=349
        self.idx_C9H3_gas=350
        self.idx_C9H3j_gas=351
        self.idx_C9H3Nj_gas=352
        self.idx_C9H4_gas=353
        self.idx_C9H4j_gas=354
        self.idx_C9H5j_gas=355
        self.idx_C9HNj_gas=356
        self.idx_C9N_gas=357
        self.idx_C9Nj_gas=358
        self.idx_C9O_gas=359
        self.idx_CCH_gas=360
        self.idx_CCD_gas=361
        self.idx_CCN_gas=362
        self.idx_CCO_gas=363
        self.idx_CCP_gas=364
        self.idx_CCPj_gas=365
        self.idx_CCS_gas=366
        self.idx_CD2ND2_gas=367
        self.idx_CD2ND2j_gas=368
        self.idx_CD2NH2_gas=369
        self.idx_CD2NH2j_gas=370
        self.idx_CD2NHD_gas=371
        self.idx_CD2NHDj_gas=372
        self.idx_CD3CN_gas=373
        self.idx_CD3CNj_gas=374
        self.idx_CD3COj_gas=375
        self.idx_CD3O2j_gas=376
        self.idx_CD3OHj_gas=377
        self.idx_CD3ODj_gas=378
        self.idx_CH2CCH_gas=379
        self.idx_CH2CCD_gas=380
        self.idx_CD2CCH_gas=381
        self.idx_CD2CCD_gas=382
        self.idx_CHDCCH_gas=383
        self.idx_CHDCCD_gas=384
        self.idx_CH2CHC2H_gas=385
        self.idx_CH2CHCHCH2_gas=386
        self.idx_CH2CHCN_gas=387
        self.idx_CH2CNj_gas=388
        self.idx_CHDCNj_gas=389
        self.idx_CD2CNj_gas=390
        self.idx_CH2DCN_gas=391
        self.idx_CH2DCNj_gas=392
        self.idx_CH2DCOj_gas=393
        self.idx_CH2DO2j_gas=394
        self.idx_CH2DODj_gas=395
        self.idx_CH2ND2_gas=396
        self.idx_CH2ND2j_gas=397
        self.idx_CH2NH_gas=398
        self.idx_CH2ND_gas=399
        self.idx_CHDNH_gas=400
        self.idx_CHDND_gas=401
        self.idx_CD2NH_gas=402
        self.idx_CD2ND_gas=403
        self.idx_CH2NH2_gas=404
        self.idx_CH2NH2j_gas=405
        self.idx_CH2NHD_gas=406
        self.idx_CH2NHDj_gas=407
        self.idx_CH3C3N_gas=408
        self.idx_CH3C4H_gas=409
        self.idx_CH3C5N_gas=410
        self.idx_CH3C6H_gas=411
        self.idx_CH3C7N_gas=412
        self.idx_CH3CCH_gas=413
        self.idx_CH3CH2OH_gas=414
        self.idx_CH3CHCH2_gas=415
        self.idx_CH3CHO_gas=416
        self.idx_CH3CHOHj_gas=417
        self.idx_CH3CN_gas=418
        self.idx_CH3CNj_gas=419
        self.idx_CH3CNHj_gas=420
        self.idx_CH3COj_gas=421
        self.idx_CH3COCH3_gas=422
        self.idx_CH3NH2_gas=423
        self.idx_CH3NH2j_gas=424
        self.idx_CH3NH3j_gas=425
        self.idx_CH3O2j_gas=426
        self.idx_CH3OCH2_gas=427
        self.idx_CH3OCH3_gas=428
        self.idx_CH3OCH3j_gas=429
        self.idx_CH3OCH4j_gas=430
        self.idx_CH3OHj_gas=431
        self.idx_CH2DOHj_gas=432
        self.idx_CH3ODj_gas=433
        self.idx_CH3OH2j_gas=434
        self.idx_CHD2CN_gas=435
        self.idx_CHD2CNj_gas=436
        self.idx_CHD2COj_gas=437
        self.idx_CHD2O2j_gas=438
        self.idx_CHD2OHj_gas=439
        self.idx_CHD2ODj_gas=440
        self.idx_CHDNH2_gas=441
        self.idx_CHDNH2j_gas=442
        self.idx_CHDNHD_gas=443
        self.idx_CHDNHDj_gas=444
        self.idx_CHDND2_gas=445
        self.idx_CHDND2j_gas=446
        self.idx_CN_gas=447
        self.idx_CNj_gas=448
        self.idx_CNk_gas=449
        self.idx_CNCj_gas=450
        self.idx_COOCH4j_gas=451
        self.idx_CS_gas=452
        self.idx_CSj_gas=453
        self.idx_D2C3Oj_gas=454
        self.idx_DC2NCHj_gas=455
        self.idx_DC2NCDj_gas=456
        self.idx_DCNDj_gas=457
        self.idx_DCNHj_gas=458
        self.idx_DCS_gas=459
        self.idx_DCSj_gas=460
        self.idx_DN2Oj_gas=461
        self.idx_DNC_gas=462
        self.idx_DNCj_gas=463
        self.idx_DNCCC_gas=464
        self.idx_DNCO_gas=465
        self.idx_DNCOj_gas=466
        self.idx_DNO_gas=467
        self.idx_DNOj_gas=468
        self.idx_DNSj_gas=469
        self.idx_DOCj_gas=470
        self.idx_DOCOj_gas=471
        self.idx_DOCSj_gas=472
        self.idx_HS_gas=473
        self.idx_HSj_gas=474
        self.idx_DS_gas=475
        self.idx_DSj_gas=476
        self.idx_FeH_gas=477
        self.idx_FeD_gas=478
        self.idx_H2C10Nj_gas=479
        self.idx_H2C3Oj_gas=480
        self.idx_H2C4Nj_gas=481
        self.idx_H2C5Nj_gas=482
        self.idx_H2C6Nj_gas=483
        self.idx_H2C8Nj_gas=484
        self.idx_H2CCN_gas=485
        self.idx_HDCCN_gas=486
        self.idx_D2CCN_gas=487
        self.idx_H2CCO_gas=488
        self.idx_H2CCOj_gas=489
        self.idx_HDCCO_gas=490
        self.idx_HDCCOj_gas=491
        self.idx_D2CCO_gas=492
        self.idx_D2CCOj_gas=493
        self.idx_H2CN_gas=494
        self.idx_H2CNj_gas=495
        self.idx_HDCN_gas=496
        self.idx_D2CN_gas=497
        self.idx_H2COj_gas=498
        self.idx_HDCOj_gas=499
        self.idx_D2COj_gas=500
        self.idx_H2COHj_gas=501
        self.idx_H2CODj_gas=502
        self.idx_HDCOHj_gas=503
        self.idx_HDCODj_gas=504
        self.idx_D2COHj_gas=505
        self.idx_D2CODj_gas=506
        self.idx_H2CS_gas=507
        self.idx_H2CSj_gas=508
        self.idx_HDCS_gas=509
        self.idx_HDCSj_gas=510
        self.idx_D2CS_gas=511
        self.idx_D2CSj_gas=512
        self.idx_H2NCj_gas=513
        self.idx_HDNCj_gas=514
        self.idx_D2NCj_gas=515
        self.idx_H2NOj_gas=516
        self.idx_HDNOj_gas=517
        self.idx_D2NOj_gas=518
        self.idx_H2S_gas=519
        self.idx_H2Sj_gas=520
        self.idx_HDS_gas=521
        self.idx_HDSj_gas=522
        self.idx_D2S_gas=523
        self.idx_D2Sj_gas=524
        self.idx_H2S2j_gas=525
        self.idx_HDS2j_gas=526
        self.idx_D2S2j_gas=527
        self.idx_H3C3Oj_gas=528
        self.idx_H3C4Nj_gas=529
        self.idx_H3C4NHj_gas=530
        self.idx_H3C6NHj_gas=531
        self.idx_H3C7Nj_gas=532
        self.idx_H3CSj_gas=533
        self.idx_H2DCSj_gas=534
        self.idx_HD2CSj_gas=535
        self.idx_D3CSj_gas=536
        self.idx_H3Sj_gas=537
        self.idx_H2DSj_gas=538
        self.idx_HD2Sj_gas=539
        self.idx_D3Sj_gas=540
        self.idx_H3S2j_gas=541
        self.idx_H2DS2j_gas=542
        self.idx_HD2S2j_gas=543
        self.idx_D3S2j_gas=544
        self.idx_H5C2O2j_gas=545
        self.idx_HC10Nj_gas=546
        self.idx_HC2N_gas=547
        self.idx_HC2Nj_gas=548
        self.idx_DC2N_gas=549
        self.idx_DC2Nj_gas=550
        self.idx_HC2NCDj_gas=551
        self.idx_HC2NCHj_gas=552
        self.idx_HC2O_gas=553
        self.idx_DC2O_gas=554
        self.idx_HC2Sj_gas=555
        self.idx_DC2Sj_gas=556
        self.idx_HC3N_gas=557
        self.idx_HC3Nj_gas=558
        self.idx_DC3N_gas=559
        self.idx_DC3Nj_gas=560
        self.idx_HC3NHj_gas=561
        self.idx_HC3NDj_gas=562
        self.idx_DC3NHj_gas=563
        self.idx_DC3NDj_gas=564
        self.idx_HC3Oj_gas=565
        self.idx_DC3Oj_gas=566
        self.idx_HC3Sj_gas=567
        self.idx_DC3Sj_gas=568
        self.idx_HC4N_gas=569
        self.idx_HC4Nj_gas=570
        self.idx_DC4N_gas=571
        self.idx_DC4Nj_gas=572
        self.idx_HC4Oj_gas=573
        self.idx_DC4Oj_gas=574
        self.idx_HC4Sj_gas=575
        self.idx_DC4Sj_gas=576
        self.idx_HC5N_gas=577
        self.idx_HC5Nj_gas=578
        self.idx_HC5Oj_gas=579
        self.idx_HC6N_gas=580
        self.idx_HC6Nj_gas=581
        self.idx_HC7N_gas=582
        self.idx_HC7Nj_gas=583
        self.idx_HC7Oj_gas=584
        self.idx_HC8N_gas=585
        self.idx_HC8Nj_gas=586
        self.idx_HC9N_gas=587
        self.idx_HC9Oj_gas=588
        self.idx_HCCNC_gas=589
        self.idx_DCCNC_gas=590
        self.idx_HCN_gas=591
        self.idx_HCNj_gas=592
        self.idx_DCN_gas=593
        self.idx_DCNj_gas=594
        self.idx_HCNCC_gas=595
        self.idx_DCNCC_gas=596
        self.idx_HCNHj_gas=597
        self.idx_HCNDj_gas=598
        self.idx_HCOOCH3_gas=599
        self.idx_HCOOHj_gas=600
        self.idx_HCOODj_gas=601
        self.idx_DCOOHj_gas=602
        self.idx_DCOODj_gas=603
        self.idx_HCS_gas=604
        self.idx_HCSj_gas=605
        self.idx_HDC3Oj_gas=606
        self.idx_HN2Oj_gas=607
        self.idx_HNC_gas=608
        self.idx_HNCj_gas=609
        self.idx_HNCCC_gas=610
        self.idx_HNCO_gas=611
        self.idx_HNCOj_gas=612
        self.idx_HNO_gas=613
        self.idx_HNOj_gas=614
        self.idx_HNSj_gas=615
        self.idx_HOCj_gas=616
        self.idx_HOCOj_gas=617
        self.idx_HOCSj_gas=618
        self.idx_HS2j_gas=619
        self.idx_DS2j_gas=620
        self.idx_HSOj_gas=621
        self.idx_DSOj_gas=622
        self.idx_HSO2j_gas=623
        self.idx_DSO2j_gas=624
        self.idx_HSS_gas=625
        self.idx_DSS_gas=626
        self.idx_HSSH_gas=627
        self.idx_HSSD_gas=628
        self.idx_DSSH_gas=629
        self.idx_DSSD_gas=630
        self.idx_N2O_gas=631
        self.idx_NC4N_gas=632
        self.idx_NC6N_gas=633
        self.idx_NC8N_gas=634
        self.idx_NCOj_gas=635
        self.idx_NH2CH2Oj_gas=636
        self.idx_NH2CHO_gas=637
        self.idx_NH2CDO_gas=638
        self.idx_NHDCHO_gas=639
        self.idx_NHDCDO_gas=640
        self.idx_ND2CHO_gas=641
        self.idx_ND2CDO_gas=642
        self.idx_NH2CN_gas=643
        self.idx_NHDCN_gas=644
        self.idx_ND2CN_gas=645
        self.idx_NH2CNDj_gas=646
        self.idx_NHDCNDj_gas=647
        self.idx_ND2CNDj_gas=648
        self.idx_NH2CNHj_gas=649
        self.idx_NHDCNHj_gas=650
        self.idx_ND2CNHj_gas=651
        self.idx_NO_gas=652
        self.idx_NOj_gas=653
        self.idx_NO2_gas=654
        self.idx_NO2j_gas=655
        self.idx_NS_gas=656
        self.idx_NSj_gas=657
        self.idx_OCN_gas=658
        self.idx_OCS_gas=659
        self.idx_OCSj_gas=660
        self.idx_S2_gas=661
        self.idx_S2j_gas=662
        self.idx_SO_gas=663
        self.idx_SOj_gas=664
        self.idx_SO2_gas=665
        self.idx_SO2j_gas=666
        self.idx_CH3OCHO_gas=667
        self.idx_Si_gas=668
        self.idx_Sij_gas=669
        self.idx_SiO_gas=670
        self.idx_SiOj_gas=671
        self.idx_SiO2_gas=672
        self.idx_SiS_gas=673
        self.idx_SiSj_gas=674
        self.idx_SiN_gas=675
        self.idx_SiNj_gas=676
        self.idx_SiC_gas=677
        self.idx_SiCj_gas=678
        self.idx_SiC2_gas=679
        self.idx_SiC2j_gas=680
        self.idx_SiH_gas=681
        self.idx_SiD_gas=682
        self.idx_SiH2_gas=683
        self.idx_SiH2j_gas=684
        self.idx_SiHD_gas=685
        self.idx_SiHDj_gas=686
        self.idx_SiD2_gas=687
        self.idx_SiD2j_gas=688
        self.idx_SiH3_gas=689
        self.idx_SiH4_gas=690
        self.idx_P_gas=691
        self.idx_Pj_gas=692
        self.idx_PH_gas=693
        self.idx_PHj_gas=694
        self.idx_PD_gas=695
        self.idx_PDj_gas=696
        self.idx_PO_gas=697
        self.idx_POj_gas=698
        self.idx_PH2_gas=699
        self.idx_PH2j_gas=700
        self.idx_PHD_gas=701
        self.idx_PHDj_gas=702
        self.idx_PD2_gas=703
        self.idx_PD2j_gas=704
        self.idx_C3P_gas=705
        self.idx_C4P_gas=706
        self.idx_C4Pj_gas=707
        self.idx_CH2PH_gas=708
        self.idx_CH2PD_gas=709
        self.idx_CHDPH_gas=710
        self.idx_CHDPD_gas=711
        self.idx_CD2PH_gas=712
        self.idx_CD2PD_gas=713
        self.idx_CH2Sij_gas=714
        self.idx_CHDSij_gas=715
        self.idx_CD2Sij_gas=716
        self.idx_CHSij_gas=717
        self.idx_CDSij_gas=718
        self.idx_CP_gas=719
        self.idx_CPj_gas=720
        self.idx_D2POj_gas=721
        self.idx_D3SiOj_gas=722
        self.idx_DCP_gas=723
        self.idx_DCPj_gas=724
        self.idx_DF_gas=725
        self.idx_DFj_gas=726
        self.idx_DPO_gas=727
        self.idx_DPOj_gas=728
        self.idx_DSiNHj_gas=729
        self.idx_DSiNDj_gas=730
        self.idx_H2CSiCH_gas=731
        self.idx_H2CSiCD_gas=732
        self.idx_HDCSiCD_gas=733
        self.idx_HDCSiCH_gas=734
        self.idx_D2CSiCH_gas=735
        self.idx_D2CSiCD_gas=736
        self.idx_H2DSiOj_gas=737
        self.idx_H2Fj_gas=738
        self.idx_HDFj_gas=739
        self.idx_D2Fj_gas=740
        self.idx_H2POj_gas=741
        self.idx_H2SiO_gas=742
        self.idx_H2SiOj_gas=743
        self.idx_HDSiO_gas=744
        self.idx_HDSiOj_gas=745
        self.idx_D2SiO_gas=746
        self.idx_D2SiOj_gas=747
        self.idx_H3SiOj_gas=748
        self.idx_HCCP_gas=749
        self.idx_DCCP_gas=750
        self.idx_HCCSi_gas=751
        self.idx_DCCSi_gas=752
        self.idx_HCP_gas=753
        self.idx_HCPj_gas=754
        self.idx_HCSi_gas=755
        self.idx_DCSi_gas=756
        self.idx_HD2SiOj_gas=757
        self.idx_HDPOj_gas=758
        self.idx_HF_gas=759
        self.idx_HFj_gas=760
        self.idx_HNSi_gas=761
        self.idx_HNSij_gas=762
        self.idx_DNSi_gas=763
        self.idx_DNSij_gas=764
        self.idx_HPNj_gas=765
        self.idx_DPNj_gas=766
        self.idx_HPO_gas=767
        self.idx_HPOj_gas=768
        self.idx_HSiNHj_gas=769
        self.idx_HSiNDj_gas=770
        self.idx_HSiOj_gas=771
        self.idx_DSiOj_gas=772
        self.idx_HSiO2j_gas=773
        self.idx_DSiO2j_gas=774
        self.idx_HSiSj_gas=775
        self.idx_DSiSj_gas=776
        self.idx_HeHj_gas=777
        self.idx_HeDj_gas=778
        self.idx_PC2Hj_gas=779
        self.idx_PC2Dj_gas=780
        self.idx_PC2H2j_gas=781
        self.idx_PC2HDj_gas=782
        self.idx_PC2D2j_gas=783
        self.idx_PC2H3j_gas=784
        self.idx_PC2H2Dj_gas=785
        self.idx_PC2HD2j_gas=786
        self.idx_PC2D3j_gas=787
        self.idx_PC2H4j_gas=788
        self.idx_PC3Hj_gas=789
        self.idx_PC3Dj_gas=790
        self.idx_PC4Hj_gas=791
        self.idx_PC4Dj_gas=792
        self.idx_PC4H2j_gas=793
        self.idx_PCH2j_gas=794
        self.idx_PCHDj_gas=795
        self.idx_PCD2j_gas=796
        self.idx_PCH3j_gas=797
        self.idx_PCH2Dj_gas=798
        self.idx_PCHD2j_gas=799
        self.idx_PCD3j_gas=800
        self.idx_PCH4j_gas=801
        self.idx_PCH3Dj_gas=802
        self.idx_PCH2D2j_gas=803
        self.idx_PCHD3j_gas=804
        self.idx_PCD4j_gas=805
        self.idx_PH3j_gas=806
        self.idx_PH2Dj_gas=807
        self.idx_PHD2j_gas=808
        self.idx_PD3j_gas=809
        self.idx_PN_gas=810
        self.idx_PNj_gas=811
        self.idx_PNH2j_gas=812
        self.idx_PNHDj_gas=813
        self.idx_PND2j_gas=814
        self.idx_PNH3j_gas=815
        self.idx_PNHD2j_gas=816
        self.idx_PNH2Dj_gas=817
        self.idx_PND3j_gas=818
        self.idx_SiC2CH3_gas=819
        self.idx_SiC2Hj_gas=820
        self.idx_SiC2Dj_gas=821
        self.idx_SiC2H2j_gas=822
        self.idx_SiC2HDj_gas=823
        self.idx_SiC2D2j_gas=824
        self.idx_SiC2H3j_gas=825
        self.idx_SiC2H2Dj_gas=826
        self.idx_SiC2HD2j_gas=827
        self.idx_SiC2D3j_gas=828
        self.idx_SiC3D_gas=829
        self.idx_SiC3Dj_gas=830
        self.idx_SiC3H_gas=831
        self.idx_SiC3Hj_gas=832
        self.idx_SiC3H2j_gas=833
        self.idx_SiC3HDj_gas=834
        self.idx_SiC3D2j_gas=835
        self.idx_SiC3H5_gas=836
        self.idx_SiC4_gas=837
        self.idx_SiC4j_gas=838
        self.idx_SiC4D_gas=839
        self.idx_SiC4Dj_gas=840
        self.idx_SiC4H_gas=841
        self.idx_SiC4Hj_gas=842
        self.idx_SiC6H_gas=843
        self.idx_SiC8H_gas=844
        self.idx_SiCD3_gas=845
        self.idx_SiCD3j_gas=846
        self.idx_SiCH2_gas=847
        self.idx_SiCHD_gas=848
        self.idx_SiCD2_gas=849
        self.idx_SiCH2D_gas=850
        self.idx_SiCH2Dj_gas=851
        self.idx_SiCH3_gas=852
        self.idx_SiCH3j_gas=853
        self.idx_SiCH4j_gas=854
        self.idx_SiCH3Dj_gas=855
        self.idx_SiCH2D2j_gas=856
        self.idx_SiCHD3j_gas=857
        self.idx_SiCD4j_gas=858
        self.idx_SiCHD2_gas=859
        self.idx_SiCHD2j_gas=860
        self.idx_SiD3_gas=861
        self.idx_SiD3j_gas=862
        self.idx_SiFj_gas=863
        self.idx_SiHj_gas=864
        self.idx_SiDj_gas=865
        self.idx_SiH2D_gas=866
        self.idx_SiH2Dj_gas=867
        self.idx_SiH2D2_gas=868
        self.idx_SiH2D2j_gas=869
        self.idx_SiH3j_gas=870
        self.idx_SiH3D_gas=871
        self.idx_SiH3Dj_gas=872
        self.idx_SiH4j_gas=873
        self.idx_SiD4_gas=874
        self.idx_SiD4j_gas=875
        self.idx_SiH5j_gas=876
        self.idx_SiH4Dj_gas=877
        self.idx_SiH3D2j_gas=878
        self.idx_SiH2D3j_gas=879
        self.idx_SiHD4j_gas=880
        self.idx_SiD5j_gas=881
        self.idx_SiHD2_gas=882
        self.idx_SiHD2j_gas=883
        self.idx_SiHD3_gas=884
        self.idx_SiHD3j_gas=885
        self.idx_SiNC_gas=886
        self.idx_SiNCj_gas=887
        self.idx_SiNCHj_gas=888
        self.idx_SiNCDj_gas=889
        self.idx_c_DCCHSi_gas=890
        self.idx_c_DCCDSi_gas=891
        self.idx_c_HCCHSi_gas=892
        self.idx_c_HCCDSi_gas=893
        self.idx_c_SiC2_gas=894
        self.idx_l_SiC3_gas=895
        self.idx_l_SiC3j_gas=896
        self.idx_l_C3H_gas=897
        self.idx_l_C3D_gas=898
        self.idx_c_C3H_gas=899
        self.idx_c_C3D_gas=900
        self.idx_l_C3H2_gas=901
        self.idx_l_C3HD_gas=902
        self.idx_l_C3D2_gas=903
        self.idx_c_C3H2_gas=904
        self.idx_c_C3HD_gas=905
        self.idx_c_C3D2_gas=906
        self.idx_l_C3H2j_gas=907
        self.idx_l_C3HDj_gas=908
        self.idx_l_C3D2j_gas=909
        self.idx_c_C3H2j_gas=910
        self.idx_c_C3HDj_gas=911
        self.idx_c_C3D2j_gas=912
        self.idx_l_C3H3j_gas=913
        self.idx_l_C3H2Dj_gas=914
        self.idx_l_C3HD2j_gas=915
        self.idx_l_C3D3j_gas=916
        self.idx_c_C3H3j_gas=917
        self.idx_c_C3H2Dj_gas=918
        self.idx_c_C3HD2j_gas=919
        self.idx_c_C3D3j_gas=920
        self.idx_Mg_gas=921
        self.idx_Mgj_gas=922
        self.idx_MgH_gas=923
        self.idx_MgD_gas=924
        self.idx_MgH2_gas=925
        self.idx_MgHD_gas=926
        self.idx_MgD2_gas=927
        self.idx_HCl_gas=928
        self.idx_HClj_gas=929
        self.idx_DCl_gas=930
        self.idx_DClj_gas=931
        self.idx_H2Cl_gas=932
        self.idx_H2Clj_gas=933
        self.idx_HDCl_gas=934
        self.idx_HDClj_gas=935
        self.idx_D2Cl_gas=936
        self.idx_D2Clj_gas=937
        self.idx_CCl_gas=938
        self.idx_CClj_gas=939
        self.idx_ClO_gas=940
        self.idx_ClOj_gas=941
        self.idx_H2CClj_gas=942
        self.idx_HDCClj_gas=943
        self.idx_D2CClj_gas=944
        self.idx_Na_gas=945
        self.idx_Naj_gas=946
        self.idx_NaH_gas=947
        self.idx_NaD_gas=948
        self.idx_NaOH_gas=949
        self.idx_NaOD_gas=950
        self.idx_NaH2j_gas=951
        self.idx_NaHDj_gas=952
        self.idx_NaD2j_gas=953
        self.idx_NaH2Oj_gas=954
        self.idx_NaHDOj_gas=955
        self.idx_NaD2Oj_gas=956
        self.idx_H_0001=957
        self.idx_H_0002=1413
        self.idx_D_0001=958
        self.idx_D_0002=1414
        self.idx_p_H2_0001=959
        self.idx_p_H2_0002=1415
        self.idx_p_D2_0001=960
        self.idx_p_D2_0002=1416
        self.idx_o_H2_0001=961
        self.idx_o_H2_0002=1417
        self.idx_o_D2_0001=962
        self.idx_o_D2_0002=1418
        self.idx_HD_0001=963
        self.idx_HD_0002=1419
        self.idx_He_0001=964
        self.idx_He_0002=1420
        self.idx_O_0001=965
        self.idx_O_0002=1421
        self.idx_O2_0001=966
        self.idx_O2_0002=1422
        self.idx_O3_0001=967
        self.idx_O3_0002=1423
        self.idx_OH_0001=968
        self.idx_OH_0002=1424
        self.idx_OD_0001=969
        self.idx_OD_0002=1425
        self.idx_H2O_0001=970
        self.idx_H2O_0002=1426
        self.idx_HDO_0001=971
        self.idx_HDO_0002=1427
        self.idx_D2O_0001=972
        self.idx_D2O_0002=1428
        self.idx_O2H_0001=973
        self.idx_O2H_0002=1429
        self.idx_O2D_0001=974
        self.idx_O2D_0002=1430
        self.idx_HOOH_0001=975
        self.idx_HOOH_0002=1431
        self.idx_HOOD_0001=976
        self.idx_HOOD_0002=1432
        self.idx_DOOD_0001=977
        self.idx_DOOD_0002=1433
        self.idx_Fe_0001=978
        self.idx_Fe_0002=1434
        self.idx_FeH_0001=979
        self.idx_FeH_0002=1435
        self.idx_N_0001=980
        self.idx_N_0002=1436
        self.idx_S_0001=981
        self.idx_S_0002=1437
        self.idx_C_0001=982
        self.idx_C_0002=1438
        self.idx_C2_0001=983
        self.idx_C2_0002=1439
        self.idx_CO_0001=984
        self.idx_CO_0002=1440
        self.idx_HCO_0001=985
        self.idx_HCO_0002=1441
        self.idx_DCO_0001=986
        self.idx_DCO_0002=1442
        self.idx_H2CO_0001=987
        self.idx_H2CO_0002=1443
        self.idx_HDCO_0001=988
        self.idx_HDCO_0002=1444
        self.idx_D2CO_0001=989
        self.idx_D2CO_0002=1445
        self.idx_CH2OH_0001=990
        self.idx_CH2OH_0002=1446
        self.idx_CD2OD_0001=991
        self.idx_CD2OD_0002=1447
        self.idx_CH2OD_0001=992
        self.idx_CH2OD_0002=1448
        self.idx_CHDOH_0001=993
        self.idx_CHDOH_0002=1449
        self.idx_CHDOD_0001=994
        self.idx_CHDOD_0002=1450
        self.idx_CD2OH_0001=995
        self.idx_CD2OH_0002=1451
        self.idx_CH3O_0001=996
        self.idx_CH3O_0002=1452
        self.idx_CHD2O_0001=997
        self.idx_CHD2O_0002=1453
        self.idx_CH2DO_0001=998
        self.idx_CH2DO_0002=1454
        self.idx_CD3O_0001=999
        self.idx_CD3O_0002=1455
        self.idx_CH3OH_0001=1000
        self.idx_CH3OH_0002=1456
        self.idx_CH3OD_0001=1001
        self.idx_CH3OD_0002=1457
        self.idx_CHD2OH_0001=1002
        self.idx_CHD2OH_0002=1458
        self.idx_CHD2OD_0001=1003
        self.idx_CHD2OD_0002=1459
        self.idx_CH2DOH_0001=1004
        self.idx_CH2DOH_0002=1460
        self.idx_CH2DOD_0001=1005
        self.idx_CH2DOD_0002=1461
        self.idx_CD3OD_0001=1006
        self.idx_CD3OD_0002=1462
        self.idx_CD3OH_0001=1007
        self.idx_CD3OH_0002=1463
        self.idx_CH_0001=1008
        self.idx_CH_0002=1464
        self.idx_CD_0001=1009
        self.idx_CD_0002=1465
        self.idx_CH2_0001=1010
        self.idx_CH2_0002=1466
        self.idx_CHD_0001=1011
        self.idx_CHD_0002=1467
        self.idx_CD2_0001=1012
        self.idx_CD2_0002=1468
        self.idx_CH3_0001=1013
        self.idx_CH3_0002=1469
        self.idx_CH2D_0001=1014
        self.idx_CH2D_0002=1470
        self.idx_CHD2_0001=1015
        self.idx_CHD2_0002=1471
        self.idx_CD3_0001=1016
        self.idx_CD3_0002=1472
        self.idx_CH4_0001=1017
        self.idx_CH4_0002=1473
        self.idx_CH3D_0001=1018
        self.idx_CH3D_0002=1474
        self.idx_CH2D2_0001=1019
        self.idx_CH2D2_0002=1475
        self.idx_CHD3_0001=1020
        self.idx_CHD3_0002=1476
        self.idx_CD4_0001=1021
        self.idx_CD4_0002=1477
        self.idx_CO2_0001=1022
        self.idx_CO2_0002=1478
        self.idx_HCOOH_0001=1023
        self.idx_HCOOH_0002=1479
        self.idx_HCOOD_0001=1024
        self.idx_HCOOD_0002=1480
        self.idx_DCOOH_0001=1025
        self.idx_DCOOH_0002=1481
        self.idx_DCOOD_0001=1026
        self.idx_DCOOD_0002=1482
        self.idx_HOCO_0001=1027
        self.idx_HOCO_0002=1483
        self.idx_DOCO_0001=1028
        self.idx_DOCO_0002=1484
        self.idx_NH_0001=1029
        self.idx_NH_0002=1485
        self.idx_ND_0001=1030
        self.idx_ND_0002=1486
        self.idx_NH2_0001=1031
        self.idx_NH2_0002=1487
        self.idx_NHD_0001=1032
        self.idx_NHD_0002=1488
        self.idx_ND2_0001=1033
        self.idx_ND2_0002=1489
        self.idx_NH3_0001=1034
        self.idx_NH3_0002=1490
        self.idx_NH2D_0001=1035
        self.idx_NH2D_0002=1491
        self.idx_NHD2_0001=1036
        self.idx_NHD2_0002=1492
        self.idx_ND3_0001=1037
        self.idx_ND3_0002=1493
        self.idx_C10_0001=1038
        self.idx_C10_0002=1494
        self.idx_C10H_0001=1039
        self.idx_C10H_0002=1495
        self.idx_C10H2_0001=1040
        self.idx_C10H2_0002=1496
        self.idx_C10N_0001=1041
        self.idx_C10N_0002=1497
        self.idx_C11_0001=1042
        self.idx_C11_0002=1498
        self.idx_C2H2_0001=1043
        self.idx_C2H2_0002=1499
        self.idx_C2HD_0001=1044
        self.idx_C2HD_0002=1500
        self.idx_C2D2_0001=1045
        self.idx_C2D2_0002=1501
        self.idx_C2H3_0001=1046
        self.idx_C2H3_0002=1502
        self.idx_C2H2D_0001=1047
        self.idx_C2H2D_0002=1503
        self.idx_C2HD2_0001=1048
        self.idx_C2HD2_0002=1504
        self.idx_C2D3_0001=1049
        self.idx_C2D3_0002=1505
        self.idx_C2H4_0001=1050
        self.idx_C2H4_0002=1506
        self.idx_C2H3D_0001=1051
        self.idx_C2H3D_0002=1507
        self.idx_C2H2D2_0001=1052
        self.idx_C2H2D2_0002=1508
        self.idx_C2HD3_0001=1053
        self.idx_C2HD3_0002=1509
        self.idx_C2D4_0001=1054
        self.idx_C2D4_0002=1510
        self.idx_C2H5_0001=1055
        self.idx_C2H5_0002=1511
        self.idx_C2H6_0001=1056
        self.idx_C2H6_0002=1512
        self.idx_C3_0001=1057
        self.idx_C3_0002=1513
        self.idx_C3N_0001=1058
        self.idx_C3N_0002=1514
        self.idx_C3O_0001=1059
        self.idx_C3O_0002=1515
        self.idx_C3S_0001=1060
        self.idx_C3S_0002=1516
        self.idx_C4_0001=1061
        self.idx_C4_0002=1517
        self.idx_C4H_0001=1062
        self.idx_C4H_0002=1518
        self.idx_C4D_0001=1063
        self.idx_C4D_0002=1519
        self.idx_C4H2_0001=1064
        self.idx_C4H2_0002=1520
        self.idx_C4HD_0001=1065
        self.idx_C4HD_0002=1521
        self.idx_C4D2_0001=1066
        self.idx_C4D2_0002=1522
        self.idx_C4H3_0001=1067
        self.idx_C4H3_0002=1523
        self.idx_C4N_0001=1068
        self.idx_C4N_0002=1524
        self.idx_C4S_0001=1069
        self.idx_C4S_0002=1525
        self.idx_C5_0001=1070
        self.idx_C5_0002=1526
        self.idx_C5H_0001=1071
        self.idx_C5H_0002=1527
        self.idx_C5D_0001=1072
        self.idx_C5D_0002=1528
        self.idx_C5H2_0001=1073
        self.idx_C5H2_0002=1529
        self.idx_C5H3_0001=1074
        self.idx_C5H3_0002=1530
        self.idx_C5H4_0001=1075
        self.idx_C5H4_0002=1531
        self.idx_C5N_0001=1076
        self.idx_C5N_0002=1532
        self.idx_C5O_0001=1077
        self.idx_C5O_0002=1533
        self.idx_C6_0001=1078
        self.idx_C6_0002=1534
        self.idx_C6H_0001=1079
        self.idx_C6H_0002=1535
        self.idx_C6H2_0001=1080
        self.idx_C6H2_0002=1536
        self.idx_C6H3_0001=1081
        self.idx_C6H3_0002=1537
        self.idx_C6H4_0001=1082
        self.idx_C6H4_0002=1538
        self.idx_C6H6_0001=1083
        self.idx_C6H6_0002=1539
        self.idx_C6N_0001=1084
        self.idx_C6N_0002=1540
        self.idx_C7_0001=1085
        self.idx_C7_0002=1541
        self.idx_C7H_0001=1086
        self.idx_C7H_0002=1542
        self.idx_C7H2_0001=1087
        self.idx_C7H2_0002=1543
        self.idx_C7H3_0001=1088
        self.idx_C7H3_0002=1544
        self.idx_C7H4_0001=1089
        self.idx_C7H4_0002=1545
        self.idx_C7N_0001=1090
        self.idx_C7N_0002=1546
        self.idx_C7O_0001=1091
        self.idx_C7O_0002=1547
        self.idx_C8_0001=1092
        self.idx_C8_0002=1548
        self.idx_C8H_0001=1093
        self.idx_C8H_0002=1549
        self.idx_C8H2_0001=1094
        self.idx_C8H2_0002=1550
        self.idx_C8H3_0001=1095
        self.idx_C8H3_0002=1551
        self.idx_C8H4_0001=1096
        self.idx_C8H4_0002=1552
        self.idx_C8N_0001=1097
        self.idx_C8N_0002=1553
        self.idx_C9_0001=1098
        self.idx_C9_0002=1554
        self.idx_C9H_0001=1099
        self.idx_C9H_0002=1555
        self.idx_C9H2_0001=1100
        self.idx_C9H2_0002=1556
        self.idx_C9H3_0001=1101
        self.idx_C9H3_0002=1557
        self.idx_C9H4_0001=1102
        self.idx_C9H4_0002=1558
        self.idx_C9N_0001=1103
        self.idx_C9N_0002=1559
        self.idx_C9O_0001=1104
        self.idx_C9O_0002=1560
        self.idx_CCH_0001=1105
        self.idx_CCH_0002=1561
        self.idx_CCD_0001=1106
        self.idx_CCD_0002=1562
        self.idx_CCN_0001=1107
        self.idx_CCN_0002=1563
        self.idx_CCO_0001=1108
        self.idx_CCO_0002=1564
        self.idx_CCS_0001=1109
        self.idx_CCS_0002=1565
        self.idx_CD3CN_0001=1110
        self.idx_CD3CN_0002=1566
        self.idx_CH2CCH_0001=1111
        self.idx_CH2CCH_0002=1567
        self.idx_CH2CCD_0001=1112
        self.idx_CH2CCD_0002=1568
        self.idx_CHDCCH_0001=1113
        self.idx_CHDCCH_0002=1569
        self.idx_CHDCCD_0001=1114
        self.idx_CHDCCD_0002=1570
        self.idx_CD2CCH_0001=1115
        self.idx_CD2CCH_0002=1571
        self.idx_CD2CCD_0001=1116
        self.idx_CD2CCD_0002=1572
        self.idx_CH2CHC2H_0001=1117
        self.idx_CH2CHC2H_0002=1573
        self.idx_CH2CHCHCH2_0001=1118
        self.idx_CH2CHCHCH2_0002=1574
        self.idx_CH2CHCN_0001=1119
        self.idx_CH2CHCN_0002=1575
        self.idx_CH2NH_0001=1120
        self.idx_CH2NH_0002=1576
        self.idx_CHDNH_0001=1121
        self.idx_CHDNH_0002=1577
        self.idx_CHDND_0001=1122
        self.idx_CHDND_0002=1578
        self.idx_CH2ND_0001=1123
        self.idx_CH2ND_0002=1579
        self.idx_CD2NH_0001=1124
        self.idx_CD2NH_0002=1580
        self.idx_CD2ND_0001=1125
        self.idx_CD2ND_0002=1581
        self.idx_CH2NH2_0001=1126
        self.idx_CH2NH2_0002=1582
        self.idx_CH2NHD_0001=1127
        self.idx_CH2NHD_0002=1583
        self.idx_CHDNH2_0001=1128
        self.idx_CHDNH2_0002=1584
        self.idx_CHDNHD_0001=1129
        self.idx_CHDNHD_0002=1585
        self.idx_CHDND2_0001=1130
        self.idx_CHDND2_0002=1586
        self.idx_CH2ND2_0001=1131
        self.idx_CH2ND2_0002=1587
        self.idx_CD2NH2_0001=1132
        self.idx_CD2NH2_0002=1588
        self.idx_CD2NHD_0001=1133
        self.idx_CD2NHD_0002=1589
        self.idx_CD2ND2_0001=1134
        self.idx_CD2ND2_0002=1590
        self.idx_CH3C3N_0001=1135
        self.idx_CH3C3N_0002=1591
        self.idx_CH3C4H_0001=1136
        self.idx_CH3C4H_0002=1592
        self.idx_CH3C5N_0001=1137
        self.idx_CH3C5N_0002=1593
        self.idx_CH3C6H_0001=1138
        self.idx_CH3C6H_0002=1594
        self.idx_CH3C7N_0001=1139
        self.idx_CH3C7N_0002=1595
        self.idx_CH3CCH_0001=1140
        self.idx_CH3CCH_0002=1596
        self.idx_CH3CH2OH_0001=1141
        self.idx_CH3CH2OH_0002=1597
        self.idx_CH3CHCH2_0001=1142
        self.idx_CH3CHCH2_0002=1598
        self.idx_CH3CHO_0001=1143
        self.idx_CH3CHO_0002=1599
        self.idx_CH3CN_0001=1144
        self.idx_CH3CN_0002=1600
        self.idx_CH2DCN_0001=1145
        self.idx_CH2DCN_0002=1601
        self.idx_CHD2CN_0001=1146
        self.idx_CHD2CN_0002=1602
        self.idx_CH3COCH3_0001=1147
        self.idx_CH3COCH3_0002=1603
        self.idx_CH3NH2_0001=1148
        self.idx_CH3NH2_0002=1604
        self.idx_CH3OCH2_0001=1149
        self.idx_CH3OCH2_0002=1605
        self.idx_CH3OCH3_0001=1150
        self.idx_CH3OCH3_0002=1606
        self.idx_CN_0001=1151
        self.idx_CN_0002=1607
        self.idx_CS_0001=1152
        self.idx_CS_0002=1608
        self.idx_H2CCN_0001=1153
        self.idx_H2CCN_0002=1609
        self.idx_HDCCN_0001=1154
        self.idx_HDCCN_0002=1610
        self.idx_D2CCN_0001=1155
        self.idx_D2CCN_0002=1611
        self.idx_H2CCO_0001=1156
        self.idx_H2CCO_0002=1612
        self.idx_HDCCO_0001=1157
        self.idx_HDCCO_0002=1613
        self.idx_D2CCO_0001=1158
        self.idx_D2CCO_0002=1614
        self.idx_H2CN_0001=1159
        self.idx_H2CN_0002=1615
        self.idx_HDCN_0001=1160
        self.idx_HDCN_0002=1616
        self.idx_D2CN_0001=1161
        self.idx_D2CN_0002=1617
        self.idx_H2CS_0001=1162
        self.idx_H2CS_0002=1618
        self.idx_HDCS_0001=1163
        self.idx_HDCS_0002=1619
        self.idx_D2CS_0001=1164
        self.idx_D2CS_0002=1620
        self.idx_H2S_0001=1165
        self.idx_H2S_0002=1621
        self.idx_HDS_0001=1166
        self.idx_HDS_0002=1622
        self.idx_D2S_0001=1167
        self.idx_D2S_0002=1623
        self.idx_HC2O_0001=1168
        self.idx_HC2O_0002=1624
        self.idx_DC2O_0001=1169
        self.idx_DC2O_0002=1625
        self.idx_HC3N_0001=1170
        self.idx_HC3N_0002=1626
        self.idx_DC3N_0001=1171
        self.idx_DC3N_0002=1627
        self.idx_HC4N_0001=1172
        self.idx_HC4N_0002=1628
        self.idx_DC4N_0001=1173
        self.idx_DC4N_0002=1629
        self.idx_HC5N_0001=1174
        self.idx_HC5N_0002=1630
        self.idx_HC6N_0001=1175
        self.idx_HC6N_0002=1631
        self.idx_HC7N_0001=1176
        self.idx_HC7N_0002=1632
        self.idx_HC8N_0001=1177
        self.idx_HC8N_0002=1633
        self.idx_HC9N_0001=1178
        self.idx_HC9N_0002=1634
        self.idx_HCCNC_0001=1179
        self.idx_HCCNC_0002=1635
        self.idx_DCCNC_0001=1180
        self.idx_DCCNC_0002=1636
        self.idx_HCN_0001=1181
        self.idx_HCN_0002=1637
        self.idx_DCN_0001=1182
        self.idx_DCN_0002=1638
        self.idx_HCNCC_0001=1183
        self.idx_HCNCC_0002=1639
        self.idx_DCNCC_0001=1184
        self.idx_DCNCC_0002=1640
        self.idx_HCOOCH3_0001=1185
        self.idx_HCOOCH3_0002=1641
        self.idx_HCS_0001=1186
        self.idx_HCS_0002=1642
        self.idx_DCS_0001=1187
        self.idx_DCS_0002=1643
        self.idx_HNC_0001=1188
        self.idx_HNC_0002=1644
        self.idx_DNC_0001=1189
        self.idx_DNC_0002=1645
        self.idx_HNCCC_0001=1190
        self.idx_HNCCC_0002=1646
        self.idx_DNCCC_0001=1191
        self.idx_DNCCC_0002=1647
        self.idx_HNCO_0001=1192
        self.idx_HNCO_0002=1648
        self.idx_DNCO_0001=1193
        self.idx_DNCO_0002=1649
        self.idx_HNO_0001=1194
        self.idx_HNO_0002=1650
        self.idx_DNO_0001=1195
        self.idx_DNO_0002=1651
        self.idx_HS_0001=1196
        self.idx_HS_0002=1652
        self.idx_DS_0001=1197
        self.idx_DS_0002=1653
        self.idx_N2_0001=1198
        self.idx_N2_0002=1654
        self.idx_N2O_0001=1199
        self.idx_N2O_0002=1655
        self.idx_NC4N_0001=1200
        self.idx_NC4N_0002=1656
        self.idx_NC6N_0001=1201
        self.idx_NC6N_0002=1657
        self.idx_NC8N_0001=1202
        self.idx_NC8N_0002=1658
        self.idx_NH2CHO_0001=1203
        self.idx_NH2CHO_0002=1659
        self.idx_NH2CDO_0001=1204
        self.idx_NH2CDO_0002=1660
        self.idx_NHDCHO_0001=1205
        self.idx_NHDCHO_0002=1661
        self.idx_NHDCDO_0001=1206
        self.idx_NHDCDO_0002=1662
        self.idx_ND2CHO_0001=1207
        self.idx_ND2CHO_0002=1663
        self.idx_ND2CDO_0001=1208
        self.idx_ND2CDO_0002=1664
        self.idx_NH2CN_0001=1209
        self.idx_NH2CN_0002=1665
        self.idx_NHDCN_0001=1210
        self.idx_NHDCN_0002=1666
        self.idx_ND2CN_0001=1211
        self.idx_ND2CN_0002=1667
        self.idx_HSO_0001=1212
        self.idx_HSO_0002=1668
        self.idx_DSO_0001=1213
        self.idx_DSO_0002=1669
        self.idx_HSS_0001=1214
        self.idx_HSS_0002=1670
        self.idx_DSS_0001=1215
        self.idx_DSS_0002=1671
        self.idx_HSSH_0001=1216
        self.idx_HSSH_0002=1672
        self.idx_HSSD_0001=1217
        self.idx_HSSD_0002=1673
        self.idx_DSSH_0001=1218
        self.idx_DSSH_0002=1674
        self.idx_DSSD_0001=1219
        self.idx_DSSD_0002=1675
        self.idx_NO_0001=1220
        self.idx_NO_0002=1676
        self.idx_NO2_0001=1221
        self.idx_NO2_0002=1677
        self.idx_NS_0001=1222
        self.idx_NS_0002=1678
        self.idx_OCN_0001=1223
        self.idx_OCN_0002=1679
        self.idx_OCS_0001=1224
        self.idx_OCS_0002=1680
        self.idx_S2_0001=1225
        self.idx_S2_0002=1681
        self.idx_SO_0001=1226
        self.idx_SO_0002=1682
        self.idx_SO2_0001=1227
        self.idx_SO2_0002=1683
        self.idx_CH3OCHO_0001=1228
        self.idx_CH3OCHO_0002=1684
        self.idx_Si_0001=1229
        self.idx_Si_0002=1685
        self.idx_SiS_0001=1230
        self.idx_SiS_0002=1686
        self.idx_SiN_0001=1231
        self.idx_SiN_0002=1687
        self.idx_SiC_0001=1232
        self.idx_SiC_0002=1688
        self.idx_SiH_0001=1233
        self.idx_SiH_0002=1689
        self.idx_SiH2_0001=1234
        self.idx_SiH2_0002=1690
        self.idx_SiH3_0001=1235
        self.idx_SiH3_0002=1691
        self.idx_SiH4_0001=1236
        self.idx_SiH4_0002=1692
        self.idx_SiC2CH3_0001=1237
        self.idx_SiC2CH3_0002=1693
        self.idx_SiC3H_0001=1238
        self.idx_SiC3H_0002=1694
        self.idx_SiC3H5_0001=1239
        self.idx_SiC3H5_0002=1695
        self.idx_SiC4_0001=1240
        self.idx_SiC4_0002=1696
        self.idx_SiC4H_0001=1241
        self.idx_SiC4H_0002=1697
        self.idx_SiC6H_0001=1242
        self.idx_SiC6H_0002=1698
        self.idx_SiC8H_0001=1243
        self.idx_SiC8H_0002=1699
        self.idx_c_HCCHSi_0001=1244
        self.idx_c_HCCHSi_0002=1700
        self.idx_c_SiC2_0001=1245
        self.idx_c_SiC2_0002=1701
        self.idx_l_C3H_0001=1246
        self.idx_l_C3H_0002=1702
        self.idx_l_C3D_0001=1247
        self.idx_l_C3D_0002=1703
        self.idx_c_C3H_0001=1248
        self.idx_c_C3H_0002=1704
        self.idx_c_C3D_0001=1249
        self.idx_c_C3D_0002=1705
        self.idx_l_SiC3_0001=1250
        self.idx_l_SiC3_0002=1706
        self.idx_l_C3H2_0001=1251
        self.idx_l_C3H2_0002=1707
        self.idx_l_C3HD_0001=1252
        self.idx_l_C3HD_0002=1708
        self.idx_l_C3D2_0001=1253
        self.idx_l_C3D2_0002=1709
        self.idx_c_C3H2_0001=1254
        self.idx_c_C3H2_0002=1710
        self.idx_c_C3HD_0001=1255
        self.idx_c_C3HD_0002=1711
        self.idx_c_C3D2_0001=1256
        self.idx_c_C3D2_0002=1712
        self.idx_Mg_0001=1257
        self.idx_Mg_0002=1713
        self.idx_MgH_0001=1258
        self.idx_MgH_0002=1714
        self.idx_MgH2_0001=1259
        self.idx_MgH2_0002=1715
        self.idx_Na_0001=1260
        self.idx_Na_0002=1716
        self.idx_NaH_0001=1261
        self.idx_NaH_0002=1717
        self.idx_F_0001=1262
        self.idx_F_0002=1718
        self.idx_HF_0001=1263
        self.idx_HF_0002=1719
        self.idx_DF_0001=1264
        self.idx_DF_0002=1720
        self.idx_MgD_0001=1265
        self.idx_MgD_0002=1721
        self.idx_MgHD_0001=1266
        self.idx_MgHD_0002=1722
        self.idx_MgD2_0001=1267
        self.idx_MgD2_0002=1723
        self.idx_NaD_0001=1268
        self.idx_NaD_0002=1724
        self.idx_SiD_0001=1269
        self.idx_SiD_0002=1725
        self.idx_Cl_0001=1270
        self.idx_Cl_0002=1726
        self.idx_HCl_0001=1271
        self.idx_HCl_0002=1727
        self.idx_DCl_0001=1272
        self.idx_DCl_0002=1728
        self.idx_H2Cl_0001=1273
        self.idx_H2Cl_0002=1729
        self.idx_HDCl_0001=1274
        self.idx_HDCl_0002=1730
        self.idx_D2Cl_0001=1275
        self.idx_D2Cl_0002=1731
        self.idx_CCl_0001=1276
        self.idx_CCl_0002=1732
        self.idx_ClO_0001=1277
        self.idx_ClO_0002=1733
        self.idx_C3H3_0001=1278
        self.idx_C3H3_0002=1734
        self.idx_C3H2D_0001=1279
        self.idx_C3H2D_0002=1735
        self.idx_C3HD2_0001=1280
        self.idx_C3HD2_0002=1736
        self.idx_C3D3_0001=1281
        self.idx_C3D3_0002=1737
        self.idx_C3H4_0001=1282
        self.idx_C3H4_0002=1738
        self.idx_C3H3D_0001=1283
        self.idx_C3H3D_0002=1739
        self.idx_C3H2D2_0001=1284
        self.idx_C3H2D2_0002=1740
        self.idx_HCCN_0001=1285
        self.idx_HCCN_0002=1741
        self.idx_DCCN_0001=1286
        self.idx_DCCN_0002=1742
        self.idx_C2H2N_0001=1287
        self.idx_C2H2N_0002=1743
        self.idx_C2H5OH_0001=1288
        self.idx_C2H5OH_0002=1744
        self.idx_C2H3D3_0001=1289
        self.idx_C2H3D3_0002=1745
        self.idx_C2H5D_0001=1290
        self.idx_C2H5D_0002=1746
        self.idx_C2H4D2_0001=1291
        self.idx_C2H4D2_0002=1747
        self.idx_C2H3N_0001=1292
        self.idx_C2H3N_0002=1748
        self.idx_C2H4O_0001=1293
        self.idx_C2H4O_0002=1749
        self.idx_CH2DCHO_0001=1294
        self.idx_CH2DCHO_0002=1750
        self.idx_CH3CDO_0001=1295
        self.idx_CH3CDO_0002=1751
        self.idx_CD3CHO_0001=1296
        self.idx_CD3CHO_0002=1752
        self.idx_CHD2CDO_0001=1297
        self.idx_CHD2CDO_0002=1753
        self.idx_CHD2CHO_0001=1298
        self.idx_CHD2CHO_0002=1754
        self.idx_CH2DCDO_0001=1299
        self.idx_CH2DCDO_0002=1755
        self.idx_C2H4D_0001=1300
        self.idx_C2H4D_0002=1756
        self.idx_C2H2D3_0001=1301
        self.idx_C2H2D3_0002=1757
        self.idx_C2H3D2_0001=1302
        self.idx_C2H3D2_0002=1758
        self.idx_C3H3N_0001=1303
        self.idx_C3H3N_0002=1759
        self.idx_H4C3N_0001=1304
        self.idx_H4C3N_0002=1760
        self.idx_C3D3N_0001=1305
        self.idx_C3D3N_0002=1761
        self.idx_HD3C3N_0001=1306
        self.idx_HD3C3N_0002=1762
        self.idx_C3H2DN_0001=1307
        self.idx_C3H2DN_0002=1763
        self.idx_H3DC3N_0001=1308
        self.idx_H3DC3N_0002=1764
        self.idx_C3HD2N_0001=1309
        self.idx_C3HD2N_0002=1765
        self.idx_H2D2C3N_0001=1310
        self.idx_H2D2C3N_0002=1766
        self.idx_HC3O_0001=1311
        self.idx_HC3O_0002=1767
        self.idx_DC3O_0001=1312
        self.idx_DC3O_0002=1768
        self.idx_C4H4_0001=1313
        self.idx_C4H4_0002=1769
        self.idx_CH5N_0001=1314
        self.idx_CH5N_0002=1770
        self.idx_CH3NH_0001=1315
        self.idx_CH3NH_0002=1771
        self.idx_FeD_0001=1316
        self.idx_FeD_0002=1772
        self.idx_H2C3N_0001=1317
        self.idx_H2C3N_0002=1773
        self.idx_H2C5N_0001=1318
        self.idx_H2C5N_0002=1774
        self.idx_H3C5N_0001=1319
        self.idx_H3C5N_0002=1775
        self.idx_H2C7N_0001=1320
        self.idx_H2C7N_0002=1776
        self.idx_H3C7N_0001=1321
        self.idx_H3C7N_0002=1777
        self.idx_H2C9N_0001=1322
        self.idx_H2C9N_0002=1778
        self.idx_H3C9N_0001=1323
        self.idx_H3C9N_0002=1779
        self.idx_H5C3N_0001=1324
        self.idx_H5C3N_0002=1780
        self.idx_C2H2O_0001=1325
        self.idx_C2H2O_0002=1781
        self.idx_C2HDO_0001=1326
        self.idx_C2HDO_0002=1782
        self.idx_C2D2O_0001=1327
        self.idx_C2D2O_0002=1783
        self.idx_HDC3N_0001=1328
        self.idx_HDC3N_0002=1784
        self.idx_D2C3N_0001=1329
        self.idx_D2C3N_0002=1785
        self.idx_H2C3O_0001=1330
        self.idx_H2C3O_0002=1786
        self.idx_HDC3O_0001=1331
        self.idx_HDC3O_0002=1787
        self.idx_D2C3O_0001=1332
        self.idx_D2C3O_0002=1788
        self.idx_C2HDN_0001=1333
        self.idx_C2HDN_0002=1789
        self.idx_C2D2N_0001=1334
        self.idx_C2D2N_0002=1790
        self.idx_CH3N_0001=1335
        self.idx_CH3N_0002=1791
        self.idx_N2H2_0001=1336
        self.idx_N2H2_0002=1792
        self.idx_N2HD_0001=1337
        self.idx_N2HD_0002=1793
        self.idx_N2D2_0001=1338
        self.idx_N2D2_0002=1794
        self.idx_NH2OH_0001=1339
        self.idx_NH2OH_0002=1795
        self.idx_NH2OD_0001=1340
        self.idx_NH2OD_0002=1796
        self.idx_N2H_0001=1341
        self.idx_N2H_0002=1797
        self.idx_N2D_0001=1342
        self.idx_N2D_0002=1798
        self.idx_CNH2_0001=1343
        self.idx_CNH2_0002=1799
        self.idx_CHNH2_0001=1344
        self.idx_CHNH2_0002=1800
        self.idx_HON_0001=1345
        self.idx_HON_0002=1801
        self.idx_DON_0001=1346
        self.idx_DON_0002=1802
        self.idx_NHNO_0001=1347
        self.idx_NHNO_0002=1803
        self.idx_CH2DN_0001=1348
        self.idx_CH2DN_0002=1804
        self.idx_CHD2N_0001=1349
        self.idx_CHD2N_0002=1805
        self.idx_CD3N_0001=1350
        self.idx_CD3N_0002=1806
        self.idx_NH2NO_0001=1351
        self.idx_NH2NO_0002=1807
        self.idx_CH3CO_0001=1352
        self.idx_CH3CO_0002=1808
        self.idx_CH2DOCH3_0001=1353
        self.idx_CH2DOCH3_0002=1809
        self.idx_CHD2OCH3_0001=1354
        self.idx_CHD2OCH3_0002=1810
        self.idx_CH2DOCH2D_0001=1355
        self.idx_CH2DOCH2D_0002=1811
        self.idx_CD3OCH3_0001=1356
        self.idx_CD3OCH3_0002=1812
        self.idx_CHD2OCH2D_0001=1357
        self.idx_CHD2OCH2D_0002=1813
        self.idx_DCOOCH3_0001=1358
        self.idx_DCOOCH3_0002=1814
        self.idx_HCOOCH2D_0001=1359
        self.idx_HCOOCH2D_0002=1815
        self.idx_DCOOCH2D_0001=1360
        self.idx_DCOOCH2D_0002=1816
        self.idx_HCOOCHD2_0001=1361
        self.idx_HCOOCHD2_0002=1817
        self.idx_DCOOCHD2_0001=1362
        self.idx_DCOOCHD2_0002=1818
        self.idx_HCOOCD3_0001=1363
        self.idx_HCOOCD3_0002=1819
        self.idx_DCOOCD3_0001=1364
        self.idx_DCOOCD3_0002=1820
        self.idx_NH2CO_0001=1365
        self.idx_NH2CO_0002=1821
        self.idx_NHDCO_0001=1366
        self.idx_NHDCO_0002=1822
        self.idx_ND2CO_0001=1367
        self.idx_ND2CO_0002=1823
        self.idx_CH2DOCHO_0001=1368
        self.idx_CH2DOCHO_0002=1824
        self.idx_CH3OCDO_0001=1369
        self.idx_CH3OCDO_0002=1825
        self.idx_CHD2OCHO_0001=1370
        self.idx_CHD2OCHO_0002=1826
        self.idx_CH2DOCDO_0001=1371
        self.idx_CH2DOCDO_0002=1827
        self.idx_CD3OCHO_0001=1372
        self.idx_CD3OCHO_0002=1828
        self.idx_CHD2OCDO_0001=1373
        self.idx_CHD2OCDO_0002=1829
        self.idx_CHDOCDO_0001=1374
        self.idx_CHDOCDO_0002=1830
        self.idx_CHDOCHDO_0001=1375
        self.idx_CHDOCHDO_0002=1831
        self.idx_CH3OCHD_0001=1376
        self.idx_CH3OCHD_0002=1832
        self.idx_CH2DOCH2_0001=1377
        self.idx_CH2DOCH2_0002=1833
        self.idx_CH2DOCHD_0001=1378
        self.idx_CH2DOCHD_0002=1834
        self.idx_CHD2OCH2_0001=1379
        self.idx_CHD2OCH2_0002=1835
        self.idx_CH3OCD2_0001=1380
        self.idx_CH3OCD2_0002=1836
        self.idx_CH2DOCD2_0001=1381
        self.idx_CH2DOCD2_0002=1837
        self.idx_CH3OCHD2_0001=1382
        self.idx_CH3OCHD2_0002=1838
        self.idx_CH3OCH2D_0001=1383
        self.idx_CH3OCH2D_0002=1839
        self.idx_CH2DCO_0001=1384
        self.idx_CH2DCO_0002=1840
        self.idx_CH3NHD_0001=1385
        self.idx_CH3NHD_0002=1841
        self.idx_CH2DNH2_0001=1386
        self.idx_CH2DNH2_0002=1842
        self.idx_P_0001=1387
        self.idx_P_0002=1843
        self.idx_PO_0001=1388
        self.idx_PO_0002=1844
        self.idx_PH_0001=1389
        self.idx_PH_0002=1845
        self.idx_PD_0001=1390
        self.idx_PD_0002=1846
        self.idx_PH2_0001=1391
        self.idx_PH2_0002=1847
        self.idx_PHD_0001=1392
        self.idx_PHD_0002=1848
        self.idx_PD2_0001=1393
        self.idx_PD2_0002=1849
        self.idx_PN_0001=1394
        self.idx_PN_0002=1850
        self.idx_CP_0001=1395
        self.idx_CP_0002=1851
        self.idx_CCP_0001=1396
        self.idx_CCP_0002=1852
        self.idx_C3P_0001=1397
        self.idx_C3P_0002=1853
        self.idx_C4P_0001=1398
        self.idx_C4P_0002=1854
        self.idx_CH2PH_0001=1399
        self.idx_CH2PH_0002=1855
        self.idx_CHD2PH_0001=1400
        self.idx_CHD2PH_0002=1856
        self.idx_CH2PD_0001=1401
        self.idx_CH2PD_0002=1857
        self.idx_CHDPD_0001=1402
        self.idx_CHDPD_0002=1858
        self.idx_CD2P2_0001=1403
        self.idx_CD2P2_0002=1859
        self.idx_HCP_0001=1404
        self.idx_HCP_0002=1860
        self.idx_DCP_0001=1405
        self.idx_DCP_0002=1861
        self.idx_HCCP_0001=1406
        self.idx_HCCP_0002=1862
        self.idx_DCCP_0001=1407
        self.idx_DCCP_0002=1863
        self.idx_HPO_0001=1408
        self.idx_HPO_0002=1864
        self.idx_DPO_0001=1409
        self.idx_DPO_0002=1865
        self.idx_surface_mask=1410
        self.idx_mantle_mask=1411
        self.idx_dummy=1412

        # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        ##END_IDXLIST

        self.lib.kemimo_dumpFluxes.restype = None
        self.lib.kemimo_dumpFluxes.argtypes = [
            array_1d_double, int_byref, dble_byref]
        self.lib.kemimo_loadVerbatim.restype = None
        self.lib.kemimo_loadVerbatim.argtypes = None
        self.lib.kemimo_computeRates.restype = None
        self.lib.kemimo_computeRates.argtypes = [array_1d_double,
            dble_byref, dble_byref, dble_byref, dble_byref]
        self.lib.kemimo_dochem.restype = None
        self.lib.kemimo_dochem.argtypes = [array_1d_double, dble_byref]
        self.lib.kemimo_theta_H2.restype = ctypes.c_double
        self.lib.kemimo_theta_H2.argtypes = [dble_byref, dble_byref]
        self.lib.kemimo_get_H2_idx.restype = None
        self.lib.kemimo_get_H2_idx.argtypes = [int_byref, int_byref, int_byref]
        # self.lib.kemimo_returnFluxes.restype = None
        # self.lib.kemimo_returnFluxes.argtypes = [
        #     array_1d_double, int_byref, array_pointer, array_pointer]
