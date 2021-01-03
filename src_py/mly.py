from utils import isNumber
import sys
import numpy as np


# molecule class, i.e. species class
class mly:
    # *******************
    # species class constructor for mly mask
    # 
    def __init__(self, seperationDistance, radius, layerThickness, name='surface'):
        self.layer = 0
        # name
        self.name = name + "_mask"
        self.dictname = name + "_mask"
        self.nameLatex =  "$f_{\mathrm{mask}}$"
        self.namebase = "mask"
        # binding energy, K
        self.Eice = self.Ebare = None
        self.isGas = False
        # enthalpy of formation, K
        self.dH = None
        # mass, amu
        self.mass = None
        # idx
        self.idx = self.idxGas = self.idxTot = None
        self.fidx = "idx_" + self.name.replace("+", "j")

        # Calculate prefactor (number of sites per grain pr monolayer)
        #self.prefactor = seperationDistance**2 / (4.0 * np.pi * radius**2)
        self.prefactor = 1e0 / (1.5e15 * 4.0 * np.pi * radius**2)

# molecule class, i.e. species class
class mly_dist:
    # *******************
    # species class constructor for mly mask WITH dust distribution
    #
    def __init__(self, seperationDistance, ngas, radius_min, radius_max):
        self.layer = 0
        # name
        self.name = "mask"
        self.dictname = "mask"
        self.nameLatex = "$f_{\mathrm{mask}}$"
        # binding energy, K
        self.Eice = self.Ebare = None
        self.isGas = False
        # enthalpy of formation, K
        self.dH = None
        # mass, amu
        self.mass = None
        # idx
        self.idx = self.idxGas = self.idxTot = None
        self.fidx = "idx_" + self.name.replace("+", "j")

        # Calculate prefactor (number of sites per cm^3)
        d2g = 1e-2  # dust/gas mass ratio
        rho0 = 3.0  # bulk density, g/cm3
        mu = 2.34 # mean molecular weight
        pmass = 1.6726219e-24 # proton mass, g
        rhod = ngas*pmass*mu*d2g # dust mass density, g/cm3
        pexp = -3.5 # size distribution exponent

        p3 = pexp + 3.0
        p4 = pexp + 4.0
        
        app = seperationDistance
        amax = radius_max
        amin = radius_min

        self.prefactor = 3.0 * rhod / app**2 / rho0 * (p4 / p3) * (amax**p3-amin**p3) / (amax**p4-amin**p4)
