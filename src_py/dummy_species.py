from utils import isNumber
import sys
import numpy as np


# molecule class, i.e. species class
class dummy_species:
    # *******************
    # species class constructor for mly mask
    # 
    def __init__(self):
        self.layer = 0
        # name
        self.name = "dummy"
        self.dictname = "dummy"
        self.namebase = "dummy"
        self.nameLatex =  "$f_{\mathrm{dummy}}$"
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

