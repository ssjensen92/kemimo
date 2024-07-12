from math import pi


# molecule class, i.e. species class
class mly:
    # *******************
    # species class constructor for mly mask
    # 
    def __init__(self, seperationDistance, radius, name='surface'):
        self.layer = 0
        self.name = name + "_mask"
        self.dictname = name + "_mask"
        self.nameLatex = "$f_{\mathrm{mask}}$"
        self.namebase = "mask"

        self.Eice = self.Ebare = None
        self.isGas = False
        self.dH = None
        self.mass = None
        self.idx = self.idxGas = self.idxTot = None
        self.fidx = "idx_" + self.name.replace("+", "j")

        # calculate the prefactor
        self.prefactor = seperationDistance**2 / (4.0 * pi * radius**2)