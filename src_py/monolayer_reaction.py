from math import exp, sqrt, pi
from utils import strF90
import sys


# reaction class for monolayer mask
class monolayer_reaction:
    # *******************
    # monolayer class constructor

    def __init__(self, mly):

        # prepare reaction verbatim
        self.verbatim = " monolayer reaction "


        # reaction type, see findType method
        self.type ="dummy"
        # reaction index, zero-based
        self.idx = None
        # dummies
        self.dH = 0.
        self.Ea = 0.
        self.yieldPD = 0.
        self.Pdelta = 0.
        self.barrierFromFile = False
        self.Tmin = 0.0
        self.Tmax = 0.0
        self.hasMultipleTranges = False

        #
        self.prefactor = mly.prefactor
        self.krateF90 = "%.16f / (ngas*xdust)" % self.prefactor
        #
        self.reactants = []
        self.products = []

        #
        self.layer = -1


    # ********************
    # prepare the RHS in F90 format and store to attribute
    def buildF90RHS(self, reactions):
        self.RHS = []
