from math import exp, sqrt, pi
from utils import strF90
import numpy as np
import sys
import re


# reaction class for surface reactions
class reaction:
    # *******************
    # reaction class constructor
    # RRs is list of reactants, list of mol objects
    # PPs is list of products, list of mol objects
    # Ealist is list of reactions to get energy barrier, K
    #  reactions are dictionaries
    # yieldPD is photodissociation yield
    def __init__(self, RRs, PPs, EaList, barrierWidths, Bratios, yieldPD=0e0, gamma=0e0, alpha=None):

        # KIDA photoreactions Av_variable name:
        self.Avvar = "variable_Av"
        self.alpha = alpha

        # list of reactants, list of species objects
        self.reactants = RRs
        # list of products, list of species objects
        self.products = PPs
        # store photodissociation yield
        self.yieldPD = yieldPD
        # KIDA gas-phase photodissocation factor gamma:
        self.gamma = gamma
        # Branching ratio:
        self.Bratio = None
        # compute enthalpy of formation (K), exothermic is positive
        try:
            self.dH = sum([x.dH for x in RRs]) - sum([x.dH for x in PPs])
        except:
            self.dH = np.nan

        # prepare reaction verbatim, e.g. OH + H2 -> H2O + H
        self.verbatim = " + ".join([x.name for x in RRs])
        self.verbatim += " -> "
        self.verbatim += " + ".join([x.name for x in PPs])

        # reaction type, see findType method
        self.type = None
        # reaction index, zero-based
        self.idx = None
        # flag if reaction has dust-phase reactants
        self.hasDustReactants = all([not x.isGas for x in RRs])
        # flag if reaction has gas-phase products
        self.hasGasProducts = any([x.isGas for x in RRs])
        # rate in F90 format
        self.krateF90 = None
        # default temperature ranges, K
        self.Tmin = -1e99
        self.Tmax = 1e99

        # reaction RHS default
        self.RHS = None

        # get barrier (K) from data, default is no barrier
        self.Ea = 0e0
        self.EaList = EaList

        # flag if barrier found from data file
        self.barrierFromFile = False

        # loop on reactions with known barrier
        for rea in EaList:
            # use namebase instead of name to include _gas species
            # condition 1: self reaction and database reaction have the same reactants
            hasReactants = (sorted(rea["reactants"])
                            == sorted([x.namebase for x in RRs]))
            # condition 2: self reaction and database reaction have the same products
            hasProducts = (sorted(rea["products"]) ==
                           sorted([x.namebase for x in PPs]))
            # if both conditions are verified, then are the same reaction
            if hasReactants and hasProducts:
                # store barrier energy
                self.Ea = rea["Ea"]
                # barrier is loaded form file
                self.barrierFromFile = True

        # flag if barrier width found from data file
        self.barrierWidthFromFile = False

        # loop on reactions with custom barrier width
        for rea in barrierWidths:
            # use namebase instead of name to include _gas species
            # condition 1: self reaction and database reaction have the same reactants
            hasReactants = (sorted(rea["reactants"])
                            == sorted([x.namebase for x in RRs]))
            # condition 2: self reaction and database reaction have the same products
            hasProducts = (sorted(rea["products"]) ==
                           sorted([x.namebase for x in PPs]))
            # if both conditions are verified, then are the same reaction
            if hasReactants and hasProducts:
                # store barrier energy
                self.barrierWidth = rea["width"]
                # barrier is loaded form file
                self.barrierWidthFromFile = True

        # flag if branching ratio from data file:
        self.bRatioFromFile = False
        # loop on reactions with known branching ratios
        if len(Bratios) > 0:
            for rea in Bratios:
                # use namebase instead of name to include _gas species
                # condition 1: self reaction and database reaction have the same reactants
                hasReactants = (sorted(rea["reactants"])
                                == sorted([x.namebase for x in RRs]))
                # condition 2: self reaction and database reaction have the same products
                hasProducts = (sorted(rea["products"]) ==
                            sorted([x.namebase for x in PPs]))
                # if both conditions are verified, then are the same reaction
                if hasReactants and hasProducts:
                    # store barrier energy
                    self.Bratio = rea["ratio"]
                    self.bRatioFromFile = True

        # prepare hash, i.e. H2_OH__H_H2O (note double underscore)
        self.hash = "_".join(sorted([x.name for x in self.reactants]))
        self.hash += "__"
        self.hash += "_".join(sorted([x.name for x in self.products]))

        # get exploded hash
        self.hashExploded = self.getHashExploded()

        # find reaction type
        self.findType()

        # if np.isnan(self.dH) and self.type.startswith('2body'):
        #     raise ValueError

        # Layer index
        try:
            if self.type == 'freezeout':
                self.layer = self.products[0].layer
            elif self.type == 'evaporation' or self.type == 'CRdesorption'\
                    or self.type == 'photodesorption' or self.type == 'photodissociation':
                self.layer = self.reactants[0].layer
            elif len(self.reactants) > 1:
                self.layer = self.reactants[0].layer
            elif self.type.startswith('chemisorption'):
                self.layer = self.reactants[0].layer
            else:
                raise ValueError(
                    'Layer number not detected for reaction: %s' % self.verbatim)
        except:
            print("reaction type: ", self.type)
            raise ValueError(
                'Layer number not detected for reaction: %s' % self.verbatim)

    # ********************
    def getHashExploded(self):
        self.hashExploded = "_".join(
            sorted(["".join(x.exploded) for x in self.reactants]))
        self.hashExploded += "__"
        self.hashExploded += "_".join(sorted(["".join(x.exploded)
                                              for x in self.products]))

        return self.hashExploded

    # ********************
    # find reaction type based on reactants/products and store as self.type attribute
    def findType(self):

        # if CRdesorption is already set ignore find type
        if self.type == "CRdesorption":
            return

        # determine type for two body reaction
        if len(self.reactants) == 2:
            # count number of gas-phase species
            nGas = sum([x.isGas for x in self.products])
            # determine type based on number of products
            # and number of gas-phase species
            if len(self.products) == 1 and nGas == 1:
                # dust+dust->gas
                self.type = "2body_gas"
            elif len(self.products) == 1 and nGas == 0:
                # dust+dust->dust
                self.type = "2body_dust"
            elif len(self.products) == 2 and nGas == 2:
                # dust+dust->gas+gas
                self.type = "2body_gas_gas"
            elif len(self.products) >= 2 and nGas == 1:
                # dust+dust->gas+dust
                self.type = "2body_gas_dust"
            elif len(self.products) >= 2 and nGas == 0:
                # dust+dust->dust+dust
                self.type = "2body_dust_dust"
            else:
                print("ERROR: unknown 2body type " + self.verbatim)
                sys.exit()

        # determine type for single body with single product
        if len(self.reactants) == 1 and len(self.products) == 1:
            if self.reactants[0].isGas:
                # gas->dust
                self.type = "freezeout"
            if self.products[0].isGas:
                # dust->gas
                self.type = "evaporation"

        # if yield is set, then is photodesorption reaction
        if self.yieldPD > 0e0:
            self.type = "photodesorption"
        if self.yieldPD < 0e0:
            self.type = "photodissociation"

        # create special treatment of yield is -42, to allow exclusion of photodesorbtion for select species:
        if int(self.yieldPD) == -42.0:
            self.type = "photodesorption"
            self.yieldPD = 0e0

        reactantNames = sorted([x.dictname for x in self.reactants])
        productNames = sorted([x.dictname for x in self.products])
        # determine chemisorption types
        if reactantNames == ["H"] and productNames == ["Hc"]:
            self.type = "chemisorption_pc"
        if reactantNames == ["Hc"] and productNames == ["H"]:
            self.type = "chemisorption_cp"
        if reactantNames == sorted(["Hc", "H"]) and productNames == ["H2"]:
            self.type = "chemisorption_cp"
        if reactantNames == sorted(["Hc", "Hc"]) and productNames == ["H2"]:
            self.type = "chemisorption_cc"

    # *********************
    # compute mobility for chemisorption, see Cazaux+2004 (Erratum 2010)
    # this routine has been tested with Iqbal+2012 plots
    def computeMobility(self, ijtype):
        from math import pi, sinh, sqrt, exp, sin, log10
        from scipy.integrate import quad
        from numpy import inf

        # some constants
        mp = 1.6726219e-24  # proton mass, g
        kb = 1.38064852e-16  # boltzmannm constant, erg/K
        hbar = 1.0545718e-27  # planck constant/2pi, erg*s
        hbar2 = hbar ** 2

        # tunnelling mobility, see Cazaux+2004 (Erratum 2010)
        def Tij1(E, Bi, Bj, Bij, Za):
            sEBij = sqrt((E - Bij) / E)
            sinh2 = sinh(sqrt(2e0 * mp * (Bi - E) * kb / hbar2) * Za) ** 2
            return 4e0 * sEBij / ((1e0 + sEBij) ** 2 + Bi * Bj * sinh2 / (Bi - E) / E)

        # hopping mobility, see Cazaux+2004 (Erratum 2010)
        def Tij2(E, Bi, Bj, Bij, Za):
            sEBij = sqrt((E - Bij) / E)
            sin2 = sin(sqrt(2e0 * mp * (E - Bi) * kb / hbar2) * Za) ** 2
            return 4e0 * sEBij / ((1e0 + sEBij) ** 2 - Bi * Bj * sin2 / (Bi - E) / E)

        # integrand, Cazaux+2004 Eqn.3
        def fint(E, Tgas, Tij, Bi, Bj, Bij, Za):
            return exp(-E / Tgas) * Tij(E, Bi, Bj, Bij, Za)

        # see Tab.2 Cazaux+2004, assuming silicate
        Za = 2.5e-8  # cm
        ZA = 2e-8  # cm
        Es = 1e2  # K
        Ep = 4e2  # K
        Ec = 1e4  # K
        Esp = 1e2  # K
        Esc = Ec / 2e0  # K

        # see Tab.1 Cazaux+2004 (Erratum 2010)
        # this change with reaction type
        dataC10 = {"pp": [Ep - Esp, Ep - Esp, 0e0, ZA, Ep],
                   "cc": [Ec - Esc, Ec - Esc, 0e0, ZA, Ec],
                   "pc": [Ep - Es, Ec - Es, Ep - Ec, Za, Ep],
                   "cp": [Ec - Es, Ep - Es, Ep - Ec, Za, Ec]}

        # get data depending on reaction type
        (Bi, Bj, Bij, Zaa, El) = dataC10[ijtype]

        Ns = 2e14  # 1/cm2, see Iqbal+2012
        nui = sqrt(2e0 * Ns * El * kb / pi ** 2 / mp)  # Cazaux+2004 Eqn.4
        xdata = []
        ydata = []
        imax = 300  # number of points to compute data
        Tmin = 1e0  # K
        Tmax = 1e3  # K

        # loop on temperature range to compute alpha_ij, Cazaux+2004 Eqn.3, 5
        for i in range(imax):
            Tgas = i * (Tmax - Tmin) / (imax - 1) + Tmin  # K
            # compute tunnelling and hopping probability
            Pij1 = quad(fint, 0e0, Bi, args=(
                Tgas, Tij1, Bi, Bj, Bij, Zaa), epsabs=0)[0]
            Pij2 = quad(fint, Bi, inf, args=(
                Tgas, Tij2, Bi, Bj, Bij, Zaa), epsabs=0)[0]
            # append data to dataset
            xdata.append(Tgas)
            ydata.append(log10(nui * (Pij1 + Pij2) / Tgas))

        return [xdata, ydata]

    # ***************
    # check if species named name is present in the self reaction
    def hasSpecies(self, name):
        return self.hasReactant(name) or self.hasProduct(name)

    # ***************
    # check if reactant named name is present in the self reaction
    def hasReactant(self, name):
        return name in [x.dictname for x in self.reactants]

    # ***************
    # check if product named name is present in the self reaction
    def hasProduct(self, name):
        return name in [x.dictname for x in self.products]

    # ***************
    # compute rate tunnelling and evaporation probability given the reaction type
    # and prepare reaction in F90 format, all stored as attribute
    def rate(self):
        # constants and parameters
        kb = 1.38064852e-16  # Boltzmann constant, erg/K
        hbar = 1.0545718e-27  # planck constant/2pi, erg*s
        mp = 1.6726219e-24  # proton mass, g
        mref = 130.  # reference mass for evaporation probability, amu
        #nu0_generic = 1e12  # Debye frequency, 1/s
        ar = 1e-8  # barrier width for tunnelling, cm
        if self.barrierWidthFromFile:
            ar = self.barrierWidth
            #pass
        else:
            self.barrierWidth = ar

        # Number of sites / cm2
        Ns = 1.5e15
        # Frequency from harmonic oscillator approximation
        nu0 = [sqrt(2e0 * Ns * x.Eice * kb / pi**2 / (x.mass*mp))
               for x in self.reactants]
        nu_max = max(nu0)
        # default evaporation probability is 1
        self.Pdelta = 1e0
        # default reaction in F90 format, None to trigger ERROR if unknown type
        self.krateF90 = None

        # Check for reactions without type:
        if self.type == None:
            print("Following reaction lacks type: ", self.verbatim)

        # special treatment:
        specials = ['p_H2', 'o_H2', 'H2']#, 'HD', 'D2']
        # generic 2body
        if self.type.startswith("2body"):
            # Garrod et al (2007) description of evaporation probability:
            if (len(self.products) == 1):
                # H2 always desorp (since we assume fixed H2 ice), THIS REACTION SHOULD NOT HAPPEN.
                if self.products[0].namebase in []: #specials:
                    if self.products[0].isGas:
                        self.Pdelta = 1.0
                    else:
                        self.Pdelta = 0.0
                else:
                    N = self.products[0].natoms
                    if N == 2:
                        s = 2.0
                    else:
                        s = 3.0*N - 5.0
                    if (self.dH == 0.0 or np.isnan(self.dH) or self.dH <= 0.0):
                        # This assumes zero evaporation probabili when enthalpy change is unknown.
                        Pevap = 0.0
                    else:
                        try:
                            P = (1e0 - self.products[0].Eice / self.dH)**(s-1.0)
                        except OverflowError:
                            P = 1e12
                            print("Error in P for: ", self.verbatim, self.dH, self.products[0].Eice)
                        ap = 0.01
                        Pevap = ap*P / (1.0 + ap*P)

                    if self.products[0].isGas:
                        self.Pdelta *= Pevap
                    else:
                        self.Pdelta *= 1e0 - Pevap
            else:
                for p in self.products:
                    #H2 = any([p.dictname in specials for p in self.products])
                    if p.isGas and p.namebase not in specials:
                        print("ERROR with reaction: ", self.verbatim, p.name, p.namebase)
                        sys.exit()
            # loop on products to compute evaporation probability
            # for species in self.products:
            #     # reduced mass for evaporation
            #     eps = (mref - species.mass) ** 2 / \
            #         (mref + species.mass) ** 2
            #     # evaporation probability
            #     try:
            #         Pevap = exp(-species.Eice /
            #                     species.dof / self.dH / eps)
            #     # In case dH is nan:
            #     except OverflowError:
            #         Pevap = 0e0

            #     if species.namebase == 'H2':
            #         self.Pdelta = 1.0

            #     # gas product has probability Pevap, while dust product 1-Pevap
            #     if species.isGas:
            #         self.Pdelta *= Pevap
            #     else:
            #         self.Pdelta *= 1e0 - Pevap
                    

            # check if probability is more than 1
            if self.Pdelta > 1.0001:
                print("ERROR: evaporation probability > 1 for " + self.verbatim)
                print("We don't trust this, setting evaporation probability to 0 (default)")
                self.Pdelta = 0.0
            elif self.Pdelta < 0:
                print("ERROR: Pdelta < 0 for " + self.verbatim)
                print(
                    "We don't trust this, setting evaporation probability to 0 (default)")
                self.Pdelta = 0.0


            # join exponential terms for reaction
            joinedSum = " + ".join(["(exp(" + strF90(-x.Ediff) + "*invTd)*" + strF90(nu0[i]) + ")"
                                    for i, x in enumerate(self.reactants)])

            # Check if homogenous reactants:
            if self.reactants[0].name == self.reactants[1].name:
                homogeneous = True
            else:
                homogeneous = False

            # Competative processes (Garrod & Pauly 2011 and Ruaud et. al. 2016) if Ea > 0.0
            if self.Ea > 0.0:
                # reduced mass
                mu = (self.reactants[0].mass * self.reactants[1].mass) \
                    / sum([x.mass for x in self.reactants])
                # reaction probability through tunneling
                Ptunnel = exp(-2e0 * ar / hbar *
                            sqrt(2e0 * mu * mp * kb * self.Ea))
                # Classical thermal crossing
                Pclassical = "exp(-" + strF90(self.Ea) + " * invTd)"
                # max of the two methods, but should it not be the sum of the probabilities?:
                if Ptunnel < 1e-30:
                    Pbarrier = Pclassical
                else:
                    Pbarrier = "max(" + strF90(Ptunnel) + ", " + Pclassical + ")"

                # In principal you should include kevap in kappa, bbut in practice, as pointed out by Garrod & Pauly it makes little difference.
                #kevap =  " + ".join([strF90(nu0[i]) + "*exp(-" + str(x.Eice) + "*invTd)" for i,x  in enumerate(self.reactants)])
                
                kappa = nu_max * self.Pdelta
                if self.Bratio != None:
                    kappa *= self.Bratio
                kappa = strF90(kappa) + "*" + Pbarrier + "&\n*(" + \
                    strF90(nu_max) + "*" + Pbarrier + " + " + \
                    joinedSum + ")**(-1d0) &\n"

                # Without reaction-diffusion competition
                # if self.Bratio == None:
                #     kappa = Pbarrier + "*" + strF90(self.Pdelta)
                # else:
                #     kappa = Pbarrier + "*" +  strF90(self.Pdelta * self.Bratio)
                
            else:
                if self.Bratio != None:
                    kappa = strF90(self.Bratio * self.Pdelta)
                else:
                    kappa = strF90(self.Pdelta)

            if homogeneous:
                kappa = "5d-1 * " + kappa
            self.krateF90 = kappa + \
                "*(" + joinedSum + ")*indns"

            # OPR ratio for H2 and D2 formation, assumed statistical for now.
            if self.products[0].name.startswith('p_H2'):
                self.krateF90 = strF90(1.0/4.0) + " * " + self.krateF90
            elif self.products[0].name.startswith('o_H2'):
                self.krateF90 = strF90(3.0/4.0) + " * " + self.krateF90
            elif self.products[0].name.startswith('p_D2'):
                self.krateF90 = strF90(1.0/3.0) + " * " + self.krateF90
            elif self.products[0].name.startswith('o_D2'):
                self.krateF90 = strF90(2.0/3.0) + " * " + self.krateF90
            else:
                pass

        # freeze-out velocity*geometric-cross-section
        if self.type == "freezeout":
            r = self.reactants[0]
            p = self.products[0]
            self.krateF90 = strF90(
                sqrt(8e0 * kb / pi / r.mass / mp)) + "*prefreezeout * kstick(%i - surface_start)" % (p.idx+2)

        # evaporation with polanyi-wigner
        if self.type == "evaporation":
            species = self.reactants[0]
            self.krateF90 = strF90(nu0[0]) + \
                "*exp(-" + str(species.Eice) + "*invTd)"

        # photodissociation, as per Furuya et al. (2013)
        if self.type == "photodissociation":
            # gamma_H2O = 1e-3
            P_H2O_isrf = 5.4e-3
            P_H2O_cr = 4.7e-3
            P_CO2 = 2.4e-3  # Same for ISRF and CR
            rateTemplate = "Fnot * Gnot * kph_factor * exp(-{:f} * {:s})"

            # LOTS OF EXCEPTIONS:
            # P-values are from Kenji Furuya. Check paper for references.
            if self.reactants[0].namebase == "H2O":
                # H2O is treated differently.
                self.krateF90 = rateTemplate.format(self.gamma, self.Avvar) + " * " + \
                    strF90(P_H2O_isrf) + " + F_cr * kph_factor * " + \
                    strF90(P_H2O_cr)
            elif self.reactants[0].namebase == "CO2":
                # CO2 is treated differently.
                self.krateF90 = "(Fnot * Gnot* kph_factor * exp(-gamma_CO2 * {:s})".format(
                    self.Avvar) + " + kph_factor * F_cr ) * " + strF90(P_CO2)
            elif self.reactants[0].namebase == 'CO':
                # Add CO self-shielding factor
                P = strF90(P_H2O_isrf) + "* ss_CO * (" + strF90(self.alpha) +\
                    "/k_H2O_ph_0) &\n * "
                P_CR = strF90(P_H2O_cr) + "* (" + strF90(self.alpha) +\
                    "/k_H2O_ph_0) &\n * "
                self.krateF90 = P + rateTemplate.format(
                    self.gamma, self.Avvar) + " + " + P_CR + "F_cr * kph_factor"
            elif self.reactants[0].namebase == 'H2':
                # Add H2 self-shielding factor
                P = strF90(P_H2O_isrf) + "* ss_H2 * (" + strF90(self.alpha) +\
                    "/k_H2O_ph_0) &\n * "
                P_CR = strF90(P_H2O_cr) + "* (" + strF90(self.alpha) +\
                    "/k_H2O_ph_0) &\n * "
                self.krateF90 = P + rateTemplate.format(
                    self.gamma, self.Avvar) + " + " + P_CR + "F_cr * kph_factor"
            elif self.reactants[0].namebase == 'HD':
                # Add HD self-shielding factor
                P = strF90(P_H2O_isrf) + "* ss_HD * (" + strF90(self.alpha) +\
                    "/k_H2O_ph_0) &\n * "
                P_CR = strF90(P_H2O_cr) + "* (" + strF90(self.alpha) +\
                    "/k_H2O_ph_0) &\n * "
                self.krateF90 = P + rateTemplate.format(
                    self.gamma, self.Avvar) + " + " + P_CR + "F_cr * kph_factor"
            elif self.reactants[0].namebase == 'N2':
                # Add N2 self-shielding factor
                P = strF90(P_H2O_isrf) + "* ss_N2 * (" + strF90(self.alpha) +\
                    "/k_H2O_ph_0) &\n * "
                P_CR = strF90(P_H2O_cr) + "* (" + strF90(self.alpha) +\
                    "/k_H2O_ph_0) &\n * "
                self.krateF90 = P + rateTemplate.format(
                    self.gamma, self.Avvar) + " + " + P_CR + "F_cr * kph_factor"
            else:
                # Adjust rate according to H2O ice. Again Furuya et al. 2013
                P = strF90(P_H2O_isrf) + "* (" + strF90(self.alpha) +\
                    "/k_H2O_ph_0) &\n * "
                P_CR = strF90(P_H2O_cr) + "* (" + strF90(self.alpha) +\
                    "/k_H2O_ph_0) &\n * "
                self.krateF90 = P + rateTemplate.format(
                    self.gamma, self.Avvar) + " + " + P_CR + "F_cr * kph_factor"
                # # Account for statistical OPR:
                # for p in self.products:
                #     if p.namebase in specials:
                #         if p.name.startswith('p_H2'):
                #             self.krateF90 = strF90(
                #                 1.0/4.0) + " * " + self.krateF90
                #         if p.name.startswith('o_H2'):
                #             self.krateF90 = strF90(
                #                 3.0/4.0) + " * " + self.krateF90
                #         if p.name.startswith('p_D2'):
                #             self.krateF90 = strF90(
                #                 1.0/3.0) + " * " + self.krateF90
                #         if p.name.startswith('o_D2'):
                #             self.krateF90 = strF90(
                #                 2.0/3.0) + " * " + self.krateF90

        # photodesorption, see Hollenbach+2009
        if self.type == "photodesorption":
            self.krateF90 = "Ffuva * " + strF90(self.yieldPD)
            if self.reactants[0].namebase == "CO":
                self.krateF90 = "Ffuva_CO * Y_CO"
            elif self.reactants[0].namebase == "H2":
                self.krateF90 = "Ffuva_H2 * " + strF90(self.yieldPD)
            elif self.reactants[0].namebase == "HD":
                self.krateF90 = "Ffuva_HD * " + strF90(self.yieldPD)
            elif self.reactants[0].namebase == "N2":
                self.krateF90 = "Ffuva_N2 * " + strF90(self.yieldPD)

        # CR desorption, see Reboussin+2014, Hasegawa+1993, Leger+1985
        if self.type == "CRdesorption":
            species = self.reactants[0]
            rate = nu_max * 3.16e-19 / 1.3e-17 * exp(-species.Eice / 7e1)
            self.krateF90 = "variable_crflux*" + strF90(rate)

        # chemisorption
        if self.type.startswith("chemisorption_"):
            # get chemisorption reaction type, pc, cp, or cc
            chemisorptionType = self.type.split("_")[1]
            # compute integral for mobility (rates)
            xdata, ydata = self.computeMobility(ijtype=chemisorptionType)
            # prepare data to fit in F90 array format
            yf90 = ", &\n".join([strF90(x) for x in ydata])
            # prepare rate using data fit function
            self.krateF90 = "1d1**fit1d(Tdust, " + str(len(xdata)) + ", " \
                            + strF90(min(xdata)) + ", " \
                            + strF90(1e0 / (max(xdata) - min(xdata))) + ", " \
                            + strF90(max(xdata) - min(xdata)) + ", " \
                            + " &\n (/" + yf90 + "/))"

        # check if krate is set
        if self.krateF90 is None:
            sys.exit("ERROR: unknown type " +
                     self.type + " when computing krate!")

        return self.krateF90

    # ********************
    # prepare the RHS in F90 format and store to attribute
    def buildF90RHS(self):
        # If layer has baseIdx, use this (layer above 1)
        try:
            idx = self.baseIdx
        except AttributeError:
            idx = self.idx+1

        # list of species number density variables in F90 format
        RHS = ["n(" + x.fidx + ")" for x in self.reactants]
        # reaction rate coefficient in F90 format
        self.RHS = "kall(" + str(idx) + ")"
        # write product of the number densities variables
        self.RHS += "*" + ("*".join(RHS))
