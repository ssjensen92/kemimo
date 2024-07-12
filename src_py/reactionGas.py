from utils import speciesToKIDA, strF90
import sys


# reaction class for gas-phase reactions
class reactionGas:

    # *******************
    # constructor for gas-phase reactions
    # row: reaction file row, here KIDA format
    # speciesList: list of species loaded from species file
    def __init__(self, row, speciesList, atomMassList, respectGasphaseLimits=True):
        # variables employed for cosmic rays and photochemistry
        self.CRvar = "variable_crflux"  # name of the CR flux variable
        self.Avvar = "variable_Av"  # name of the Av variable
        self.Tgasvar = "variable_Tgas"  # name of the Tgas variable

        # list of products and reactants
        self.reactants = []
        self.products = []

        # Layer = 0 for gas reactions
        self.layer = 0

        # this flag is false when unknown species are present
        self.hasSpeciesFlag = True

        # index and reaction rate in F90 format
        self.idx = self.krateF90 = None

        # true if there are multiple temperature ranges
        self.hasMultipleTranges = False

        # gas-phase reaction type
        self.type = "gasphase"

        # flag if barrier found from data file
        self.barrierFromFile = False

        # unknowns species of this reactions wrt to speciesList
        # this is not empty when unknown species are found
        self.unknownSpecies = []

        # default limits
        self.Tmin = -1e99
        self.Tmax = 1e99

        # respect temperature limits?
        self.respectGasphaseLimits = respectGasphaseLimits

        # parse a KIDA file line
        self.parseFormatKIDA(row, speciesList, atomMassList)

        # prepare verbatim, e.g. H + O -> OH
        self.verbatim = " + ".join([x.name for x in self.reactants])
        self.verbatim += " -> "
        self.verbatim += " + ".join([x.name for x in self.products])

        # prepare hash, e.g. H_O__OH, note double underscore
        self.hash = "_".join(sorted([x.name for x in self.reactants]))
        self.hash += "__"
        self.hash += "_".join(sorted([x.name for x in self.products]))

        # exploded hash, e.g. CH3 + H -> CH4 becomes CHHH_H__CHHHH, note double undersocre
        self.hashExploded = self.getHashExploded()

    # ********************
    # get reaction exploded hash, e.g. CH3 + H -> CH4 becomes CHHH_H__CHHHH
    def getHashExploded(self):
        self.hashExploded = "_".join(sorted(["".join(x.exploded) for x in self.reactants]))
        self.hashExploded += "__"
        self.hashExploded += "_".join(sorted(["".join(x.exploded) for x in self.products]))

        return self.hashExploded

    # ********************
    # get missing species with the given combination of atoms
    def getMissingSpecies(self, atoms):

        relevantSpecies = []
        numbers = [str(x) for x in range(10)]
        for species in self.unknownSpecies:
            rep = species + "#"
            for atom in atoms + ["+", "-", "E"] + numbers:
                rep = rep.replace(atom, "")
            if rep == "#":
                relevantSpecies.append(species)

        return relevantSpecies

    # ***************
    # check if species is present
    def hasSpecies(self, name):
        return self.hasReactant(name) or self.hasProduct(name)

    # ***************
    # check if reactant is present
    def hasReactant(self, name):
        return name in [x.name for x in self.reactants]

    # ***************
    # check if product is present
    def hasProduct(self, name):
        return name in [x.name for x in self.products]

    # ********************
    # parse KIDA database file row
    def parseFormatKIDA(self, row, speciesList, atomMassList):

        # names that will be ignored, note upper case
        specials = ["", "G", "PHOTON", "CR", "CRP"]

        CRvar = self.CRvar  # name of the CR flux variable
        Avvar = self.Avvar  # name of the Av variable
        Tgasvar = self.Tgasvar  # name of the Tgas variable

        # number of reactants and products expected in KIDA file
        maxReactants = 3
        maxProducts = 5

        # spacing format
        fmt = [11] * 3 + [1] + 5 * [11] + [1] + 3 * [11] + [8, 9] + [1, 4, 3] \
            + 2 * [7] + [3, 6, 3, 2]
        # keys names
        keys = ["R" + str(i) for i in range(maxReactants)] + ["x"] \
            + ["P" + str(i) for i in range(maxProducts)] + ["x"] \
            + ["a", "b", "c"] + ["F", "g"] + ["x", "unc", "type"] \
            + ["tmin", "tmax"] + ["formula", "num", "subnum", "recom"]

        srow = row.strip()

        dataRow = dict()
        position = 0
        # loop on format to get data from the row as a dictionary
        for i in range(len(fmt)):
            # fill dictionary using format keys
            dataRow[keys[i]] = srow[position:position + fmt[i]].strip()
            # determine if spaces are present to check correct format
            startSpace = dataRow[keys[i]].startswith(" ")
            endSpace = dataRow[keys[i]].endswith(" ")
            hasSpace = (" " in dataRow[keys[i]])
            # error message if problem with spaces
            if (not startSpace) and (not endSpace) and hasSpace:
                print("ERROR: in KIDA network row element has spaces in the middle!")
                print(" Probably format problems:", dataRow[keys[i]])
                print(" Line here below triggered the ERROR")
                print(srow)
                sys.exit()
            # increase format position
            position += fmt[i]

        # fill species dictionary
        #speciesDict = {"".join(x.exploded): x for x in speciesList if x.isGas}
        speciesDict = {x.name.replace("_gas", ""): x for x in speciesList if x.isGas}

        # loop on expected reactants
        for i in range(maxReactants):
            specName = dataRow["R" + str(i)].strip()
            specName = speciesToKIDA(specName)
            # ignore specials
            if specName.upper() in specials:
                continue

            # add to unknown species if not auto-generated from species file
            if specName not in speciesDict:
                self.hasSpeciesFlag = False
                self.unknownSpecies.append(specName)
                continue
            self.reactants.append(speciesDict[specName])

        # loop on expected products
        for i in range(maxProducts):
            specName = dataRow["P" + str(i)].strip()
            specName = speciesToKIDA(specName)
            if specName.upper() in specials:
                continue

            if specName not in speciesDict:
                self.hasSpeciesFlag = False
                self.unknownSpecies.append(specName)
                continue
            self.products.append(speciesDict[specName])

        # read temperature ranges
        if self.respectGasphaseLimits:
            self.Tmin = float(dataRow["tmin"])
        else:
            if float(dataRow["tmin"]) == 10.0:
                self.Tmin = 5.0 # Set to 5 for dense/prestellar core conditions.
            else:
                self.Tmin = float(dataRow["tmin"]) 
        self.Tmax = float(dataRow["tmax"])

        # Formula is a number that refers to the formula needed to compute the rate coefficient of the reaction
        # see http://kida.obs.u-bordeaux1.fr/help
        # 1: Cosmic-ray ionization (direct and secondary)
        # 2: Photo-dissociation (Draine style)
        # 3: Kooij
        # 4: ionpol1
        # 5: ionpol2
        # 6: Troe fall-off (NOT SUPPORTED!)
        # 7: grain recombination
        arow = dataRow
        arow["formula"] = int(arow["formula"])
        self.formula = arow["formula"]
        self.kidatype = int(arow["type"])
        if arow["formula"] == 0:
            KK = arow["a"]
            mass = 0.0
            isElectron = False
            for r in self.reactants:
                if r.name == 'E_gas':
                    isElectron = True
                mass += r.mass
            KK += "* 1d0/sqrt(" + strF90(mass) + ")"
            if isElectron:
                KK += " * kgr_neutral"
            else:
                KK += " * kgr_ion"
        elif arow["formula"] == 1:
            if self.kidatype == 1 or self.kidatype == 0:
                self.type = "gasphase_CR"
                KK = arow["a"] + "*" + CRvar
            elif self.kidatype == 2 or self.kidatype == 3:
                self.type = "gasphase_CRP"
                if "o_H2" not in speciesDict:
                    KK = arow["a"] + "*0.5*" + "n(idx_H2_gas)/(n(idx_H_gas)+2*n(idx_H2_gas))"+ "*" + CRvar 
                else:
                    KK = arow["a"] + "*2.0*" + "(n(idx_o_H2_gas)+n(idx_p_H2_gas))/(n(idx_H_gas)+2*(n(idx_o_H2_gas)+n(idx_p_H2_gas)))"+ "*" + CRvar 
            else:
                print('error in: ', self.kidatype, ' ;  ', print(srow))
                sys.exit()
        elif arow["formula"] == 2:
            self.type = "gasphase_Av"
            self.alpha = float(arow["a"])
            self.gamma = float(arow["c"])
            KK = arow["a"]
            if float(arow["c"]) != 0e0:
                KK += "*Gnot*exp(-" + arow["c"] + "*" + Avvar + ")"
            # Check for self-shielding functions:
            if len(self.reactants) > 0:
                if self.reactants[0].name == 'CO_gas':
                    KK = "ss_CO*" + KK
                elif self.reactants[0].name == 'p_H2_gas' or self.reactants[0].name == 'o_H2_gas' or self.reactants[0].name == 'H2_gas':
                    KK = "ss_H2*" + KK
                elif self.reactants[0].name == 'HD_gas':
                    KK = "ss_HD*" + KK
                elif self.reactants[0].name == 'N2_gas':
                    KK = "ss_N2*" + KK
        elif arow["formula"] == 3 or str(arow["formula"]).startswith('3'):
            KK = arow["a"]
            if float(arow["b"]) != 0e0:
                KK += "*(" + Tgasvar + "/3d2)**(" + arow["b"] + ")"
            if float(arow["c"]) != 0e0:
                KK += "*exp(-" + arow["c"] + "/" + Tgasvar + ")"
        elif arow["formula"] == 4:
            KK = arow["a"]
            if float(arow["b"]) != 1e0:
                KK += "*" + arow["b"]
            gpart = ""
            if float(arow["c"]) != 0e0:
                gpart = "+ 0.4767d0*(" + arow["c"] + \
                    ")*sqrt(3d2/" + Tgasvar + ")"
            KK += "*(0.62e0 " + gpart + ")"
        elif arow["formula"] == 5:
            KK = arow["a"]
            if float(arow["b"]) != 1e0:
                KK += "*" + arow["b"]
            gpart = ""
            if float(arow["c"]) != 0e0:
                gpart = "+ 0.0967e0*(" + arow["c"] + \
                    ")*sqrt(3d2/" + Tgasvar + ") + ("
                gpart += arow["c"] + ")**2*28.501d0/" + Tgasvar
                KK += "*(1d0 " + gpart + ")"
        else:
            print(srow)
            print("ERROR: KIDA formula " + \
                str(arow["formula"]) + " not supported!")
            sys.exit()

        # replace double signs using standard algebra rules
        KK = KK.replace("--", "+").replace("++", "+")
        self.krateF90 = KK.replace("-+", "-").replace("+-", "-")

    # ********************
    # prepare the RHS in F90 format and store to attribute
    def buildF90RHS(self):

        RHS = ["n(" + x.fidx + ")" for x in self.reactants]
        self.RHS = "kall(" + str(self.idx + 1) + ")"
        self.RHS += "*" + ("*".join(RHS))

    # *****************
    def getJacPD(self, fidx):
        RHS = ["n(" + x.fidx + ")" for x in self.reactants]
        RHS.remove("n(" + fidx + ")")
        JPD = "kall(" + str(self.idx + 1) + ")"
        if len(RHS) == 0:
            return JPD
        return JPD + "*" + ("*".join(RHS))
