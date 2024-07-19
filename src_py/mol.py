from utils import isNumber
import sys
import copy


# molecule class, i.e. species class
class mol:
    # *******************
    # species class constructor
    # name is species name as string
    # atomsMass is a dictionary, key=atom name string, value=mass in amu
    def __init__(self, name, atomsMass, layer = 1):
        # only allow layer in [0, 1, 2]:
        if layer > 2:
            raise ValueError("Layer > 2 for %s. This is not allowed."% name)
        self.atomsMass = atomsMass
        # species name, gas phase species ends with _gas
        self.dictname = name

        # ions are always gas phase species
        if ("+" in name) and ("_gas" not in name):
            self.name = name + "_gas"
        elif ("–" in name) and ("_gas" not in name):
            self.name = name + "_gas"
        elif ("j" in name) and ("_gas" not in name):
            self.name = name + "_gas"
        elif ("k" in name) and ("_gas" not in name):
            self.name = name + "_gas"
        else:
            if layer == 1:
                self.name = name + "_surface"
            elif layer == 2:
                self.name = name + "_mantle"
            else:
                raise ValueError("ERROR in the molecular input. Check layer in [1,2] or gas species")
        # binding energy, K
        self.Eice = self.Ebare = self.Ediff = None
        # enthalpy of formation, K
        self.dH = None
        # mass, amu
        self.mass = None
        # index, zero-based
        self.idx = self.idxGas = self.idxTot = None
        # base name (i.e. without ending _gas)
        self.namebase = self.dictname.replace("_gas", "")

        # treat ortho, para, etc. the same as H2 when reading Ebind, deltaH etc.
        if "p_" in self.namebase:
            self.namebase = self.namebase.replace("p_", "")
        if "o_" in self.namebase:
            self.namebase = self.namebase.replace("o_", "")
        if "l_" in self.namebase:
            self.namebase = self.namebase.replace("l_", "")
        if "c_" in self.namebase:
            self.namebase = self.namebase.replace("c_", "")
            
        # flag if species is gas phase
        self.isGas = (self.dictname.endswith("_gas"))

        # Avoid layer extension to gas-phase species
        if self.isGas:
            self.name = self.dictname
            self.fidx = "idx_" + self.dictname.replace("+", "j").replace("-", "k")
        else:
            if "+" in self.dictname or "-" in self.dictname:
                print("WARNING: surface species %s has a +/- sign. This is not allowed." % self.dictname)
                self.fidx = "idx_" + self.dictname.replace("+", "j").replace("-", "k")
            # F90 index variable, + is replaced with j
            if layer == 1:
                self.fidx = "idx_" + self.dictname + "_surface"
            else:
                self.fidx = "idx_" + self.dictname + "_mantle"
        
        # charge (count positive)
        self.charge = self.name.count("+")
        # charge (count negative)
        self.charge -= self.name.count("-")

        # electron has negative charge
        if self.dictname == "E_gas":
            self.charge = -1

        # gas phase species always layer 0
        if self.isGas:
            self.layer = 0
        else:
            # species layer:
            self.layer = layer

        # get latex name
        self.nameLatex = self.getNameLatex(name)

        # exploded is alphabetically sorted atom list, e.g. C2H -> ["C","C","H"]
        self.exploded = self.getExploded(self.namebase, atomsMass)

        # determine mass from exploded, NOTE: amu
        electronMass = 9.10938356e-28  # g
        protonMass = 1.6726219e-24  # g
        electronAMU = electronMass / (electronMass + protonMass)
        self.mass = sum([atomsMass[x] for x in self.exploded])
        self.mass += (name.count("-") - name.count("+")) * electronAMU
        if self.mass < 0e0:
            sys.exit("ERROR: species " + name + " has mass < 0e0!")

        # count number of atoms
        self.natoms = len([x for x in self.exploded if (x != "+")])
        # degrees of freedom, 3*N
        self.dof = 3 * self.natoms

    # ********************
    # convert name string into latex species name
    def getNameLatex(self, name):

        # electrons are special
        if name == "E_gas":
            return "e$^-$"
            
        latexName = ""
        # loop on name characters
        for char in list(name):
            if char in ["+", "-"]:
                # signs are superscripts
                latexName += "$^" + char + "$"
            elif isNumber(char):
                # numbers are subscripts
                latexName += "$_" + char + "$"
            else:
                # standard characters
                latexName += char

        # replace _gas with dedicated latex command
        latexName = latexName.replace("_gas", "$\\gas$")

        # return latex name
        return latexName

    # ********************
    # get species exploded name using namebase
    def getExploded(self, namebase, atomsMass):
        import itertools
        if "p_" in self.name:
            namebase = "p_" + namebase
        if "o_" in self.name:
            namebase = "o_" + namebase
        if "i_" in self.name:
            namebase = "i_" + namebase
        if "c_" in self.name:
            namebase = "c_" + namebase
        if "l_" in self.name:
            namebase = "l_" + namebase
        if "m_" in self.name:
            namebase = "m_" + namebase
        # copy namebase to replace after
        specName = namebase + ""
        # store keys sorted by inverse length
        atoms = sorted(list(atomsMass.keys()), key=lambda xx: len(xx), reverse=True)
        # produce unique character combinations
        alpha = ["".join(x) for x in list(itertools.product("XYZ", repeat=4))]
        # check to have enough combinations
        if len(atoms) > len(alpha):
            sys.exit("ERROR: in species parser alpha needs to be extended!")

        # replace atoms slash-separated
        for i in range(len(atoms)):
            specName = specName.replace(atoms[i], "/" + alpha[i] + "/")
        # replace double slashes
        while "//" in specName:
            specName = specName.replace("//", "/")
        # split at slashes
        aspec = [x for x in specName.split("/") if x != ""]

        # search for number and when found multiply previous non-number
        exploded = []
        for a in aspec:
            if isNumber(a):
                for j in range(int(a) - 1):
                    exploded.append(aold)
            else:
                exploded.append(a)
            aold = a

        # store exploded with real atom names
        try:
            exploded = [atoms[alpha.index(x)] for x in exploded]
        except:
            print("ERROR: wanted to parse ", namebase)
            print(" but something went wrong with ", exploded)
            print(" Available atoms are:", atoms)
            print(" Add to atom list file if needed.")
            sys.exit()

        return sorted(exploded)


