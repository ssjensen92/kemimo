
# molecule class, i.e. species class
class dummy_species():
    # *******************
    # species class constructor for mly mask
    # 
    def __init__(self):
        self.layer = 0
        self.name = "dummy"
        self.dictname = "dummy"
        self.namebase = "dummy"
        self.nameLatex = "$f_{\mathrm{dummy}}$"
        self.Eice = self.Ebare = None
        self.isGas = False
        self.dH = None
        self.mass = None
        self.idx = self.idxGas = self.idxTot = None
        self.fidx = "idx_" + self.name.replace("+", "j")
        # Step 4: Adjust attributes and methods as needed, leveraging inheritance