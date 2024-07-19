from utils import tabrow, strF90, doPP, doPP_datfile, doPP_Python, latexExp, latexInt
from utils import speciesToKIDA, getReactionType
from mol import mol
from mly import mly
from dummy_species import dummy_species
from reaction import reaction
from monolayer_reaction import monolayer_reaction
from reactionGas import reactionGas
from subprocess import check_output
from math import sqrt, pi
import sys
import re
import os
import numpy as np
from shutil import copyfile
import time
import copy
from multiprocessing import Pool
import multiprocessing as mp

import itertools

# ***************************
# database class to load data and do the main operations
class database:
    # *****************
    # database class constructor
    def __init__(self, datadir="./standard_data/", nlayers=2, layerThickness=4.0,
                 ignoreMissingH=True, ignoreMissingEbind=True, 
                 allSpeciesToDust=False, include_H2=True, limit2body=True, respectGasphaseLimits=True,
                 multiprocessing=True, doSwap=False, H2spin=True, encounterDesorption=False):
        '''
            nlayers is the number of ice layers to consider, which will be hardcoded in the ODE system.
            layerThickness is the thickness of one individual layer in mly.
            reactionType is for benchmarking against various models.

            nlayers: the number of ice layer (2 in a three-phase model)
            layer thickness: the thickness of the *surface* layer in mly
        '''
        self.layerThickness = layerThickness
        # do (physical) swapping between surface and mantle:
        self.doSwap = doSwap
        # Include H2 formation in network or keep static?
        self.include_H2 = include_H2
        # include H2 spin (this flag is necessary for kemimo_ode where o/p H2 is hardcoded)
        self.H2spin = H2spin
        # Limit 2body reactions to Ea only or use combinatorics to look for reactions?
        self.limit2body = limit2body
        # allSpeciesToDust: Add all species present in the gasphase to the dust phase.
        self.allSpeciesToDust = allSpeciesToDust
        # Respect (lower) gasphase limits. Many reactions in KIDA limited to > 10 K. This will lower the limit to 5 K.
        self.respectGasphaseLimits = respectGasphaseLimits  #  USE AT OWN RISK
        
        self.datadir = datadir

        # multiprocessing
        self.multiprocessing = multiprocessing
        # *****************
        # fileSpecies: name of the file with list of species
        # fileEbin: name of the file species binding energies
        # fileDeltaH: name of the file with species enthalpy of formations
        # fileAtoms: name of the file with atoms and their masses
        # fileEa: name of the file with reactions activation barriers
        # fileYields: name of the file with photodissociation yields
        fileSpecies = datadir + "species.dat"
        fileDustSpecies = datadir + "speciesDust.dat"
        fileEbind = datadir + "Ebind.dat"
        fileDeltaH = datadir + "deltaH.dat"
        fileAtoms = datadir + "atoms.dat"
        fileEa = datadir + "Ea.dat"
        fileBarrierWidth = datadir + "barrierWidths.dat"
        fileBratio = datadir + "branchingRatios.dat"
        fileYields = datadir + "yieldPD.dat"
        fileCheckBarriers = datadir + "check_barrier.dat"
        # Start timing:
        start = time.time()
        # load data from files, see details in the corresponding methods
        self.loadSpecies(fileSpecies, fileDustSpecies)
        self.loadEbinds(fileEbind)
        self.loadDeltaHs(fileDeltaH)
        self.loadAtoms(fileAtoms)
        self.loadBarriers(fileEa)
        self.loadBarrierWidths(fileBarrierWidth)
        self.loadBranchingRatios(fileBratio)
        self.loadYieldsPD(fileYields)
        self.loadCheckBarriers(fileCheckBarriers)

        # prepare species based on files data (including gas)
        self.species = [mol(x, self.mass, layer=1) for x in self.speciesNames]

        # get H2 ortho/para idx:
        for idx, s in enumerate(self.species):
            if s.name == 'p_H2_gas':
                self.idx_p_H2_gas = idx
            elif s.name == 'o_H2_gas':
                self.idx_o_H2_gas = idx
            else:
                continue

        # compute number of species
        self.nSpecies = len(self.species)
        nGasSpecies = 0
        for s in self.species:
            if s.isGas:
                nGasSpecies += 1
        self.nDustSpecies = self.nSpecies - nGasSpecies
        self.nGasSpecies = nGasSpecies

        # store names of the missing species (found in KIDA)
        self.missingSpecies = []

        self.icount = 0
        # set the properties to the species found
        if ignoreMissingH:
            print("We ignore species with missing enthalpy from the surface network.")

        # Book keeping
        n_missing_H = 0
        n_missing_Ebind = 0
        missing_H_species = []
        missing_Ebind_species = []
        # Read species data
        for species in self.species:
            species.idx = self.icount
            # only for neutral species
            if species.charge == 0 and not species.name == "GRAIN0":
                # check if species has no binding energy
                if species.namebase not in self.Eice:
                    # Exception for species starting with "l_", "c_", etc.
                    if not species.isGas and species.name.strip('_0001') in self.Eice: 
                        self.Eice[species.namebase] = self.Eice[species.name.strip(
                            '_0001')]
                        self.Ebare[species.namebase] = self.Ebare[species.name.strip(
                            '_0001')]
                    elif species.isGas and species.name.strip('_gas') in self.Eice:
                        self.Eice[species.namebase] = self.Eice[species.name.strip(
                            '_gas')]
                        self.Ebare[species.namebase] = self.Ebare[species.name.strip(
                            '_gas')]
                    # If the above exception does not work:
                    # usual approach: ignore missing binding energy, set to H2O.
                    elif ignoreMissingEbind:
                        n_missing_Ebind += 1
                        missing_Ebind_species.append(species.namebase)
                        self.Eice[species.namebase] = self.Eice['H2O']
                        self.Ebare[species.namebase] = self.Ebare['H2O']
                    else:
                        print(
                            ("ERROR: binding energy missing for " + species.namebase))
                        print((species.charge, species.dictname))
                        sys.exit()
                # check if species has no enthalpy of formation
                if species.namebase not in self.deltaH:
                    # Exception for species starting with "l_", "c_", etc.
                    if not species.isGas and species.name.strip('_0001') in self.deltaH:
                        self.deltaH[species.namebase] = self.deltaH[species.name.strip(
                            '_0001')]
                        self.deltaH[species.namebase] = self.deltaH[species.name.strip(
                            '_0001')]
                    elif species.isGas and species.name.strip('_gas') in self.deltaH:
                        self.deltaH[species.namebase] = self.deltaH[species.name.strip(
                            '_gas')]
                        self.deltaH[species.namebase] = self.deltaH[species.name.strip(
                            '_gas')]
                    elif ignoreMissingH:
                        n_missing_H += 1
                        missing_H_species.append(species.namebase)
                        self.deltaH[species.namebase] = np.nan
                    else:
                        print(("ERROR: enthalpy missing for " + species.namebase))
                        sys.exit()
                # get properties form the dictionaries
                # binding energy ice, K
                species.Eice = self.Eice[species.namebase]
                # binding energy bare, K
                species.Ebare = self.Ebare[species.namebase]
                # diffusion energy:
                species.Ediff = 0.5 * self.Eice[species.namebase]
                # enthalpy of formation, K
                species.dH = self.deltaH[species.namebase]
            self.icount += 1

        if ignoreMissingH and len(missing_H_species) > 0:
            print("Missing enthalpy information on %i species!" % n_missing_H)
            print(np.array(missing_H_species))

        if ignoreMissingEbind and len(missing_Ebind_species) > 0:
            print("Ignoring missing Binding energies on %i species!" %
                  n_missing_Ebind)
            print(np.array(missing_Ebind_species))
            print("Assuming H2O binding energy for all species without information!")

        print("****************************************************")
        if nlayers == 6:
            print("Running 7-phase model. This is in alpha stage.")
            self.nlayers = nlayers
        elif nlayers == 2:
            print("Running three phase model")
            # create layers above 1:
            self.nlayers = nlayers
        else:
            self.nlayers = 1
            print("Running with single surface layer (bulk ice)")


        # use combinatorics to find reactions
        findReactionsStart = time.time()
        if self.limit2body == False:
            self.findReactions()
        else:
            self.loadReactions()

        # add Hincelin+2015 H2 encounter-desorption?
        if encounterDesorption:
            self.addEncounterDesorption()



        # add photorates using reactions from yield file
        self.addPhotoRates()

        # add CR desorption rates
        self.addCRdesorptionRates()

        # load gas phase form KIDA, given the species from dictionary
        self.loadKIDA()

        findReactionsTime = time.time() - findReactionsStart
        print(("Time to find/load network reactions: %1.3f seconds" %
               findReactionsTime))

        # if there are missing species print a warning
        self.checkMissingSpecies()

        # save chemical network to a file
        self.saveNetworkToFile()

        # prepare rates
        self.computeRates()

        # Create surface and mantle masks
        mly_surface = mly(seperationDistance=3e-8, radius=1e-5, name='surface')
        mly_mantle = mly(seperationDistance=3e-8, radius=1e-5, name='mantle')
        # Adjust indexing
        mly_surface.idx = self.icount
        mly_mantle.idx = self.icount+1

        # Add masks:
        self.species.append(mly_surface)
        self.species.append(mly_mantle)
        # Add mask reaction, only used for rate of mask change:
        self.maskReaction = [monolayer_reaction(mly_surface)]

        # Add dummy species used in reaction loop in kemimo_ode::
        dummy = dummy_species()
        dummy.idx = self.icount+2
        self.dummy_idx = dummy.idx
        self.species.append(dummy)

        # Finished
        elapsed = time.time() - start
        print(("Time to initiate network database: %1.3f seconds" % elapsed))

    # *****************
    # check if reactions of a list are unique (for Tgas limit)
    def checkReactionsUnique(self, reactions):
        for rea1 in reactions:
            foundEqual = 0
            for rea2 in reactions:
                # skip gasphase_CR and gasphase_Av
                if rea2.type != "gasphase":
                    continue
                # count equal reactions
                if rea1.hash == rea2.hash:
                    foundEqual += 1
                # if more than 1 (includes self) set attribute to multiple rate
                if foundEqual >= 2:
                    # print "WARNING: reaction "+rea1.verbatim+" has multiple temperature ranges!"
                    rea1.hasMultipleTranges = True
                    break

    # *****************
    # loop over photodesorption data to build rates
    def addPhotoRates(self):
        desorptionReactants = []
        # get a species dictionary where key=name, value=species object
        speciesDict = {x.dictname: x for x in self.species}
        layeredSpeciesDict = {x.name: x for x in self.species}
        RR_names = []
        # loop on data to search for yields
        for data in self.yieldPD:
            #i = 1
            knownSpecies = True
            # check if all the species in the photodissociation rate are known
            for speciesName in data["reactants"] + data["products"]:
                if speciesName not in speciesDict:
                    knownSpecies = False
            # skip reaction if not all the species are known
            if not knownSpecies:
                continue
            # get reactants and products from the species dictionary
            RRs = []
            PPs = []
            for x in data["reactants"]:
                if x.endswith('_gas'):
                    RRs.append(speciesDict[x])
                else:
                    RRs.append(layeredSpeciesDict[x + '_%0.4i' % 1])
                RR_names.append(RRs[-1].name)
            for x in data["products"]:
                if x.endswith('_gas'):
                    PPs.append(speciesDict[x])
                else:
                    PPs.append(layeredSpeciesDict[x + '_%0.4i' % 1])
            # create a new reaction object with the given yield
            rea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios,
                           yieldPD=data["yield"])
            desorptionReactants.append(RRs[0].name)

            # append the reaction to the list of dust reactions
            self.reactions.append(rea)

        # #######################################################
        # Check which species are missing photodesorption-reactions and add these with generic yield:
        # Specials for which we do not consider photodesorption
        for s in layeredSpeciesDict.keys():
            if s not in desorptionReactants and not layeredSpeciesDict[s].isGas:
                RRs = [layeredSpeciesDict[s]]
                try:
                    PPs = [speciesDict[s.replace("_0001", "_gas")]]
                except KeyError:
                    # Species has not gas equivalant. Print warning
                    #print("Species is surface-only: ", s)
                    continue
                srea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios,
                                yieldPD=1e-3)
                srea.layer = 1
                self.reactions.append(srea)
            else:
                continue



    # *****************
    # loop over photodesorption data to build rates
    def addCRdesorptionRates(self):

        # get a species dictionary where key=name, value=species object
        speciesDict = {x.dictname: x for x in self.species}

        # loop on data to search for yields
        for species in [x for x in self.species if not x.isGas]:
            #if species.namebase == 'H2' or species.namebase == 'Hc':
            #    continue
            # get reactants and products from the species dictionary
            RRs = [species]
            try:
                PPs = [speciesDict[species.dictname + "_gas"]]
            except KeyError:
                # Species has not gas equivalant. Print warning if generating from combinations. If limit2body, then species are printed in loadReactions.
                if not self.limit2body:
                    print("Species is surface-only: ", species.dictname)
                continue
            # create a new reaction object with the given yield
            try:
                rea = reaction(RRs, PPs, self.Ea,
                               self.barrierWidths, self.Bratios)
            except ValueError:
                # Skip
                continue
            rea.type = "CRdesorption"
            # append the reaction to the list of dust reactions
            self.reactions.append(rea)

    # *****************
    # print a warning when missing species are present
    def checkMissingSpecies(self):

        # this attribute is filled in loadKIDA() method
        if self.missingSpecies:
            print(
                "WARNING: These species are missing from the gas-phase (but present in KIDA):")
            print((sorted(list(set(self.missingSpecies)))))
            print("Add them to species file to include their gas-phase reactions.")

    # *******************
    # set reaction index and build RHS in F90 format
    def updateReactions(self):
        self.icount = 0
        # set reaction index and build RHS in F90 format
        for rea in self.reactionsAll:
            rea.idx = self.icount
            rea.buildF90RHS()
            self.icount += 1
            rea.baseIdx = rea.idx + 1  # Adjust to F90 indexing
            if rea.reactants[0].name == 'CO_0001' and rea.type == "photodesorption":
                self.CO_photodesorption = rea.baseIdx

    def updateMaskReaction(self):
        self.maskReaction[0].idx = len(self.reactionsAll)
        #self.maskReaction[0].buildF90RHS(self.reactions)
        self.maskReaction[0].baseIdx = len(self.reactionsAll) + 1  # F90 format

    # ********************
    # store reaction verbatim to file (loaded by F90 at runtime)
    def verbatimToFile(self, fname="verbatim.dat"):

        # open file to write
        fout = open(fname, "w")
        for rea in self.reactionsAll:
            # each line is a different reaction
            if isinstance(rea, reaction):
                #for ilayer in range(1, self.nlayers+1):
                ilayer = 1
                iverbatim = rea.verbatim.replace(
                    "_0001", "_%0.4i" % ilayer)
                fout.write(iverbatim + "\n")
            else:
                iverbatim = rea.verbatim
                fout.write(iverbatim + "\n")

        fout.close()

    # ********************
    # preprocess and run the F90 code with subprocess
    # clean: do make clean before compiling
    # debug: compile in debug mode, default is optimized
    # dry: don't preprocess nor compile
    # doDNF: run getFormODE creating network without gas phase formation.
    #  (i.e. reactions are there, but gas-phase species are not updated)
    def run(self, clean=True, debug=False, dry=False, run=False, library=False):

        # check unique reactions gas
        self.checkReactionsUnique(self.reactionsGas)

        # get a purged list of surface reactions using gas phase from KIDA as reference
        # (this is to avoid creating unknown surface reactions)
        # NB: Also keeps reactions with Ea from file!
        self.reactions = self.purgeUsingKIDA(self.reactions, self.mass)

        # merge dust and gas reactions
        # if self.reactionType == 'grainoble1':
        #     self.reactionsAll = self.reactions
        # else:
        self.reactionsAll = self.reactions + self.reactionsGas

        self.nReactions = len(self.reactions) + \
            len(self.reactionsGas)

        # update reactions objects (compute index and RHS)
        self.updateReactions()

        # print a recap with info on the network
        self.printReacap()

        # Adjust the mask reaction, add index and compute RHS
        self.updateMaskReaction()
        self.reactionsAll += self.maskReaction

        # dry run skips F90 preprocessing, compilation, and run
        if not dry:
            # copy verbatim rates to file
            self.verbatimToFile()

            # preprocess F90 files to replace pragmas
            self.preproc()

            if run:
                print("\ncompiling...")
                start = time.time()
                # make clean
                if clean:
                    process_output = check_output(["make", "clean"])
                    print(process_output.decode('utf-8'))

                # compile
                if debug:
                    process_output = check_output(["make", "debug"])
                    print(process_output.decode('utf-8'))
                else:
                    process_output = check_output(["make", "-j"])
                    print(process_output.decode('utf-8'))

                elapsed = time.time() - start
                print(("Compilation time: %1.3f seconds" % elapsed))

                # run executable
                print("running...")
                start = time.time()
                process_output = check_output(["./main"])
                print(process_output.decode('utf-8'))
                elapsed = time.time() - start
                print(("Run time: %1.3f seconds" % elapsed))

                # say goodbye
                print("executable done!")
            elif library:
                print("\ncompiling shared library...")
                # make clean
                if clean:
                    process_output = check_output(["make", "clean"])
                    print(process_output.decode('utf-8'))

                # compile
                process_output = check_output(["make", "sharedlib"])
                print(process_output.decode('utf-8'))
            else:
                print("WARNING: not compiling or running, only preprocessing!")
        else:
            # say goodbye
            print("WARNING: dry mode, no preprocessing and no running!")

    # *****************************
    def printReacap(self):
        nGas = len([x for x in self.species if x.isGas])
        nSurface = self.nSpecies - nGas
        message = "*************** RECAP ***************\n"
        message += "Number of species: " + str(self.nSpecies) + "\n"
        message += "Number of gas species: " \
                   + str(nGas) + "\n"
        message += "Number of surface species: " \
                   + str(nSurface) + "\n"
        message += "Number of surface layers: " \
                   + str(self.nlayers) + "\n"
        message += "Number of reactions: " + str(self.nReactions) + "\n"
        message += "Number of gas reactions: " + \
            str(len(self.reactionsGas)) + "\n"
        message += "Number of surface reactions: " + \
            str(len(self.reactions)) + "\n"

        print(message)

    # *****************
    # load reaction from KIDA database
    # fileName: KIDA database file
    def loadKIDA(self, fileName="gasNetwork.dat"):
        print(("Reading gas-phase network " + fileName))
        i = 1  # Relict layer index.
        # get a species dictionary where key=name, value=species object
        speciesDict = {x.dictname: x for x in self.species}
        layeredSpeciesDict = {x.name: x for x in self.species}

        # get list of atoms from species exploded
        availableAtoms = []
        # loop on species
        for species in self.species:
            availableAtoms += species.exploded
        # unique list, avoid atom repetitions
        availableAtoms = list(set(availableAtoms))

        print("Adding photochemistry for surface species present in KIDA network where applicable.")
        self.reactionsGas = []

        # #######################################################
        # loop on file lines
        for row in open(self.datadir + fileName, "rb"):
            srow = row.strip().decode('ascii')
            if srow == "":
                continue
            if srow.startswith("#"):
                continue
            # skip KIDA comment
            if srow.startswith("!"):
                continue

            # parse row line for reaction
            # self.species is needed to get only reactions between
            # known species
            rea = reactionGas(srow, self.species, self.mass, self.respectGasphaseLimits)

            # skip reactions with species not included in the list of species
            # but also check if reaction has some species that could be missing
            # and store them into missingSpecies class attribute
            if not rea.hasSpeciesFlag:
                self.missingSpecies += rea.getMissingSpecies(availableAtoms)
                continue


            # Check if gas-phase photodissociation:
            if rea.type == "gasphase_Av" and rea.formula == 2 and rea.kidatype in [1, 2, 3]:
                addReaction = True
                # First add photodissociation for all relevant reactions present in gas-phase (KIDA)

                # Ignore CH3OH (special treatment)
                if rea.reactants[0].name in ['CH3OH_gas', 'CH2DOH_gas', 'CH3OD_gas']: #, 'H2O_gas', 'HDO_gas', 'D2O_gas']:
                    addReaction = False
                    #continue

                # Avoid ions since these are not considered in grain surface chemistry
                pnames = [x.name for x in rea.products]
                for p in pnames:
                    if ('+' in p) or ('-' in p):
                        addReaction = False

                # Avoid ions since these are not considered in grain surface chemistry
                rnames = [x.name for x in rea.reactants]
                for r in rnames:
                    if ('+' in r) or ('-' in r):
                        addReaction = False

                # Don't add photodissociation into three species for now:
                # if len(rea.products) > 2:
                #     continue

                if addReaction:
                    # Create surface reaction
                    # get reactants and products from the species dictionary
                    RRs = []
                    PPs = []
                    for x in rea.reactants:
                        name = x.dictname.replace('_gas', '')
                        # if keyerror then the species is not in the surface chemistry network and the reaction is ignored
                        try:
                            RRs.append(layeredSpeciesDict[name + '_%0.4i' % i])
                        except KeyError:
                            addReaction = False
                    for x in rea.products:
                        name = x.dictname.replace('_gas', '')
                        # if keyerror then the species is not in the surface chemistry network and the reaction is ignored
                        try:
                            PPs.append(
                                layeredSpeciesDict[name + '_%0.4i' % i])
                        except KeyError:
                            addReaction = False
                    # create a new reaction object with the given gamma. Yield negative to signal dissociation, not desorption.
                    if addReaction:
                        try:
                            srea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios,
                                            yieldPD=-1.0, gamma=rea.gamma, alpha=rea.alpha)
                            srea.layer = 1
                            self.reactions.append(srea)
                        except ValueError:
                            continue


            # skip reactions with unknown formula
            if rea.formula not in [0, 1, 2, 3, 4, 5]:
                continue
            # append reaction found
            self.reactionsGas.append(rea)

        # #######################################################
        # NOTE: exception for CH3OH, following Karin Oberg et al (2009) branching ratios 5:1:1
        if 'CH3OH_gas' in speciesDict.keys():
            k_CH3OH_total = 1.69e-9  # KIDA total gasphase alpha.
            k_H2O_total = 8.01e-10
            #P_H2O_isrf = 5.4e-3
            #P_H2O_cr = 4.7e-3
            #P_CH3OH_isrf = P_H2O_isrf * (k_CH3OH_total/k_H2O_total)
            #P_CH3OH_cr = P_H2O_cr * (k_CH3OH_total/k_H2O_total)
            
            CH3OH_gamma = 2.76  # Mix of KIDA values

            # CH3OH -> CH2OH + H
            RRs = [layeredSpeciesDict['CH3OH_{:04d}'.format(i)]]
            PPs = [layeredSpeciesDict['CH2OH_{:04d}'.format(
                i)], layeredSpeciesDict['H_{:04d}'.format(i)]]
            srea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios, yieldPD=-
                            1.0, gamma=CH3OH_gamma, alpha=(5.0/7.0)*k_CH3OH_total)
            srea.layer = i
            self.reactions.append(srea)

            # CH3OH -> CH3O + H
            PPs = [layeredSpeciesDict['CH3O_{:04d}'.format(
                i)], layeredSpeciesDict['H_{:04d}'.format(i)]]
            srea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios,
                            yieldPD=-1.0, gamma=CH3OH_gamma, alpha=(1.0/7.0)*k_CH3OH_total)
            srea.layer = i
            self.reactions.append(srea)

            # CH3OH -> CH3 + OH
            PPs = [layeredSpeciesDict['CH3_{:04d}'.format(
                i)], layeredSpeciesDict['OH_{:04d}'.format(i)]]
            srea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios,
                            yieldPD=-1.0, gamma=CH3OH_gamma, alpha=(1.0/7.0)*k_CH3OH_total)
            srea.layer = i
            self.reactions.append(srea)

        if 'CH3OD_gas' in speciesDict.keys():
            ####
            # CH3OD -> CH2OD + H
            RRs = [layeredSpeciesDict['CH3OD_{:04d}'.format(i)]]
            PPs = [layeredSpeciesDict['CH2OD_{:04d}'.format(
                i)], layeredSpeciesDict['H_{:04d}'.format(i)]]
            srea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios, yieldPD=-
                            1.0, gamma=CH3OH_gamma, alpha=(5.0/7.0)*k_CH3OH_total)
            srea.layer = i
            self.reactions.append(srea)

            # CH3OD -> CH3O + D
            PPs = [layeredSpeciesDict['CH3O_{:04d}'.format(
                i)], layeredSpeciesDict['D_{:04d}'.format(i)]]
            srea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios,
                            yieldPD=-1.0, gamma=CH3OH_gamma, alpha=(1.0/7.0)*k_CH3OH_total)
            srea.layer = i
            self.reactions.append(srea)

            # CH3OD -> CH3 + OD
            PPs = [layeredSpeciesDict['CH3_{:04d}'.format(
                i)], layeredSpeciesDict['OD_{:04d}'.format(i)]]
            srea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios,
                            yieldPD=-1.0, gamma=CH3OH_gamma, alpha=(1.0/7.0)*k_CH3OH_total)
            srea.layer = i
            self.reactions.append(srea)
        if 'CH2DOH_gas' in speciesDict.keys():   
            ####
            # CH2DOH -> CH2OH + D
            RRs = [layeredSpeciesDict['CH2DOH_{:04d}'.format(i)]]
            PPs = [layeredSpeciesDict['CH2OH_{:04d}'.format(
                i)], layeredSpeciesDict['D_{:04d}'.format(i)]]
            srea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios, yieldPD=-
                            1.0, gamma=CH3OH_gamma, alpha=(1.0/3.0 * 5.0/7.0)*k_CH3OH_total)
            srea.layer = i

            self.reactions.append(srea)
            # CH2DOH -> CHDOH + H
            PPs = [layeredSpeciesDict['CHDOH_{:04d}'.format(i)], layeredSpeciesDict['H_{:04d}'.format(i)]]
            srea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios, yieldPD=-
                            1.0, gamma=CH3OH_gamma, alpha=(2.0/3.0  * 5.0/7.0)*k_CH3OH_total)
            srea.layer = i
            self.reactions.append(srea)

            # CH2DOH -> CH2DO + H
            PPs = [layeredSpeciesDict['CH2DO_{:04d}'.format(
                i)], layeredSpeciesDict['H_{:04d}'.format(i)]]
            srea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios,
                            yieldPD=-1.0, gamma=CH3OH_gamma, alpha=(1.0/7.0)*k_CH3OH_total)
            srea.layer = i
            self.reactions.append(srea)

            # CH2DOH -> CH2D + OH
            PPs = [layeredSpeciesDict['CH2D_{:04d}'.format(
                i)], layeredSpeciesDict['OH_{:04d}'.format(i)]]
            srea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios,
                            yieldPD=-1.0, gamma=CH3OH_gamma, alpha=(1.0/7.0)*k_CH3OH_total)
            srea.layer = i
            self.reactions.append(srea)

    # *****************
    # load reaction from KIDA database and compare with a set of reactions
    def purgeUsingKIDA(self, reactions, atomsMass, fileName="gasNetwork.dat"):
        print(("Comparing gas-phase with surface network " + fileName))

        # dict of surface reactions, key=exploded hash, value = known in the gas phase
        # the condition in the first list element will set as true (i.e. known) reactions that
        # are not in the gas phase but they have known barrier. it will also keep true
        # reactions that are not 2body
        hashList = {x.hashExploded: x.barrierFromFile or not x.type.startswith("2body")
                    for x in reactions}

        # loop on file lines
        for row in open(self.datadir + fileName, "rb"):
            srow = row.strip().decode('ascii')
            if srow == "":
                continue
            if srow.startswith("#"):
                continue
            # skip KIDA comment
            if srow.startswith("!"):
                continue

            # parse row line for reaction
            # self.species is needed to get only reactions between
            # known species, atomMass is to parse species
            rea = reactionGas(srow, self.species, atomsMass)

            # if gas reaction exploded hash is in the surface reaction exploded hashes
            # sets to known
            if rea.hashExploded in hashList:
                hashList[rea.hashExploded] = True

        updatedReactions = []
        # loop on reactions to get found reactions
        for rea in reactions:
            # check if reaction found
            if hashList[rea.hashExploded]:
                updatedReactions.append(rea)
            # if Ea from file, also keep:
            elif rea.barrierFromFile:
                updatedReactions.append(rea)
            else:
                print("removing:", rea.verbatim)
                #print "MAYBE ADD TEXTFILE OF REMOVED REACTIONS HERE"
                continue

        # return updated list
        return updatedReactions

    # ******************
    # do pre-processing on F90 files
    def preproc(self):
        # prepare ODE
        if self.nlayers == 2:
            if self.include_H2 and self.H2spin:
                copyfile('./f90templates/kemimo_ode_include_H2.f90', './kemimo_ode.f90')
                copyfile('./f90templates/kemimo_include_H2.f90', './kemimo.f90')
            elif self.include_H2 and not self.H2spin:
                copyfile('./f90templates/kemimo_ode_include_H2_nospin.f90', './kemimo_ode.f90')
                copyfile('./f90templates/kemimo_include_H2.f90', './kemimo.f90')
            else:
                copyfile('./f90templates/kemimo_ode_fixed_H2.f90', './kemimo_ode.f90')
                copyfile('./f90templates/kemimo_fixed_H2.f90', './kemimo.f90')
            copyfile('./f90templates/kemimo_flux_threephase.f90',
                     './kemimo_flux.f90')
        else:
            copyfile('./f90templates/kemimo_ode_twophase.f90',
                     './kemimo_ode.f90')
            copyfile('./f90templates/kemimo_twophase.f90',
                     './kemimo.f90')
            copyfile('./f90templates/kemimo_flux_twophase.f90',
                     './kemimo_flux.f90')

                
        # copy remaining f90 files if missing:
        files = ["kemimo_commons.f90", "kemimo_sticking.f90", \
            "kemimo_gas_rates.f90", "kemimo_dust_rates.f90", "kemimo_reactionarray.f90",\
                "kemimo_rates.f90", "kemimo_swappingrates.f90"]

        for fname in files:
            if not os.path.exists(fname):
                copyfile('./f90templates/%s' % fname, './%s' % fname)            

        # prepare COMMONS
        doPP("kemimo_commons.f90", {"ARRAYSIZE": self.getArraySizes(),
                                    "IDXLIST": self.getIdxList(),
                                    "SPECIESNAMES": self.getSpeciesNamesArray()})
        # prepare REACTIONARRAY
        doPP_datfile("reactionarray.dat", self.getReactionArray())
        doPP_datfile("jacarray.dat", self.getJacArray())


        # prepare SWAPPINGRATES
        if self.doSwap:
           doPP("kemimo_swappingrates.f90", {
             "SWAPPINGRATES": self.getSwappingRates()})
        
        # prepare STICKINGRATES
        doPP("kemimo_sticking.f90", {
             "STICKING": self.getSticking()})


        # prepare RATES
        doPP("kemimo_dust_rates.f90", {
             "RATES": self.getRatesF90(dustOnly=True)})
        doPP("kemimo_gas_rates.f90", {"RATES": self.getRatesF90(gasOnly=True)})

        # prepare Python interface
        doPP_Python("pykemimo.py", {"ARRAYSIZE": self.getArraySizes(python=True),
                                    "IDXLIST": self.getIdxList(),
                                    "SPECIESNAMES": self.getSpeciesNamesArray(python=True)})

    # ************************
    # load species from file. Species are names space- or tab-separated
    # one or more per line. A line with ENDIFLE stops reading.
    # store data into a class attribute list
    def loadSpecies(self, fname, fname2):
        self.speciesNames = []
        # Load fname:
        # loop on file lines
        for row in open(fname, "rb"):
            # strip and replace tabs with spaces
            srow = row.strip().decode('utf-8').replace("\t", " ")
            # ENDFILE line stops reading
            if srow == "ENDFILE":
                break
            if srow == "":
                continue
            if srow.startswith("#"):
                continue

            # read species as gas-phase
            self.speciesNames += [x +
                                  "_gas" for x in srow.split(" ") if (x != "")]

        # If not allSpeciesToDust, then load fname2
        if not self.allSpeciesToDust:
            # create dust-phase names form gas-phase names, skip ions and electrons
            dustSpeciesNames = []
            for row in open(fname2, "rb"):
                # strip and replace tabs with spaces
                srow = row.strip().decode('utf-8').replace("\t", " ")
                # ENDFILE line stops reading
                if srow == "ENDFILE":
                    break
                if srow == "":
                    continue
                if srow.startswith("#"):
                    continue

                # read species as dust-phase
                for x in srow.split(" "):
                    if (x != ""):
                        if (not "+" in x and x != "E_gas"):
                            if (not "-" in x):
                                dustSpeciesNames.append(x)

            self.speciesNames += dustSpeciesNames
            print(("Species names loaded from " + fname + ", and " + fname2))

        # Else we add all species to dust network
        else:
            # create dust-phase names form gas-phase names, skip ions and electrons
            dustSpeciesNames = []
            for x in self.speciesNames:
                if (not "+" in x and x != "E_gas"):
                    if (not "-" in x):
                        dustSpeciesNames.append(x.replace("_gas", ""))

            self.speciesNames += dustSpeciesNames

            print(("Species names loaded from " + fname))

    def _read_file_lines(self, fname):
        """Generator that yields processed lines from a file."""
        with open(fname, "rb") as file:
            for row in file:
                srow = row.strip().decode('ascii').replace("\t", " ")
                if srow == "ENDFILE":
                    break
                if srow == "":
                    continue
                if srow.startswith("#"):
                    continue
                yield srow.split("#")[0]

    def _parse_ascii_data(self, line, data_type):
        """Parse ascii data from a line based on the data type."""
        arow = [x for x in line.split(" ") if x != ""]
        print(arow)
        if data_type == 'Ea':
            return {"reactants": arow[:2], "products": arow[2:-1], "Ea": float(arow[-1])}
        elif data_type == 'width':
            return {"reactants": arow[:2], "products": arow[2:-1], "width": float(arow[-1])*1e-8}
        elif data_type == 'ratio':
            return {"reactants": arow[:2], "products": arow[2:-1], "ratio": float(arow[-1])}
        elif data_type == 'yield':
            return {"reactants": arow[:1], "products": arow[1:-1], "yield": float(arow[-1])}
        elif data_type == 'gamma':
            return {"reactants": arow[:1], "products": arow[1:-1], "gamma": float(arow[-1])}
        elif data_type == 'alpha':
            return {"reactants": arow[:1], "products": arow[1:-1], "alpha": float(arow[-1])}
        elif data_type == 'Ebind':
            return {"name": arow[0], "Ebare": float(arow[1]), "Eice": float(arow[2])}
        elif data_type == 'atom':
            return {"name": arow[0], "mass": float(arow[1])}
        elif data_type == 'dH':
            return {"name": arow[0], "dH": float(arow[1]) * 120.274} # Convert from KJ/mol to K
        elif data_type == 'species':
            return arow[0]
        else:
            raise ValueError(f"Unknown data type {data_type}")
        

    def loadBarriers(self, fname):
        self.Ea = []
        """Load 2body reaction barriers from file."""
        for line in self._read_file_lines(fname):
            self.Ea.append(self._parse_ascii_data(line, 'Ea'))
        print(f"Barriers loaded from {fname}")
    
    def loadEbinds(self, fname):
        self.Eice = dict()
        self.Ebare = dict()
        """Load binding energies from file."""
        for line in self._read_file_lines(fname):
            name, Eb, Ei = self._parse_ascii_data(line, 'Ebind').values()
            self.Eice[name] = Ei
            self.Ebare[name] = Eb
        print(f"Binding energies loaded from {fname}")

    def loadAtoms(self, fname):
        self.mass = dict()
        """Load species masses from file."""
        for line in self._read_file_lines(fname):
            name, mass = self._parse_ascii_data(line, 'atom').values()
            self.mass[name] = mass
        print(f"Species masses loaded from {fname}")

    def loadDeltaHs(self, fname):
        self.deltaH = dict()
        """Load enthalpies from file."""
        for line in self._read_file_lines(fname):
            name, dH = self._parse_ascii_data(line, 'dH').values()
            self.deltaH[name] = dH
        print(f"Enthalpies loaded from {fname}")

    def loadCheckBarriers(self, fname):
        self.speciesCheckBarrier = []
        """Load 2body reaction barriers from file."""
        for line in self._read_file_lines(fname):
            self.speciesCheckBarrier.append(self._parse_ascii_data(line, 'species'))
        print(f"Check barriers loaded from {fname}")

    def loadBarrierWidths(self, fname):
        self.barrierWidths = []
        """Load 2body reaction barrier width from file."""
        for line in self._read_file_lines(fname):
            self.barrierWidths.append(self._parse_ascii_data(line, 'width'))
        print(f"Barrier widths loaded from {fname}")

    def loadBranchingRatios(self, fname):
        self.Bratios = []
        """Load 2body branching ratios from file."""
        for line in self._read_file_lines(fname):
            self.Bratios.append(self._parse_ascii_data(line, 'ratio'))
        print(f"Branching ratios loaded from {fname}")

    def loadYieldsPD(self, fname):
        self.yieldPD = []
        """Load photodesorption yields from file."""
        for line in self._read_file_lines(fname):
            self.yieldPD.append(self._parse_ascii_data(line, 'yield'))
        print(f"Photodesorption yields loaded from {fname}")


    # ***********************
    # return a reaction object with the given index idx (zero-based)
    def getReactionByIndex(self, idx):

        # loop on reactions to find reaction with index idx
        for rea in self.reactionsAll:
            if rea.idx == idx:
                return rea

        # if no reaction found rises ERROR
        print(("ERROR: cannot get reaction with index " + idx))
        sys.exit()

    # ****************
    # show table (see header below) with species properties, if fileName is set save to file
    # fileName: name of the file where to save the table, if None print to stdout
    # sortby: sorting criterion based on species attribute
    def showSpecies(self, fileName=None, sortby="name"):

        # table header
        message = "----------------------\n"
        message += tabrow(["species", "Eice/K", "Ebare/K",
                           "dH/K", "mass/amu"]).strip() + "\n"
        message += "----------------------\n"
        # loop on species sorted by sortby criterion
        for species in sorted(self.species, key=lambda x: getattr(x, sortby)):
            message += tabrow([species.name, species.Eice, species.Ebare,
                               species.dH, species.mass]).strip() + "\n"
        message += "\n\n"

        # write list of species for gnuplot users
        message += "List of species: " + \
            " ".join([x.name for x in self.species]) + "\n\n"

        # write list of species in LaTeX, mass sorted
        sortedSpecies = sorted(
            self.species, key=lambda xx: xx.mass if isinstance(xx.mass, float) else -1.0)
        message += "\\newcommand{\\gas}{_{(g)}}\n"
        message += "List of species (LaTeX): " + \
            ", ".join([x.nameLatex for x in sortedSpecies]) + "\n\n"

        # write to stdout or to file depending on the filename option
        if fileName is None:
            # print message to screen
            print(message)
        else:
            # save message to file
            fout = open(fileName, "w")
            fout.write(message)
            fout.close()
            print(("Table with species saved to " + fileName))

    # ****************
    # save a file with a latex table of the species employed
    # fileName: name of the file where to save the table, if None print to stdout
    # sortby: sorting criterion based on species attribute
    def saveSpeciesLatexTable(self, fileName="species_latex.out"):

        # table header
        message = "#############\n"
        message += "\\begin{table}\n"
        message += "\\begin{tabular}{llll}\n"
        message += "\\hline\n"
        message += " & ".join(["Species (gas)", "Species (dust)",
                               "$E_{bind}$/K", "$\\Delta H$/K"]) + "\\\\\n"
        message += "\\hline\n"

        # species dictionary, key=name, value=species
        speciesNames = {x.name: x for x in self.species}

        # sort species
        sortedSpecies = sorted(
            self.species, key=lambda x: x.mass if isinstance(x.mass, float) else -1.0)

        # loop on gas species sorted by sortby criterion
        for species in [x for x in sortedSpecies if x.isGas]:
            # get the surface name
            surfaceName = species.name.replace("_gas", "")
            surfaceSpecies = ""
            cationName = species.name.replace("_gas", "+_gas")
            cationSpecies = ""
            # if surface is in list get its latex name
            if surfaceName in speciesNames:
                surfaceSpecies = speciesNames[surfaceName].nameLatex
            if cationName in speciesNames:
                cationSpecies = speciesNames[cationName].nameLatex
            message += " & ".join([species.nameLatex, surfaceSpecies, cationSpecies,
                                   latexInt(species.Eice), latexExp(species.dH)]) + "\\\\\n"
        message += "\\hline\n"
        message += "\\end{tabular}\\caption{Your caption here}\\label{your_label_here}\n"
        message += "\\end{table}\n"

        # save message to file
        fout = open(fileName, "w")
        fout.write(message)
        fout.close()
        print(("LaTeX table with species saved to " + fileName))

    # ****************
    # show table (see header below) with DUST-phase reactions, if fileName is set save to file
    # fileName: name of the file where to save the table, if None print to stdout
    # sortby: sorting criterion based on reaction attribute
    def showReactions(self, fileName=None):

        # header, where lmax=[40, 20] means 40 characters to verbatim and 20 to other columns
        message = "----------------------\n"
        message += tabrow(["#", "type", "dH/K", "Ea/K", "barrier width /Å", "yield", "Pdelta", "Ea from file", "Branching ratio", "Branching ratio from file"],
                          lmax=[40, 20]).strip() + "\n"
        message += "----------------------\n"

        # loop on reactions sorted by sortby criterion
        for rea in sorted(self.reactions, key=lambda x: x.dH if isinstance(x.dH, float) else 0.0):
            message += tabrow([rea.verbatim, rea.type, rea.dH, rea.Ea, rea.barrierWidth*1e8,
                               rea.yieldPD, rea.Pdelta, rea.barrierFromFile, rea.Bratio, rea.bRatioFromFile],
                              lmax=[40, 20]).strip() + "\n"

        # write to stdout or to file depending on the filename option
        if fileName is None:
            print(message)
        else:
            fout = open(fileName, "w")
            fout.write(message)
            fout.close()
            print(("Table with reactions saved to " + fileName))

    # ****************
    # show table (see header below) with GAS-phase reactions, if fileName is set save to file
    # fileName: name of the file where to save the table, if None print to stdout
    # sortby: sorting criterion based on reaction attribute
    def showReactionsGas(self, fileName=None):

        # print a table with reactions
        message = "----------------------\n"
        message += tabrow(["#", "type", "Tmin", "Tmax",
                           "rateF90"], lmax=[50, 20]).strip() + "\n"
        message += "----------------------\n"

        # loop on gas-phase reaction sorted by sortby
        for rea in sorted(self.reactionsGas, key=lambda x: x.type if isinstance(x.type, str) else ""):
            message += tabrow([rea.verbatim, rea.type, rea.Tmin, rea.Tmax, rea.krateF90],
                              lmax=[50, 20]).strip() + "\n"

        # write to stdout or to file depending on the filename option
        if fileName is None:
            print(message)
        else:
            fout = open(fileName, "w")
            fout.write(message)
            fout.close()
            print(("Table with gas reactions saved to " + fileName))

    # ****************
    # get 2body reactions from Ea, create adsorption/desorption from dust species.
    def loadReactions(self):
        speciesDict = {x.dictname: x for x in self.species}
        temp = []
        combinations = []
        combinationNames = []
        reactions = []
        print("Loading reactions.")
        
        # ------------------------------------------
        # Add every species to combinations list, to create adsorption/desorption reactions for these:
        for sp in self.species:
            temp.append({"mols": [sp], "exploded": sp.exploded})

        # clean
        for comb in temp:
            spsp = [comb["mols"][0].dictname]
            if spsp in combinationNames:
                print("Already found, ", spsp)
                continue
            combinationNames.append(spsp)
            combinations.append(comb)

        print("Found adsorption/desorption combinations.")
        # ------------------------------------------
        ncpu=int(mp.cpu_count())
        # We limit our subprocess count to 4. May not be worth it at all in this routine (inherited from findReactions)
        ncpu=np.min([ncpu, 4])

        if self.multiprocessing:
            p=Pool(ncpu)
            # loop on combinations to find reactions
            result = p.map(findReactions_parser, list(zip(
                combinations, itertools.repeat(combinations), itertools.repeat(self.Ea), itertools.repeat(self.barrierWidths), itertools.repeat(self.Bratios), itertools.repeat(self.include_H2))))
            p.close()
            p.join()
        else:
            result = []
            for i in range(len(combinations)):
                res = findReactions_parser(list([combinations[i], combinations, self.Ea, self.barrierWidths, self.Bratios, self.include_H2]))
                result.append(res)

        for r in result:
            if r is None:
                continue
            for rr in r:
                reactions.append(rr)


        print("Created reactions from combinations")
        # ------------------------------------------
        # create spin variants for species where this is included in KIDA, as we need this information when loading from Ea.dat:
        spinSpecies = {}
        for s in self.species:
            # Skip gas spin species
            if s.isGas:
                continue
            # check the name starts with "*_", where * is *lower* character
            if s.name[0].isalpha() and s.name[0].islower() and s.name[1] == '_':
                name = s.name.replace('_0001', '')
                try:
                    spinSpecies[name[2:]].append(name)
                except KeyError:
                    spinSpecies[name[2:]] = [name]

        print("Detected the follow spin-variants in dust network:", spinSpecies)
        # ------------------------------------------
        combinationNames = []
        # Store a list for the species that are present in Ea.dat, but not found in speciesDust.dat
        missingSpecies = []
        missingGasSpecies = []
        # Add Ea reactions to combinations
        for ni, entry in enumerate(self.Ea):
            isSpinSpecies = False
            nSpinIsomers = 0
            rs = list(entry["reactants"])
            ps = list(entry["products"])
            entry_spinSpecies = []
            spinDesignator1, spinDesignator2 = [], []
            for i, species in enumerate(rs+ps):
                if species in spinSpecies.keys():
                    isSpinSpecies = True
                
                    nSpinIsomers = len(spinSpecies[species])
                    # If one spinSpecies, simply change name and move on. 
                    if nSpinIsomers == 1:
                        if i < 2:
                            # i < 2: update reactants:
                            rs[i] = spinSpecies[species][0]
                        else:
                            # update products:
                            ps[i-2] = spinSpecies[species][0]
                    elif nSpinIsomers == 2:
                        if len(entry_spinSpecies) == 0:
                            rs2, ps2 = list(entry["reactants"]), list(entry["products"])
                        entry_spinSpecies.append(species)
                        if i < 2:
                            # i < 2: update reactants:
                            rs[i] = spinSpecies[species][0]
                            rs2[i] = spinSpecies[species][1]
                            spinDesignator1.append(rs[i].strip("_"+species))
                            spinDesignator2.append(rs2[i].strip("_"+species))
                        else:
                            # update products:
                            ps[i-2] = spinSpecies[species][0]
                            ps2[i-2] = spinSpecies[species][1]
                            spinDesignator1.append(ps[i-2].strip("_"+species))
                            spinDesignator2.append(ps2[i-2].strip("_"+species))
                    else:
                        print("More than two spin variants of species. Not supported for now")
                else:
                    continue

            # Perform check that we have no spin conversion:
            if len(spinDesignator1) > 1:
                if spinDesignator1[0] != spinDesignator1[1]:
                    print("Error. A spin conversion occured. Check that this is not problematic.")
                    print(entry, spinDesignator1)
            if len(spinDesignator2) > 1:
                if spinDesignator2[0] != spinDesignator2[1]:
                    print("Error. A spin conversion occured. Check that this is not problematic.")
                    print(entry, spinDesignator2)


            # Perform check to see if reaction is already present
            spsp = "__".join(sorted([s for s in rs])+sorted([s for s in ps]))
            if spsp not in combinationNames:
                combinationNames.append(spsp)
            else:
                print("Reaction already in network. Check Ea.dat for duplicate reactions. Skipping: ", ni, entry)
                continue

            # Get mol instances for the species:
            try:
                reactants = np.array([speciesDict[i] for i in rs])
                products = np.array([speciesDict[i] for i in ps])
                # check if H2 in products, due to special treatment:
                if self.include_H2 or len(ps) < 2:
                    pass
                else:
                    # If H2 in products and not including H2 explicitly, then release to gasphase:
                    if any([x == 'o_H2' for x in ps]) or any([x == 'p_H2' for x in ps]):
                        products = []
                        for x in ps:
                            if x == 'o_H2' or x == 'p_H2':
                                products.append(speciesDict[x+'_gas'])
                            else:
                                products.append(speciesDict[x])
                        products = np.array(products)



                # sort by mass:
                reactantsMass = [x.mass for x in reactants]
                productsMass = [x.mass for x in products]
                reactants = reactants[np.argsort(reactantsMass)]
                products = products[np.argsort(productsMass)]
            except KeyError:
                # species not in network. Skipping this reaction
                for i in rs+ps:
                    if i not in speciesDict.keys():
                        missingSpecies.append(i)
                continue
            
            try:
                gasProducts = [speciesDict[i+'_gas'] for i in ps]
            except KeyError:
                # We move on if a species is missing from the gasphase, but store the name.
                for i in ps:
                    if i+"_gas" not in speciesDict.keys():
                        missingGasSpecies.append(i)
                gasProducts = []
                pass

            # create a reaction object
            try:
                rea = reaction(reactants, products, self.Ea, self.barrierWidths, self.Bratios)
                rea.barrierFromFile = True

                # add reactions of type 2body_gas (but not 2body_gas_gas)
                if len(gasProducts) == 1:
                    reaGas = reaction(reactants, gasProducts,
                                      self.Ea, self.barrierWidths, self.Bratios)
                    reaGas.barrierFromFile = True
            except ValueError:
                # Skip
                print("ValueError in reaction creation")
                for i in reactants+products:
                    print(i.dictname)
                continue
            except TypeError:
                print("TypeError in reaction creation")
                print((reactants[0].dictname,
                        products[0].dictname))
                print((reactants[0].charge,
                        products[0].charge))

            # Add to list
            reactions.append(rea)
            if 'reaGas' in locals():
                reactions.append(reaGas)
                del reaGas
            
            # ------------------------------------------
            # If another spin variant exists:
            if isSpinSpecies and nSpinIsomers == 2:
                # Perform check to see if reaction is already present
                spsp = "__".join(
                    sorted([s for s in rs2])+sorted([s for s in ps2]))
                if spsp not in combinationNames:
                    combinationNames.append(spsp)
                else:
                    print("Reaction already in network. Check Ea.dat for duplicate reactions. Skipping: ", spsp,  ni, entry)
                    continue

                # Get mol instances for the species:
                try:
                    reactants = np.array([speciesDict[i] for i in rs2])
                    products = np.array([speciesDict[i] for i in ps2])

                    # check if H2 in products, due to special treatment:
                    if self.include_H2 or len(ps2) < 2:
                        pass
                    else:
                        # If H2 in products and not including H2 explicitly, then release to gasphase:
                        if any([x == 'o_H2' for x in ps2]) or any([x == 'p_H2' for x in ps2]):
                            products = []
                            for x in ps2:
                                if x == 'o_H2' or x == 'p_H2':
                                    products.append(speciesDict[x+'_gas'])
                                else:
                                    products.append(speciesDict[x])
                            products = np.array(products)


                    # sort by mass:
                    reactantsMass = [x.mass for x in reactants]
                    productsMass = [x.mass for x in products]
                    reactants = reactants[np.argsort(reactantsMass)]
                    products = products[np.argsort(productsMass)]
                except KeyError:
                    for i in rs+ps:
                        if i not in speciesDict.keys():
                            missingSpecies.append(i)
                    continue
                try:
                    gasProducts = np.array([speciesDict[i+'_gas'] for i in ps2])
                    # sort by mass:
                    productsMass = [x.mass for x in gasProducts]
                    gasProducts = gasProducts[np.argsort(productsMass)]
                except KeyError:
                    # We move on if a species is missing from the gasphase, but store the name.
                    for i in ps2:
                        if i+"_gas" not in speciesDict.keys():
                            missingGasSpecies.append(i)
                    gasProducts = []
                    pass

                # create a reaction object
                try:
                    rea = reaction(reactants, products, self.Ea,
                                self.barrierWidths, self.Bratios)
                    rea.barrierFromFile = True

                    # add reactions of type 2body_gas (but not 2body_gas_gas)
                    if len(gasProducts) == 1:
                        reaGas = reaction(reactants, gasProducts,
                                        self.Ea, self.barrierWidths, self.Bratios)
                        reaGas.barrierFromFile = True
                except ValueError:
                    # Skip
                    print("ValueError in reaction creation")
                    for i in reactants+products:
                        print(i.dictname)
                    continue
                except TypeError:
                    print("TypeError in reaction creation")
                    print((reactants[0].dictname,
                        products[0].dictname))
                    print((reactants[0].charge,
                        products[0].charge))

                # Add to list
                reactions.append(rea)
                if 'reaGas' in locals():
                    reactions.append(reaGas)
                    del reaGas

            # ------------------------------------------

        print("Added 2body reactions")
        print("Following species are present in Ea.dat, but missing from either species.dat or speciesDust.dat: ")
        print(np.array(missingSpecies))
        print(np.array(missingGasSpecies))

        # remove barrierless reactions (with some reactants, default:[H2, H2O, O2, CH4, ...])
        reactions = self.checkBarriers(reactions, self.speciesCheckBarrier)
        # Results:
        self.reactions = reactions

    # ****************
    # use combinatorics and endothermicity to create reactions for the DUST-PHASE
    def findReactions(self):
        combinations = []
        combinationNames = []
        print("Finding reactions.")
        # prepare [species] and [species, species] combinations (with exploded)
        # loop on species
        tempspecies = copy.deepcopy(self.species)
        ncpu = int(mp.cpu_count())
        # We limit our subprocess count to 36.
        if self.multiprocessing:
            ncpu = np.min([ncpu, 20])
            print("Multiprocessing with %i threads" % ncpu)

            p = Pool(ncpu)
            result = p.map(findCombinations_parser, list(zip(
                tempspecies, itertools.repeat(tempspecies))))
            print("Found combinations")
        else:
            result = []
            for i in range(len(tempspecies)):
                res = findCombinations_parser(tempspecies[i], tempspecies)
                result.append(res)

        for r in result:
            if r is None:
                continue
            for rr in r:
                if len(rr["mols"]) > 1:
                    spsp = sorted(
                        [rr["mols"][0].dictname, rr["mols"][1].dictname])
                    if spsp in combinationNames:
                        continue
                    combinationNames.append(spsp)
                    combinations.append(rr)
                else:
                    spsp = [rr["mols"][0].dictname]
                    if spsp in combinationNames:
                        continue
                    combinationNames.append(spsp)
                    combinations.append(rr)

        print("Matched combinations")
        reactions = []
        # loop on combinations to find reactions
        if self.multiprocessing:
            result = p.map(findReactions_parser, list(zip(
                combinations, itertools.repeat(combinations), itertools.repeat(self.Ea), itertools.repeat(self.barrierWidths), itertools.repeat(self.Bratios), itertools.repeat(self.include_H2))))
            p.close()
            p.join()
        else:
            result = []
            for i in range(len(combinations)):
                res = findCombinations_parser(combinations[i], combinations, self.Ea, self.barrierWidths, self.Bratios, self.include_H2)
                result.append(res)

        for r in result:
            if r is None:
                continue
            for rr in r:
                reactions.append(rr)

        print("Created reactions from combinations")

        speciesDict = {x.dictname: x for x in self.species}

        # remove barrierless reactions (with some reactants, default:[H2, H2O, O2, CH4, ...])
        reactions = self.checkBarriers(reactions, self.speciesCheckBarrier)

        # Now check branching ratios:
        #reactions = self.checkBranchingRatios(reactions)
        # Results:
        self.reactions = reactions

    # ******************************
    # compute reaction rates
    def computeRates(self):

        # loop on reaction to compute rates evaporation and tunnelling probabilities
        for rea in self.reactions:
            rea.rate()

    # ************************
    # get array size parameters to be replaced in commons.f90
    def getArraySizes(self, python=False):
        # ---------------------------------------
        # Find mantle and surface boundaries:
        ilower_surface = 1e6
        iupper_surface = 0

        for s in self.species:
            if s.isGas:
                continue
            if s.name == 'dummy' or s.namebase == 'mask':
                continue
            else:
                if s.layer == 1:
                    if int(s.idx) < ilower_surface:
                        ilower_surface = int(s.idx)+1
                    if int(s.idx) >= iupper_surface:
                        iupper_surface = int(s.idx)+1
        offset = self.nDustSpecies*(self.nlayers-1)+3

        if self.nlayers < 3:
            ilower_mantle = ilower_surface + offset
            iupper_mantle = iupper_surface + offset

            # 2 for mlys and dummy.
            nmols = self.nGasSpecies + (self.nlayers*self.nDustSpecies) + 3
            nrea = self.nReactions
            # get number of species and reactions as commons
            arraySize = "!number of species\n"
            arraySize += "integer,parameter::nmols=" + str(nmols) + "\n"
            arraySize += "!number of unique species\n"
            arraySize += "integer,parameter::nmolsu=" + \
                str(self.nGasSpecies + self.nDustSpecies + 3) + "\n"
            arraySize += "!number of dust species\n"
            arraySize += "integer,parameter::nmols_dust=" + \
                str((self.nlayers*self.nDustSpecies)) + "\n"
            arraySize += "!number of reactions (nlayer*dust+gas)\n"
            arraySize += "integer,parameter::nrea=" + str(nrea+1) + "\n"
            arraySize += "!number of unique reactions (dust+gas)\n"
            arraySize += "!number of dust-phase reactions\n"
            arraySize += "integer,parameter::nreadust=" + \
                str(len(self.reactions)) + "\n"
            arraySize += "!number of gas-phase reactions\n"
            if python:
                arraySize += "integer,parameter::nreagas=self.nrea-self.nreadust\n\n"
            else:
                arraySize += "integer,parameter::nreagas=nrea-nreadust\n\n"
            arraySize += "!monolayer thickness of each layer in model: \n"
            arraySize += "integer,parameter::layerThickness=" + \
                str(int(self.layerThickness))+"\n"
            arraySize += "!idx for mantle and surface species: \n"
            arraySize += "integer,parameter::surface_start=" + \
                str(int(ilower_surface))+"\n"
            arraySize += "integer,parameter::surface_end=" + \
                str(int(iupper_surface))+"\n"
            arraySize += "integer,parameter::mantle_start=" + \
                str(int(ilower_mantle))+"\n"
            arraySize += "integer,parameter::mantle_end=" + \
                str(int(iupper_mantle))+"\n"
            arraySize += "!do swapping? \n"
            if self.doSwap:
                arraySize += "logical,parameter::doSwap=.true.\n"
            else:
                arraySize += "logical,parameter::doSwap=.false.\n"

        else:
            if not self.nlayers == 7:
                print("Error, number of layers not supported. Should be: 2,3, or 6")
            offset = self.nDustSpecies
            # Create arrays for 7-phase model
            mstarts, mends = [], []
            for i in range(1, self.nlayers):
                mstarts.append(ilower_surface + offset*i + 3)
                mends.append(iupper_surface + offset*i + 3)

            ilower_mantle = ilower_surface + offset + 3
            iupper_mantle = iupper_surface + offset*(self.nlayers-1) + 3

            # Thickness of each mantle layer
            mantleThickness = "integer, dimension(5) :: Mthickness = (/ 10, 10, 10, 10, 10 /)\n"
            # Array for mstarts, mends:
            mstarts_str = "integer, dimension(5) :: mstarts = (/ {:d}, {:d}, {:d}, {:d}, {:d} /)\n".format(
                *mstarts)
            mends_str = "integer, dimension(5) :: mends = (/ {:d}, {:d}, {:d}, {:d}, {:d} /)\n".format(
                *mends)

            # 2 for mlys and dummy.
            nmols = self.nGasSpecies + (self.nlayers*self.nDustSpecies) + 3
            nrea = self.nReactions
            # get number of species and reactions as commons
            arraySize = "!number of species\n"
            arraySize += "integer,parameter::nmols=" + str(nmols) + "\n"
            arraySize += "!number of unique species\n"
            arraySize += "integer,parameter::nmolsu=" + \
                str(self.nGasSpecies + self.nDustSpecies + 3) + "\n"
            arraySize += "!number of dust species\n"
            arraySize += "integer,parameter::nmols_dust=" + \
                str(self.nlayers*self.nDustSpecies) + "\n"
            arraySize += "!number of reactions (nlayer*dust+gas)\n"
            arraySize += "integer,parameter::nrea=" + str(nrea+1) + "\n"
            arraySize += "!number of unique reactions (dust+gas)\n"
            arraySize += "integer,parameter::nreau=" + \
                str(nrea) + "\n"
            arraySize += "!number of dust-phase reactions\n"
            arraySize += "integer,parameter::nreadust=" + \
                str(len(self.reactions)) + "\n"
            arraySize += "!number of gas-phase reactions\n"
            if python:
                arraySize += "integer,parameter::nreagas=self.nrea-self.nreadust\n\n"
            else:
                arraySize += "integer,parameter::nreagas=nrea-nreadust\n\n"
            arraySize += "!monolayer thickness of each layer in model: \n"
            arraySize += "integer,parameter::layerThickness=" + \
                str(int(self.layerThickness))+"\n"
            arraySize += "!idx for mantle and surface species: \n"
            arraySize += "integer,parameter::surface_start=" + \
                str(int(ilower_surface))+"\n"
            arraySize += "integer,parameter::surface_end=" + \
                str(int(iupper_surface))+"\n"
            arraySize += "integer,parameter::mantle_start=" + \
                str(int(ilower_mantle))+"\n"
            arraySize += "integer,parameter::mantle_end=" + \
                str(int(iupper_mantle))+"\n"
            arraySize += "! 7-phase model arrays: \n"
            arraySize += mantleThickness
            arraySize += mstarts_str
            arraySize += mends_str

        try:
            arraySize += "integer,parameter:: CO_desorption_idx=" + str(self.CO_photodesorption) + "\n"
        except AttributeError:
            # exception in case CO photodesorption is not present.
            pass

        return arraySize

    # ************************
    # return a F90 array with names of the species
    def getSpeciesNamesArray(self, python=False):

        # find largest string length, since all the strings in a parameter array
        # must have the same length when gfortran is employed
        maxlen = max([len(x.name) for x in self.species])

        # function that add missing spaces at the end of the arguments to reach
        # maxlen characters
        def fillSpaces(arg, maxlen):
            return arg + " " * (maxlen - len(arg))

        if not python:
            return "character(len=maxVerbatimSize),parameter::speciesNames(nmolsu) = (/" \
                   + ", &\n".join(["\"" + fillSpaces(x.name, maxlen)
                                   + "\"" for x in self.species]) + "/)"
        else:
            return "self.speciesNames = (" \
                   + ", &\n".join(["\"" + fillSpaces(x.name, maxlen)
                                   + "\"" for x in self.species]) + ")"


    # ************************
    # get a matrix with reaction idx for each species to be used in the JAC loop
    def getJacArray(self):
        print(" Running getJacArray ")
        dummy = str(-1)  # dummy species idx + 1 for F90
        # Size of array (fill with dummies)
        arraySize = len(self.reactionsAll) - 1 # Remove mask reactions removed
        # nSpecies:
        nSpecies = self.nGasSpecies + self.nlayers*self.nDustSpecies
        # Create array 
        allIdx = [[] for _ in range(nSpecies)]
        # Loop in reactions to check reactants and products 
        for i, rea in enumerate(self.reactionsAll):
            for rp in rea.reactants:
                if "dummy" in rp.name or "mask" in rp.name:
                    continue
                # Append reaction to all rp.idx for the reactants and products
                try:
                    allIdx[rp.idx].append(str(rea.idx+1))
                except IndexError:
                    print("not in Jac array:", rp.idx, rp.name)
        
        # fill list, use dummy as fill value (to be ignored in F90)
        out = []
        for i in allIdx:
            [i.append(str(dummy)) for _ in range(arraySize - len(i))] 
            out.append(" ".join(i))
        
        return out

    # ************************
    # get a matrix with reactants and products indexes to be replaced in kemimo_reactionarray.dat
    def getReactionArray(self):
        print(" Running getReactionArray ")
        maxR = 4  # max number of reactants
        maxP = 4  # max number of products
        dummy = str(self.dummy_idx+1)  # dummy species idx + 1 for F90
        allIdx = []
        offset = self.nDustSpecies*(self.nlayers-1)+3
        # loop on reactions to store indexes
        for rea in self.reactionsAll:
            if rea.layer < 0:
                continue
            if isinstance(rea, reaction):
                for ilayer in range(1, self.nlayers+1):
                    # Skip everything but 2body in layers below surface
                    if ilayer > 1:
                        continue
                    # get indexes in F90 format
                    idxR = [str(rea.baseIdx), str(ilayer)]
                    for x in rea.reactants: 
                        if x.isGas:
                            idxR.append(str(x.idx + 1))
                        else:
                            idxR.append(str(int(x.idx + 1) + (ilayer-1)*offset))
                    idxP = []
                    for x in rea.products:
                        if x.isGas:
                            idxP.append(str(x.idx + 1))
                        else:
                            #if x.name in ["o_H2_0001", "p_H2_0001"]:
                            #    idxP.append(dummy)
                            #else:
                            idxP.append(str(int(x.idx + 1) + (ilayer-1)*offset))
                    # fill missing species with dummies
                    idxR += [dummy] * (maxR - len(idxR))
                    idxP += [dummy] * (maxP - len(idxP))

                    # ----------------------------------------------
                    # add mask reaction integer for the reaction
                    rtype = getReactionType(rea)
                    idxP += [str(rtype)]
                    # store the joined row
                    allIdx.append(" ".join(idxR + idxP))
            else:
                # get indexes in F90 format
                idxR = [str(rea.baseIdx), str(0)]
                for x in rea.reactants:
                    idxR.append(str(x.idx + 1))
                idxP = []
                for x in rea.products:
                    idxP.append(str(x.idx + 1))
                # fill missing species with dummies
                idxR += [dummy] * (maxR - len(idxR))
                idxP += [dummy] * (maxP - len(idxP))

                # ----------------------------------------------
                # add mask reaction integer for the reaction
                rtype = getReactionType(rea)
                idxP += [str(rtype)]
                # store the joined row
                allIdx.append(" ".join(idxR + idxP))

        # join all rows
        print(" Number of reactions in getReactionArray: ", len(allIdx))

        return allIdx

    def getReactionArrayIdx(self):
        print(" Running getReactionArray ")
        maxR = 4  # max number of reactants
        maxP = 4  # max number of products

        # prepare F90 statement
        rpIdx = "integer,parameter::maxRP2=" + str(maxR + maxP + 1) + "\n"
        rpIdx += "integer, dimension(nrea-1, maxRP2):: reactionArray"
        return rpIdx

    # ************************
    # get commons to be replaced in commons.f90
    def getIdxList(self):
        idxList = ""
        # loop on species to get indexes
        for species in self.species:
            if species.isGas or species.namebase in ['dummy', 'mask']:
                fidx = species.fidx
                idx = species.idx
                idxList += "integer,parameter::" + fidx + \
                    "=" + str(idx + 1) + "\n"
            else:
                for ilayer in range(1, self.nlayers+1):
                    offset = self.nDustSpecies*(ilayer-1)
                    if ilayer > 1:
                        offset += 3  # Dummy + mask
                    fidx = species.fidx.replace("_0001", "_%0.4i" % ilayer)
                    idx = species.idx + offset
                    idxList += "integer,parameter::" + fidx + \
                        "=" + str(idx + 1) + "\n"

        return idxList

    # ********************
    # get reaction fluxes to be replaced in flux.f90
    def getFluxes(self):
        fluxes = ""
        # loop on all rates
        for rea in self.reactionsAll:
            # Skip layers above one for freezeout as these will be duplicates and skip dummy
            if (rea.layer > 1 and rea.type == 'freezeout') or rea.type == 'dummy':
                continue
            fluxes += "flux(" + str(rea.idx + 1) + ") = " + rea.RHS + "\n"

        return fluxes

    # ********************
    # get reaction fluxes to be replaced in flux.f90, F90 format
    def getFluxesForm(self):
        fluxes = ""
        # loop on all rates
        for rea in self.reactions:
            # Skip layers above one for freezeout as these will be duplicates and skip dummy
            if (rea.layer > 1 and rea.type == 'freezeout') or rea.type == 'dummy':
                continue
            fluxes += "flux(" + str(rea.idx + 1) + ") = " + rea.RHS + "\n"

        return fluxes

    # ***************************
    # get Jacobian to be replaced in ode.f90
    def getJacobian(self):
        # Jacobian is a dictionary where df(i)/dx(j) = JPD[j][i]
        JPD = dict()

        # loop on reactions
        for rea in self.reactionsAll:
            # loop on reactants (you need derivative only for them)
            for species in rea.reactants:
                # create dictionary if species not in matrix
                if species.fidx not in JPD:
                    JPD[species.fidx] = dict()
                # loop on reactants
                for RR in rea.reactants:
                    # create list if matrix element not present
                    if RR.fidx not in JPD[species.fidx]:
                        JPD[species.fidx][RR.fidx] = []
                    # append jacobian expression
                    JPD[species.fidx][RR.fidx].append(
                        "-" + rea.getJacPD(species.fidx))
                # loop on products
                for PP in rea.products:
                    # create list if matrix element not present
                    if PP.fidx not in JPD[species.fidx]:
                        JPD[species.fidx][PP.fidx] = []
                    # append Jacobian expression
                    JPD[species.fidx][PP.fidx].append(
                        "+" + rea.getJacPD(species.fidx))

        jacobianFull = ""
        idxSpecies = [x.fidx for x in sorted(
            self.species, key=lambda xx: xx.idx)]
        # loop on indexes
        for jdx in JPD:
            # IF wrapper (see DLSODES documentation)
            jacobianFull += "if(j==" + jdx + ") then\n"

            # loop on derived RHS elements
            for idx in idxSpecies:
                if idx in JPD[jdx]:
                    jacobianFull += "pdj(" + idx + ") = " + \
                        (" &\n".join(JPD[jdx][idx])) + "\n"
                else:
                    jacobianFull += "pdj(" + idx + ") = 0d0\n"

            jacobianFull += "return\n"
            jacobianFull += "endif\n\n"

        return jacobianFull


    # ********************
    # get reaction rates to be replaced in rates.f90, F90 format
    def getRatesF90(self, gasOnly=False, dustOnly=False):
        if gasOnly and dustOnly:
            print("ERROR: getRatesF90: cannot be both gasOnly and dustOnly")
            sys.exit()

        rates = ""
        # loop on gas rates
        if gasOnly:
            reactions = self.reactionsGas
        # loop in dust rates
        elif dustOnly:
            reactions = self.reactions
            reactions.append(self.maskReaction[0])

        # loop on all rates(dust+gas)
        else:
            reactions = self.reactionsAll

        for rea in reactions:
            # reaction comment
            rates += "!" + rea.verbatim + " (" + rea.type + ")\n"
            # check if limits are present (NOTE: dust application so Tgas<800K is OK)
            hasLimits = (rea.Tmin > 1e0 or rea.Tmax <
                            8e2) and rea.hasMultipleTranges
            # if temperature limits are present wrap rate with IF
            if hasLimits:
                rates += "if(variable_Tgas>=" + strF90(rea.Tmin) \
                    + " .and. variable_Tgas<" + \
                    strF90(rea.Tmax) + ") then\n"
            rates += "kall(" + str(rea.idx + 1) + ") = " + \
                str(rea.krateF90) + "\n"
            # end IF wrapper
            if hasLimits:
                rates += "end if\n"
            rates += "\n"

        return rates

    # ********************
    # check barrierless reactions, reactants to check reactivity are listed in molscheck
    # (i.e. check only reactions with at least one of these reactants, if barrier is zero
    # and it's autogenerated, then exclude reaction, since it is unlikely to be barrierless)
    # returns reactions list without suspicious reactions
    def checkBarriers(self, reactions, molscheck):

        noBarriers = []
        okReactions = []
        # loop on reactants to check
        for rea in reactions:
            # check if at least one of the reactants is present
            hasReactant = any([rea.hasReactant(x) for x in molscheck])
            # check if reactant is present and barrier is zero
            if (not rea.barrierFromFile) and hasReactant and (len(rea.reactants) == 2):
                noBarriers.append(rea)
            else:
                okReactions.append(rea)

        # write message if barrierless are present
        if noBarriers:
            print("\nReactions unknown barrier and these reactants will be removed:")
            print((" " + (", ".join(molscheck))))
            print((" (i.e. from " + str(len(reactions)) + " to "
                   + str(len(okReactions)) + " reactions)"))

        # wait for press enter
        # raw_input("Press enter to continue:")

        return okReactions

    # ********************
    # check branching ratios for unspecified surface reactions
    def checkBranchingRatios(self, reactions):
        indices, count, reactants_list = [], [], []

        print("Estimating branching ratios for 2body reactions without information from the input file. May want to check these to avoid errors.")
        # loop on reactants to check
        for i, rea in enumerate(reactions):
            if not rea.type.startswith('2body'):
                continue
            if rea.bRatioFromFile:
                continue

            reactants = sorted([r.name for r in rea.reactants])
            reactants = "_".join(reactants)
            # Check if reactants already in list
            if reactants in reactants_list:
                idx = reactants_list.index(reactants)
                indices[idx].append(i)
                # 2body reactions that end with gas-phase products are not counted in the bratios, since the branching between chemical desorption calculated from dH
                if rea.type.endswith('_gas'):
                    continue
                else:
                    # H2/D2 ortho/para are done based on statistical OPR ratio. Therefore we ignore the ortho-species when 'counting' the Bratio in this function.
                    if any(pp.name.startswith('o_') for pp in rea.products):
                        count[idx] += 0
                    else:
                        count[idx] += 1
            else:
                reactants_list.append(reactants)
                indices.append([i])
                # 2body reactions that end with gas-phase products are not counted in the bratios, since the branching between chemical desorption calculated from dH
                if rea.type.endswith('_gas'):
                    count.append(0)
                else:
                    # H2/D2 ortho/para are done based on statistical OPR ratio. Therefore we ignore the ortho-species when 'counting' the Bratio in this function.
                    if any(pp.name.startswith('o_') for pp in rea.products):
                        count.append(0)
                    else:
                        count.append(1)

                if len(count) != len(indices) or len(indices) != len(reactants_list):
                    print("Error!!")
                    sys.exit()

        # now add Bratios to reactions.
        for i, idx in enumerate(indices):
            if count[i] > 1:
                ratio = 1.0 / count[i]
                for ii in idx:
                    reactions[ii].Bratio = ratio

        return reactions

    # *********************
    # save network to a file, see header below
    def saveNetworkToFile(self, fileName="network.ntw"):

        # open file to write
        fout = open(fileName, "w")
        fout.write("############\n")
        fout.write("#this network is autogenerated\n")
        # header, where lmax is the column size
        fout.write(tabrow(["#type", "Ea/K", "R", "R",
                           "P", "P"], lmax=18).strip() + "\n")
        # loop on reactions
        for rea in self.reactions:
            reactantsName = [x.name for x in rea.reactants]
            productsName = [x.name for x in rea.products]
            # add an empty species if less reactants
            if len(reactantsName) == 1:
                reactantsName += [""]
            # write table row to file
            fout.write(tabrow([rea.type, rea.Ea] + reactantsName
                              + productsName, lmax=18).strip() + "\n")

        fout.close()

        print(("Network saved to " + fileName))

    # ***********************************
    def getSwappingRates(self, ratio = 0.8):
        # We assume gas-phase comes first in the species list, need ofset for species
        offset = self.nGasSpecies
        speciesDict = {x.dictname: x for x in self.species}

        # constants and parameters
        kb = 1.38064852e-16  # Boltzmann constant, erg/K
        mp = 1.6726219e-24  # proton mass, g
        Ns = 1.5e15  # Number of sites / cm2
        species_specific = ['H', 'H2', 'o_H2', 'p_H2', 'C', 'N', 'O', 'CO', 'CO2', 'DCN', 'DNC', 'HCN', 'HNC', 'HCO', 'DCO']
        species_H2O = speciesDict['H2O']
        rates = ""
        for i, s in enumerate(self.species):
            if s.isGas:
                continue
            if s.Eice == None:
                continue
            # Frequency from harmonic oscillator approximation
            if s.namebase in species_specific:
                nu0 = sqrt(2e0 * Ns * s.Eice * kb / pi**2 / (s.mass*mp))
                #rate = strF90(nu0) + "*exp(-" + strF90(s.Eice*ratio) + "*invTd)"
                rate = strF90(nu0) + "*exp(-" + strF90(max([550.0, s.Eice*0.8])) + "*invTd)"
            else:
                nu0 = sqrt(2e0 * Ns * species_H2O.Eice * kb / pi**2 / (s.mass*mp))
                rate = strF90(nu0) + "*exp(-" + strF90(species_H2O.Eice*ratio) + "*invTd)"
                

            rates += "kswap(%i) = " % (s.idx + 1 - offset) + rate + "\n"

        return rates

    # ***********************************
    def getSticking(self):
        # We assume gas-phase comes first in the species list, need ofset for species
        offset = self.nGasSpecies
        # constants and parameters, from He+2016
        alpha = 0.5
        # average:
        beta = 0.11
        gamma = 0.042
        # specific:
        betas = {'H2':0.059, 'HD': 0.072, 'D2':0.072, 'N2':0.12, 'O2':0.17, 'CO':0.08,\
            'CH4':0.18, 'CO2':0.082}
        gammas = {'H2':0.051, 'HD':0.029, 'D2':0.029, 'N2':0.043, 'O2':0.042, 'CO':0.04,\
            'CH4':0.045, 'CO2':0.044}

        E_lc = {'H2': 315.0, 'HD': 650.0, 'D2': 650.0}

        lines = ""
        for i, s in enumerate(self.species):
            if s.isGas:
                continue
            if s.Eice == None:
                continue

            # Check if species has specific value:
            try:
                b = betas[s.namebase.strip('_0001')]
                g = gammas[s.namebase.strip('_0001')]
            except KeyError:
                b = beta
                g = gamma

            try:
                e = E_lc[s.namebase.strip('_0001')]
            except KeyError:
                e = s.Eice / 1.5

            stick = strF90(alpha) + " * (1d0 - tanh(" + strF90(b) + "*(Tdust - " + strF90(g) + "*" + strF90(e) + ")))" 
            lines += "kstick(%i) = " % (s.idx + 1 - offset) + stick + "\n"

        return lines


    # ***********************************
    def getEbind(self):
        # We assume gas-phase comes first in the species list, need ofset for species
        offset = self.nGasSpecies
        # following furuya+2015
        Edes_H2 = 23.0
        Edes_H2_water = 450.0
        #

        lines = ""
        for i, s in enumerate(self.species):
            if s.isGas:
                continue
            if s.Eice == None:
                continue

            stick = strF90(s.Eice) + " * (1d0 - H2_coverage) + H2_coverage * (" +  strF90(s.Eice) + " * %.1f / %.1f )" %(Edes_H2, Edes_H2_water)
            lines += "ebind(%i) = " % (s.idx + 1 - offset) + stick + "\n"

        return lines

    # ***********************************
    # function to add special encounter desorption for H2. Hincelin+2015. 
    def addEncounterDesorption(self):
        if self.H2spin:
            RRs = []
            PPs = []
            for s in self.species:
                if s.name == 'o_H2_0001':
                    RRs.append(s)
                    RRs.append(s)
                    PPs.append(s)
                if s.name == 'o_H2_gas':
                    PPs.append(s)
            if len(RRs) == 2 and len(PPs) == 2:
                rea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios, yieldPD=-43)
                self.reactions.append(rea)
            RRs = []
            PPs = []
            for s in self.species:
                if s.name == 'p_H2_0001':
                    RRs.append(s)
                    RRs.append(s)
                    PPs.append(s)
                if s.name == 'p_H2_gas':
                    PPs.append(s)
            if len(RRs) == 2 and len(PPs) == 2:
                rea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios, yieldPD=-43)
                self.reactions.append(rea)
        else:
            RRs = []
            PPs = []
            for s in self.species:
                if s.name == 'H2_0001':
                    RRs.append(s)
                    RRs.append(s)
                    PPs.append(s)
                if s.name == 'H2_gas':
                    PPs.append(s)
            if len(RRs) == 2 and len(PPs) == 2:
                rea = reaction(RRs, PPs, self.Ea, self.barrierWidths, self.Bratios, yieldPD=-43)
                self.reactions.append(rea)

# ***********************************
def findReactions_parser(args):
    res = findReactions_func(*args)
    return res

# ****************
# Loop over combinations to find reactions
def findReactions_func(comb1, combinations, Ea, barrierWidths, Bratios, include_H2):
    reactions = []

    # loop on combinations
    for comb2 in combinations:
        # different combination exploded is not a reaction, e.g. [C,C,H] and [C,H,H]
        if comb1["exploded"] != comb2["exploded"]:
            continue

        # check if all species in combinations are gas-phase
        allGas1 = all([x.isGas for x in comb1["mols"]])
        allGas2 = all([x.isGas for x in comb2["mols"]])
        # check if there is at least a gas-phase species
        anyGas1 = any([x.isGas for x in comb1["mols"]])
        anyGas2 = any([x.isGas for x in comb2["mols"]])

        # if all the species are gas-phase skip
        if allGas1 and allGas2:
            continue

        # avoid 2body_gas_gas, 2body_gas_dust, 2body_dust_gas for now, except for H2.
        # Garrod et al. chemical desorption only for single product reactions
        specials_gas = ['p_H2_gas', 'o_H2_gas']
        #specials_gas =['H2_gas', 'HD_gas', 'D2_gas']
        specials_dust = ['p_H2_0001', 'o_H2_0001']

        if anyGas2:
            if len(comb2["mols"]) == 1:
                pass
            else:
                if include_H2:
                    continue
                else:
                    continueMainLoop = False
                    # We allow H2 gas products (for 2body with products) only:
                    for i in comb2["mols"]:
                        if i.isGas and i.name not in specials_gas:
                            continueMainLoop = True
                            break
                        else:
                            pass
                    if continueMainLoop:
                        continue

        else:
            if len(comb2["mols"]) == 1:
                pass
            else:
                if include_H2:
                    pass
                else:
                    continueMainLoop = False
                    for i in comb2["mols"]:
                        # We don't want H2_ice formation, as we want to release it straight into the gas phase.
                        if i.name in specials_dust:
                            continueMainLoop = True
                            break
                        else:
                            pass
                    if continueMainLoop:
                        continue
                    

        # ignore 2 body with gas-phase species
        if anyGas1 and (len(comb1["mols"]) == 2):
            continue

        # avoid destruction and direct sublimation for now
        # i.e. H2CO -> H2_gas + CO or H2CO -> H2_gas + CO_gas raising errors
        if anyGas2 and (len(comb2["mols"]) >= 2):
            if include_H2:
                continue
            else:
                # Check if species in specials (H2)
                nn = sorted([x.dictname for x in comb2["mols"]])
                if any([i.startswith('p_H2_gas') for i in nn]):
                    pass
                elif any([i.startswith('o_H2_gas') for i in nn]):
                    pass
                else:
                    continue

        # avoid freezout into multiple species
        if anyGas1 and (len(comb2["mols"]) == 2):
            continue

        # avoid (photo)dissociation at the moment:
        if (len(comb1["mols"]) == 1) and (len(comb2["mols"]) >= 2):
            continue

        # avoid isomeric change reactions
        if len(comb1["mols"]) == 1 and len(comb2["mols"]) == 1 \
                and not allGas1 and not allGas2:
            continue
        if len(comb1["mols"]) == 1 and len(comb2["mols"]) == 1:
            if comb1["mols"][0].namebase != comb2["mols"][0].namebase:
                continue
            if comb1["mols"][0].name.startswith('p_') and comb2["mols"][0].name.startswith('o_'):
                continue
            if comb1["mols"][0].name.startswith('o_') and comb2["mols"][0].name.startswith('p_'):
                continue

        # get exploded species
        expl1 = sorted([x.exploded for x in comb1["mols"]])
        expl2 = sorted([x.exploded for x in comb2["mols"]])

        # get species names
        names1 = sorted([x.dictname for x in comb1["mols"]])
        names2 = sorted([x.dictname for x in comb2["mols"]])

        # skip same names since A+B->A+B (this includes A+B->A+B_gas)
        # only for 2body, since A->A_gas or A_gas->A is allowed
        if (expl1 == expl2) and (len(comb1["mols"]) == 2):
            continue

        # skip same names since A->A (this NOT includes A->A_gas and
        # A_gas->A, because use names instead of exploded)
        if names1 == names2:
            continue

        # create a reaction object
        try:
            rea = reaction(comb1["mols"], comb2["mols"],
                           Ea, barrierWidths, Bratios)
        except ValueError:
            # Skip
            print("Value error in reaction creation")
            for i in comb1["mols"]+comb2["mols"]:
                print(i.dictname)
            continue
        except TypeError:
            print((comb1["mols"][0].dictname, comb2["mols"][0].dictname))
            print((comb1["mols"][0].charge, comb2["mols"][0].charge))

        # skip endothermic (note that sign convention is flipped)
        if rea.dH < 0e0 and rea.barrierFromFile == False:
            continue
            #print "reaction with negative dH", rea.verbatim

        if not include_H2:
            # Ignore H2 evap/freezeout, as we treat this seperately.
            if rea.type == 'evaporation' and rea.products[0].namebase == 'H2':
               continue
            if rea.type == 'freezeout' and rea.reactants[0].namebase == 'H2':
               continue

        # append reaction object to reaction list
        reactions.append(rea)

    return reactions


def findCombinations_parser(args):
    combinations = findCombinations_func(*args)
    return combinations

# ****************
# Find combinations:
def findCombinations_func(sp1, species):
    combinations = []
    combinationNames = []
    dictname1 = sp1.dictname
    tempspecies = copy.deepcopy(species)
    # ignore chemisorbed hydrogen (added later)
    if sp1.dictname == "Hc":
        return

    # add single species to combinations
    combinations.append({"mols": [sp1], "exploded": sp1.exploded})
    # ignore charged species (dust-phase)
    if sp1.charge != 0:
        return

    # remove p_ and o_ from exploded
    ex1 = sp1.exploded
    while "p_" in ex1:
        ex1.remove("p_")
    while "o_" in ex1:
        ex1.remove("o_")
    # loop on species to get the partner to sp1
    for sp2 in tempspecies:
        dictname2 = sp2.dictname
        # ignore chemisorbed hydrogen (added later)
        if sp2.dictname == "Hc":
            continue
        # ignore charged species (dust-phase)
        if sp2.charge != 0:
            continue
        # put species in sorted list
        spsp = sorted([dictname1, dictname2])

        # if species combination is known then skip
        if spsp in combinationNames:
            continue
        # append combination name to combination list
        combinationNames.append(spsp)
        # remove p_ and o_ from exploded
        ex2 = sp2.exploded
        while "p_" in ex2:
            ex2.remove("p_")
        while "o_" in ex2:
            ex2.remove("o_")
        # append dictionary to combination list
        combinations.append({"mols": [sp1, sp2],
                             "exploded": sorted(sp1.exploded + sp2.exploded)})
    return combinations



def LoopCombinations(results, combinationNames, combinations, i):
    r = results[i]
    if r is None:
        return
    for rr in r:
        if len(rr["mols"]) > 1:
            spsp = sorted(
                [rr["mols"][0].dictname, rr["mols"][1].dictname])
            if spsp in combinationNames:
                continue
            combinationNames.append(spsp)
            combinations.append(rr)
        else:
            combinations.append(rr)
    return
