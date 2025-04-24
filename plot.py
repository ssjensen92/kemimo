# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import numpy as np
import re

try:
    plt.style.use('mydefault')
except IOError:
    print("Using default mpl style")

def all_ice(dictionary, species_in):
    species = species_in + "_"
    first = True
    for key in dictionary.keys():
        if (key.startswith(species_in)) and (species in key) and (key != (species + "gas")):
            if first:
                ice = np.array(species_dict[key])
                first = False
            else:
                ice = np.add(ice, species_dict[key])
        else:
            continue

    return ice


def read_commons(fname="kemimo_commons.f90"):
    species = []
    sortkey = []
    idxListBegin = False
    idxListEnd = False
    for row in open(fname, "r"):
        srow = row.strip()
        if "BEGIN_IDXLIST" in srow:
            idxListBegin = True
        if "END_IDXLIST" in srow:
            idxListEnd = True
        if "integer,parameter::idx_" in srow and idxListBegin and not idxListEnd:
            match = re.compile('idx_(.*)=').findall(srow)
            match2 = re.compile('=\d+').findall(srow)
            species.append(match[0])
            sortkey.append(int(match2[0][1:]))
    species = np.array(species)
    return species[np.argsort(sortkey)]





if __name__ == "__main__":
    ## ----------------------------------------------------------
    ## LOAD DATA:
    fname = 'output.dat'
    try:
        os.path.isfile(fname)
    except:
        "%s not found" % fname
        pass

    try:
        data = np.genfromtxt(fname)
    except:
        print("Error reading output in file: ", fname)
        sys.exit()

    time = data[:,0]
    Tgas = data[:,1]
    ntot = data[:,2]
    Av = data[:,3]

    # Get species from commons:
    species = read_commons()
    # Establish surface species
    ice_species = []
    for s in species:
        if not s.endswith('_gas') and not s == 'mask' and not s == 'dummy':
            ice_species.append(s)

    # Create python dictionary for each species. 
    # e.g. species_dict['CO_gas'] returns CO gas abundances.
    species_dict = {species[i] : data[:,4+i] for i in range(len(species))}

    # Some sorting for testing etc.
    sortkey = np.argsort(data[-1, 4:])
    species = np.array(species)
    sorted_species = species[sortkey]
    sorted_abundances = data[:, 4:][:, sortkey]

    # Get the mask
    mask = species_dict['surface_mask']
    ## ----------------------------------------------------------
    # Get total surface abundances, using all_ice function.
    H2O = all_ice(species_dict, 'H2O')
    HDO = all_ice(species_dict, 'HDO')
    D2O = all_ice(species_dict, 'D2O')
    CO = all_ice(species_dict, 'CO')
    CO2 = all_ice(species_dict, 'CO2')
    CH3OH = all_ice(species_dict, 'CH3OH')
    H2CO = all_ice(species_dict, 'H2CO')
    HCN = all_ice(species_dict, 'HCN')
    HNC = all_ice(species_dict, 'HNC')

    HDCO = all_ice(species_dict, 'HDCO')
    CH3OD = all_ice(species_dict, 'CH3OD')
    CH2DOH = all_ice(species_dict, 'CH2DOH')

    ## ----------------------------------------------------------
    ## PLOT DATA:
    fig1 = plt.figure(figsize=(12, 7))
    ax1 = fig1.add_subplot(111)
    ax1.plot(time, species_dict['H_gas'], c='k', ls=':', label=r'H gas')
    ax1.plot(time, H2O, c='r', ls='-', label=r'H$_2$O')
    ax1.plot(time, HDO, c='r', ls=':', label=r'HDO')
    ax1.plot(time, D2O, c='r', ls='-.', label=r'D$_2$O')
    ax1.plot(time, CO, c='c', ls= '-', label = r'CO')
    ax1.plot(time, species_dict['CO_gas'], c='c', ls='-.', label=r'CO gas')
    ax1.plot(time, H2CO, c='b', ls='-', label=r'H$_2$CO')
    ax1.plot(time, HDCO, c='b', ls=':', label=r'HDCO')
    ax1.plot(time, CO2, c='g', ls= '-', label = r'CO$_2$')
    ax1.plot(time, HCN, c='y', ls= ':', label = r'HCN')
    ax1.plot(time, HNC, c='y', ls= '-.', label = r'HNC')
    ax1.plot(time, CH3OH, c='m', ls='-', label=r'CH$_3$OH')
    ax1.plot(time, CH3OD, c='m', ls='-.', label=r'CH$_3$OD')
    ax1.plot(time, CH2DOH, c='m', ls=':', label=r'CH$_2$DOH')

    plt.legend(fontsize=10)
    ax1.set_ylim([0.9e-10, 5e-4])
    ax1.set_xlim([1000, 2e7])
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Abundance')
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    plt.tight_layout()
    plt.show()

