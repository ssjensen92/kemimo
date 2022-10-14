# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import numpy as np
import re
from glob import glob

def readRow(row):
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
            if row.startswith('#'):
                continue
            print("ERROR: in KIDA network row element has spaces in the middle!")
            print(" Probably format problems:", dataRow[keys[i]])
            print(" Line here below triggered the ERROR")
            print(srow)
            #sys.exit()
        # increase format position
        position += fmt[i]
    return dataRow


def readRow2(row):
    # number of reactants and products expected in KIDA file
    maxReactants = 2
    maxProducts = 2

    # spacing format
    fmt = [5] + [7] * maxReactants + maxProducts * [7] + [1] + 3 * [12] + [10, 10] + [6] \
        + [5,4] + [7, 2, 2, 2]
    # keys names
    keys = ["n"] + ["R" + str(i) for i in range(maxReactants)] \
        + ["P" + str(i) for i in range(maxProducts)] + ['x'] \
        + ["a", "b", "c"] + ["F", "g"] + ["unc", "type"] \
        + ["tmin", "tmax"] + ["formula", "num", "subnum", "recom"]

    keys = ["n"] + ["R" + str(i) for i in range(maxReactants)] \
        + ["P" + str(i) for i in range(maxProducts)] \
        + ["a", "b", "c"] + ["F", "g"] + ["unc", "type"] \
        + ["tmin", "tmax"] + ["formula", "num", "subnum", "recom"]

    srow = row.split()

    dataRow = dict()
    for i, entry in enumerate(srow):
        dataRow[keys[i]] = entry.strip()
    return dataRow



def read_file(fname='gasNetwork.dat'):
    with open(fname, 'r') as f:
        lines = f.readlines()
    return lines

## read:        
old = read_file()
new = read_file()
sipila = read_file('2017_ground_state_to_species_arrhenius_1.dat')

rows_old = []
rows_sipila = []
for s in old:
    rows_old.append(readRow(s))
for s in sipila:
    rows_sipila.append(readRow2(s))

sipila_found = np.zeros_like(rows_sipila)

n_match = 0  
n_pmatch = 0  
match, no_match = [], []
imatch = []
for i, row in enumerate(rows_old):
    R0, R1, P0, P1 = row['R0'], row['R1'], row['P0'], row['P1']
    for j, row2 in enumerate(rows_sipila):
        R0_2, R1_2, P0_2, P1_2 = row2['R0'], row2['R1'], row2['P0'], row2['P1']
        if (R0 == R0_2) and (R1 == R1_2) and (P0 == P0_2) and (P1 == P1_2):
            if all([(abs(float(row[i])-float(row2[i])) < 1e-1) for i in ['a', 'b', 'c']]):
                n_pmatch += 1
            else:
                print(row)
                print(row2)
            n_match += 1
            match.append(row2)
            s = "{:>10} {:>10} {:>10} ".format("%.3e" % float(match[-1]['a']), "%0.3e" % float(match[-1]['b']), "%0.3e" % float(match[-1]['c']))
            sold = new[i]
            snew = sold[:90] + s + sold[123:153] + "{:<3}".format(int(match[-1]['tmin'])) + sold[156:]
            new[i] = snew
            sipila_found[j] = 1
        elif (R1 == R0_2) and (R0 == R1_2) and (P0 == P0_2) and (P1 == P1_2):
            if all([(abs(float(row[i])-float(row2[i])) < 1e-1) for i in ['a', 'b', 'c']]):
                n_pmatch += 1
            else:
                print(row)
                print(row2)
            n_match += 1
            match.append(row2)
            s = "{:>10} {:>10} {:>10} ".format("%.3e" % float(match[-1]['a']), "%0.3e" % float(match[-1]['b']), "%0.3e" % float(match[-1]['c']))
            sold = new[i]
            snew = sold[:90] + s + sold[123:153] + "{:<3d}".format(int(match[-1]['tmin'])) + sold[156:]
            new[i] = snew
            sipila_found[j] = 1
        elif (R0 == R0_2) and (R1 == R1_2) and (P1 == P0_2) and (P0 == P1_2):
            if all([(abs(float(row[i])-float(row2[i])) < 1e-1) for i in ['a', 'b', 'c']]):
                n_pmatch += 1
            else:
                print(row)
                print(row2)
            n_match += 1
            match.append(row2)
            s = "{:>10} {:>10} {:>10} ".format("%.3e" % float(match[-1]['a']), "%0.3e" % float(match[-1]['b']), "%0.3e" % float(match[-1]['c']))
            sold = new[i]
            snew = sold[:90] + s + sold[123:153] + "{:<3}".format(int(match[-1]['tmin'])) + sold[156:]
            new[i] = snew
            sipila_found[j] = 1
        elif (R1 == R0_2) and (R0 == R1_2) and (P1 == P0_2) and (P0 == P1_2):
            if all([(abs(float(row[i])-float(row2[i])) < 1e-1) for i in ['a', 'b', 'c']]):
                n_pmatch += 1
            else:
                print(row)
                print(row2)
            n_match += 1
            match.append(row2)
            s = "{:>10} {:>10} {:>10} ".format("%.3e" % float(match[-1]['a']), "%0.3e" % float(match[-1]['b']), "%0.3e" % float(match[-1]['c']))
            sold = new[i]
            snew = sold[:90] + s + sold[123:153] + "{:<3}".format(int(match[-1]['tmin'])) + sold[156:]
            new[i] = snew
            sipila_found[j] = 1
            

for row in rows_sipila:
    if row not in match:
        no_match.append(row)

new2 = []
n_identity = 0
for line in new:
    row = readRow(line)
    R0, R1, P0, P1 = row['R0'], row['R1'], row['P0'], row['P1']
    if sorted([R0, R1]) == sorted([P0, P1]):
        n_identity += 1
        continue # identity reactions
    else:
        new2.append(line)

new = new2
del new2
        

fmt = [11] * 3 + [1] + 5 * [11] + [1] + 3 * [11] + [8, 9] + [1, 4, 3] \
        + 2 * [7] + [3, 6, 3, 2]
t = ''
i = 0
left = [0,1,2,3,4,5,6,7,8]
right = [13, 14, 15, 16, 17, 18, 19, 20, 21, 23]
for iii in fmt:
    if i in left:
        t += "{:<%i}" % iii
    elif i in right:
        t += "{:>%i}" % iii
    else:
        t += "{:^%i}" % iii
    i += 1

new[-1] += '\n'
n_notfound = 0
n = int(new[-1].split()[-3])
for i in range(len(sipila_found)):
    if sipila_found[i] == 0:
        srow = rows_sipila[i]
        if sorted([srow['R0'], srow['R1']]) == sorted([srow['P0'], srow['P1']]):
            #print(srow)
            continue
        if float(srow['a']) == 0.0:
            continue
        n += 1
        #print(rows_sipila[i])
        n_notfound += 1
        line = t.format(srow['R0'], srow['R1'], '', '', srow['P0'], srow['P1'], '', '', '', '', "%0.3e" % float(srow['a']), "%0.3e" % float(srow['b']), "%0.3e" % float(srow['c']), srow['F'], srow['g'], '', srow['unc'], srow['type'], srow['tmin'], srow['tmax'], str(3), str(n), srow['subnum'], srow['recom'])
        new.append(line + '\n')


                 
print("n_match: ", n_match)
with open('gasNetwork_fixed.dat', 'w') as f:
    f.writelines(new)
        

