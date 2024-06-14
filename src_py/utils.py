import shutil
import re
import glob
import os
import sys
from math import log10
import numpy as np


# *****************
# return a list with the models
def getModelsList(modelFolder="models/"):
    return [x[0].replace(modelFolder, "") for x in os.walk(modelFolder) if x[0] != modelFolder]


# ***************
# copy model files to main folder
def makeModel(name):
    # model folder
    sourceFolder = "models/" + name

    # check if model folder exists
    if not os.path.exists(sourceFolder):
        print("ERROR: model " + name + " does not exist!")
        print(" models are:")
        print(getModelsList())
        sys.exit()

    # copy all files
    for filename in glob.glob(os.path.join(sourceFolder, "*.*")):
        print("copying " + filename + " from model " + name)
        shutil.copy(filename, ".")


# ***************
# convert the numbers in a F90 expression (arg) to Python readable numbers
# e.g. exp(-1d2*x)+3d-2 -> exp(-1e2*x)+3e-2
def F90toPy(arg):
    for match in re.compile('[0-9][d][-+0-9]').findall(arg):
        repmatch = match.replace("d", "e")
        arg = arg.replace(match, repmatch)
    return arg


# ******************
# return true if arg is a number
def isNumber(arg):
    try:
        float(arg)
        return True
    except ValueError:
        return False


# ******************
# pretty-print a list vlist as a table (column size in characters is lmax)
def tabrow(vlist, lmax=40):
    # if lmax is a list just fills the missing values with the last value of lmax
    # e.g. if vlist has 5 elements lmax becomes [20,10] -> [20,10,10,10,10]
    if isinstance(lmax, list):
        lamax = lmax + [lmax[-1]] * (len(vlist) - len(lmax))
    else:
        # if lmax is not a list just replicate lmax for len(vlist) times
        lamax = [lmax] * len(vlist)

    # convert to string
    slist = [str(x) for x in vlist]
    # fill with missing spaces
    jlist = [slist[i] + (" " * (lamax[i] - len(slist[i]))) for i in range(len(slist))]

    return "".join(jlist)


# *************************
# F90 style for floating values, numbers smaller than lowerLimit are flushed to zero
def strF90(arg, fmt='%17.8E', lowerLimit=1e-302):
    if abs(arg) < lowerLimit:
        return "0d0"
    erg = fmt % arg
    # if negative number, add parenthesis
    if arg < 0:
        erg =  "(" + erg + ")"
    return (erg.lower().replace("e", "d")).strip()


# *************************
# LaTeX style for floating values
def latexExp(arg, rnd=2):
    if (arg == "None") or (arg is None) or (np.isnan(arg)):
        return "-"
    if arg == 0e0:
        return "0.0"
    larg = int(log10(abs(arg)))
    return "$" + str(round(arg / 1e1 ** larg, rnd)) + "\\times10^" + str(larg) + "$"


# *************************
# LaTeX style for integer values
def latexInt(arg):
    if (arg == "None") or (arg is None) or (np.isnan(arg)):
        return "-"
    return "$" + str(int(arg)) + "$"


# ********************
# replace names from KIDA style to our format
def speciesToKIDA(speciesName):
    if speciesName == "e-":
        speciesName = "E"
    if speciesName.startswith("p"):
        speciesName = "p_" + speciesName[1:]
    if speciesName.startswith("o"):
        speciesName = "o_" + speciesName[1:]
    if speciesName.startswith("m"):
        speciesName = "m_" + speciesName[1:]
    if speciesName.startswith("l"):
        speciesName = "l_" + speciesName[2:]
    if speciesName.startswith("c"):
        speciesName = "c_" + speciesName[2:]

    return speciesName


# **********************
# return human-friendly formatted time
def fuzzyTime(time):
    if time < 1:
        return str(round(time, 2)) + " s"
    if time < 60:
        return str(round(time, 1)) + " s"
    if time < 3600:
        return str(round(time / 60., 1)) + " min"

    return str(round(time / 3600., 2)) + " hr"


# **********************
# generic preprocessor operation (handles pragmas)
def doPP(fname, pragmaDict):
    # temp file
    ftmp = "temp.f90"

    # dict with blocks, key=block name, value=True inside block
    blocks = dict()

    # pragama begin and end block
    pragmaBegin = "!!BEGIN_"
    pragmaEnd = "!!END_"

    # open to write
    fout = open(ftmp, "w")

    # open to read
    for row in open(fname, "r"):
        srow = row.strip()

        # loop on pragmas
        for (pragma, rep) in pragmaDict.items():
            # check for BEGIN block
            if srow == pragmaBegin + pragma:
                fout.write(row)
                blocks[pragma] = True
                comment = "! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
                comment += "! NOTE: This block is auto-generated\n"
                comment += "! WHEN: " + getCurrentTime() + "\n"
                comment += "! CHANGESET: " + getChangeset() + "\n"
                comment += "! BY: " + getUser() + "\n\n"
                endcomment = "\n! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                repx = comment + rep + endcomment
                if blocks[pragma]:
                    fout.write(indentF90(repx))
            # check for END block
            if srow == pragmaEnd + pragma: blocks[pragma] = False

        # if at least one block is true: skip content
        if any(blocks.values()):
            continue

        fout.write(row)
    fout.close()

    # check if all pragmas are properly closed
    if any(blocks.values()):
        print({k: v for (k, v) in blocks.items() if v})
        sys.exit("ERROR: pragma still open in " + fname)

    # store original as a BAK file
    fback = fname.replace(".f90", ".bak")
    shutil.copyfile(fname, fback)

    # copy temporary file to original
    shutil.copyfile(ftmp, fname)

    print("PRE-PROCESSOR: " + fname + " done (original stored in " + fback + ")!")


# **********************
# generic preprocessor operation (handles pragmas)
def doPP_datfile(fname, data):
    # temp file
    ftmp = "temp.dat"

    # open to write
    fout = open(ftmp, "w")

    # loop on data
    for entry in data:
        fout.write(entry + "\n")
    fout.close()

    # store original as a BAK file
    fback = fname.replace(".dat", ".bak")
    try:
        shutil.copyfile(fname, fback)
    except FileNotFoundError:
        pass

    # copy temporary file to original
    shutil.copyfile(ftmp, fname)

    print("PRE-PROCESSOR: " + fname + " done (original stored in " + fback + ")!")

# **********************
def doPP_Python(fname, pragmaDict):
    # temp file
    ftmp = "temp.py"

    # dict with blocks, key=block name, value=True inside block
    blocks = dict()

    # pragama begin and end block
    pragmaBegin = "##BEGIN_"
    pragmaEnd = "##END_"

    # open to write
    fout = open(ftmp, "w")

    # open to read
    for row in open(fname, "r"):
        srow = row.strip()

        # loop on pragmas
        for (pragma, rep) in pragmaDict.items():
            # check for BEGIN block
            if srow == pragmaBegin + pragma:
                fout.write(row)
                blocks[pragma] = True
                comment = "# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
                comment += "# NOTE: This block is auto-generated\n"
                comment += "# WHEN: " + getCurrentTime() + "\n"
                comment += "# CHANGESET: " + getChangeset() + "\n"
                comment += "# BY: " + getUser() + "\n\n"
                endcomment = "\n# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                repx = comment + rep + endcomment
                repx = repx.replace("!", "#")
                if blocks[pragma]:
                    fout.write(indentPy(repx, depth=2))
            # check for END block
            if srow == pragmaEnd + pragma:
                row = row.replace("!!", "##")
                blocks[pragma] = False

        # if at least one block is true: skip content
        if any(blocks.values()): continue

        fout.write(row)
    fout.close()

    # check if all pragmas are properly closed
    if any(blocks.values()):
        print({k: v for (k, v) in blocks.items() if v})
        sys.exit("ERROR: pragma still open in " + fname)

    # store original as a BAK file
    fback = fname.replace(".py", ".bak")
    shutil.copyfile(fname, fback)

    # copy temporary file to original
    shutil.copyfile(ftmp, fname)

    print("PRE-PROCESSOR: " + fname + " done (original stored in " + fback + ")!")



# **********************
# get current git changeset, if fails returns "xxxxxxx"
def getChangeset():
    import os
    # name of the git master file
    masterfile = "./.git/refs/heads/master"
    changeset = ("x" * 7)  # default unknown changeset
    # if git master file exists grep the changeset
    if os.path.isfile(masterfile):
        changeset = open(masterfile, "r").read()
    return changeset.strip()


# ***********************
# get current time with format fmt
def getCurrentTime(fmt="%Y-%m-%d %H:%M:%S"):
    import datetime
    return datetime.datetime.now().strftime(fmt)


# **********************
# get user@hostname, if fails returns unknown@unknown
def getUser():
    import getpass
    import socket
    try:
        return getpass.getuser().decode('utf-8') + "@" + socket.gethostname().decode('utf-8')
    except:
        return "unknown@unknown"


# *******************************
# indent F90 text and remove multiple blank lines, depth is the root indenting depth
def indentF90(text, depth=1):
    # opening and closing tokens
    tokenopen = ["do ", "function", "subroutine", "contains", "else", "else if", "elseif", "module",
                 "program", "type,", "interface", "module procedure"]
    tokenclose = ["end do", "end if", "end function", "end subroutine", "else if", "elseif", "else",
                  "enddo", "end module", "endif", "end type", "endtype", "contains", "endfunction",
                  "endsubroutine", "endmodule", "end program", "endprogram", "end interface",
                  "endinterface", "module procedure"]

    # sort tokens by size
    tokenopen = sorted(tokenopen, key=lambda x: len(x), reverse=True)
    tokenclose = sorted(tokenclose, key=lambda x: len(x), reverse=True)

    # size of offset
    offset = (" " * 4)

    # ampersand block flag, for F90 multiline
    inamper = False

    indentedText = ""
    # loop on text lines
    for row in text.split("\n"):
        srow = row.strip()
        # check closing tokens (skip comments)
        if not srow.startswith("!"):
            # check for closing tokens
            for token in tokenclose:
                if srow.startswith(token):
                    depth -= 1
                    break

        # indent text (rstrip to avoid empty lines)
        indentedText += ((offset * depth) + srow).rstrip() + "\n"

        # check opening tokens (skip comments)
        if not srow.startswith("!"):
            # check for opening tokens
            for token in tokenopen:
                if (srow.startswith(token)):
                    depth += 1
                    break
            if srow.startswith("if") and srow.endswith("then"):
                depth += 1

        # check for ampersand at line end (skip comments)
        if not srow.startswith("!"):
            # increase depth if ampersand found (only if not already in ampersand block)
            if srow.endswith("&") and not inamper:
                depth += 1
                inamper = True
            # decrease depth if no more ampersand found (only if in ampersand block)
            if not srow.endswith("&") and inamper:
                depth -= 1
                inamper = False

    # remove multiple blank lines
    while "\n\n\n" in indentedText:
        indentedText = indentedText.replace("\n\n\n", "\n\n")

    return indentedText


# *******************************
# indent Python text and remove F90 tokens, F90 multiline characters;
# depth is the root indenting depth
def indentPy(text, depth=1):
    # F90 tokens
    F90tokens = ["integer", "character"]

    # size of offset
    offset = (" " * 4)

    indentedText = ""
    print("Entered `indentPy`.")
    inamper = False
    # loop on text lines
    for row in text.split("\n"):
        srow = row.strip()
        # check closing tokens (skip comments)
        if not srow.startswith("#"):
            # check for multiline ampersand character
            if srow.endswith("&") and not inamper:
                srow = srow.strip("&")
                inamper = True
                depth += 1
            if not srow.endswith("&") and inamper:
                depth -= 1
                inamper = False
            # check for F90 tokens
            print(srow)
            for token in F90tokens:
                if srow.startswith(token):
                    idx = srow.find("::")
                    if idx != -1:
                        srow = 'self.' + srow[idx + 2:]

        # indent text (rstrip to avoid empty lines)
        indentedText += ((offset * depth) + srow).rstrip() + "\n"

    # remove multiple blank lines
    while "\n\n\n" in indentedText:
        indentedText = indentedText.replace("\n\n\n", "\n\n")

    return indentedText


def getReactionType(rea):
    '''
        Reaction type for reactionarray.
    '''
    if rea.type.startswith('freezeout'):
        rtype = 0
    elif rea.type.startswith('evaporation'):
        rtype = 1
    elif rea.type.startswith('CRdesorption'):
        rtype = 2
    elif rea.type.startswith('photodesorption'):
        rtype = 3
    elif rea.type.startswith('photodissociation'):
        rtype = 4
    elif rea.type.startswith('2body'):
        rtype = 5
    elif rea.type.startswith('"gasphase"'):
        rtype = 6
    else:
        rtype = -1
    return rtype
