import sys
import os

# move to the python source folder
sys.path.insert(0, "./src_py/")
from utils import makeModel, getModelsList

# check usage
if len(sys.argv) != 2:
    print("ERROR: usage is")
    print((" python " + sys.argv[0] + " model"))
    print("where models are:")
    print((getModelsList()))
    sys.exit()

# copy model files to root folder
makeModel(sys.argv[1])

# run the model
os.system("python pathway.py")
