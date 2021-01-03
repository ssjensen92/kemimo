# This is a template model

import sys, os
import datetime as dt

#move to the python source folder
sys.path.insert(0, "./src_py/")
from database import database

#create a database from data files
db = database(datadir="./data_deuterated_total/", nlayers=2, layerThickness=4.0)

#save species data to a file
#db.showSpecies(fileName="species.out")
#db.saveSpeciesLatexTable()

#save reactions data to files
#db.showReactions(fileName="reactions.out")
#db.showReactionsGas(fileName="reactionsGas.out")

#preprocess F90 files, compile, and run
db.run(run=True, library=False)
