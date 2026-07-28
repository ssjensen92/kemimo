# This is a template model

import sys

#move to the python source folder
sys.path.insert(0, "./src_py/")
from database import database

#create a database from data files
db = database()

#save species data to a file
#db.showSpecies(fileName="species.out")
#db.saveSpeciesLatexTable()

#save reactions data to files
#db.showReactions(fileName="reactions.out")
#db.showReactionsGas(fileName="reactionsGas.out")

#preprocess F90 files, compile, and run
db.run(run=True, library=False)
