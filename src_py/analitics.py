from utils import tabrow, fuzzyTime

class analitics:

	#***************
	#class constructor requires reactions database
	def __init__(self, db):
		self.db = db

	#****************
	#load fluxes from file
	def loadFlux(self, fname):
		import time

		data = dict()
		print("Loading data from "+fname)
		with open(fname) as fh:
			fileLines = sum(1 for _ in fh)
		print("File lines:", fileLines)

		timeStart = time.time()

		#loop on rows
		icount = 0
		for row in open(fname, "rb"):
			srow = row.strip()
			if(srow==""): continue
			#read data
			(xvar, idx, flux) = [float(x) for x in srow.split(" ") if(x!="")]
			#index is integer
			idx = int(idx)
			flux = max(flux, 0e0)

			#get reaction from database, idx-1 because F90 counts form one
			rea = self.db.getReactionByIndex(idx-1)
			if(rea.type.startswith("gasphase")):
				continue

			#create dictionary element if not present
			if(not xvar in data):
				data[xvar] = {"idx":[], "flux":[], "reaction":[]}
			#store data in dictionary
			data[xvar]["idx"].append(idx)
			data[xvar]["flux"].append(flux)
			data[xvar]["reaction"].append(rea)
			icount += 1
			if(icount%10000==0):
				dtime = (time.time() - timeStart)*(fileLines-icount)/icount
				print("LOADING:", icount*1e2/fileLines,"%", fuzzyTime(dtime))

		print("done!")
		#copy data to attribute
		self.fluxData = data

	#*******************
	#recursively find the path to the final product
	def findPath(self, finalProduct, epsProb=1e-4, epsCount=1e-1, itermax=1000):
		from random import random as rand
		from random import shuffle
		import time, subprocess

		#species to find in the reactions products

		print("**********************")
		print("Finding reaction paths to "+finalProduct)
		print(" probability threshold:", epsProb)
		print(" count threshold:", epsCount)
		print(" max iterations:", itermax)

		allReactions = []
		maxAbsCount = dict()

		timeStart = time.time()

		icount = 0
		#loop on time dumps
		for xvar, data in self.fluxData.items():
			icount += 1
			if(icount%10==0):
				dtime = (time.time() - timeStart)*(len(self.fluxData)-icount)/icount
				print("MC:", round(icount*1e2/len(self.fluxData), 2), "%", fuzzyTime(dtime))

			#get indexes from data
			idxs = list(range(len(data["reaction"])))
			#store indexes in order of flux
			idxs = sorted(idxs, key=lambda x:data["flux"][x], reverse=True)


			#precompute some quantities to speed-up
			#key=species.name, value=[flux/totFlux]
			fluxDict = dict()
			#key=species.name, value=[reactions]
			reactionsDict = dict()
			#key=species.name, value=[sumFlux/totFlux]
			sumFluxDict = dict()
			#loop on species names to prepare dictionaries
			for speciesName in [x.name for x in self.db.species]:
				totFlux = 0e0
				#compute total flux for the reactions involving speciesName
				for ii in idxs:
					rea = data["reaction"][ii]
					#if reaction has species add to total flux
					if(rea.hasSpecies(speciesName)):
						totFlux += data["flux"][ii]

				#skip species with no flux
				if(totFlux==0e0):
					continue

				#init lists in dictionaries for the given species
				sumFluxDict[speciesName] = []
				reactionsDict[speciesName] = []
				fluxDict[speciesName] = []
				sumFlux = 0e0
				#loop on reactions to compute sum progression of fluxes
				for ii in idxs:
					rea = data["reaction"][ii]
					if(rea.hasSpecies(speciesName)):
						sumFlux += data["flux"][ii]
						#store sum progression as list
						sumFluxDict[speciesName].append(sumFlux/totFlux)
						#store flux of corresponding progression
						fluxDict[speciesName].append(data["flux"][ii]/totFlux)
						#store reaction object as corresponding progression
						reactionsDict[speciesName].append(rea)

			toFind = []
			foundReactions = []
			totalCount = dict()

			pathOK = pathKO = 0
			#iterates to do MC
			for i in range(itermax):
				#initialize path if no species to find
				if(toFind==[]):
					istart = i
					prob = 1e0 #total path probability
					path = [] #list of reactions in path
					mcCount = dict() #number of time iteration pass by reaction, key=reaction.hash
					toFind = [finalProduct] #initially the species to find is the final product

				shuffle(toFind)
				speciesName = toFind[0]

				#throw a dice to determine path
				randFlux = rand()
				#loop on partial sums to determine reaction chosen by dice
				for ii in range(len(reactionsDict[speciesName])):
					#when reaction found
					if(randFlux<=sumFluxDict[speciesName][ii]):
						#compute total path probability
						prob *= fluxDict[speciesName][ii]
						#get reaction
						rea = reactionsDict[speciesName][ii]
						break

				#if probability drops below limit, ignore path and unset toFind to restart
				if(prob<epsProb):
					toFind = []
					pathKO += 1
					continue

				#append reaction to path if not already found
				if(not rea in path):
					path.append(rea)

				#initialize iteration count
				if(not rea.hash in mcCount):
					mcCount[rea.hash] = 0
				#add iteration count to the chosen reaction
				mcCount[rea.hash] += 1

				#add products or reactants to path (gas species are excluded)
				if(rea.hasReactant(speciesName)):
					toFind = [x.name for x in rea.products if not x.isGas]
				else:
					toFind = [x.name for x in rea.reactants if not x.isGas]

				#if path reached a dead end and probability is still high
				# stores path into full reaction list
				if(toFind==[]):
					pathOK += 1
					#loop on path reactions
					for rea in path:
						#init total count if not present
						if(not rea.hash in totalCount):
							totalCount[rea.hash] = 0
						#sum path counts into total count
						totalCount[rea.hash] += mcCount[rea.hash]
						#if reaction is not known append to all reactions
						if(not rea in foundReactions):
							foundReactions.append(rea)

			if(float(pathOK)/(pathKO+pathOK)<.5):
				print("WARNING: more than 50% of paths ended because of low probability!")

			maxCount = max(totalCount.values())
			#check if counts are enough to store reaction
			for rea in sorted(foundReactions, key=lambda x:totalCount[x.hash], reverse=True):
				#if counts are not many just break loop (that is ascending count sorted)
				if(totalCount[rea.hash]<maxCount*epsCount): break
				#add reaction to full list if not already there
				if(not rea in allReactions):
					allReactions.append(rea)
				#init a total counter to be displayed later
				if(not rea.hash in maxAbsCount):
					maxAbsCount[rea.hash] = 0e0
				#store maxium count/max_count to be displayed later
				maxAbsCount[rea.hash] = max(maxAbsCount[rea.hash], float(totalCount[rea.hash])/maxCount)
				#print tabrow([rea.verbatim, mcCount[rea.hash]])

		#search recursively for sinks and sources and remove reactions
		while(True):
			#get list of sink/source names
			skrList = self.checkSinksSources(allReactions, onlySurface=True)
			#empty list is OK
			if(skrList==[]):
				break
			print("Removing sink/sources:",skrList)
			reaUpdate = []
			#loop on reactions
			for rea in allReactions:
				#if any source/sink found remove reaction
				if(any([rea.hasSpecies(x) for x in skrList])):
					continue
				#add reactions to list
				reaUpdate.append(rea)
			#update list of reactions
			allReactions = reaUpdate[:]


		#write reactions found
		print("")
		print("Reduced network found is:")
		for rea in sorted(allReactions, key=lambda x:maxAbsCount[x.hash], reverse=True):
			#write in a table
			print(tabrow(["ktmp("+str(rea.idx+1)+")", \
				"= kall("+str(rea.idx+1)+")", \
				"!"+rea.verbatim, \
				maxAbsCount[rea.hash]], \
				lmax=[12,15,30]).strip())

		print("reactions found:", len(allReactions))

		try:
			subprocess.Popen(['notify-send', "KEMIKO: done with analysis!"])
		except:
			print("WARINING: notify-send not working!")


	#**********************
	#given a list of reactions looks for sink and source species
	def checkSinksSources(self, reactionList, onlySurface=False):

		allReactants = []
		allProducts = []
		#store reactant and product names
		for rea in reactionList:
			allReactants += [x.name for x in rea.reactants if(onlySurface and not x.isGas)]
			allProducts += [x.name for x in rea.products if(onlySurface and not x.isGas)]

		#unique list of reactants and products
		allReactants = list(set(allReactants))
		allProducts = list(set(allProducts))


		skrList = []
		#loop on reactants for sources
		for species in allReactants:
			if(not species in allProducts):
				skrList.append(species)

		#loop on products for sinks
		for species in allProducts:
			if(not species in allReactants):
				skrList.append(species)

		return skrList

	#********************
	#find reactions that affects the ODE of a given species, up to eps fraction of the species ODE
	def findChannels(self, speciesName, eps = 1e-3):

		relevantReactions = []
		#loop on stored flux data
		for xvar, data in self.fluxData.items():
			#loop on reactions to determine flux if species is present
			for ii in range(len(data["reaction"])):
				if(data["reaction"][ii].hasSpecies(speciesName)):
					data["reaction"][ii].flux = data["flux"][ii]
				else:
					data["reaction"][ii].flux = 0e0

			totFlux = 1e-40
			#loop on reactions to compute the effect on the flux
			# starting from the most fluxy and decreasing
			for rea in sorted(data["reaction"], key=lambda x:x.flux, reverse=True):
				totFluxOld = totFlux
				#sum or subtract of species is reactant or product
				if(rea.hasProduct(speciesName)):
					totFlux += rea.flux
				else:
					totFlux -= rea.flux
				#if relative variation is smaller than eps than quit
				if(abs(totFluxOld-totFlux)/totFluxOld<eps):
					break
				#store reactions skipping known
				if(not rea in relevantReactions):
					relevantReactions.append(rea)

		#print results
		print("Reactions including "+speciesName+" with relative threshold ", eps)
		for rea in relevantReactions:
			print(rea.idx+1, rea.verbatim)

