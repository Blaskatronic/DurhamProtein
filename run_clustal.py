import numpy as np
import pylab as P
import os
import cPickle as pickle
import operator
import csv
import matplotlib.gridspec as gridspec
import json
from matplotlib.ticker import MultipleLocator
from PIL import Image, ImageTk
#from Tkinter import Tk, Label, Text, LEFT, RIGHT, BOTH, RAISED, W, N, E, S, Toplevel
#from ttk import Frame, Style, Button

wtvalue = 1.138069

class clustalw:
	def __init__(self, directory):
		self.directory = directory	#	Directory containing sequences
		self.listSeq = []		#	List of sequences
		self.listAln = []		#	List of alignments
		self.orderedList = []		#	Output list of aligned sequences
		self.overwriteFlag = 0		#	Selects whether to re-calculate alignments

	def runSetup(self):
		self.importSequences()
		self.runClustal()
		for element in self.listAln:	
			self.orderedList.append(self.reordering(element))
		return self.orderedList, len(self.listSeq)

	def importSequences(self):
		"""Finds txt files in the input directory that contain the sequences"""
		listFiles = os.listdir(self.directory)
		for n in range(0,len(listFiles)):
			if listFiles[n][-3:] == 'txt':
				self.listSeq.append(self.directory+'/'+listFiles[n])
		print 'Number of sequences found:',len(self.listSeq),'\n'
	#		if len(list_seq) == 0:
	#			directory = str(raw_input('Error: No sequences found. Please enter a new directory: '))
	#		else:
	#			list_seq.insert(0,list_seq[-1])
	#			del list_seq[-1]
	#			break

	def checkAlign(self, name):
		"""Checks to see if the corresponding binary sequence alignment for each txt file already exists"""
		nameAln = name[:-3]+'aln'
		check = os.system('ls '+nameAln+' >/dev/null 2>/dev/null')
		return check, nameAln

	def runClustal(self):
		"""Runs Clustalw if the sequence alignment doesn't already exist"""
		for element in self.listSeq:
			print 'Checking for binary alignment of '+str(element)+'...'
			check, nameAln = self.checkAlign(element)
			if check == 0 and self.overwriteFlag == 0:
				print 'Alignment file already present!\n'
			else:
				print 'Alignment file not present. Executing CLUSTAL...'
				os.system('clustalw ' + str(element))
				check, nameAln = self.checkAlign(element)
			self.listAln.append(nameAln)

	def reordering(self, aln):
		"""Converts the aln files into something that Python can interpret"""
		sequences = []
		alnFileHandle = open(aln, 'r')
		for line in alnFileHandle:
			if '|' in line:
				breakdown = line.split(' ')
				sequences.append(breakdown[0])
				sequences.append(breakdown[-1][:-1])
		alnFileHandle.close()
		sequencesMistress = []
		sequencesMaster = []
		notin = 0
		for element in range(0, len(sequences), 2):
			if len(sequencesMaster) == 0:
				sequencesMaster.append([sequences[element]])
			else:
				for n in range(len(sequencesMaster)):
					if sequences[element] not in sequencesMaster[n]:
	#					print sequencesMaster[n]
	#					print sequences[element], 'not in sequencesMaster!'
						notin = 1
					else:
						notin = 0
						break
				if notin == 1:
					sequencesMaster.append([sequences[element]])				
		for n in range(len(sequencesMaster)):
			sequencesMistress.append([])
		for element in range(len(sequences)):
			for n in range(len(sequencesMaster)):		
				if sequences[element] == sequencesMaster[n][0]:
					sequencesMistress[n].append(sequences[element+1])
		for element in range(len(sequencesMistress)):
			splitseq = list(''.join(sequencesMistress[element]))
			for acid in splitseq:
	#			if acid != '-':					# REMOVE GAPS - MAYBE DO GAP ANALYSIS FIRST THEN REMOVE THEM FOR RESIDUE NUMBERS/STRUCTURE ANALYSIS?
				sequencesMaster[element].append(acid)
		return sequencesMaster
	#	for element in range(0, len(sequences), 2):


class sequences:	# Seperate class instance for each binary alignment
	def __init__(self, alignment, groupUserDefined=None):
		self.groupAlphaHelix = [[8, 16], [99, 109], [110, 137], [138, 150], [151, 153], [168, 177], [179, 193]]
		self.groupBetaSheet = [[19, 23], [92, 98], [38, 44], [69, 70], [28, 30], [85, 88], [47, 52], [58, 64], [164, 166], [202, 205], [196, 199]]
		self.groupLigand = [[71], [72], [82], [83], [123], [127]]	# Also 128 from monomer B (omitted)
		self.groupUserDefined = groupUserDefined
		self.alignment = alignment
		self.matrix = np.asarray(self.alignment)
		self.residueList = []	# List of residues for each sequence (no names or gaps)
		self.params = []	# Parameters of the alignment: [Name, W/ Gaps Ident, W/O N-Terminus, W/O C-Terminus, Alpha, Beta, Ligand, Custom]
		self.nTerminusExtension = 0
		self.cTerminusExtension = 0
		self.warnedOnce = 0
		self.yvals = []		# Frequency list of residue number mutations
		
	def runAnalysis(self):
		#### Identify any C/N terminus and then display alternative identicalities.
		if self.matrix[0,1] == '-':	# If gaps present before E-coli M, then sequence has extended N terminus
			self.nTerminusExtension = 1
			print 'Compared sequence has an N-Terminus Extension...'
			print '-=Taking into account N-Terminus Extension: =-'
		if self.matrix[0,-1] == '-':
			self.cTerminusExtension = 1
			print 'Compared sequence has a C-Terminus Extension...'
			print '-=Taking into account C-Terminus Extension: =-'
		identityPercentage, identity = self.identityMatrix(self.matrix)
		print 'Sequence', self.matrix[1,0], 'is %.2f%% identical to sequence' % identityPercentage , str(self.matrix[0,0])+'.\n'

		location = findIndex(self.matrix[1,0], '|')
		proteinName = self.matrix[1,0][location[2]+1:]
		if proteinName[-1] == '|':
			proteinName = proteinName[:-1]
		if proteinName[-1] == '.':
			proteinName += '1'
		self.params.append(proteinName)
		self.params.append(identityPercentage)
		
		truncatedMatrix = np.asarray(self.matrix)

		if self.nTerminusExtension == 1 and self.cTerminusExtension == 1:
			print '-= Removing both N- and C-Terminus Extensions: =-'
			while True:
				if truncatedMatrix[0][1] == '-':
					truncatedMatrix = np.delete(truncatedMatrix, 1, 1)
				else:
					break
			while True:
				if truncatedMatrix[0][-1] == '-':
					truncatedMatrix = np.delete(truncatedMatrix, -1, 1)
				else:
					break
			identityPercentage, identity = self.identityMatrix(truncatedMatrix)
			print 'Without Terminus Extensions, sequence', self.matrix[1,0], 'is %.2f%% identical to sequence' % identityPercentage , str(self.matrix[0,0])+'.\n'
			self.params.append(identityPercentage)

		elif self.nTerminusExtension == 1:
			print '-= Removing N-Terminus Extension: =-'
			while True:
				if truncatedMatrix[0][1] == '-':
					truncatedMatrix = np.delete(truncatedMatrix, 1, 1)
				else:
					break
			identityPercentage, identity = self.identityMatrix(truncatedMatrix)
			print 'Without N-Terminus Extension, sequence', self.matrix[1,0], 'is %.2f%% identical to sequence' % identityPercentage , str(self.matrix[0,0])+'.\n'
			self.params.append(identityPercentage)
		
		elif self.cTerminusExtension == 1:
			print '-= Removing C-Terminus Extension: =-'
			while True:
				if truncatedMatrix[0][-1] == '-':
					truncatedMatrix = np.delete(truncatedMatrix, -1, 1)
				else:
					break
			identityPercentage, identity = self.identityMatrix(truncatedMatrix)
			print 'Without C-Terminus Extension, sequence', self.matrix[1,0], 'is %.2f%% identical to sequence' % identityPercentage , str(self.matrix[0,0])+'.\n'
			self.params.append(identityPercentage)
		
		else:
			self.params.append(None)

#		print truncatedMatrix
		#### This truncated matrix can include gaps, so for now just going to delete the gaps and the residues that correspond to the gaps in the second sequence - not sure if this is allowed! ####

#		print identity
		colsToBeDeleted = []
		for residue in range(len(truncatedMatrix[0])):
			if truncatedMatrix[0][residue] == '-':
				colsToBeDeleted.append(residue)
		colsToBeDeleted.reverse()
#		print identity
		for column in colsToBeDeleted:
			truncatedMatrix = np.delete(truncatedMatrix, column, 1)
			identity.pop(column)
#		print identity
#		print len(identity)
#		print truncatedMatrix
		# Obtain the gapless sequence alignment for residue numbers
		for sequence in self.alignment:
			self.residueList.append([])
			for residue in sequence[1:]:
				if residue != '-':
					self.residueList[-1].append(residue)
		self.groupAnalysis()
		paramsListElement, loopIndexElement, alpha, beta, custom = self.output()
		return paramsListElement, identity, loopIndexElement, alpha, beta, custom

	def identityMatrix(self, matrix):
#		print self.matrix
		### Calculate identity ###
		identity = []
		height, width = matrix.shape
		for i in range(1, width):			#### This will only work for binary alignments
#			print "Comparing", matrix[0,i], "with", matrix[1,i]
			if matrix[0,i] == matrix[1,i]:
				identity.append(1)
			else:
				identity.append(0)
		identityPercentage = float(100*sum(identity))/float(len(identity))		# Don't include the name in the identity calculation!
		test = 100
#		print "\n Identity Percentage =", identityPercentage, "Identity =", identity
#		print identity
		return identityPercentage, identity


	def identityList(self, structureList):			## Only works for binary alignments
#		print structureList
		breakOut = 0
		identityMaster = []
		comparison = []
		for group in range(len(structureList[0])):
			comparison.append([])
			for sequence in structureList:
				comparison[-1].append(sequence[group])
#		print comparison					# comparison = [[Seq1Group1, Seq2Group1], [Seq1Group2, Seq2Group2]...]
		for group in comparison:
#			print "group +1"
			groupidentity = 0
			for element in range(len(group[0])):
#				print "element +1"
				if len(group[0]) != len(group[1]):
					print group
					print "***WARNING*** RESIDUES NOT PRESENT FOR GROUP (IS SEQUENCE TOO SHORT?)"
#					print "SELF.WARNEDONCE =", self.warnedOnce
					if self.warnedOnce == 0:
						raw_input("Press return to continue...")
						self.warnedOnce = 1
					breakOut = 1
					break
				if group[0][element] == group[1][element]:
					groupidentity += 1
			if breakOut == 1:
				break
			identityMaster.append(float(100*groupidentity)/float(len(group[0])))
#		print identityMaster		# identityMaster = [Conservation of Alpha[20,30], Conservation of Alpha[34, 37]]
		return identityMaster, np.average(identityMaster)				

	def groupAnalysis(self):
		# Structure parts listed as: [sequence1[group1[residues], group2[residues]...], sequence2,...]
		alphaBit = []
		betaBit = []
		ligandBit = []
		for sequence in self.residueList:
			alphaBit.append([])
			betaBit.append([])
			ligandBit.append([])
			for section in self.groupAlphaHelix:
				alphaBit[-1].append(sequence[section[0]-1:section[1]])			## The -1 is to keep the residues inclusive (e.g. 20-30 inclusive for alpha)
			for section in self.groupBetaSheet:						## whilst still taking into account the python list[0].
				betaBit[-1].append(sequence[section[0]-1:section[1]])
			for section in self.groupLigand:						# This is individual residues hence difference!
				ligandBit[-1].append([sequence[section[0]-1]])
		if self.groupUserDefined != None:
			self.groupUserDefined = json.loads(self.groupUserDefined)
			customBit = []
#			print self.groupAlphaHelix
#			print self.groupUserDefined
			for sequence in self.residueList:
				customBit.append([])
				for section in self.groupUserDefined:
#					print "Custom:", section
					customBit[-1].append(sequence[section[0]-1:section[1]])
			customidentity, customConservation = self.identityList(customBit)
#		print self.residueList
#		print self.groupAlphaHelix
#		print alphaBit
		alphaidentity, alphaConservation = self.identityList(alphaBit)
		betaidentity, betaConservation = self.identityList(betaBit)
#		print "\n\n\n"
		ligandidentity, ligandConservation = self.identityList(ligandBit)
		print "\nBased on the inputted groups of:"
		print "Alpha helices at residues", self.groupAlphaHelix, "inclusive,"
		print "Beta sheets at residues", self.groupBetaSheet, "inclusive,"
		if self.groupUserDefined != None:
			print "Ligand binding sites at residues", self.groupLigand, "inclusive,"
			print "And user defined residue group of", self.groupUserDefined, "inclusive..."
		else:
			print "And ligand binding sites at residues", self.groupLigand, "inclusive..."
		print "\nThe alpha helix identity percentages are:", alphaidentity, "- a conservation factor of %.2f%%," % (alphaConservation)
		print "The beta sheet identity percentages are:", betaidentity, "- a conservation factor of %.2f%%," % (betaConservation)
		self.params.append(alphaConservation)
		self.params.append(betaConservation)
		if self.groupUserDefined != None:
			print "The Ligand binding site identity percentages are:", ligandidentity, "- a conservation factor of %.2f%%," % (ligandConservation)
			print "And the user defined group identity percentages are:", customidentity, "- a conservation factor of %.2f%%.\n" % (customConservation)
			self.params.append(ligandConservation)
			self.params.append(customConservation)
		else:
			print "And the ligand binding site identity percentages are:", ligandidentity, "- a conservation factor of %.2f%%.\n" % (ligandConservation)
			self.params.append(ligandConservation)
			self.params.append(None)
		self.identifyLoops()

	def identifyLoops(self):
		loops = []
		residueList = self.residueList[:]
		for sequence in residueList:
			loops.append([])
			for section in self.groupAlphaHelix:
				for i in range(section[0]-1, section[1]):
#					print 'Removing element', i
					sequence[i] = None
			for section in self.groupBetaSheet:
				for i in range(section[0]-1, section[1]):
					sequence[i] = None
			loops[-1].append([])
#			print 'Pre-residue loop', loops
			for residue in sequence:
#				print 'Residue =', residue
				if residue != None:
#					print 'Residue != None, noneFound = 0'
					noneFound = 0
					loops[-1][-1].append(residue)
#					print loops
				elif noneFound == 0:
#					print 'Residue == None, noneFound now = 1'
					noneFound = 1
					loops[-1].append([])
#					print loops
				else:			# elif residue == None and noneFound == 1
#					print 'Residue == None, noneFound == 1,  therefore skipping'
					pass
#		print residueList
#		print loops
		loopIndices = []
		for sequence in range(len(residueList)):
			loopIndices.append([])
			for loop in loops[sequence]:
#				print 'Looking at loop', loop
				loopIndices[-1].append([])
				for residue in loop:
#					print 'Looking at residue', residue
					index = residueList[sequence].index(residue)
					loopIndices[-1][-1].append(index)
					residueList[sequence][index] = None
		del(residueList)
#		print loops
#		print loopIndices			# THESE ARE LIST INDICES NOT RESIDUE NUMBERS!!!!!
#		print self.groupAlphaHelix
#		print self.groupBetaSheet
		self.loopIndices = loopIndices
			
		

	def output(self):	# Output the required data back to the main program (maybe pickle it and then have a seperate program for analysis comparison?)
		return self.params, self.loopIndices, self.groupAlphaHelix, self.groupBetaSheet, self.groupUserDefined

class analysis:
	def __init__(self, paramsList, sortFactor1, mutationalityList, loopIndexElement, alpha, beta, custom):
		self.paramsList = paramsList
		self.sortFactor1 = sortFactor1
		self.mutationalityList = mutationalityList
		self.loopxvals, self.loopyvals, self.loopLengths = self.getLoopData(loopIndexElement)
		self.alpha = alpha
		self.beta = beta
		if custom != None:
			self.custom = json.loads(custom)
		else:
			self.custom = custom

	def getLoopData(self, loopIndexElement):
#		print loopIndexElement
		yvals = []
		loopLengths = []
		xvals = np.arange(0, 210)
		for i in xvals:
			yvals.append(0)
		for element in loopIndexElement[0]:
			if len(element) > 2:
#				print 'Length of loop > 2'
#				print loopIndexElement
				for i in element:
#					print 'Increasing yvals', i, 'by 1'
					yvals[i] = 1
					loopLengths.append(len(element))
		return xvals, yvals, loopLengths

	def plotter(self, mode='noloop'):
		xvals = list(np.arange(0, 210))
		P.figure()
		P.ion()
		P.xlabel('Residue Number')
#		if mode == 'Ident':
#			P.ylabel('Conservation Factor in inputted Sequences')
		if mode == 'loop':
			P.bar(self.loopxvals, self.loopyvals, color=self.makeColourMap(self.loopLengths, mode='yg'), linewidth=0)
		P.ylabel('Mutation Factor in inputted Sequences')
		P.bar(xvals, self.mutationalityList)
		P.xlim(0, 210)
		P.ylim(0, 1)
		P.draw()
	#	P.show()

	def structure(self):
		structureMatrix = np.empty((11, 210, 3), dtype=np.uint8)
		structureMatrix.fill([255, 255, 255])
#		print self.alpha
#		print self.beta
		for section in self.alpha:
			structureMatrix[0:11, section[0]:section[1]] = [0, 156, 46]
			width = section[1]-section[0]
			if width >= 8:		# Draw an Alpha
				indent = section[0]+int((width-8)/2)
				structureMatrix[5, indent+1] = [255, 255, 0]
				structureMatrix[4, indent+1] = [255, 255, 0]
				structureMatrix[6, indent+1] = [255, 255, 0]
				structureMatrix[3, indent+2] = [255, 255, 0]
				structureMatrix[7, indent+2] = [255, 255, 0]
				structureMatrix[3, indent+3] = [255, 255, 0]
				structureMatrix[7, indent+3] = [255, 255, 0]
				structureMatrix[4, indent+4] = [255, 255, 0]
				structureMatrix[6, indent+4] = [255, 255, 0]
				structureMatrix[5, indent+5] = [255, 255, 0]
				structureMatrix[4, indent+6] = [255, 255, 0]
				structureMatrix[6, indent+6] = [255, 255, 0]
				structureMatrix[3, indent+7] = [255, 255, 0]
				structureMatrix[7, indent+7] = [255, 255, 0]
		for section in self.beta:
			structureMatrix[0:11, section[0]:section[1]] = [0, 192, 237]
			width = section[1]-section[0]
			if width >= 6:		# Draw a Beta
				indent = section[0]+int((width-5)/2)
				structureMatrix[2:10, indent+1] = [255, 255, 0]
				structureMatrix[8, indent+2] = [255, 255, 0]
				structureMatrix[5, indent+2] = [255, 255, 0]
				structureMatrix[2, indent+2] = [255, 255, 0]
				structureMatrix[8, indent+3] = [255, 255, 0]
				structureMatrix[5, indent+3] = [255, 255, 0]
				structureMatrix[4, indent+3] = [255, 255, 0]
				structureMatrix[2, indent+3] = [255, 255, 0]
				structureMatrix[7, indent+4] = [255, 255, 0]
				structureMatrix[6, indent+4] = [255, 255, 0]
				structureMatrix[3, indent+4] = [255, 255, 0]
		if self.custom != None:
			for section in self.custom:
	#			print self.custom
				structureMatrix[0:11, section[0]:section[1]] = [192, 46, 46]
				width = section[1]-section[0]
				if width >= 8:		# Draw a 'U'
					indent = section[0]+int((width-5)/2)
					structureMatrix[2:7, indent] = [255, 255, 0]
					structureMatrix[7, indent+1] = [255, 255, 0]
					structureMatrix[8, indent+2] = [255, 255, 0]
					structureMatrix[8, indent+3] = [255, 255, 0]
					structureMatrix[7, indent+4] = [255, 255, 0]
					structureMatrix[2:7, indent+5] = [255, 255, 0]
#		print structureMatrix
		im = Image.fromarray(structureMatrix)
		return im

	def allostery(self, mode='noLoop', barlinewidth=None):
		allosteryReader = csv.reader(open('allval.csv', 'rb'))
		allosteryResidues = []
		allosteryValues = []
		wtallostery = []
		for row in allosteryReader:
			allosteryResidues.append(int(row[0]))
			allosteryValues.append(float(row[1])-wtvalue)	# To get bars centered around wtvalue
	#		wtallostery.append(wtvalue)
		xvals = list(np.arange(1, 211))
		for xval in range(len(xvals)):
			if xvals[xval] not in allosteryResidues:
				allosteryResidues.insert(xval, xvals[xval])
				allosteryValues.insert(xval, 0)		# To get Allostery and Residue lists the same length
	#			wtallostery.append(wtvalue)
		if mode == 'valsOnly':
			return allosteryResidues, allosteryValues, self.mutationalityList[:]
		RGBHex, colourBar = self.makeColourMap(allosteryValues, mode='rb')
		P.ion()
		fig = P.figure(figsize = (18, 10))
#		gs = gridspec.GridSpec(2, 2, width_ratios=[17, 1], height_ratios=[15, 1])
#		ax0 = P.subplot(gs[0])
#		ax1 = P.subplot(gs[1])
#		ax2 = P.subplot(gs[2])
		gs = gridspec.GridSpec(17, 15)
		axtop = P.subplot(gs[0,1:-2])
		ax0 = P.subplot(gs[2:-2,1:-2])	# Main plot
		ax1 = P.subplot(gs[2:-2,-1])	# ColourBar
		ax2 = P.subplot(gs[-1,1:-2])	# ResidueBar

		ax0.set_xlabel('Residue Number')
		ax0.set_ylabel('Magnitude of allostery difference')
#	Plot bars up to 1.5 that are yellow or something to highlight the loop residues
		if mode == 'loop':
#			print 'Plotting the Loops'
#			print self.loopxvals
#			print self.loopyvals
			ax0.bar(self.loopxvals, self.loopyvals, bottom=0.5, color=self.makeColourMap(self.loopLengths, mode='yg'), linewidth=0)
#		print allosteryValues
		ax0.bar(allosteryResidues, allosteryValues, bottom=wtvalue, color=RGBHex, linewidth=barlinewidth)
		ax0.axhline(wtvalue, color='black', lw=2)
		ax0.set_xlim(0, 210)
		ax0.set_ylim(0.776138, 1.5)	# Ensures wtVal is central

		ax0b = ax0.twinx()
		mutationVals = self.mutationalityList[:]
		xvals = list(np.arange(1, 211))
#		for element in self.mutationalityList:
#			mutationNeg.append(-element)
#		ax0b.plot(xvals, mutationNeg, '-k', linewidth=0.8)
#		for element in range(len(allosteryValues)):
#			if allosteryValues[element] < 0:
##				print "Before", mutationVals[element]
#				mutationVals[element] = -mutationVals[element]
##				print "After", mutationVals[element], "\n"
		ax0b.plot(xvals, mutationVals, 'xk', linewidth=0.8)
		ax0b.set_xlim(0, 210)
#		mutationYlim = max(self.mutationalityList)
		ax0b.set_ylim(-max(self.mutationalityList), max(self.mutationalityList))
		ax0b.set_ylabel('Mutation Factor', horizontalalignment='right', rotation=270)
#		ax0.set_ylim(-self.maximum+wtvalue, self.maximum+wtvalue)


		ax1x = []
		ax1y = np.arange(-self.maximum+wtvalue, self.maximum+wtvalue+self.increment, self.increment)
		for element in ax1y:
			ax1x.append(1)
		ax1.set_xticklabels([])
#		print colourBar
		ax1.barh(ax1y, ax1x, color=colourBar, linewidth=0)
#		print self.maximum
#		print ax1y
#		print self.colourBar
		ax1.set_xlim(0, 1)
		ax1.set_ylim(-self.maximum+wtvalue, self.maximum+wtvalue)

		ax2y = []
		for element in allosteryResidues:
			ax2y.append(1)
		ax2.set_yticklabels([])
		ax2.set_yticks((-1, 2))
		ax2.set_xticks([0, 50, 100, 150, 200], minor=True)
		ax2.xaxis.set_minor_locator(MultipleLocator(1))
	#	ax2.set_xticks(np.arange(0, 211, 1))
	#	ax2.set_xticklabels([0, 50, 100, 150, 200])
		ax2.bar(allosteryResidues, ax2y, color=RGBHex, linewidth=0)
		ax2.set_xlim(0, 210)
		ax2.set_ylim(0, 1)
	#		P.plot(allosteryResidues, wtallostery, 'r--')

		structure = self.structure()
#		print structure.shape
		axtop.set_xlim(0, 210)
		axtop.set_ylim(0, 10)
		axtop.set_yticklabels([])
#		axtop.set_xticklabels([])
#		axtop.set_xticks((-1, 211))
#		axtop.set_yticks((-1, 2))
		axtop.imshow(structure, aspect='auto', interpolation='nearest', cmap='jet')
		P.draw()
		fig.subplots_adjust(wspace=0.1,hspace=0.1,top=0.95,bottom=0.05,left=0.1,right=0.9)
		P.savefig('./temp.png')
#		P.show()
	#	print allosteryResidues
	#	print allosteryValues

	def makeColourMap(self, yvals, mode='rb'):
		colours = []
		RGBHex = []
		colourBar = []
		yvals2 = []
		if mode == 'rb':
			for element in yvals:
				yvals2.append(abs(element))
			self.maximum = np.max(yvals2)
			self.increment = (2*self.maximum)/510.	# Each increment modifies colour by '1', from (255, 0, 0) to (0, 0, 255)
			self.comparator = np.arange(-self.maximum, self.maximum+(2*self.increment), self.increment)
			for element in yvals:
				for value in range(len(self.comparator)):
					if self.comparator[value] > element:
						colours.append(value-1)
						break
			colourBarY = np.arange(0, 511)
			for element in colourBarY:
				if element <= 255:	# This is blue
					elementHex = "%.2x" % (element)
					colourBar.append('#'+(elementHex*2)+'ff')
				else:			# This is red
					elementHex = "%.2x" % (510-element)
					colourBar.append('#ff'+(elementHex*2))
			for element in colours:
				if element <= 255:	# This is blue
					elementHex = "%.2x" % (element)
					RGBHex.append('#'+(elementHex*2)+'ff')
				else:			# This is red
					elementHex = "%.2x" % (510-element)
					RGBHex.append('#ff'+(elementHex*2))
			return RGBHex, colourBar

		else:
			hexTemp = []
			maximum = np.max(yvals)
			minimum = np.min(yvals)
			colourIncrement = int(255/(maximum-minimum))
			for element in yvals:
				hexTemp.append(255-((element-minimum)*colourIncrement))
#			print yvals
#			print hexTemp
			for element in hexTemp:
				elementHex = "%.2x" % (element)
				RGBHex.append('#'+elementHex+'ff00')
#			print RGBHex
			return RGBHex

	def displayOutput(self):
		### Sort lists by SortFactor column
		if self.sortFactor1 != 0:
			reverseSort = True
		else:
			reverseSort = False
		self.sortedParamsList = sorted(self.paramsList, key=operator.itemgetter(self.sortFactor1), reverse=reverseSort)
	#	if sortFactor2 != None:
	#		sortedParamsList = sorted(sortedParamsList, key=operator.itemgetter(sortFactor2))
	#	if sortFactor3 != None:	
	#		sortedParamsList = sorted(sortedParamsList, key=operator.itemgetter(sortFactor3))
		for element in self.sortedParamsList:
			print element

def findIndex(string, logical):
	'''This function returns the locations of an inputted character (logical) in an inputted string'''
	index = 0
	locations = []
	while index < len(string):
		if string[index] == logical:
			locations.append(index)
		index += 1
	return locations

def sortParams(paramsList, identityList, loopIndexElement, alpha, beta, cutOff, sortFactor1, totalSequences):
	truncatedParams = paramsList[:]
	truncatedIdentityList = identityList[:]
	newParams = []
	newIdentity = []
	newIdentityList = []
	mutationalityList = []


#	for element in range(len(paramsList)):
#		print paramsList[element]
#		print identityList[element]


	for element in range(len(paramsList)):
#		print element
#		print element[1]
#		print paramsList[element][1], 'Cut Off =', cutOff
		if paramsList[element][1] <= cutOff:
#			print paramsList[element][1], '<', cutOff
			truncatedParams[element] = None
			truncatedIdentityList[element] = None
			totalSequences -= 1
	for element in truncatedParams:
		if element != None:
			newParams.append(element)
	for element in truncatedIdentityList:
		if element != None:
			newIdentity.append(element)
	for element in newIdentity:
		if len(newIdentityList) == 0:
			newIdentityList += element
		else:
			for ident in range(len(element)):
				newIdentityList[ident] += element[ident]



	for element in range(len(newIdentityList)):
		newIdentityList[element] /= float(totalSequences)
	for element in newIdentityList:
		mutationalityList.append(1-element)
#	print mutationalityList
	return newParams, mutationalityList, loopIndexElement, alpha, beta

def CLI(paramsList, identityList, loopIndexElement, alpha, beta, custom, totalSequences):
	print "\n---------------===============***************===============---------------"
	print "\nRESULTS"
	print "\n---------------===============***************===============---------------\n"
	sortFactor1 = 0
	cutOff = 0
	newParams, mutationalityList, loopIndexElement, alpha, beta = sortParams(paramsList, identityList, loopIndexElement, alpha, beta, cutOff, sortFactor1, totalSequences)
	analyse = analysis(newParams, sortFactor1, mutationalityList, loopIndexElement, alpha, beta, custom)
############################################################################
#		if len(identityList) == 0:
#			identityList += identity
#		else:
#			for element in range(len(identity)):
#				identityList[element] += identity[element]
#	for element in range(len(identityList)):
#		identityList[element] /= float(totalSequences)
############################################################################
#	for element in identityList:
#		mutationalityList.append(1-element)
#	print "Number of sequences with Extended N-Terminus:", len(NTerms)
#	print "Number of sequences with Extended C-Terminus:", len(CTerms)
#	print "Making new Analyse"
#	for element in range(len(paramsList)):
#		print paramsList[element][1]
#		print identityList[element]
	while True:
#		print analyse
		print '\nMAIN MENU:\n'
		print '1) List Protein Statistics'
		print '2) Plot Residue Allostery'
		print '3) Plot Residue Allostery with loop overlay'
		print '4) Settings'
		print '5) Quit\n'
		choice = raw_input('Please select an option: ')
		try:
			choice = int(choice)
		except:
			pass
		if choice == 1:
			analyse.displayOutput()
		elif choice == 2:
			analyse.allostery(barlinewidth=0)
		elif choice == 3:
			analyse.allostery(mode='loop', barlinewidth=0)
		elif choice == 4:
			while choice == 4:
				print '\nElement Choices:'
				print '\n1) Name'
				print '2) Residue identity (includes all terminus extensions and gaps)'
				print '3) Residue identity (omitting N- and C- terminus extensions)'
				print '4) Residue identity in Alpha Helices'
				print '5) Residue identity in Beta Sheets'
				print '6) Residue identity in Ligand Binding Sites'
				print '7) Residue identity in User-Defined Custom Group\n'
				sortFactor1 = int(raw_input('Please select a primary sorting option: '))-1
				cutOff = float(raw_input('Please input a residue-identity cutoff (e.g. "70.0"): '))
				newParams, sortFactor1, mutationalityList, loopIndexElement, alpha, beta = sortParams(paramsList, identityList, loopIndexElement, alpha, beta, cutOff, sortFactor1, totalSequences)
				try:
					delete(analyse)
				except:
					pass
				analyse = analysis(newParams, sortFactor1, mutationalityList, loopIndexElement, alpha, beta)
				print 'Returning to main menu.\n'
				choice = 0
		elif choice == 5:
			break
	print 'Thank you for using this program. Have a nice day!'
	

if __name__ == "__main__":
	print "\n---------------===============***************===============---------------"
	print "\nMatty Jones' Wonderful Protein Sequence Alignment Analysis Tool"
	print "\n---------------===============***************===============---------------\n"

	orderedList, totalSequences = clustalw(directory='./Sequences').runSetup()
	paramsList = []
	identityList = []
	for element in orderedList:
		paramsListElement, identity, loopIndexElement, alpha, beta, custom = sequences(element).runAnalysis()
		paramsList.append(paramsListElement)
		identityList.append(identity)

	CLI(paramsList, identityList, loopIndexElement, alpha, beta, custom, totalSequences)
#	print paramsList




