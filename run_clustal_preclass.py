import numpy as np
import os
import cPickle as pickle

def importSequences(directory='./Sequences'):
	list_files = os.listdir(directory)
	list_seq = []
	for n in range(0,len(list_files)):
		if list_files[n][-3:] == 'txt':
			list_seq.append(directory+'/'+list_files[n])
	print 'Number of sequences found:',len(list_seq),'\n'
#		if len(list_seq) == 0:
#			directory = str(raw_input('Error: No sequences found. Please enter a new directory: '))
#		else:
#			list_seq.insert(0,list_seq[-1])
#			del list_seq[-1]
#			break
	return list_seq

def checkAlign(name):
	nameAln = name[:-3]+'aln'
	check = os.system('ls '+nameAln+' >/dev/null 2>/dev/null')
	return check, nameAln

def runClustal(list_seq):
	overwriteFlag = 1
	list_aln = []
	for element in list_seq:
		print 'Checking for binary alignment of '+str(element)+'...'
		check, nameAln = checkAlign(element)
		if check == 0 and overwriteFlag == 0:
			print 'Alignment file already present!'
		else:
			print 'Alignment file not present. Executing CLUSTAL...'
			os.system('clustalw ' + str(element))
			check, nameAln = checkAlign(element)
		list_aln.append(nameAln)
	return list_aln

def reordering(aln):
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

def analysis(alignment, groupUserDefined):
	groupAlphaHelix = [20, 30]
	groupBetaSheet = [40, 50]
	groupLigand = [60, 70]
	### Obtain residue numbers ###
	residueNo1 = []
	residueNo2 = []
	##############################

	matrix = np.asarray(alignment)
	print matrix
	### Calculate identicality ###
	identicality = [0]
	height, width = matrix.shape
	for i in range(1, width):			#### This will only work for binary alignments
#		print "Comparing", matrix[0,i], "with", matrix[1,i]
		if matrix[0,i] == matrix[1,i]:
			identicality.append(1)
		else:
			identicality.append(0)
	identicalityPercentage = float(100*sum(identicality))/float(len(identicality)-1)		# Don't include the name in the identicality calculation!
	test = 100
#	print identicalityPercentage
	print 'Sequence', matrix[1,0], 'is %.2f%% identical to sequence' % identicalityPercentage , str(matrix[0,0])+'.'

		

if __name__ == "__main__":
	print "\n---------------===============***************===============---------------"
	print "\nTitle Placeholder"
	print "\n---------------===============***************===============---------------\n"
	list_seq = importSequences()
	list_aln = runClustal(list_seq)
	orderedList = []
	for element in list_aln:	
		orderedList.append(reordering(element))
	groupUserDefined = [80, 90]
	for element in orderedList:
		analysis(element, groupUserDefined)






