import numpy as np
import os

def importAlignments(directory='./Sequences'):
	list_files = os.listdir(directory)
	list_aln = []
	for n in range(0,len(list_files)):
		if list_files[n][-3:] == 'aln':
			list_aln.append(directory+'/'+list_files[n])
	print 'Number of alignments found:',len(list_aln),'\n'
#		if len(list_seq) == 0:
#			directory = str(raw_input('Error: No sequences found. Please enter a new directory: '))
#		else:
#			list_seq.insert(0,list_seq[-1])
#			del list_seq[-1]
#			break
	return list_aln


if __name__ == "__main__":
	print "\n---------------===============***************===============---------------"
	print "\nTitle Placeholder"
	print "\n---------------===============***************===============---------------\n"
	list_aln = importAlignments()
	print list_aln
