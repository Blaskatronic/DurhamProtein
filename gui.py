#import ttk
#import Tkinter
import run_clustal
import subprocess as sub
import tkMessageBox as box
import os
import numpy as np
import pylab as P
import matplotlib.gridspec as gridspec
import csv
from ttk import *
from Tkinter import *
from PIL import Image, ImageTk

class mainScreen(Frame):
	def __init__(self, parent):
		Frame.__init__(self, parent)
		self.parent = parent
#		self.centreWindow(parent)
		self.initUI()
		self.groupUserDefined = None
#		self.output = 0

	def centreWindow(self, win):
		win.update_idletasks()
		self.width = win.winfo_width()
		self.height = win.winfo_height()
		screenWidth = self.parent.winfo_screenwidth()
		screenHeight = self.parent.winfo_screenheight()
		self.x = (screenWidth - self.width)/2
		self.y = (screenHeight - self.height)/2
#		print self.parent.winfo_geometry()
		win.geometry('%dx%d+%d+%d' % (self.width, self.height, self.x, self.y))

	def inputScreen(self):
		self.win = Toplevel(self.parent)
		self.win.title("Input Sequence")
		self.pack(fill=BOTH, expand=1)
		Style().configure("TFrame")#, background="#FFFFFF")
#		win.geometry('%dx%d+%d+%d' % (self.width, self.height, self.x, self.y))
		frame = Frame(self.win, relief=RAISED, borderwidth=1)
		frame.pack(fill=BOTH)

		lbl = Label(self.win, text="Please input base sequence in FASTA format (or leave blank for Ecoli CAP)")
		lbl.pack()

		area = Text(self.win, width=16, height=8)
		area.pack(fill=X, expand=1, padx=5, pady=5)

		lbl2 = Label(self.win, text="Please input protein sequences to compare to the base in FASTA format")
		lbl2.pack()

		area2 = Text(self.win)
		area2.pack(fill=X, expand=1, padx=5, pady=5)

		obtn = Button(self.win, text="Input", command= lambda:self.inputSeq(area,area2))
		obtn.pack(side=RIGHT, padx=5, pady=5)

		cbtn = Button(self.win, text="Exit", command=self.win.destroy)
		cbtn.pack(side=LEFT, padx=5, pady=5)
		
		hbtn = Button(self.win, text="Help")
#		hbtn.grid(row=5, column=1)
		hbtn.pack(side=RIGHT, padx=5, pady=5)
	
		self.centreWindow(self.win)

	def inputSeq(self, textBoxBase, textBoxSeqList):
		contentsBase = textBoxBase.get(1.0, END)
		if len(contentsBase) == 1:
			ecoliFile = open('./ecoli.txt', 'r')
			baseSequence = ecoliFile.read()
			ecoliFile.close()
		else:
			baseSequence = contentsBase
		contents = textBoxSeqList.get(1.0, END)
		separatedContents = contents.split('>')[1:]
		try:
			for sequence in separatedContents:
				baseSequencePipes = run_clustal.findIndex(baseSequence, '|')
				sequencePipes = run_clustal.findIndex(sequence, '|')
				sequence = '>'+sequence
				if len(contentsBase) == 1:
					name = sequence[sequencePipes[-2]+2:sequencePipes[-1]+1]+'_with_K12CAP'

				else:
					name = sequence[sequencePipes[-2]+2:sequencePipes[-1]+1]+'_with_'+baseSequence[baseSequencePipes[-2]+2:baseSequencePipes[-1]+1]
				print name
				testFile = open('./test.txt', 'w+')
				testFile.write(baseSequence)
				testFile.write('\n')
				testFile.write(sequence)
				testFile.close()
				os.rename('./test.txt', './Sequences/'+str(name)+'.txt')
			if len(separatedContents) >= 2:
				self.onInfo(str(len(separatedContents))+' sequences Inputted into System')
			else:
				self.onInfo('Sequence Inputted into System')
			textBox.delete(1.0, END)
		except:
			self.onInfo("Invalid Format, please input a sequence in a valid FASTA format (including pipes '|')")
		self.win.destroy()


	def initUI(self):
		self.parent.title("Main Windows")
		self.pack(fill=BOTH, expand=1)
		Style().configure("TFrame")#, background="#FFFFFF")

		self.columnconfigure(0, weight=1)

		self.rowconfigure(0, weight=1)	# Image
		self.rowconfigure(1, weight=1)	# Run
		self.rowconfigure(2, weight=1)	# Input
		self.rowconfigure(3, weight=1)	# Blank
		self.rowconfigure(4, weight=1)	# Settings
		self.rowconfigure(5, weight=1)	# Quit

		BSILogoOpen = Image.open("./UIFiles/BSILogo.jpg")
		BSILogo = ImageTk.PhotoImage(BSILogoOpen)
		label1 = Label(self, image=BSILogo)
		label1.image = BSILogo
		label1.grid(row=0, column=0, sticky=N+E+W, pady=5, padx=5)


		try:
			self.paramsList[:]
			runButton = Button(self, text="Run Analysis", command=self.executeScreen)
			runButton.grid(row=1, column=0, padx=5, pady=4)
		except:
			self.runButton = Button(self, text="Load Sequences", command=self.load)
			self.runButton.grid(row=1, column=0, padx=5, pady=4)

		inputButton = Button(self, text="Input New Sequences", command=self.inputScreen)
		inputButton.grid(row=2, column=0, padx=5, pady=4)
		
		settingsButton = Button(self, text="Settings", command=self.settingsScreen)
		settingsButton.grid(row=4, column=0, padx=5, pady=4)
		
		quitButton = Button(self, text="Exit", command=self.quit)
		quitButton.grid(row=5, column=0, padx=5, pady=4)

		self.centreWindow(self.parent)


#		self.parent.title("Main Window")
#		self.pack(fill=Tkinter.BOTH, expand=1)
#		ttk.Style().configure("TFrame", background="#FFFFFF")

#		frame = ttk.Frame(self, relief=Tkinter.RAISED, borderwidth=1)
#		frame.pack(fill=Tkinter.BOTH, expand=1)

#		BSILogoOpen = Image.open("./UIFiles/BSILogo.jpg")
#		BSILogo = ImageTk.PhotoImage(BSILogoOpen)
#		label1 = Tkinter.Label(self, image=BSILogo)
#		label1.image = BSILogo
#		label1.place(x=(self.width/2)-82, y=5)

#		quitButton = ttk.Button(self, text="Exit", command=self.quit)
#		quitButton.pack(side=Tkinter.RIGHT, padx=5, pady=4)

#		inputButton = ttk.Button(self, text="Input Sequence", command=self.inputScreen)
#		inputButton.pack(side=Tkinter.RIGHT)

#		runButton = ttk.Button(self, text="Run CLI", command=self.test)
#		runButton.pack(side=Tkinter.LEFT, padx=5, pady=4)

	def onInfo(self, text):
		box.showinfo("Information", text)

	def onScale(self, val):
		val = int(float(val))
		self.cutOff.set(val)
		self.onUpdate(val)		

	def settingsScreen(self):
		self.win4 = Toplevel(self.parent)
		self.win4.title("Settings")
		self.pack(fill=BOTH, expand=1)
		Style().configure("TFrame")

		customGroupInfo = Label(self.win4, text="Input custom group residues:")
		customGroupInfo.grid(row=0, column=0, padx=5, pady=5)

		self.customGroup = Entry(self.win4)
		self.customGroup.grid(row=1, column=0, padx=5, pady=5)

		acceptButton = Button(self.win4, text="Accept", command=self.closeSettings)
		acceptButton.grid(row=2, column=0, padx=5, pady=20)

	def closeSettings(self):
		self.groupUserDefined = self.customGroup.get()
		self.win4.destroy()

	def executeScreen(self):
		win2 = Toplevel(self.parent)
		win2.title("Execute")
		self.pack(fill=BOTH, expand=1)
		Style().configure("TFrame")#, background="#FFFFFF")
#		win.geometry('%dx%d+%d+%d' % (self.width, self.height, self.x, self.y))
#		frame = Frame(win2, relief=RAISED, borderwidth=1)
#		frame.pack(fill=BOTH)

#		area = Text(win2)
#		area.pack(fill=BOTH, expand=1, padx=5, pady=5)

		cutOffInfo = Label(win2, text="Select residue identity cut-off:")
		cutOffInfo.grid(row=0, columnspan=2)

		self.cutOff = IntVar()
		self.loopOverlay = IntVar()
		self.splittage = IntVar()

		cutOffScale = Scale(win2, from_=60, to=99, length=400, orient=HORIZONTAL, command=self.onScale)
		cutOffScale.set(60)
		cutOffScale.grid(row=1, columnspan=2, sticky=W)

		loopButton = Checkbutton(win2, text="Show Loop Overlay", variable=self.loopOverlay)
		loopButton.grid(row=2, padx=5, columnspan=2)

#		cutOffOutput = Label(win2, text=0, textvariable=self.cutOff, width=3)
#		cutOffOutput.grid(row=1, column=2)

		plotButton = Button(win2, text="Plot residue mutation and allostery", command=self.makePlot)
		plotButton.grid(row=4, padx=5, pady=15, columnspan=2)

		distroButton = Button(win2, text="Plot Mutation Distribution", command=self.makeDistro)
		distroButton.grid(row=5, padx=5, pady=5, columnspan=2)

		corrButton = Button(win2, text="Plot Correlation Function", command=self.makeCorre)
		corrButton.grid(row=6, padx=5, pady=5)

		corrModeButton = Checkbutton(win2, text="Split plots", variable=self.splittage)
		corrModeButton.grid(row=6, column=1, padx=5)

		quitButton = Button(win2, text="Exit", command=win2.destroy)
		quitButton.grid(row=7, column=1, padx=5, pady=5, sticky=SE)

		self.centreWindow(win2)

	def makeDistro(self):
		P.ion()
		fig = P.figure(figsize = (18, 10))
#		print self.mutationOccurence
#		print np.sum(self.mutationOccurence.values())
		xvals = self.mutationOccurence.keys()
		yvals = self.mutationOccurence.values()
		P.plot(xvals, yvals, 'xb')
#		P.xlim(0, 210)
		P.xlabel('Number of mutations per sequence')
		P.ylabel('Frequency of mutations')
		P.draw()
		
	def makeCorre(self):
		mode = self.splittage.get()
		if mode == 1:
			allosteryResidues, allosteryValues, mutationalityList = self.analyse.allostery(mode='valsOnly')
	#		print allosteryResidues, allosteryValues
			allosteryNegative = []
			allosteryPositive = []
			mutNegative = []
			mutPositive = []
			for element in range(len(allosteryResidues)):
				if allosteryValues[element] <= 0:
					allosteryNegative.append(-allosteryValues[element])
					allosteryPositive.append(0)
					mutNegative.append(mutationalityList[element])
					mutPositive.append(0)
				else:
					allosteryNegative.append(0)
					allosteryPositive.append(allosteryValues[element])
					mutNegative.append(0)
					mutPositive.append(mutationalityList[element])

			allo1, allo2, allo3 = self.threeMax(allosteryNegative)
			mut1, mut2, mut3 = self.threeMax(mutNegative)

			allo4, allo5, allo6 = self.threeMax(allosteryPositive)
			mut4, mut5, mut6 = self.threeMax(mutPositive)

			residueAllo1 = [i for i,x in enumerate(allosteryNegative) if x == allo1][0]
			residueAllo2 = [i for i,x in enumerate(allosteryNegative) if x == allo2][0]
			residueAllo3 = [i for i,x in enumerate(allosteryNegative) if x == allo3][0]

			residueMut1 = [i for i,x in enumerate(mutNegative) if x == mut1][0]
			residueMut2 = [i for i,x in enumerate(mutNegative) if x == mut2][0]
			residueMut3 = [i for i,x in enumerate(mutNegative) if x == mut3][0]

			residueAllo4 = [i for i,x in enumerate(allosteryPositive) if x == allo4][0]
			residueAllo5 = [i for i,x in enumerate(allosteryPositive) if x == allo5][0]
			residueAllo6 = [i for i,x in enumerate(allosteryPositive) if x == allo6][0]

			residueMut4 = [i for i,x in enumerate(mutPositive) if x == mut4][0]
			residueMut5 = [i for i,x in enumerate(mutPositive) if x == mut5][0]
			residueMut6 = [i for i,x in enumerate(mutPositive) if x == mut6][0]

			P.ion()
			fig = P.figure(figsize = (18, 10))
			gs = gridspec.GridSpec(1, 2)
			ax1 = P.subplot(gs[0])
			ax2 = P.subplot(gs[1])
	#		print mutationalityList
	#		print mutationalityList
			ax1.plot(allosteryNegative, mutNegative, 'xb')
			ax1.set_ylabel('Mutation Frequency')
			ax1.set_xlabel('$\mathrm{\Delta}$k [Negative Values]')
			ax1.annotate(residueAllo1, xy=(allo1, mutNegative[residueAllo1]), xytext=(allo1-(0.007*(len(str(residueAllo1))/3)), mutNegative[residueAllo1]+0.01))
			ax1.annotate(residueAllo2, xy=(allo2, mutNegative[residueAllo2]), xytext=(allo2-(0.007*(len(str(residueAllo2))/3)), mutNegative[residueAllo2]+0.01))
			ax1.annotate(residueAllo3, xy=(allo3, mutNegative[residueAllo3]), xytext=(allo3-(0.007*(len(str(residueAllo3))/3)), mutNegative[residueAllo3]+0.01))
			ax1.annotate(residueMut1, xy=(allosteryNegative[residueMut1], mut1), xytext=(allosteryNegative[residueMut1]-(0.007*(len(str(residueAllo1))/3)), mut1+0.01))
			ax1.annotate(residueMut2, xy=(allosteryNegative[residueMut2], mut2), xytext=(allosteryNegative[residueMut2]-(0.007*(len(str(residueAllo2))/3)), mut2+0.01))
			ax1.annotate(residueMut3, xy=(allosteryNegative[residueMut3], mut3), xytext=(allosteryNegative[residueMut3]-(0.007*(len(str(residueAllo3))/3)), mut3+0.01))
			ax1.set_ylim(-0.001, 0.75)
			ax1.set_xlim(0, 0.35)

			ax2.plot(allosteryPositive, mutPositive, 'xr')
			ax2.set_ylabel('Mutation Frequency')
			ax2.set_xlabel('$\mathrm{\Delta}$k [Positive Values]')
			ax2.annotate(residueAllo4, xy=(allo4, mutPositive[residueAllo4]), xytext=(allo4-(0.007*(len(str(residueAllo4))/3)), mutPositive[residueAllo4]+0.01))
			ax2.annotate(residueAllo5, xy=(allo5, mutPositive[residueAllo5]), xytext=(allo5-(0.007*(len(str(residueAllo5))/3)), mutPositive[residueAllo5]+0.01))
			ax2.annotate(residueAllo6, xy=(allo6, mutPositive[residueAllo6]), xytext=(allo6-(0.007*(len(str(residueAllo6))/3)), mutPositive[residueAllo6]+0.01))
			ax2.annotate(residueMut4, xy=(allosteryNegative[residueMut4], mut4), xytext=(allosteryNegative[residueMut4]-(0.007*(len(str(residueAllo4))/3)), mut4+0.01))
			ax2.annotate(residueMut5, xy=(allosteryNegative[residueMut5], mut5), xytext=(allosteryNegative[residueMut5]-(0.007*(len(str(residueAllo5))/3)), mut5+0.01))
			ax2.annotate(residueMut6, xy=(allosteryNegative[residueMut6], mut6), xytext=(allosteryNegative[residueMut6]-(0.007*(len(str(residueAllo6))/3)), mut6+0.01))
			ax2.set_ylim(-0.001, 0.75)

		else:
			allosteryResidues, allosteryValues, mutationalityList = self.analyse.allostery(mode='valsOnly')
	#		print allosteryResidues, allosteryValues
			for element in range(len(allosteryValues)):
				if allosteryValues[element] <= 0:
					allosteryValues[element] *= -1

			averageAllostery = np.average(allosteryValues)
			averageMutation = np.average(mutationalityList)

			allo1, allo2, allo3 = self.threeMax(allosteryValues)
			mut1, mut2, mut3 = self.threeMax(mutationalityList)

			residueAllo1 = [i for i,x in enumerate(allosteryValues) if x == allo1][0]
			residueAllo2 = [i for i,x in enumerate(allosteryValues) if x == allo2][0]
			residueAllo3 = [i for i,x in enumerate(allosteryValues) if x == allo3][0]

			residueMut1 = [i for i,x in enumerate(mutationalityList) if x == mut1][0]
			residueMut2 = [i for i,x in enumerate(mutationalityList) if x == mut2][0]
			residueMut3 = [i for i,x in enumerate(mutationalityList) if x == mut3][0]

	#		print residueAllo1
	#		print residueAllo2
	#		print residueAllo3
		
			P.ion()
			fig = P.figure(figsize = (18, 10))

			P.plot(allosteryValues, mutationalityList, 'xb')
			P.fill([0, averageAllostery, averageAllostery, 0], [0, 0, averageMutation, averageMutation], 'r', alpha=0.5) 
			P.ylabel('Mutation Frequency')
			P.xlabel('$\mathrm{\Delta}$k')
			P.annotate(residueAllo1, xy=(allo1, mutationalityList[residueAllo1]), xytext=(allo1+0.001, mutationalityList[residueAllo1]+0.01))
			P.annotate(residueAllo2, xy=(allo2, mutationalityList[residueAllo2]), xytext=(allo2+0.001, mutationalityList[residueAllo2]+0.01))
			P.annotate(residueAllo3, xy=(allo3, mutationalityList[residueAllo3]), xytext=(allo3+0.001, mutationalityList[residueAllo3]+0.01))
			P.annotate(residueMut1, xy=(allosteryValues[residueMut1], mut1), xytext=(allosteryValues[residueMut1]+0.001, mut1+0.01))
			P.annotate(residueMut2, xy=(allosteryValues[residueMut2], mut2), xytext=(allosteryValues[residueMut2]+0.001, mut2+0.01))
			P.annotate(residueMut3, xy=(allosteryValues[residueMut3], mut3), xytext=(allosteryValues[residueMut3]+0.001, mut3+0.01))
			P.ylim(-0.001, 0.75)
	#		print allo1, mutationalityList[residueAllo1]
	#		print len(allosteryValues), len(mutationalityList)
			P.draw()
			
			aboveKaboveM = 0
			aboveKbelowM = 0
			belowKaboveM = 0
			belowKbelowM = 0
			for i in range(len(allosteryValues)):
				if allosteryValues[i] >= averageAllostery:
					if mutationalityList[i] >= averageMutation:
						aboveKaboveM += 1
					else:
						aboveKbelowM += 1
				else:
					if mutationalityList[i] >= averageMutation:
						belowKaboveM += 1
					else:
						belowKbelowM += 1
			print "Given that the allostery is below the average:"
			print "There are", belowKaboveM, "residues above the average mutation factor, and", belowKbelowM, "below it."
			print "Given that the allostery is above the average:"
			print "There are", aboveKaboveM, "residues above the average mutation factor, and", aboveKbelowM, "below it."


			csvFile = open('./dataFile.csv', 'w+')
			dataWriter = csv.writer(csvFile)
			dataWriter.writerow(['Residue Number', 'Allostery', 'Mutation Factor'])
			for i in range(len(allosteryValues)):
				dataWriter.writerow([i, allosteryValues[i], mutationalityList[i]])
			#print len(allosteryValues)
			#print len(mutationalityList)

	def threeMax(self, listname):
		listName = listname[:]
		max1 = np.max(listName)
		listName.remove(max1)
		max2 = np.max(listName)
		listName.remove(max2)
		max3 = np.max(listName)
		listName.remove(max3)
		return max1, max2, max3		

	def makePlot(self):
		mode = self.loopOverlay.get()
		if mode == 1:
			mode = "loop"
		else:
			mode = "noLoop"
		self.analyse.allostery(barlinewidth=0, mode=mode)

#		self.win3 = Toplevel(self.parent)
#		self.win3.title("PUT CUTOFF HERE")
#		Style().configure("TFrame")
#		screenWidth = self.parent.winfo_screenwidth()
#		screenHeight = self.parent.winfo_screenheight()
#		figWidth = screenWidth*0.7
#		figHeight = screenHeight*0.7
##		self.x = (screenWidth - figWidth)/2
##		self.y = (screenHeight - figHeight)/2
###		print self.parent.winfo_geometry()
##		win3.geometry('%dx%d+%d+%d' % (figWidth, figHeight, self.x, self.y))
##		graphOpen = Image.open("./temp.png")
##		graph = ImageTk.PhotoImage(graphOpen)
#		Grid.columnconfigure(self.win3, 0, weight=1)
#		Grid.rowconfigure(self.win3, 0, weight=1)	
##		canvas = Canvas(win3, width = figWidth, height = figHeight)
##		canvas.pack(expand = YES, fill = BOTH)
#		self.updateImage(figWidth, figHeight)

#	def reCalculate(self, figWidth, figHeight):
#		self.analyse.allostery(mode='loop', barlinewidth=0)
#		self.label1.destroy()
##		graph.destroy()
#		self.updateImage(figWidth, figHeight)
#		
#	def updateImage(self, figWidth, figHeight):
##		self.analyse.allostery(mode='loop', barlinewidth=0)
#		graphOpen = Image.open("./temp.png")
#		graphOpen = graphOpen.resize((int(figWidth), int(figHeight)), Image.ANTIALIAS)
#		global graph
#		graph = ImageTk.PhotoImage(graphOpen)
##		graphOpen.close()
##		canvas.create_image(50, 10, image = graph, anchor = NW)
#		self.label1 = Label(self.win3, image=graph)
#		self.label1.image = graph
#		self.label1.grid(sticky=N+E+S+W, pady=5, padx=5)
##		self.win3.after(3000, lambda: self.reCalculate(figWidth, figHeight))

# Every TIME_TO_PLOT seconds check if bar has changed, if has then replot, otherwise sleep.



	def load(self):
		orderedList, totalSequences = run_clustal.clustalw(directory='./Sequences').runSetup()
		paramsList = []
		identityList = []
		mutationDistribution = []
		for element in orderedList:
			paramsListElement, identity, loopIndexElement, alpha, beta, custom = run_clustal.sequences(element, self.groupUserDefined).runAnalysis()
			paramsList.append(paramsListElement)
			identityList.append(identity)


		mutationOccurence = {}
		for sequence in identityList:
			numberOfMutations = len(sequence)-np.sum(sequence)
#			print "This sequence has", numberOfMutations
#			print mutationOccurence
			if numberOfMutations in mutationOccurence:
				mutationOccurence[numberOfMutations] += 1
			else:
				mutationOccurence[numberOfMutations] = 1
#		print mutationOccurence

###########################################################################################
		for residues in range(len(identityList[0])):
			mutationDistribution.append(0)
		for sequence in range(len(identityList)):
			for residue in range(len(identityList[sequence])):
				if identityList[sequence][residue] == 0:
#					print "Mutation found in", sequence, residue
					mutationDistribution[residue] += 1
############################################################################################

		self.paramsList = paramsList
		self.identityList = identityList
		self.loopIndexElement = loopIndexElement
		self.alpha = alpha
		self.beta = beta
		self.totalSequences = totalSequences
		self.mutationOccurence = mutationOccurence
		self.mutationDistribution = mutationDistribution

		cutOff = 60
		newParams, mutationalityList, loopIndexElement, alpha, beta = run_clustal.sortParams(self.paramsList, self.identityList, self.loopIndexElement, self.alpha, self.beta, cutOff, 1, self.totalSequences)
		try:
			delete(analyse)
		except:
			pass
		self.analyse = run_clustal.analysis(newParams, 1, mutationalityList, loopIndexElement, alpha, beta, self.groupUserDefined)
		self.onInfo("Sequences Loaded.")
		self.runButton.destroy()
		self.initUI()
#		CLI(paramsList, identityList, loopIndexElement, alpha, beta, totalSequences)

	def onUpdate(self, cutOff):
		newParams, mutationalityList, loopIndexElement, alpha, beta = run_clustal.sortParams(self.paramsList, self.identityList, self.loopIndexElement, self.alpha, self.beta, cutOff, 1, self.totalSequences)
		try:
			delete(analyse)
		except:
			pass
		self.analyse = run_clustal.analysis(newParams, 1, mutationalityList, loopIndexElement, alpha, beta, self.groupUserDefined)
		print len(newParams), "sequences to be included into the plot."

def GUI():
	root = Tk()
	app = mainScreen(root)
	root.mainloop()

if __name__ == "__main__":
	print "\n---------------===============***************===============---------------"
	print "\nInitialising GUI"
	print "\n---------------===============***************===============---------------\n"
	GUI()














