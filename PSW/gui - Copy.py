#import ttk
#import Tkinter
import run_clustal
import subprocess as sub
import tkMessageBox as box
import os
import glob
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
#		try:
		for sequence in separatedContents:
		         	baseSequencePipes = run_clustal.findIndex(baseSequence, '|')
				sequencePipes = run_clustal.findIndex(sequence, '|')
				sequence = '>'+sequence
				if len(contentsBase) == 1:
					name = sequence[sequencePipes[-2]+2:sequencePipes[-1]+1]+'_with_K12CAP'

				else:
					name = sequence[sequencePipes[-2]+2:sequencePipes[-1]+1]+'_with_'+baseSequence[baseSequencePipes[-2]+1:baseSequencePipes[-1]]
				print name, 'added to the system.'
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
#			textBoxSeqList.delete(1.0, END)
#		except:
#          		self.onInfo("Invalid Format, please input a sequence in a valid FASTA format (including pipes '|')")
		self.win.destroy()
		self.initUI()


	def initUI(self):
	        self.numberOfSeqs = 0
	        self.numberOfAlns = 0
		self.parent.title("Protein Sequencer")
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

                listFiles = os.listdir("./Sequences")
		for n in range(0,len(listFiles)):
			if listFiles[n][-3:] == 'txt':
				self.numberOfSeqs += 1
		        elif listFiles[n][-3:] == 'aln':
		                self.numberOfAlns += 1
				
		try:
			self.paramsList[:]
			runButton = Button(self, text="Run Analysis ("+str(self.numberOfAlns)+" alignments)", command=self.executeScreen)
			runButton.grid(row=1, column=0, padx=5, pady=4)
		except:
			self.runButton = Button(self, text="Load Sequences ("+str(self.numberOfSeqs)+" sequences)", command=self.load)
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
		
		self.purgeSeq = IntVar()
		self.purgeAln = IntVar()

		customGroupInfo = Label(self.win4, text="Input custom group residues:")
		customGroupInfo.grid(row=0, column=0, padx=5, pady=5)

		self.customGroup = Entry(self.win4)
		self.customGroup.insert(0, "None")
		self.customGroup.grid(row=1, column=0, padx=5, pady=5)
		
		purgeSeqModeButton = Checkbutton(self.win4, text="Purge All Sequences", variable=self.purgeSeq)
		purgeSeqModeButton.grid(row=2, column=0, padx=5, pady=20)
		
		purgeAlnModeButton = Checkbutton(self.win4, text="Purge All Alignments", variable=self.purgeAln)
		purgeAlnModeButton.grid(row=3, column=0, padx=5)

		acceptButton = Button(self.win4, text="Accept", command=self.closeSettings)
		acceptButton.grid(row=4, column=0, padx=5, pady=40)
		
	def onSure(self, text):
	       box.askokcancel("Are you sure?", text)
		
        def purge(self, items):
                noOfFiles = 0
                if items == 'seq':
                    self.onSure("About to delete "+str(self.numberOfSeqs)+" sequences...")
                    for file in glob.glob(".\\Sequences\\*.txt"):
                        os.remove(file)
                        noOfFiles += 1
                    self.onInfo("Purge complete.")
                elif items == 'aln':
                    self.onSure("About to delete "+str(self.numberOfAlns)+" alignments...")
                    for file in glob.glob(".\\Sequences\\*.aln"):
                        os.remove(file)
                        noOfFiles += 1
                    for file in glob.glob(".\\Sequences\\*.dnd"):
                        os.remove(file)
                        noOfFiles += 1
                    self.onInfo("Purge complete.")
                elif items == 'both':
                    self.onSure("About to delete "+str(self.numberOfSeqs)+" sequences, and their alignments...")
                    for file in glob.glob(".\\Sequences\\*.txt"):
                        os.remove(file)
                        noOfFiles += 1
                    for file in glob.glob(".\\Sequences\\*.aln"):
                        os.remove(file)
                        noOfFiles += 1
                    for file in glob.glob(".\\Sequences\\*.dnd"):
                        os.remove(file)
                        noOfFiles += 1
                    self.onInfo("Purge complete.")
                print "Successfully deleted", noOfFiles, "files."
                self.initUI()
                    

	def closeSettings(self):
	        seq = self.purgeSeq.get()
	        aln = self.purgeAln.get()
	        if seq == 1 and aln == 1:
	               self.purge('both')
	        elif seq == 1:
	               self.purge('seq')
	        elif aln == 1:
	               self.purge('aln')
		self.groupUserDefined = self.customGroup.get()
		if self.groupUserDefined == "None":
         	      self.groupUserDefined = None
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

                self.fullAllosteryResidues, self.fullAllosteryValues, self.fullMutationalityList = self.analyse.allostery(mode='valsOnly')
                
                ######## Truncate Allostery Values to length of Mutation list (Allostery data not gathered yet) ########
                
                del self.fullAllosteryResidues[len(self.fullMutationalityList):]
                del self.fullAllosteryValues[len(self.fullMutationalityList):]
                
#                print len(self.fullMutationalityList), "MUT LIST"
#                print len(self.fullAllosteryValues), "ALL LIST"
#                print len(self.fullAllosteryResidues), "ALL RES LIST"
                
                ########################################################################################################
                
		cutOffInfo = Label(win2, text="Select residue identity cut-off:")
		cutOffInfo.grid(row=0, columnspan=3)

		self.cutOff = IntVar()
		self.loopOverlay = IntVar()
		self.splittage = IntVar()
		self.annotate = IntVar()
#		self.corrXminVal = IntVar()
#		self.corrXmaxVal = IntVar()
#		self.corrYminVal = IntVar()
#		self.corrYmaxVal = IntVar()

		cutOffScale = Scale(win2, from_=30, to=99, length=900, orient=HORIZONTAL, command=self.onScale)
		cutOffScale.set(30)
		cutOffScale.grid(row=1, columnspan=3, sticky=W)

#		loopButton = Checkbutton(win2, text="Show Loop Overlay", variable=self.loopOverlay)
#		loopButton.grid(row=2, padx=5, columnspan=3)

#		cutOffOutput = Label(win2, text=0, textvariable=self.cutOff, width=3)
#		cutOffOutput.grid(row=1, column=2)

		plotButton = Button(win2, text="Plot Mutation Factor for each Residue", command=self.makePlot)
		plotButton.grid(row=4, padx=5, pady=15, columnspan=3)

		distroButton = Button(win2, text="Plot Mutation Distribution", command=self.makeDistro)
		distroButton.grid(row=5, padx=5, pady=5, columnspan=3)

		corrButton = Button(win2, text="Plot Correlation Function", command=self.makeCorre)
		corrButton.grid(row=6, padx=5, pady=5)

		corrModeButton = Checkbutton(win2, text="Split plots", variable=self.splittage)
		corrModeButton.grid(row=6, column=1, padx=5)
		
		annotateModeButton = Checkbutton(win2, text="Display residues", variable=self.annotate)
		annotateModeButton.grid(row=6, column=2, padx=5)
		
           	corrXMinInfo = Label(win2, text="Select minimum allostery cut-off (percentage of max allostery value):")
		corrXMinInfo.grid(row=7, pady = 5)

           	corrXMaxInfo = Label(win2, text="Select maximum allostery cut-off (percentage of max allostery value):")
		corrXMaxInfo.grid(row=7, column=2, pady = 5)	
						
		corrXmin = Scale(win2, from_=0, to=100, length=400, orient=HORIZONTAL, command=self.xmin)
		corrXmin.set(0)
		corrXmin.grid(row=8, sticky=W)
		
		corrXmax = Scale(win2, from_=0, to=100, length=400, orient=HORIZONTAL, command=self.xmax)
		corrXmax.set(100)
		corrXmax.grid(row=8, column=2, sticky=E)
		
           	corrYMinInfo = Label(win2, text="Select minimum mutation cut-off (percentage of min mutation value):")
		corrYMinInfo.grid(row=9, pady = 5)

           	corrYMaxInfo = Label(win2, text="Select maximum mutation cut-off (percentage of max mutation value):")
		corrYMaxInfo.grid(row=9, column=2, pady = 5)

       		corrYmin = Scale(win2, from_=0, to=100, length=400, orient=HORIZONTAL, command=self.ymin)
		corrYmin.set(0)
		corrYmin.grid(row=10, sticky=W)
		
		corrYmax = Scale(win2, from_=0, to=100, length=400, orient=HORIZONTAL, command=self.ymax)
		corrYmax.set(100)
		corrYmax.grid(row=10, column=2, sticky=E)

		quitButton = Button(win2, text="Exit", command=win2.destroy)
		quitButton.grid(row=11, column=2, padx=5, pady=5, sticky=SE)

		self.centreWindow(win2)
		
        def xmin(self, val):
                val = int(float(val))
                self.corrXminVal = val

        def xmax(self, val):
                val = int(float(val))
                self.corrXmaxVal = val
                
        def ymin(self, val):
                val = int(float(val))
                self.corrYminVal = val

        def ymax(self, val):
                val = int(float(val))
                self.corrYmaxVal = val

	def makeDistro(self):
   	        P.close()
#		P.ion()
		fig = P.figure(figsize = (18, 10))
#		print self.mutationOccurence
#		print np.sum(self.mutationOccurence.values())
		xvals = self.mutationOccurence.keys()
		yvals = self.mutationOccurence.values()
		P.plot(xvals, yvals, 'xb')
#		P.xlim(0, 210)
		P.xlabel('Number of mutations per sequence')
		P.ylabel('Frequency of mutations')
		self.makeCSVs()
		P.show()
		
	def makeCorre(self):
	        P.close()
		mode = self.splittage.get()
		annotate = self.annotate.get()
		allosteryResidues = self.fullAllosteryResidues
		allosteryValues = self.fullAllosteryValues
		mutationalityList = self.fullMutationalityList
		
		#### If we only want to output the plotted values, need to sort that out here
		
		if mode == 1:
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
#			P.ion()
			fig = P.figure(figsize = (18, 10))
			gs = gridspec.GridSpec(1, 2)
			ax1 = P.subplot(gs[0])
			ax2 = P.subplot(gs[1])
	#		print mutationalityList
	#		print mutationalityList
			ax1.plot(allosteryNegative, mutNegative, 'xb')
			ax1.set_ylabel('Mutation Frequency')
			ax1.set_xlabel('$\mathrm{\Delta}$(T$\mathrm{\Delta}$S) [Negative Values]')
			ax1.set_ylim(-0.001, 0.75)
			ax1.set_xlim(0, 0.35)

			if annotate == 1:
         			for i in range(len(allosteryNegative)):
         			    if mutNegative != 0:
                 			    P.annotate(i, xy=(allosteryNegative[i], mutNegative[i]), xytext=(allosteryNegative[i], mutNegative[i]))
			ax2.plot(allosteryPositive, mutPositive, 'xr')
			ax2.set_ylabel('Mutation Frequency')
			ax2.set_xlabel('$\mathrm{\Delta}$(T$\mathrm{\Delta}$S) [Positive Values]')
			if annotate == 1:
         			for i in range(len(allosteryPositive)):
         			    if mutPositive != 0:
                 			    P.annotate(i, xy=(allosteryPositive[i], mutPositive[i]), xytext=(allosteryPositive[i], mutPositive[i]))
			ax2.set_ylim(-0.001, 0.75)
			P.show()

		else:
	#		print allosteryResidues, allosteryValues
			for element in range(len(allosteryValues)):
				if allosteryValues[element] <= 0:
					allosteryValues[element] *= -1

			averageAllostery = np.average(allosteryValues)
			averageMutation = np.average(mutationalityList)
		
#			P.ion()
			fig = P.figure(figsize = (18, 10))
#                        print "\n\n", len(allosteryValues), len(mutationalityList), "\n\n"
			P.plot(allosteryValues, mutationalityList, 'xb')
			P.fill([0, averageAllostery, averageAllostery, 0], [0, 0, averageMutation, averageMutation], 'r', alpha=0.5) 
			P.ylabel('Mutation Frequency')
			P.xlabel('$\mathrm{\Delta}$(T$\mathrm{\Delta}$S)')
			if annotate == 1:
         			for i in range(len(allosteryValues)):
         			    P.annotate(i, xy=(allosteryValues[i], mutationalityList[i]), xytext=(allosteryValues[i], mutationalityList[i]))
         		if self.corrYminVal > self.corrYmaxVal:
         		         print "Lower bound for Y axis set higher than upper bound, ignoring upper bound"
         		         self.corrYmaxVal = 100
         		if self.corrXminVal > self.corrXmaxVal:
         		         print "Lower bound for X axis set higher than upper bound, ignoring upper bound"
         		         self.corrXmaxVal = 100			
			P.ylim(0.01*self.corrYminVal*np.max(mutationalityList), 0.01*self.corrYmaxVal*np.max(mutationalityList))
			P.xlim(0.01*self.corrXminVal*np.max(allosteryValues), 0.01*self.corrXmaxVal*np.max(allosteryValues))
			P.show()
#			P.draw()
#			P.ioff()
			P.savefig('./CorrelationPlot.png')
			
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
			
			self.makeCSVs()

			#print len(allosteryValues)
			#print len(mutationalityList)

        def makeCSVs(self):
            	csvFile = open('./CorrelationDataFull.csv', 'w+')
		dataWriter = csv.writer(csvFile)
		dataWriter.writerow(['Residue Number', 'Allostery', 'Mutation Factor'])
#		print len(self.fullAllosteryValues)
#		print len(self.fullMutationalityList)
		for i in range(len(self.fullAllosteryValues)):
			dataWriter.writerow([i+1, self.fullAllosteryValues[i], self.fullMutationalityList[i]])
				
		csvFile = open('./CorrelationDataCut.csv', 'w+')
		dataWriter = csv.writer(csvFile)
		dataWriter.writerow(['Residue Number', 'Allostery', 'Mutation Factor'])
		for i in range(len(self.fullAllosteryValues)):
			dataWriter.writerow([i+1, self.fullAllosteryValues[i], self.mutationalityListCut[i]])
	
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
#		mode = self.loopOverlay.get()
#		if mode == 1:
#			mode = "loop"
#		else:
#			mode = "noLoop"
                self.makeCSVs()
		self.analyse.allostery(barlinewidth=0, mode='noLoop')

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
#			paramsListElement, identity, loopIndexElement, alpha, beta, custom = run_clustal.sequences(element, self.groupUserDefined).runAnalysis() # CAP Stuff
			paramsListElement, identity = run_clustal.sequences(element, self.groupUserDefined).runAnalysis()
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
#		self.loopIndexElement = loopIndexElement
#		self.alpha = alpha
#		self.beta = beta
		self.totalSequences = totalSequences
		self.mutationOccurence = mutationOccurence
		self.mutationDistribution = mutationDistribution

		cutOff = 60
		newParams, mutationalityList = run_clustal.sortParams(self.paramsList, self.identityList, cutOff, 1, self.totalSequences)
		try:
			delete(analyse)
		except:
			pass
#		self.analyse = run_clustal.analysis(newParams, 1, mutationalityList, loopIndexElement, alpha, beta, self.groupUserDefined)   # CAP stuff
                self.analyse = run_clustal.analysis(newParams, 1, mutationalityList, None, None, None, None)
		self.onInfo("Sequences Loaded.")
		self.runButton.destroy()
		self.initUI()
#		CLI(paramsList, identityList, loopIndexElement, alpha, beta, totalSequences)

	def onUpdate(self, cutOff):
		newParams, mutationalityList = run_clustal.sortParams(self.paramsList, self.identityList, cutOff, 1, self.totalSequences)
		try:
			delete(analyse)
		except:
			pass
		self.analyse = run_clustal.analysis(newParams, 1, mutationalityList, None, None, None, None)
		self.mutationalityListCut = mutationalityList
		print len(newParams), "sequences to be included into the plot."

def GUI():
	root = Tk()
	app = mainScreen(root)
	root.mainloop()

if __name__ == "__main__":
	print "\n---------------===============***************===============---------------"
	print "\nMatty Jones' Protein Sequence Alignment Analysis Tool"
	print "\n---------------===============***************===============---------------\n"
	GUI()














