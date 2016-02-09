#import ttk
#import Tkinter
import run_clustal
import subprocess as sub
import tkMessageBox as box
import os
from ttk import *
from Tkinter import *
from PIL import Image, ImageTk

class mainScreen(Frame):
	def __init__(self, parent):
		Frame.__init__(self, parent)
		self.parent = parent
#		self.centreWindow(parent)
		self.initUI()
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
		win = Toplevel(self.parent)
		win.title("Input Sequence")
		self.pack(fill=BOTH, expand=1)
		Style().configure("TFrame")#, background="#FFFFFF")
#		win.geometry('%dx%d+%d+%d' % (self.width, self.height, self.x, self.y))
		frame = Frame(win, relief=RAISED, borderwidth=1)
		frame.pack(fill=BOTH)

		lbl = Label(win, text="Please input Protein sequence in FASTA format")
		lbl.pack()

		area = Text(win)
		area.pack(fill=BOTH, expand=1, padx=5, pady=5)

		obtn = Button(win, text="Input", command= lambda:self.inputSeq(area))
		obtn.pack(side=RIGHT, padx=5, pady=5)

		cbtn = Button(win, text="Exit", command=win.destroy)
		cbtn.pack(side=LEFT, padx=5, pady=5)
		
		hbtn = Button(win, text="Help")
#		hbtn.grid(row=5, column=1)
		hbtn.pack(side=RIGHT, padx=5, pady=5)

		self.centreWindow(win)

	def inputSeq(self, textBox):
		contents = textBox.get(1.0, END)
		ecoliFile = open('./ecoli.txt', 'r')
		ecoli = ecoliFile.read()
		ecoliFile.close()
		sequencePipes = run_clustal.findIndex(contents, '|')
		name = contents[sequencePipes[-2]+1:sequencePipes[-1]]
		testFile = open('./test.txt', 'w+')
		testFile.write(ecoli)
		testFile.write('\n')
		testFile.write(contents)
		testFile.close()
		os.rename('./test.txt', './Sequences/'+str(name)+'.txt')
		self.onInfo('Sequence Inputted into System')
		textBox.delete(1.0, END)

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
		
		settingsButton = Button(self, text="Settings")#, command=self.settingsScreen)
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

		cutOffScale = Scale(win2, from_=40, to=100, length=400, orient=HORIZONTAL, command=self.onScale)
		cutOffScale.set(40)
		cutOffScale.grid(row=1, sticky=W)

#		cutOffOutput = Label(win2, text=0, textvariable=self.cutOff, width=3)
#		cutOffOutput.grid(row=1, column=2)

		plotButton = Button(win2, text="Plot residue mutation and allostery with loop overlay", command=self.makePlot)
		plotButton.grid(row=2, padx=5, pady=15)

		quitButton = Button(win2, text="Exit", command=win2.destroy)
		quitButton.grid(row=3, column=1, padx=5, pady=5, sticky=SE)

		self.centreWindow(win2)

	def makePlot(self):
		self.analyse.allostery(mode='loop', barlinewidth=0)
		self.win3 = Toplevel(self.parent)
		self.win3.title("PUT CUTOFF HERE")
		Style().configure("TFrame")
		screenWidth = self.parent.winfo_screenwidth()
		screenHeight = self.parent.winfo_screenheight()
		figWidth = screenWidth*0.7
		figHeight = screenHeight*0.7
#		self.x = (screenWidth - figWidth)/2
#		self.y = (screenHeight - figHeight)/2
##		print self.parent.winfo_geometry()
#		win3.geometry('%dx%d+%d+%d' % (figWidth, figHeight, self.x, self.y))
#		graphOpen = Image.open("./temp.png")
#		graph = ImageTk.PhotoImage(graphOpen)
		Grid.columnconfigure(self.win3, 0, weight=1)
		Grid.rowconfigure(self.win3, 0, weight=1)	
#		canvas = Canvas(win3, width = figWidth, height = figHeight)
#		canvas.pack(expand = YES, fill = BOTH)
		self.updateImage(figWidth, figHeight)

	def reCalculate(self, figWidth, figHeight):
		self.analyse.allostery(mode='loop', barlinewidth=0)
		self.label1.destroy()
#		graph.destroy()
		self.updateImage(figWidth, figHeight)
		
	def updateImage(self, figWidth, figHeight):
#		self.analyse.allostery(mode='loop', barlinewidth=0)
		graphOpen = Image.open("./temp.png")
		graphOpen = graphOpen.resize((int(figWidth), int(figHeight)), Image.ANTIALIAS)
		global graph
		graph = ImageTk.PhotoImage(graphOpen)
#		graphOpen.close()
#		canvas.create_image(50, 10, image = graph, anchor = NW)
		self.label1 = Label(self.win3, image=graph)
		self.label1.image = graph
		self.label1.grid(sticky=N+E+S+W, pady=5, padx=5)
		self.win3.after(3000, lambda: self.reCalculate(figWidth, figHeight))

# Every TIME_TO_PLOT seconds check if bar has changed, if has then replot, otherwise sleep.



	def load(self):
		orderedList, totalSequences = run_clustal.clustalw(directory='./Sequences').runSetup()
		paramsList = []
		identityList = []
		for element in orderedList:
			paramsListElement, identity, loopIndexElement, alpha, beta = run_clustal.sequences(element).runAnalysis()
			paramsList.append(paramsListElement)
			identityList.append(identity)
		self.paramsList = paramsList
		self.identityList = identityList
		self.loopIndexElement = loopIndexElement
		self.alpha = alpha
		self.beta = beta
		self.totalSequences = totalSequences
		cutOff = 40
		newParams, mutationalityList, loopIndexElement, alpha, beta = run_clustal.sortParams(self.paramsList, self.identityList, self.loopIndexElement, self.alpha, self.beta, cutOff, 1, self.totalSequences)
		try:
			delete(analyse)
		except:
			pass
		self.analyse = run_clustal.analysis(newParams, 1, mutationalityList, loopIndexElement, alpha, beta)
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
		self.analyse = run_clustal.analysis(newParams, 1, mutationalityList, loopIndexElement, alpha, beta)
		print len(newParams)

def GUI():
	root = Tk()
	app = mainScreen(root)
	root.mainloop()

if __name__ == "__main__":
	print "\n---------------===============***************===============---------------"
	print "\nInitialising GUI"
	print "\n---------------===============***************===============---------------\n"
	GUI()














