import ttk
import Tkinter
import run_clustal
import subprocess as sub
import tkMessageBox as box
import os
from PIL import Image, ImageTk

class mainScreen(ttk.Frame):
	def __init__(self, parent):
		ttk.Frame.__init__(self, parent)
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
		win = Tkinter.Toplevel(self.parent)
		win.title("Input Sequence")
		self.pack(fill=Tkinter.BOTH, expand=1)
		ttk.Style().configure("TFrame", background="#FFFFFF")
#		win.geometry('%dx%d+%d+%d' % (self.width, self.height, self.x, self.y))
#		frame = Tkinter.Frame(self, relief=Tkinter.RAISED, borderwidth=1)
#		frame.pack(fill=Tkinter.BOTH, expand=1)

		self.columnconfigure(1, weight=1)
		self.columnconfigure(3, pad=5)
		self.rowconfigure(3, weight=1)
		self.rowconfigure(5, pad=5)

		lbl = Tkinter.Label(win, text="Please input Protein sequence in FASTA format")
		lbl.grid(row=0, column=0, pady=5, padx=5, sticky=Tkinter.W)

		area = Tkinter.Text(win)
		area.grid(row=1, column=0, columnspan=2, rowspan=4, padx=5, sticky=Tkinter.E+Tkinter.S)

		cbtn = ttk.Button(win, text="Cancel", command=win.destroy)
		cbtn.grid(row=2, column=3, padx=5, pady=4, sticky=Tkinter.E)
		
		hbtn = ttk.Button(win, text="Help")
#		hbtn.grid(row=5, column=1)
		hbtn.grid(row=5, column=0, padx=5, pady=4, sticky=Tkinter.W)

		obtn = ttk.Button(win, text="Input", command= lambda:self.inputSeq(area))
		obtn.grid(row=5, column=3, padx=5, pady=4)

		self.centreWindow(win)

	def inputSeq(self, textBox):
		contents = textBox.get(1.0, Tkinter.END)
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
		textBox.delete(1.0, Tkinter.END)

	def initUI(self):
		self.parent.title("Main Windows")
		self.pack(fill=Tkinter.BOTH, expand=1)
		ttk.Style().configure("TFrame", background="#FFFFFF")

		self.columnconfigure(0, weight=1)

		self.rowconfigure(0, weight=1)	# Image
		self.rowconfigure(1, weight=1)	# Run
		self.rowconfigure(2, weight=1)	# Input
		self.rowconfigure(3, weight=1)	# Blank
		self.rowconfigure(4, weight=1)	# Settings
		self.rowconfigure(5, weight=1)	# Quit

		BSILogoOpen = Image.open("./UIFiles/BSILogo.jpg")
		BSILogo = ImageTk.PhotoImage(BSILogoOpen)
		label1 = Tkinter.Label(self, image=BSILogo)
		label1.image = BSILogo
		label1.grid(row=0, column=0, sticky=Tkinter.N+Tkinter.E+Tkinter.W, pady=5, padx=5)

		runButton = ttk.Button(self, text="Load Sequences", command=self.test)
		runButton.grid(row=1, column=0, padx=5, pady=4)

		inputButton = ttk.Button(self, text="Input New Sequences", command=self.inputScreen)
		inputButton.grid(row=2, column=0, padx=5, pady=4)
		
		settingsButton = ttk.Button(self, text="Settings")#, command=self.settingsScreen)
		settingsButton.grid(row=4, column=0, padx=5, pady=4)
		
		quitButton = ttk.Button(self, text="Exit", command=self.quit)
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


	def test(self):
		orderedList, totalSequences = run_clustal.clustalw(directory='./Sequences').runSetup()
		paramsList = []
		identityList = []
		for element in orderedList:
			paramsListElement, identity, loopIndexElement, alpha, beta = run_clustal.sequences(element).runAnalysis()
			paramsList.append(paramsListElement)
			identityList.append(identity)
		self.onInfo("Sequences Loaded.")
#		CLI(paramsList, identityList, loopIndexElement, alpha, beta, totalSequences)

def GUI():
	root = Tkinter.Tk()
	app = mainScreen(root)
	root.mainloop()

if __name__ == "__main__":
	print "\n---------------===============***************===============---------------"
	print "\nInitialising GUI"
	print "\n---------------===============***************===============---------------\n"
	GUI()














