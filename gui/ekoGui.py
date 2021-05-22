from Tkinter import *
import ttk
import tkFileDialog
import tkMessageBox
import os
import re
import subprocess
import shutil
import ntpath
from sepFiles import createConformers
from genhtml  import genHtmlPage 

def importDirFrame():
   dirFrame = ttk.Frame(root, padding = "3 3 12 12")
   dirFrame.grid(column=0, row=0, sticky=(N, W, E, S))
   dirFrame.columnconfigure(0, weight=1)
   dirFrame.rowconfigure(2, weight=1)
   ttk.Label(dirFrame, text="Select a directory where the output files(*_lt3.pdb) for QMD are present").grid(column=2, row=1, sticky=(W,E), rowspan=2)
   ttk.Entry(dirFrame, width=7, textvariable=dirName).grid(column=2, row=3, sticky=(W,E))
   ttk.Button(dirFrame, text="SelectDirectory", command= lambda: selDirectory(dirFrame)).grid(column=6, row=3, sticky=W)
   nxtBtn = ttk.Button(dirFrame, text="Next",  command= lambda: loadFrame('m', 1))
   nxtBtn.grid(column=2, row=10, sticky=W)
   dirName.set("")
   for child in dirFrame.winfo_children(): child.grid_configure(padx=5, pady=5)


def importMenuFrame():
   menuFrame = ttk.Frame(root, padding = "3 3 12 12")
   menuFrame.grid(column=0, row=0, sticky=(N, W, E, S))
   menuFrame.columnconfigure(0, weight=1)
   menuFrame.rowconfigure(2, weight=1)
   
   ttk.Label(menuFrame, text="Select Target Protein",  font = "Helvetica 16 bold", foreground="blue").grid(column=1, row=1, sticky=W)
   pdbLbl = ttk.Label(menuFrame, text='<No PDB>')
   pdbLbl.grid(column=2, row=1, sticky=W, columnspan=2)
   chain1Lbl = ttk.Label(menuFrame, text='Chain1')
   chain2Lbl = ttk.Label(menuFrame, text='Chain2')
   chain1Drop = ttk.Combobox(menuFrame, width=2)
   chain2Drop = ttk.Combobox(menuFrame, width=2)
   chain1Drop.bind('<<ComboboxSelected>>', chain1Sel)
   chain2Drop.bind('<<ComboboxSelected>>', chain2Sel)
   chain1Nick.set("")
   chain2Nick.set("")
   chain1Lbl.grid(column=2, row=2, sticky=W, padx=5, pady=5)
   chain1Drop.grid(column=3, row=2, sticky=W, padx=5, pady=5)
   ttk.Label(menuFrame, text="Chain1 Name").grid(column=4, row=2, sticky=(W,E), padx=5, pady=5)
   ttk.Entry(menuFrame,  width=7, textvariable=chain1Nick).grid(column=5, row=2, sticky=W, padx=5, pady=5)

   chain2Lbl.grid(column=2, row=3, sticky=W, padx=5, pady=5)
   chain2Drop.grid(column=3, row=3, sticky=W, padx=5, pady=5)
   ttk.Label(menuFrame, text="Chain2 Name").grid(column=4, row=3, sticky=(W,E), padx=5, pady=5)
   ttk.Entry(menuFrame,  width=7, textvariable=chain2Nick).grid(column=5, row=3, sticky=W, padx=5, pady=5)

   ttk.Button(menuFrame, text="Browse", command= lambda: getPdbFile(pdbLbl, chain1Drop, chain2Drop )).grid(column=4, row=1, sticky=W)
   ttk.Label(menuFrame, text="OR", font = "Times 12 bold").grid(column=5, row=1, sticky=E)
   ttk.Label(menuFrame, text="Enter PDB ID:").grid(column=6, row=1, sticky=E)
   ttk.Entry(menuFrame, width=7, textvariable=pdbID).grid(column=7, row=1, sticky=W)
   ttk.Button(menuFrame, text="Download", command= lambda: downloadPdb(pdbLbl, chain1Drop, chain2Drop)).grid(column=8, row=1, sticky=W)
   
   ttk.Label(menuFrame, text="Select Conformer(s)", font = "Helvetica 16 bold", foreground="blue").grid(column=1, row=4, sticky=W, columnspan=2) 
   var1 = StringVar()
   conformerList = getItems()
   selConformers = []
   conformerAtoms = []
   
   validCmd = menuFrame.register(isCommaSeparated)
   ttk.Label(menuFrame, text="Comformer", font = "Times 12 bold").grid(column=1, row=5, sticky=(W,E), columnspan=3)
   ttk.Label(menuFrame, text="atoms4rms (Comma Separated)", font = "Times 12 bold").grid(column=4, row=5, sticky=W, columnspan=20) 
   r = 6
   for conformer in conformerList:
      chkd = IntVar()
      c = ttk.Checkbutton(menuFrame, text=conformer, variable=chkd).grid(column=1, row=r, sticky=(W,E), columnspan=3)
      selConformers.append(chkd)
      atoms = StringVar()
      inFilePath = dirName.get() + "/" + conformer.replace("_lt3.pdb", ".in")
      if os.path.isfile(inFilePath) and os.access(inFilePath, os.R_OK):
         inFile = open(inFilePath, "r")
         pat = r'atoms4rms="(.*)"'
	 for line in inFile:
            if re.match(pat, line):
               atoms.set(re.match(pat, line).group(1).replace(" ",","))
               break
         inFile.close()

      ttk.Entry(menuFrame, width=15, textvariable=atoms, validate="key", validatecommand=(validCmd, '%P')).grid(column=4, row=r,sticky=W,  columnspan=20)
      conformerAtoms.append(atoms)
      r = r+1
   
   r = r+2
   ttk.Label(menuFrame, text="Parameters", font = "Helvetica 16 bold", foreground="blue").grid(column=1, row=r, sticky=W, columnspan=2)
   
   r = r+1
   vcmd = menuFrame.register(isDouble)
   ttk.Label(menuFrame, text="Distance").grid(column=1, row=r, sticky=W, columnspan=2)
   ttk.Entry(menuFrame, width=2, textvariable=distVal, validate="key", validatecommand=(vcmd, '%P')).grid(column=3, row=r, sticky=(W,E)) 
   
   r = r+1
   ttk.Label(menuFrame, text="Angle").grid(column=1, row=r, sticky=W, columnspan=2)
   ttk.Entry(menuFrame, width=2, textvariable=anglVal, validate="key", validatecommand=(vcmd, '%P')).grid(column=3, row=r, sticky=(W,E))
   
   r = r+1
   ttk.Label(menuFrame, text="Filter").grid(column=1, row=r, sticky=W, columnspan=2)
   ttk.Entry(menuFrame, width=2, textvariable=fltrVal, validate="key", validatecommand=(vcmd, '%P')).grid(column=3, row=r, sticky=(W,E))
   
   r = r+1
   ttk.Label(menuFrame, text="Pointing").grid(column=1, row=r, sticky=W, columnspan=2)
   ttk.Entry(menuFrame, width=2, textvariable=pntVal, validate="key", validatecommand=(vcmd, '%P')).grid(column=3, row=r, sticky=(W,E))

   r = r+1
   ttk.Label(menuFrame, text="Contact").grid(column=1, row=r, sticky=W, columnspan=2)
   ttk.Entry(menuFrame, width=2, textvariable=contVal, validate="key", validatecommand=(vcmd, '%P')).grid(column=3, row=r, sticky=(W,E))
   
   r = r+3
   ttk.Button(menuFrame, text="Previous", command= lambda: loadFrame('d')).grid(column=1, row=r, sticky=E)
   ttk.Button(menuFrame, text="Run", command= lambda: runEkoMatching(conformerList, selConformers, conformerAtoms)).grid(column=6, row=r, sticky=W)
   
   for child in menuFrame.winfo_children(): child.grid_configure(padx=5, pady=5)


def downloadPdb( pdbLbl,chain1Drop, chain2Drop):
   if pdbID.get() == "":
      return
   cwd = os.getcwd()
   os.chdir(dirName.get())
   cmd  = "wget http://www.rcsb.org/pdb/files/" + pdbID.get() + ".pdb"
   process = subprocess.Popen(cmd, shell=True)
   process.wait()
   
   global pdbFilePath
   pdbFilePath = os.path.join(dirName.get(), pdbID.get() + ".pdb")
   if not os.path.exists(pdbFilePath):
      tkMessageBox.showwarning(title="Warning", message="No pdb was found. Please make sure that the pdb Id is correct.")
   else:
      pdbLbl.config(text = pdbFilePath) 
      chains = getChains(pdbFilePath)
      chain1Drop['values'] = chains
      chain2Drop['values'] = chains 
   os.chdir(cwd) 

def getPdbFile(pdbLbl,chain1Drop, chain2Drop):
    filename = tkFileDialog.askopenfilename(initialdir=dirName.get(), title= "Select PDB File",filetypes = (("PDB Files","*.pdb"),("all files","*.*")))
    global pdbFilePath
    global chain1, chain2
    if filename:
       pdbLbl.config(text=filename)
       pdbFilePath =  filename
       chains = getChains(pdbFilePath)
       #print chains
       #TODO check if there are no chains
       chain1Drop['values'] = chains
       chain2Drop['values'] = chains
       #chain1Drop.current(0)
       #chain2Drop.current(0)
       #chain1=chain1Drop.get()
       #chain2=chain2Drop.get()

def getChains(pdbFilePath):
   f = open(pdbFilePath, "r")
   chain = set()
   for line in f:
      if line.startswith("ATOM"):
         chain.add(filter(None,re.split(' |\t', line))[4])
   f.close()
   return list(chain)

def chain1Sel(event):
   global chain1
   chain1 = event.widget.get()
   #print 'Chain1:' + chain1

def chain2Sel(event):
   global chain2
   chain2 = event.widget.get()
   #print chain2

def isCommaSeparated(P):
   #print P
   for c in P:
      #print c
      if not ((c>='0' and c<='9') or (c==',')):
         return False
   return True

def isDouble(P):
   try:
      float(P)
      return True
   except ValueError:
      return False

def getItems():
   lst = []
   for lt3File in os.listdir(dirName.get()):
      if lt3File.endswith('_lt3.pdb'):
         lst.append(lt3File)
   return lst

def insertItems(itemsListBox):
   for lt3File in os.listdir(dirName.get()):
      if lt3File.endswith('_lt3.pdb'):
         itemsListBox.insert(END, lt3File)
   
def selDirectory(dirFrame):
   options = {}
   options['initialdir'] = '.'
   options['mustexist'] = True
   options['parent'] = dirFrame
   options['title'] = 'Choose a directory'
   dirName.set(tkFileDialog.askdirectory(**options))

def loadFrame(opt, nxtFlg = 0):
   if nxtFlg ==1 :
      if (not (os.path.isdir(dirName.get()) and os.access(dirName.get(), os.R_OK))) or (not any(fname.endswith('_lt3.pdb') for fname in os.listdir(dirName.get()))):
          tkMessageBox.showwarning(title="Warning", message="Please make sure that the directory exists and it has files ending in _lt3.pdb") 
          return
   deleteAllFrames()
   if opt=='d':
      importDirFrame()
   elif opt=='m':
      importMenuFrame()

def deleteAllFrames():
   for child in root.winfo_children():
      child.destroy()

def findNumConf(filePath):
   pattern = r'ENDMDL.*'
   f =  open(filePath, "r")
   count  = 0
   for line in f:
      if re.search(pattern,line):
         count = count + 1
   f.close()
   return count 

def validateBeforeMatching(conformersList, selConformers, conformerAtoms):
   if pdbFilePath == "":
      tkMessageBox.showwarning(title="Warning", message="Please select a target protein or download one.")
      return False

   chains = getChains(pdbFilePath)
   if (chain1 not in chains) or (chain2 not in chains):
      tkMessageBox.showwarning(title="Warning", message="Please select a valid chain.")
      return False
   
   selectedConformers = []
   selectedConfAtoms = []
   for isConfChkd in selConformers:
      if isConfChkd.get() ==1 :
         idx = selConformers.index(isConfChkd)
         selectedConformers.append(conformersList[idx])
         selectedConfAtoms.append(conformerAtoms[idx].get())
         if conformersList[idx] == '' or conformerAtoms[idx].get() == '':
            tkMessageBox.showwarning(title="Warning", message="Please select atoms for the selected conformer.")
            return False
   if len(selectedConformers) == 0 or len(selectedConfAtoms) == 0:
      tkMessageBox.showwarning(title="Warning", message="Please select at least one conformer.")
      return False

   return True

def load3D(progressLbl, chainDict):
   global browserFilesCreated
   cwd = os.getcwd()
   os.chdir(dirName.get())
   conformersPath = os.path.join(dirName.get(), "conformers")
   if not browserFilesCreated:
      progressLbl.configure(text="Generating all Conformers...")
      root.update()
      
      dirPath = os.path.join(dirName.get(), "conformers")
      if os.path.isdir(conformersPath):
         decision = tkMessageBox.askyesno(title="Overwrite Directory", message="The directory " + conformersPath + " already exists. Do you want to overwrite it?")
         root.update()
         if decision:
            shutil.rmtree(conformersPath)
         else:
            root.destroy()
            return

      createConformers()
      progressLbl.configure(text="Configuring html page for 3D view")
      root.update()
      genHtmlPage(os.environ.get('EKOGUIHOME'), chainDict)
      browserFilesCreated =  True
   
   os.chdir(os.environ.get('EKOGUIHOME'))
   progressLbl.configure(text="Opening Browser")
   root.update()
   cmd = 'chromium-browser --args --allow-file-access-from-files ' + os.path.join(conformersPath, "viewconformers.html") + " &"
   process = subprocess.Popen(cmd, shell=True)
   progressLbl.configure(text="Browser Opened")
   root.update()
   os.chdir(cwd)

def runEkoMatching(conformersList, selConformers, conformerAtoms):
   isValid = validateBeforeMatching(conformersList, selConformers, conformerAtoms)
   if not isValid:
       return 

   deleteAllFrames()
   
   resultFrame = ttk.Frame(root, padding = "3 3 12 12")
   resultFrame.grid(column=0, row=0, sticky=(N, W, E, S))
   resultFrame.columnconfigure(0, weight=1)
   resultFrame.rowconfigure(2, weight=1)

   resultsWind = Text(resultFrame, height=15)
   resultsWind.grid(column=1, row=1, sticky=E )
   progressLbl = ttk.Label(resultFrame, text="Running...", font = "Times 16 bold", foreground="blue")
   root.update()
   progressLbl.grid(column=1, row=26,columnspan=2)   
   

   selectedConformers = []
   selectedConfAtoms = []
   for isConfChkd in selConformers:
      if isConfChkd.get() ==1 :
         idx = selConformers.index(isConfChkd)
         selectedConformers.append(conformersList[idx])
         selectedConfAtoms.append(conformerAtoms[idx].get())
   global pdbFilePath 
   #PPIHOME="/apps/burgess/bin"
   for i in range(0, len(selectedConformers)):
      filename = selectedConformers[i]
      #print filename
      resultsWind.insert(END, filename+"\n")
      filePath = dirName.get() + '/' + filename
      outPath = dirName.get() + '/' + filename.replace(".pdb", "_overlay.out")
      numConformations = findNumConf(filePath)
      resultsWind.insert(END, str(numConformations)+"\n")
      root.update()
      if not pdbFilePath.startswith(dirName.get()): 
         shutil.copy(pdbFilePath, dirName.get())
      pdbFilePath = ntpath.basename(pdbFilePath)
      cwd = os.getcwd()
      os.chdir(dirName.get())
      cmdPrefix = os.environ.get('PPIHOME') + "/triplet_search --Wdist=" + str(distVal.get()) + " --Wang=" + str(anglVal.get()) + " --filter=" + str(fltrVal.get()) +  " --contact=" + str(contVal.get()) + " --pointing=" + str(pntVal.get()) + " " + filename  + " " + selectedConfAtoms[i] + " " + pdbFilePath +" "
      cmd1 = cmdPrefix + chain1 + " " + chain2 + " > " + outPath
      cmd2 = cmdPrefix + chain2 + " " + chain1 + " >> " + outPath
     
      resultsWind.insert(END, cmd1+"\n")
      root.update()
      process = subprocess.Popen(cmd1, shell=True)
      process.wait()
      resultsWind.insert(END, cmd2+"\n")
      root.update()
      process = subprocess.Popen(cmd2, shell=True)
      process.wait()
      os.chdir(cwd)
      createSummaryFile()
      progressLbl.configure(text="Completed")
      root.update()

      chainDict = {}
      if chain1Nick.get().strip(' \t\n\r') != "":
         chainDict[chain1] = chain1Nick.get().strip(' \t\n\r')
      if chain2Nick.get().strip(' \t\n\r') != "":
         chainDict[chain2] = chain2Nick.get().strip(' \t\n\r')

      ttk.Button(resultFrame, text="Show 3D Conformers", command= lambda: load3D(progressLbl, chainDict)).grid(column=1, row=28)


def createSummaryFile():
   summaryFile = dirName.get() + "/score_summary.txt"   
   f = open(summaryFile, "w")
   f.write("#line frame energy chain   resids     score     rms     angle   pdbname mimic_name\n")
   f.write("# 1     2     3       4  5   6   7      8        9        10       11      12\n")
   mimicFiles = []
   for mimicFile in os.listdir(dirName.get()):
      if mimicFile.endswith("_lt3_overlay.out"):
         mimicFiles.append(mimicFile)
   
   #resultsWind.insert(END, "Summary\n")
   count  = 1
   for mimic in mimicFiles:
     filePath = dirName.get() + "/" +mimic
     mimicName = mimic.replace("_lt3_overlay.out", "") 
     #resultsWind.insert(END, mimicName+"\n")
     pattern  = r'.*Energy.*'
     mimicFile = open(filePath, "r")
     line  = mimicFile.readline()
     while line:
        if re.search(pattern, line):
          origFrame = filter(None,re.split(' |\t', line))[2]
          energy = filter(None,re.split(' |\t', line))[4]
          while not (re.search(r'REMARK.*EKO.*', line)):
             line = mimicFile.readline()
          if re.search(r'REMARK.*EKO.*', line):
             params = filter(None,re.split(' |\t', line))
             entry = str(count).ljust(5)+ " " + origFrame + " " + energy + " " + params[5] + " " + params[9] + " " +  params[10] + " " + params[11] + " "+ params[17].strip()+ " " + params[13] + " " + params[15] + " " + params[3] + " " + mimicName + "\n"
             f.write(entry)
             #resultsWind.insert(END, entry)
             count =  count + 1 
        
        line =  mimicFile.readline()
     
     mimicFile.close()
     
   f.close()  


root = Tk()
root.title("EKO GUI")
root.geometry("900x600+30+30")
dirName = StringVar()
distVal = DoubleVar()
distVal.set(0.1)
anglVal = DoubleVar()
anglVal.set(1.0)
fltrVal = DoubleVar()
fltrVal.set(3.0)
pntVal = DoubleVar()
pntVal.set(90.0)
contVal = DoubleVar()
contVal.set(4.0)
chain1 = ''
chain2 = ''
chain1Nick =  StringVar()
chain2Nick =  StringVar()
pdbFilePath = ''
pdbID = StringVar()
browserFilesCreated  = False
if os.environ.get('PPIHOME') is None:
   os.environ['PPIHOME'] = "/apps/burgess/bin"
if os.environ.get('EKOGUIHOME') is None:
   os.environ['EKOGUIHOME'] = "/apps/burgess/eko_qmd_gui"

importDirFrame()
root.mainloop()

