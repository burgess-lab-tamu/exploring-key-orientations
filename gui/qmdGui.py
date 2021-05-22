from Tkinter import *
import ttk
import tkFileDialog
import tkMessageBox
#from tkinter.filedialog
import os
import re
import subprocess
import shutil
import ntpath
import time

def importDirFrame():
   dirFrame = ttk.Frame(root, padding = "3 3 12 12")
   dirFrame.grid(column=0, row=0, sticky=(N, W, E, S))
   dirFrame.columnconfigure(0, weight=1)
   dirFrame.rowconfigure(2, weight=1)
   ttk.Label(dirFrame, text="Select a directory where you want to run the QMD").grid(column=2, row=1, sticky=(W,E), rowspan=2)
   ttk.Entry(dirFrame, width=7, textvariable=dirName).grid(column=2, row=3, sticky=(W,E))
   ttk.Button(dirFrame, text="SelectDirectory", command= lambda: selDirectory(dirFrame)).grid(column=6, row=3, sticky=W)
   ttk.Button(dirFrame, text="Next",  command= lambda: loadFrame('e', 1)).grid(column=2, row=10, sticky=W)
   dirName.set("")
   for child in dirFrame.winfo_children(): child.grid_configure(padx=5, pady=5)


def importEntryFrame():
   entryFrame = ttk.Frame(root, padding = "3 3 12 12")
   entryFrame.grid(column=0, row=0, sticky=(N, W, E, S))
   entryFrame.columnconfigure(0, weight=5)
   entryFrame.rowconfigure(0, weight=5)
    
   
   validCmd = entryFrame.register(isValidId)
   ttk.Label(entryFrame, text="Mimic Id:").grid(column=2, row=1, sticky=E)
   ttk.Entry(entryFrame, width=15, textvariable=mimicId, validate="key", validatecommand=(validCmd, '%P')).grid(column=5, row=1, sticky=W)

   ttk.Label(entryFrame, text="Smiles String:").grid(column=2, row=3, sticky=E)
   ttk.Entry(entryFrame, width=30, textvariable=smilesString).grid(column=5, row=3, sticky=W, columnspan=2)

   vcmd = entryFrame.register(isCommaSeparated)
   ttk.Label(entryFrame, text="Ca Cb Serial Number:").grid(column=2, row=4, sticky=E)
   ttk.Entry(entryFrame, width=15, textvariable=cacbSerial, validate="key", validatecommand=(vcmd, '%P')).grid(column=5, row=4, sticky=W)
   ttk.Label(entryFrame, text="OR", font = "Times 12 bold").grid(column=6, row=4, sticky=W)
   ttk.Button(entryFrame, text="Select Atoms", command= lambda: openPage()).grid(column=8, row=4, sticky=W)
   
   processingLabel.set("") 
   ttk.Label(entryFrame, textvariable=processingLabel).grid(column=6, row=6)

   ttk.Button(entryFrame, text="Next",  command= lambda: loadFrame('m',2)).grid(column=10, row=20, sticky=W)
   ttk.Button(entryFrame, text="Previous", command= lambda: loadFrame('d')).grid(column=2, row=20, sticky=E)
   smilesString.set("")
   mimicId.set("")
   cacbSerial.set("")
   for child in entryFrame.winfo_children(): child.grid_configure(padx=5, pady=5)


def importMenuFrame():
   menuFrame = ttk.Frame(root, padding = "3 3 12 12")
   menuFrame.grid(column=0, row=0, sticky=(N, W, E, S))
   menuFrame.columnconfigure(0, weight=1)
   menuFrame.rowconfigure(2, weight=1)

   ttk.Label(menuFrame, text="Mimic Info",  font = "Helvetica 16 bold", foreground="blue").grid(column=1, row=1, sticky=W)
   ttk.Label(menuFrame, text="Mimic Id:", font = "Times 12 bold").grid(column=1, row=3, sticky=(W,E), columnspan=3)
   ttk.Label(menuFrame, text=mimicId.get()).grid(column=4, row=3, sticky=W, columnspan=20)
   ttk.Label(menuFrame, text="Smiles String:", font = "Times 12 bold").grid(column=1, row=4, sticky=(W,E), columnspan=3)
   ttk.Label(menuFrame, text=smilesString.get()).grid(column=4, row=4, sticky=W, columnspan=20)
   ttk.Label(menuFrame, text="Ca Cb Serial Number:", font = "Times 12 bold").grid(column=1, row=5, sticky=(W,E), columnspan=3)
   ttk.Label(menuFrame, text=cacbSerial.get()).grid(column=4, row=5, sticky=W, columnspan=20)
   
   ttk.Label(menuFrame, text="Run Parameters", font = "Helvetica 16 bold", foreground="blue").grid(column=1, row=7, sticky=W, columnspan=2) 
   vcmd = menuFrame.register(isDouble)
   ttk.Label(menuFrame, text="Dielectric Constant").grid(column=1, row=9, sticky=W, columnspan=2)
   ttk.Entry(menuFrame, width=4, textvariable=diel, validate="key", validatecommand=(vcmd, '%P')).grid(column=3, row=9, sticky=(W,E)) 
   ttk.Label(menuFrame, text="Energy Cutoff").grid(column=1, row=10, sticky=W, columnspan=2)
   ttk.Entry(menuFrame, width=2, textvariable=ecutoff, validate="key", validatecommand=(vcmd, '%P')).grid(column=3, row=10, sticky=(W,E))
   ttk.Label(menuFrame, text="RMSD Cutoff").grid(column=1, row=11, sticky=W, columnspan=2)
   ttk.Entry(menuFrame, width=2, textvariable=rmsdcutoff, validate="key", validatecommand=(vcmd, '%P')).grid(column=3, row=11, sticky=(W,E))
   
   #TODO Check if this previous button should be removed as the mol2 files would have been generated 
   #ttk.Button(menuFrame, text="Previous", command= lambda: loadFrame('e')).grid(column=1, row=16, sticky=E)
   ttk.Button(menuFrame, text="Run", command= lambda: runQMD()).grid(column=6, row=16, sticky=W)
   
   for child in menuFrame.winfo_children(): child.grid_configure(padx=5, pady=5)


def isValidId(P):
   for c in P:
      if not ((c>='0' and c<='9') or (c=='_') or (c>='a' and c<='z') or (c>='A' and c<='Z')):
         return False
   return True

def isCommaSeparated(P):
   for c in P:
      if not ((c>='0' and c<='9') or (c==',')):
         return False
   return True


def isDouble(P):
   try:
      float(P)
      return True
   except ValueError:
      return False

   
def selDirectory(dirFrame):
   options = {}
   options['initialdir'] = '.'
   options['mustexist'] = True
   options['parent'] = dirFrame
   options['title'] = 'Choose a directory'
   dirName.set(tkFileDialog.askdirectory(**options))

def loadFrame(opt, nxtFlg = 0):
   global browserOpened
   if nxtFlg ==1 :
      if not (os.path.isdir(dirName.get()) and os.access(dirName.get(), os.R_OK)):
          return
   elif nxtFlg == 2:
      if mimicId.get() == "" or smilesString.get() == "" or cacbSerial.get() == "":
         tkMessageBox.showwarning(title="Warning", message="Please fill in all the fields.")
         return
      cacb = cacbSerial.get().split(",")
      if len(cacb) != 6:
         tkMessageBox.showwarning(title="Warning", message="Please enter 6 atom numbers as Ca Cb Serial Number.")
         return

      dirPath = os.path.join(dirName.get(), mimicId.get())
      if (os.path.exists(dirPath) and browserOpened == False) or (not os.path.exists(dirPath)):
          createMol2Files()

   deleteAllFrames()
   if opt=='d':
      importDirFrame()
   elif opt=='e':
      importEntryFrame()
   elif opt=='m':
      importMenuFrame()

def deleteAllFrames():
   for child in root.winfo_children():
      child.destroy()


def createMol2Files():
   dirPath = os.path.join(dirName.get(), mimicId.get())
   smilesPath = os.path.join(dirPath, mimicId.get()+".smiles")
   nohMolPath = os.path.join(dirPath, mimicId.get()+"_noh.mol2")
   molPath = os.path.join(dirPath, mimicId.get()+".mol2")
   tempPath = os.path.join(dirPath, mimicId.get()+"_temp.mol2")
   if os.path.exists(dirPath):
      decision = tkMessageBox.askyesno(title="Overwrite Directory", message="The directory " + dirPath + " already exists. Do you want to overwrite it?")
      root.update()
      if decision:
         processingLabel.set("Overwriting Directory")
         root.update()
         shutil.rmtree(dirPath)
      else:
         root.destroy()
   
   processingLabel.set("Processing...Generating Mol2 Files")
   root.update()
   os.makedirs(dirPath)
   f = open(smilesPath, "w")
   f.write(smilesString.get())
   f.close()
   cmd = "babel -ismiles " + smilesPath + " -omol2 " + nohMolPath + " -d --gen3D"
   process = subprocess.Popen(cmd, shell=True)
   process.wait()
   cmd = "babel -imol2 " + nohMolPath + " -omol2 " + tempPath + " -h"
   process = subprocess.Popen(cmd, shell=True)
   process.wait()
   cmd = "babel --partialcharge none -imol2 " + tempPath +" -omol2 "+ molPath
   process = subprocess.Popen(cmd, shell=True)
   process.wait()
   processingLabel.set("")
   root.update()
   os.remove(tempPath)

   return molPath

def openPage():
   global browserOpened
   if mimicId.get() == "" or smilesString.get() == "":
        tkMessageBox.showwarning(title="Warning", message="Please fill in the Mimic Id and Smiles String first")
        return;
   molPath  = createMol2Files()
   processingLabel.set("Opening Browser")
   root.update()
   browserOpened = True
   cmd = "chromium-browser --args --allow-file-access-from-files file://" + os.environ.get('EKOGUIHOME') + "/selectatoms.html?path=" + molPath
   #cmd = "chromium-browser --args --allow-file-access-from-files file:///home/rhljain08/gui/test.html?path=/home/rhljain08/gui/temp/KB/KB.mol2"
   process = subprocess.Popen(cmd, shell=True)
   process.wait()
   processingLabel.set("")
   root.update()
   filePath  = os.path.join("Downloads", "cacbserial.txt")
   serialFile = os.path.join(os.path.expanduser("~"), filePath)
   if  os.path.isfile(serialFile):
      f = open(serialFile,"r")
      cacbSerial.set(f.readline().strip())
      f.close()
      os.remove(serialFile)
   else:
       tkMessageBox.showwarning(title="Warning", message="You have not downloaded the file from the browser. Make sure you download the text file or fill in the Ca Cb Serial Number manually.")


def runQMD():
   global browserOpened
   deleteAllFrames()
   
   resultFrame = ttk.Frame(root, padding = "3 3 12 12")
   resultFrame.grid(column=20, row=20, sticky=(N, W, E, S))
   resultFrame.columnconfigure(0, weight=1)
   resultFrame.rowconfigure(2, weight=1)
   ttk.Label(resultFrame, text="").grid(column=1, row=1)
   progressLbl = ttk.Label(resultFrame, text="Running...", font = "Times 16 bold", foreground="blue")
   progressLbl.grid(column=20, row=26,columnspan=2)   
   
   
   dirPath = os.path.join(dirName.get(), mimicId.get())
   progressLbl.config(text='Running...Creating .in File')
   root.update()
   createInFile(dirPath) 
   
   progressLbl.config(text='Running...Creating job File')
   root.update()
   createJobFile(dirPath)
   browserOpened = False
   #Change the directory as the qmd is run from the containing directory.
   cwd = os.getcwd()
   os.chdir(dirPath) 
   cmd = "qsub " + mimicId.get() + ".job"

   progressLbl.config(text='Running...')
   root.update()

   
   process = subprocess.Popen(cmd, shell=True)
   process.wait()
   
   progressLbl.config(text='Submitted to Queue')
   root.update()
   os.chdir(cwd)

def createInFile(dirPath):
   inFile = os.path.join(dirPath , mimicId.get() + ".in")   
   f = open(inFile, "w")
   #f.write(mimicId.get() + " " + smilesString.get() + " " + cacbSerial.get())
   f.write('atoms4rms="'+cacbSerial.get().replace(",", " ") + '"\n')
   f.write("overwritetop=1\ndiel="+str(diel.get())+"\necutoff="+str(ecutoff.get())+"\nrmsdcutoff="+str(rmsdcutoff.get())+"\nrun_qmd=1")
   f.close()

def createJobFile(dirPath):
   jobFile = os.path.join(dirPath, mimicId.get() + ".job")
   f = open(jobFile, "w")
   f.write("""#!/bin/bash
#
#PBS -l walltime=240:00:00,mem=2gb,nodes=1:ppn=1

cd $TMPDIR

echo -n Start Date: 
date

export LD_LIBRARY_PATH=/apps/burgess/gromacs/lib:$LD_LIBRARY_PATH

. /apps/intel/2016_update3/bin/compilervars.sh intel64

export AMBERHOME=/apps/burgess/amber10

export PATH=$AMBERHOME/bin:/apps/burgess/bin:$PATH

cp $PBS_O_WORKDIR/""" + mimicId.get()+".mol2 .\ncp $PBS_O_WORKDIR/" + mimicId.get()+".in .\n\n" +
"""which gromacs_cluster

gromacs_cluster """ + mimicId.get()+".in >& " + mimicId.get() + """_cluster.log 

echo -n End Date: 
date

tar -zcvf """ + mimicId.get()+ """_results_$PBS_JOBID.tgz *
mv """ + mimicId.get() + """_results_$PBS_JOBID.tgz $PBS_O_WORKDIR
#cd $PBS_O_WORKDIR
#tar -zxvf """ + mimicId.get() + """_results_$PBS_JOBID.tgz

exit""")
   f.close()


root = Tk()
root.title("QMD GUI")
root.geometry("900x600+30+30")
dirName = StringVar()
smilesString = StringVar()
mimicId = StringVar()
cacbSerial = StringVar()
processingLabel = StringVar()
diel = DoubleVar()
diel.set(80.0)
ecutoff = DoubleVar()
ecutoff.set(3.0)
rmsdcutoff = DoubleVar()
rmsdcutoff.set(0.5)
browserOpened = False

if os.environ.get('EKOGUIHOME') is None:
   os.environ['EKOGUIHOME'] = "/apps/burgess/eko_qmd_gui"

importDirFrame()
root.mainloop()

