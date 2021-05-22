#!/usr/bin/python
import sys
import os
import re
from sets import Set

def createConformers():
   filesProcessed = Set()
   processedLines = Set()
   dateCount = 0
   pdbName = ''
   commandLines = []

   if not os.path.isdir("./conformers"):
      os.makedirs("conformers")
      #os.chdir("./conformers")
      allFiles = []
      allFiles += [each for each in os.listdir("./") if each.endswith('_overlay.out')]
      cnt = 1
      infoFile = open("./conformers/info.pdb", "a+")
      lines= []
      pdbLines = []
      newConformer =  False
      for File in allFiles:
         f = open(File, "r+")
         for line in f:
           #line =  line.strip('\n')
           if newConformer:
              pdbLines.append(line)
              line =  line.strip('\n')
              if line == "END":
                 newConformer = False
                 fname = File.replace('lt3_overlay.out','File'+str(cnt)+'.pdb')
                 fpath = './conformers/'+fname
                 pdbFile =  open(fpath, 'w+')
                 pdbFile.writelines(pdbLines)
                 pdbFile.close()
                 pdbLines = []
                 cnt += 1
           else:
	      if re.match(r'REMARK.*origframe.*', line):
                 newConformer = True
                 pdbLines.append(line)
                 continue
              elif re.match(r'REMARK.*command.*', line) and File not in filesProcessed:
                 filesProcessed.add(File)
              
                 if pdbName == '':
                    regx = r'REMARK.*'+File.replace("_overlay.out", ".pdb")+'.* (.*).pdb.*'
                    pdbName = re.match(regx, line).group(1) + '.pdb'
		    infoFile.write('<h1 align="center"> Matching Results </h1>\n')
                    commandLines.append('<font size="5"><b> Run Parameters </b></font>\n<br/>')
	      
	         #infoFile.write(line)
                 commandLines.append(line+'\n<br/>\n')
	      line =  line.strip('\n')
              if (re.match(r'REMARK.*version.*', line) or re.match(r'REMARK.*run date.*', line) or re.match(r'REMARK.*params.*', line) or re.match(r'REMARK.*interface.*', line)) and line not in processedLines:
                 if re.match(r'REMARK.*run date.*', line):
                    if dateCount != 0:
                       continue
                    dateCount += 1
                 processedLines.add(line)
                 #lines.append(line+'\n')
	         lines.append(line+'\n<br/>\n')
         f.close()   
      lines.append( '<br/>\n')
      commandLines =  commandLines + lines
      #infoFile.writelines(lines)
      #infoFile.write('<br/>\n')
      lines = []
      infoFile.write('<font size="5"><b> PDB Info: ' + pdbName+' </b></font>\n<br/>\n')
      with open(pdbName) as pdbFile:
         head = [next(pdbFile) for x in xrange(50)]
      for line in head:
         if re.match(r'TITLE.*', line):
            #lines.append(line)
            lines.append(line.strip('\n')+'\n<br/>\n')
         elif re.match(r'COMPND.*MOLECULE.*', line):
            #lines.append(re.match(r'COMPND.*(MOLECULE.*);.*', line).group(1) + '\n')
            molecule = re.match(r'COMPND.*(MOLECULE.*);.*', line).group(1)
            lines.append('<br/>\n'+molecule+'\n<br/>\n')
         elif re.match(r'COMPND.*CHAIN.*', line):
            #lines.append(re.match(r'COMPND.*(CHAIN.*);.*', line).group(1) + '\n')
            chain = re.match(r'COMPND.*(CHAIN.*);.*', line).group(1)
            lines.append(chain+'\n<br/>\n')
      infoFile.writelines(lines)
      infoFile.write('<br/>\n')
      infoFile.writelines(commandLines)
      infoFile.write('<br/>\n')
      commandLines = []
      lines = []
   
      infoFile.close()                   
