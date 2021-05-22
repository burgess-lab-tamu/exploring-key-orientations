#!/usr/bin/python
import sys
import os
import re
from itertools import islice

def genHtmlPage(jsmolPath, chainNames):
   htmlText = '''
<!DOCTYPE html>
<html>
<head>
<script type="text/javascript" src="http://www.kryogenix.org/code/browser/sorttable/sorttable.js"></script>
<script type="text/javascript">
function sortByRms() {
  var myTH = document.getElementsByTagName("th")[5];
  sorttable.innerSortFunction.apply(myTH, []);
}
</script>
<style type="text/css">
table {
    font-family: arial, sans-serif;
    border-collapse: collapse;
    width: 100%;
}

td, th {
    border: 1px solid #dddddd;
    text-align: center;
    padding: 3px !important;
}

tr:nth-child(even) {
    background-color: #dddddd;
}

/* Sortable tables */
table.sortable thead {
    background-color:orange;
    color:black;
    font-weight: bold;
    cursor: default;
}

</style>
</head>
<body onload=sortByRms()>
'''

   #print htmlText

   folderPath = os.getcwd()
   os.chdir(os.path.join(folderPath, "conformers"))
   confPath = os.getcwd()
   infoFile = open(os.path.join(confPath, "info.pdb"), "r")
   summary = open(os.path.join(folderPath, "score_summary.txt"), "r")
   htmlFile =  open(os.path.join(confPath, "viewconformers.html"), "w")
  
   htmlFile.write(htmlText)

   minEnergy = {}
   for line in summary:
      if line.startswith("#"):
         continue
      results = filter(None,re.split(' |\t', line))
      results[11] = results[11].strip()
      results[2] = results[2].strip()
      if results[11] not in minEnergy:
         minEnergy[results[11]] = float(results[2])
      else:
         if float(results[2]) < minEnergy[results[11]] :
            minEnergy[results[11]] = float(results[2])
   summary.close()
   #print minEnergy

   for line in infoFile:
      htmlFile.write(line)
      #print line
   infoFile.close()
   #os.remove(confPath + "/info.pdb")

   htmlText = '''
<table border=1 class=sortable>
  <tr>
    <th title="Scaffold type and the file name containing the atomic coordinates of mimic superimposed on the target structure.">File</th>
    <th title="Identifier of the conformer in the QMD output.">Origframe</th>
    <th title="Energy of conformer from QMD.">&Delta; Energy (kcal/mol)</th>
    <th title="Member of the complex on which the conformer is superimposed.">Chain</th>
    <th title="Residues in the interface on which the conformer is superimposed.">Residues</th>
    <th title="Root Mean Square distance between 6 C&alpha; and C&beta; atomsof matched residues with corresponding atoms in mimic.">RMS(&Aring;)</th>
    <th title="Mean angle in degrees between 3 C&alpha; - C&beta; vectors and the corresponding vectors in the mimic.">Angle</th>
    <th title="Weighted sum of RMS and angle.">Score</th>
    <th title="Click on the link to view superimposed conformer and target protein.">Jsmol Link</th>
  </tr>
'''

   #print htmlText
   htmlFile.write(htmlText)

   files = [f for f in os.listdir('.') if os.path.isfile(f)]
   for File in files:
      if File == "info.pdb" or File == "viewconformers.html":
         continue
      with open(File) as myfile:
         head = [next(myfile) for x in xrange(4)]
         #head = list(islice(File, 4))
      #print head
      energy = 0.0
      orig = ''
      rmsd = 0.0
      chain =''
      rsd = ''
      score = 0.0
      rms =0.0
      angle =0.0
      prot = ''
      deltaEnergy = 0.0
      for line in head:
         if re.match(r'.*origframe ([0-9]+).*Energy: ([0-9.]+).*rmsd: *([0-9.]+).*', line):
            pat = r'.*origframe ([0-9]+).*Energy: ([0-9.]+).*rmsd: *([0-9.]+).*'
            orig = re.match(pat, line).group(1)
            energy = float(re.match(pat, line).group(2))
            mimic = File.split("_")[0] + "_" + File.split("_")[1]
            deltaEnergy = energy - minEnergy[mimic]
            #rmsd = float(re.match(pat, line).group(3))
         elif re.match(r'.*chain=([^ ]+).*resids=\(([^)]+).*', line):
            pat = r'.*chain=([^ ]+).*resids=\(([^)]+).*'
            chain = re.match(pat, line).group(1)
            rsd = re.match(pat, line).group(2).replace(", ", ",")
         elif re.match(r'.*score=([0-9.]+).*rms=([0-9.]+).*angle=([0-9.]+).*', line):
            pat = r'.*score=([0-9.]+).*rms=([0-9.]+).*angle=([0-9.]+).*'
            score = float(re.match(pat, line).group(1))
            rms = float(re.match(pat, line).group(2))
            angle = float(re.match(pat, line).group(3))
         elif re.match(r'.*EKO pdb: *([^ ]+).*', line):
            prot = re.match(r'.*EKO pdb: *([^ ]+).*', line).group(1)
      if rms >= 0.5:
         continue
      #print "   <tr>"
      htmlFile.write("   <tr>\n")
      #print "      <td>"+File+"</td>"
      htmlFile.write("      <td>"+File+"</td>\n")
      #print "      <td>"+orig+"</td>"
      htmlFile.write( "      <td>"+orig+"</td>\n")
      #print "      <td>"+str(round(deltaEnergy, 2))+"</td>"
      htmlFile.write("      <td>"+str(round(deltaEnergy, 2))+"</td>\n")
      #print "      <td>"+chain+"</td>"
      if chain in chainNames:
         htmlFile.write("      <td>"+chain+ " (" + chainNames[chain] + ") " + "</td>\n")
      else:
         htmlFile.write("      <td>"+chain+"</td>\n")
      #print "      <td>"+rsd+"</td>"
      htmlFile.write("      <td>"+rsd+"</td>\n")
      #print "      <td>"+str(round(rms, 2))+"</td>"
      htmlFile.write("      <td>"+str(round(rms, 2))+"</td>\n")
      #print "      <td>"+str(round(angle, 2))+"</td>"
      htmlFile.write("      <td>"+str(round(angle, 2))+"</td>\n")
      #print "      <td>"+str(round(score, 2))+"</td>"
      htmlFile.write("      <td>"+str(round(score, 2))+"</td>\n")
      filePath = confPath+"/"+File
      protPath = folderPath+"/"+prot
      #print '      <td><a href="'+sys.argv[1]+'/JSMOLScripting.html?conf='+filePath+'&prot='+protPath+'&chainRsd='+chain+' '+rsd+'">View in 3D</a></td>'
      htmlFile.write('      <td><a href="'+jsmolPath+'/JSMOLScripting.html?conf='+filePath+'&prot='+protPath+'&chainRsd='+chain+' '+rsd+'">View in 3D</a></td>\n')
      #print '   </tr>'
      htmlFile.write('   </tr>\n')
   os.chdir(folderPath)
   #print '</table>'
   htmlFile.write('</table>\n')
   #print '</body>'
   htmlFile.write('</body>\n')
   #print '</html>'
   htmlFile.write('</html>\n')
