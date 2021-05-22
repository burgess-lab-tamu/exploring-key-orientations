import sys
import re
import datetime
from dateutil.tz import tzlocal

now = datetime.datetime.now(tzlocal())
idDoc = open(sys.argv[1],'r')
dbFile = open(sys.argv[3], 'a')
idSet =  set()

for line in idDoc:
   if 'structureId' in line:
      id = line[line.find('"')+1: line.rfind('"')]
      idSet.add(id)
idDoc.close()

print '<entityInfo>'
reg = r'.*structureId="([^"]*)".*release_date="([^"]*).*'
doc = open(sys.argv[2], 'r')

line  = doc.readline()
while line:
   if 'structureId' in line:
      id = re.match(reg, line).group(1)
      releaseDate = re.match(reg, line).group(2)
      newPDB =  False
      if id not in idSet:
         dbFile.write(id)
         dbFile.write("\t\t\t")
         dbFile.write(releaseDate.ljust(30))
         dbFile.write("\t\t\t")
         dbFile.write(now.strftime("%m/%d/%Y %H:%M %Z"))
         dbFile.write("\n")
         newPDB = True
      while '</PDB>' not in line:
         if newPDB:
            print line,
         line = doc.readline()
      if '</PDB>' in line and newPDB:
         print line,
   line = doc.readline()  
doc.close()         
print '</entityInfo>'

for i in range(0,100):
   dbFile.write("-")
dbFile.write("\n")
dbFile.close()
