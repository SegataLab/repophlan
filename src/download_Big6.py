import urllib 
import tarfile
import os
import subprocess
import argparse
import csv


def parseArgs():
	parser = argparse.ArgumentParser()
	parser.add_argument("--settings", help="Settings file name", default='settings.txt')
	parser.add_argument("--rsversion", help="Version of RefSeq catalog to download", default='00')
	parser.add_argument("--prok_file", help="Prokaryotes file name", default="prokaryotes.txt")
	parser.add_argument("--euk_file", help="Eukaryotes file name", default="eukaryotes.txt")
	parser.add_argument("--vir_file", help="Viruses file name", default="viruses.txt")
	parser.add_argument("--ncbi_nodes", help="NCBI nodes file name", default="nodes.txt")
	parser.add_argument("--ncbi_names", help="NCBI names file name", default="names.txt")
	parser.add_argument("--logfile", help="Log file name", default='log.txt')
	args = parser.parse_args()
	return args

def readSettings():
	settings = []
	with open("settings.txt") as tsv:
		for line in csv.reader(tsv,dialect="excel-tab"):
			settings.append(line)
	return settings


def makeNewDirectory(newpath): 
	if not os.path.exists(newpath): 
  		os.makedirs(newpath)

args = parseArgs()
makeNewDirectory(r'Big6')
settings = readSettings()
f = open('Big6/'+args.logfile,'w')
f.write('Log file for download_Big6.py:\n')

#Download all files
taxdumpurl = settings[0][1]
taxdump = settings[1][1]
refsequrl = settings[2][1]
rf = settings[3][1]
textfilesurl = settings[4][1]
prok_file = settings[5][1]
euk_file = settings[6][1]
vir_file = settings[7][1]
delfiles = settings[8][1:]

try:
	os.system('wget -P Big6/'+taxdump+' '+taxdumpurl+taxdump)
	f.write("Downloaded",taxdump)
except:
	f.write("Error downloading "+taxdump)

refseq = rf+args.rsversion+'.catalog'
refseqgz = refseq+'.gz'
if not os.path.isfile(refseq):
  try:
	os.system('wget -P Big6/'+refseqgz+' '+refsequrl+refseqgz)
    	f.write("Downloaded latest version of RefSeq catalog")
  except:
    	f.write("Error downloading specified version of RefSeq catalog. Please ensure you specify latest RefSeq catalog version.")
else:
	f.write("Latest version of RefSeq catalog already on file")

tfs = [[prok_file,args.prok_file],[euk_file,args.euk_file], [vir_file,args.vir_file]]
for tf in tfs:
	try:
		os.system('wget -P Big6/'+tf[1]+' '+textfilesurl+tf[0])
  		f.write("Downloaded "+tf[0]+" and saved to local file "+tf[1])
	except:
		f.write("Error downloading "+tf[0])
  
#Untar/gunzip
os.chmod('Big6')
tfile = tarfile.open(taxdump)
if tarfile.is_tarfile(taxdump):
  # extract all contents
  tfile.extractall('.')
else:
  f.write(taxdump + " is not a tarfile.")
    
bashCommand = "gunzip "+refseqgz
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output = process.communicate()[0]    
f.write("Unzipped refseq catalog")    

for df in delfiles:
  os.remove(df)
os.rename('names.dmp',args.ncbi_names)
os.rename('nodes.dmp',args.ncbi_nodes)

f.write("Deleted extra files")
