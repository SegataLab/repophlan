import urllib 
import tarfile
import os
import subprocess
import argparse
import csv
import time

def parseArgs():
	parser = argparse.ArgumentParser()
	parser.add_argument("--settings", help="Settings file name", default='settings.txt')
	parser.add_argument("--rsversion", help="Version of RefSeq catalog to download", default='00')
	parser.add_argument("--prok_file", help="Prokaryotes file name", default="prokaryotes.txt")
	parser.add_argument("--euk_file", help="Eukaryotes file name", default="eukaryotes.txt")
	parser.add_argument("--vir_file", help="Viruses file name", default="viruses.txt")
	parser.add_argument("--ncbi_nodes", help="NCBI nodes file name", default="nodes.dmp")
	parser.add_argument("--ncbi_names", help="NCBI names file name", default="names.dmp")
	parser.add_argument("--logfile", help="Log file name", default="log.txt")
	parser.add_argument("--outputdir", help="Output directory name", default="Big6")
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
	os.chdir(newpath)

def download(url,logfile):
	try:
        	os.system('wget -q -N '+url)
        	f.write("Downloaded file from "+url+"\n")
	except:
        	f.write("Error downloading file from "+url+"\n")

args = parseArgs()
settings = readSettings()
makeNewDirectory(args.outputdir)
f = open(args.logfile,'a')
f.write(time.strftime("\nRun of download_Big6 on %Y-%m-%d %H:%M:%S\n", time.gmtime()))

#Initialize settings variables
taxdumpurl = settings[0][1]
taxdump = settings[1][1]
refsequrl = settings[2][1]
rf = settings[3][1]
textfilesurl = settings[4][1]
prok_file = settings[5][1]
euk_file = settings[6][1]
vir_file = settings[7][1]
delfiles = settings[8][1:]

download(taxdumpurl+taxdump,args.logfile)

refseq = rf+args.rsversion+'.catalog'
refseqgz = refseq+'.gz'

if not os.path.isfile(refseq):
	download(refsequrl+refseqgz,args.logfile)
else:
	f.write("Latest version of RefSeq catalog already on file\n")

tfs = [[prok_file,args.prok_file],[euk_file,args.euk_file], [vir_file,args.vir_file]]
for tf in tfs:
	download(textfilesurl+tf[0],args.logfile)
	os.rename(tf[0],tf[1])
  
#Untar/gunzip
tfile = tarfile.open(taxdump)
if tarfile.is_tarfile(taxdump):
  # extract all contents
  tfile.extractall('.')
else:
  f.write(taxdump + " is not a tarfile.\n")

if os.path.isfile(refseqgz):    
	bashCommand = "gunzip "+refseqgz
	process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
	output = process.communicate()[0]    
	f.write("Unzipped refseq catalog\n")    

for df in delfiles:
  os.remove(df)
os.rename('names.dmp',args.ncbi_names)
os.rename('nodes.dmp',args.ncbi_nodes)

f.write("Deleted extra files\n")
