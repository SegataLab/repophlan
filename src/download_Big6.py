# ==============================================================================
# download_Big6.py
#
# Authors: Roman Stolyarov (r.m.stolyarov@gmail.com)
#
# Downloads latest versions of prokaryotes.txt, eukaryotes.txt, and viruses.txt from NCBI.
# Downloads specified version of refseq catalog (if it is latest).
# Downloads latest taxdump.tar.gz. Extracts names.dmp and nodes.dmp.
# ==============================================================================

__author__ = 'Roman Stolyarov (r.m.stolyarov@gmail.com)'
__version__ = '1.1.1'
__date__ = '9 July 2013'


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

def readSettings(args):
	settings = {}
	f = open(args.settings)
	data = f.readlines()
	for line in data:
    		key, value = line.split("\t")
    		settings[key] = value.split()
	return settings


def makeNewDirectory(newpath): 
	if not os.path.exists(newpath): 
  		os.makedirs(newpath)
	os.chdir(newpath)

def download(url,logfile):
	try:
		bashCommand = 'wget -q -N '+url
                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
		output = process.communicate()[0]
		f.write("Downloaded file from "+url+"\n")
	except:
        	f.write("Error downloading file from "+url+"\n")

args = parseArgs()
settings = readSettings(args)
makeNewDirectory(args.outputdir)
f = open(args.logfile,'a')
f.write(time.strftime("\nRun of download_Big6 on %Y-%m-%d %H:%M:%S\n", time.gmtime()))

#Initialize settings variables
taxdumpurl = settings['taxdumpurl'][0]
taxdump = settings['taxdump'][0]
refsequrl = settings['refsequrl'][0]
rf = settings['refseq'][0]
textfilesurl = settings['textfilesurl'][0]
prok_file = settings['prok_file'][0]
euk_file = settings['euk_file'][0]
vir_file = settings['vir_file'][0]
delfiles = settings['delete_files']

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
