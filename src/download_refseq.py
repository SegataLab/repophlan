import urllib
import subprocess
from ftplib import FTP
import csv
import os

def readSettings():
	settings = {}
	f = open("settings.txt")
	data = f.readlines()
	for line in data:
		key, value = line.split("\t")
		settings[key] = value.split()
	return settings

#Parse settings file, set relevant variables
settings = readSettings()
refseqdbftp = settings['refseqdbftp'][0]
refseqdbpath = settings['refseqdbpath'][0]
refseqdbloc = settings['refseqdbloc'][0]
refseqdburl = 'ftp://'+refseqdbftp+'/'+refseqdbpath

#Access the FTP, download files with names *.genomic.fna.gz or *.genomic.gbff.gz or *.mstr.gbff.gz
ftp = FTP(refseqdbftp)
ftp.login()
ftp.cwd(refseqdbpath)
all_files = []
ftp.dir(all_files.append)
print 'Got list of files'
for file in all_files:
	fn = file.split()[-1]
	if fn[:8] == 'complete' and (fn[-14:] == 'genomic.fna.gz' or fn[-15:] == 'genomic.gbff.gz' or fn[-12] == 'mstr.gbff.gz'):
		print 'Downloading',fn,'....'
		try:
  			bashCommand = 'wget -N -P '+refseqdbloc+' '+refseqdburl+fn
          		process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
           		output = process.communicate()[0]
           		print "Downloaded "+fn
    		except:
           		print "Error downloading "+fn

#Uncompress *.fna.gz, but keep the compressed files to allow timestamping of future downloads
bashCommand = 'for FILE in '+refseqdbloc+'*.fna.gz; do gunzip -c $FILE > ${FILE%.*}; done'
os.system(bashCommand)
