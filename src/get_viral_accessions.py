# ==============================================================================
# get_viral_accessions.py

#
# Authors: Roman Stolyarov (r.m.stolyarov@gmail.com)
#
# #Goes through viral.1.1.genomic.fna.gz and constructs table relating 
#full name of virus to NCBI accession. Prints table to file.

# ==============================================================================
__author__ = 'Roman Stolyarov (r.m.stolyarov@gmail.com)'
__version__ = '1.1.1'
__date__ = '23 July 2013'

from collections import deque

def readSettings():
	settings = {}
	f = open("settings.txt")
	data = f.readlines()
	for line in data:
		key, value = line.split("\t")
		settings[key] = value.split()
	return settings


settings = readSettings()
vf = settings['refseqdbloc'][0]+'viral.1.1.genomic.fna'
virus_file = open(vf,'r')
data = virus_file.readlines()
viral_accessions = {}
skipped = 0
skippedlist = []
for line in data:
	if line[0] != ">":
		continue
	header = deque(line.split())
	if ("complete" not in header and "genomic" not in header and "genome" not in header):
		skipped = skipped + 1
		skippedlist.append(header)
		continue
	accession = header.popleft()
	ncbi_accession = accession.split("|")[-2]
	str = None
	name = ""
	i = 0
	while str == None or (header[i] != "complete" and header[i] != "genomic" and header[i] != "genome"):
		str = header[i]
		if str[-1] == ",":
			name = name + str[:-1] + " "
			break
		else:
			name = name + str + " "
		i = i+1
	viral_accessions[name] = ncbi_accession

print "skipped:",skipped
for el in skippedlist:
	print el
with open( "virus_accessions.txt", "w" ) as out:
	for key in viral_accessions:
        	out.write( key+"\t"+viral_accessions[key]+"\n" )