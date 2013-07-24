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
import argparse

def parseArgs():
        parser = argparse.ArgumentParser()
		parser.add_argument("--settings", help="Input settings file name", default='settings.txt')
        parser.add_argument("--output", help="Output file name", default='virus_accessions.txt')
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

def categorize(virus_accessions):
	categories = {}
	for full_name in virus_accessions:
		n = len(full_name[0])
		for i in range(n):
			el = full_name[0][i]
			if i+1<n:
				el1 = full_name[0][i+1]
			else:
				el1 = None
			if "virus" in el or "phage" in el:
				if len(el1) == 2:
					scientific_name = full_name[0][:full_name[0].index(el)+2]
				else:
					scientific_name = full_name[0][:full_name[0].index(el)+1]
				if scientific_name in categories:
					categories[scientific_name].append([full_name[0],full_name[1]])
				else:
					categories[scientific_name] = [[full_name[0],full_name[1]]]
	return categories

settings = readSettings(args)
args = parseArgs()
vf = settings['refseqdbloc'][0]+settings['virusfna'][0]
virus_file = open(vf,'r')
data = virus_file.readlines()
viral_accessions_string = {}
viral_accessions_list = []
skipped = 0
skippedlist = []
virus_names = []
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
	virus_names.append(list(header))
	while str == None or (header[i] != "complete" and header[i] != "genomic" and header[i] != "genome"):
		str = header[i]
		if str[-1] == ",":
			name = name + str[:-1] + " "
			break
		else:
			name = name + str + " "
		i = i+1
	viral_accessions_string[name] = ncbi_accession
	viral_accessions_list.append([tuple(header),ncbi_accession])

with open( args.output, "w" ) as out:
	for key in viral_accessions:
        	out.write( key+"\t"+viral_accessions[key]+"\n" )
