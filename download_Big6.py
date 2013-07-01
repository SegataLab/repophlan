import urllib 
import tarfile
import os
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--rsversion", help="Version of RefSeq catalog to download")
parser.add_argument("--prok_file", help="Prokaryotes file")
parser.add_argument("--euk_file", help="Eukaryotes file")
parser.add_argument("--vir_file", help="Viruses file")
parser.add_argument("--ncbi_nodes", help="NCBI nodes file")
parser.add_argument("--ncbi_names", help="NCBI names file")
args = parser.parse_args()

newpath = r'Big6' 
if not os.path.exists(newpath): 
  os.makedirs(newpath)

os.chdir(newpath)
#Download all files
taxdump = 'taxdump.tar.gz'
urllib.urlretrieve('ftp://ftp.ncbi.nih.gov/pub/taxonomy/'+taxdump, taxdump)
print "Downloaded",taxdump

refseq = 'RefSeq-release'+args.rsversion+'.catalog'
refseqgz = refseq+'.gz'
if not os.path.isfile(refseq):
  try: 
    urllib.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/'+refseqgz, refseqgz)
    print "Downloaded latest version of RefSeq catalog"
  except:
    print "Specified version of RefSeq NOT downloaded. Please specify latest RefSeq catalog version."
else:
  print "Latest version of RefSeq catalog already on file"
    
textfiles = [args.prok_file, args.euk_file, args.vir_file]
for tf in textfiles:
  urllib.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/'+tf, tf)
  print "Downloaded",tf
  
#Untar/gunzip
tfile = tarfile.open(taxdump)
if tarfile.is_tarfile(taxdump):
  # extract all contents
  tfile.extractall('.')
else:
  print taxdump + " is not a tarfile."
    
bashCommand = "gunzip "+refseqgz
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output = process.communicate()[0]    
print "Unzipped refseq catalog"    
delfiles = ['taxdump.tar.gz','delnodes.dmp','merged.dmp','gencode.dmp','gc.prt','division.dmp','citations.dmp']
for df in delfiles:
  os.remove(df)

print "Deleted extra files"