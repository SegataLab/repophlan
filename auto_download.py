import urllib 
import tarfile
import os

newpath = r'Big6' 
if not os.path.exists(newpath): os.makedirs(newpath)
os.chdir(newpath)
#Download all files
taxdump = 'taxdump.tar.gz'
refseq = 'RefSeq-release59.catalog.gz'
urllib.urlretrieve('ftp://ftp.ncbi.nih.gov/pub/taxonomy/'+taxdump, taxdump)
print "Downloaded",taxdump
urllib.urlretrieve(' ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/'+refseq, refseq)
print "Downloaded RefSeq catalog"
textfiles = ['prokaryotes.txt', 'eukaryotes.txt', 'viruses.txt']
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
    
bashCommand = "gunzip "+refseq
import subprocess
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output = process.communicate()[0]    

print "Unzipped refseq catalog"    

delfiles = ['taxdump.tar.gz','delnodes.dmp','merged.dmp','gencode.dmp','gc.prt','division.dmp','citations.dmp']
for df in delfiles:
  os.remove(df)

print "Deleted extra files"
  

