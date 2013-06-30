import urllib
import subprocess

urllib.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/complete*.genomic.fna.gz')
urllib.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/complete*.genomic.gbff.gz')
#urllib.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/complete*.protein.gpff.gz')
bashCommand = "gunzip *.fna.gz"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output = process.communicate()[0]  