import urllib
import subprocess
from ftplib import FTP

ftp = FTP('ftp.ncbi.nlm.nih.gov')
ftp.login()
ftp.cwd('refseq/release/complete/')
all_files = []
ftp.dir(all_files.append)
print 'Got list of files'
for file in all_files:
  fn = file.split()[-1]
  if fn[:8] == 'complete' and fn[-14:] == 'genomic.fna.gz':
    print 'Downloading',fn,'....'
    urllib.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/'+fn,'/banche_dati/sharedCM/genomicDB/'+fn)
    print 'Downloaded',fn







bashCommand = "gunzip /banche_dati/sharedCM/genomicDB/*.gz"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output = process.communicate()[0]  
