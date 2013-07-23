#!/usr/bin/env python
import os
import sys
import glob
import gzip
import argparse

def read_params(args):
    parser = argparse.ArgumentParser(description='')
    arg = parser.add_argument
    arg( '--refseq_dir', metavar='refseq_dir', required = True, type=str,
         help="The refseq main folder")
    arg( '--output', metavar='output', required = True, type=str,
         help="The output file")
    return vars(parser.parse_args())

def get_scaffolds( inp ):
    scfs = []
    #Go through lines of master record
    for l in gzip.open( inp ):
        #Get scaffold master accession number
        if l.startswith( "ACCESSION" ):
            acc = l.strip().split()[1]
        #Get individual contig accession numbers
        elif l.startswith( "WGS_SCAFLD" ): 
            #There will either be a single accession or a range. 
            #If range, get all accessions in the range. Append to scfs.
            scf = l.strip().split()[1].split("-")
            if len(scf) == 1:
                scfs.append( scf[0] )
            else:
                fr,to = scf
                fri,toi = int(fr[7:]),int(to[7:])
                for r in range(fri,toi+1):
                    scfs.append( fr[:7]+str(r)+".1" )
    return (acc,",".join(scfs))

def run(pars):
    scafs = []
    
    #For all draft genome master records, if the annotation file 
    #is not available, it means they are actually stored in the 
    #same files where "final" genomes are. Only, they would have 
    #more than one contig per chromosome, since they are still 
    # draft genomes. For draft genomes stored in this way, we need 
    #to create scaffs file which gives accessions of every contig
    #pertaining to that project.
    for mstr in glob.glob(pars['refseq_dir']+"/completeNZ_*.mstr.gbff.gz" ):
        fna_in = pars['refseq_dir']+".".join(os.path.basename( mstr ).split(".")[:1]) + ".genomic.gbff.gz" 
        if not os.path.exists( fna_in ):
            scaf = get_scaffolds( mstr )
            if scaf:
                scafs.append(scaf)
                
    #Output all hash table entries as separate lines to file.
    with open( pars['output'], "w" ) as out:
        for f,t in scafs:
            out.write( f+"\t"+t+"\n" )


if __name__ == '__main__':
    par = read_params(sys.argv)
    run( par )


