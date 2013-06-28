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
    for l in gzip.open( inp ):
        if l.startswith( "ACCESSION" ):
            acc = l.strip().split()[1]
        elif l.startswith( "WGS_SCAFLD" ): 
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
    for mstr in glob.glob(pars['refseq_dir']+"/completeNZ_*.mstr.gbff.gz" ):
        fna_in = pars['refseq_dir']+".".join(os.path.basename( mstr ).split(".")[:1]) + ".genomic.gbff.gz" 
        if not os.path.exists( fna_in ):
            scaf = get_scaffolds( mstr )
            if scaf:
                scafs.append(scaf)
    with open( pars['output'], "w" ) as out:
        for f,t in scafs:
            out.write( f+"\t"+t+"\n" )


if __name__ == '__main__':
    par = read_params(sys.argv)
    run( par )


