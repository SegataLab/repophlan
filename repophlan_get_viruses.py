#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import urllib2
import sys
import argparse
from StringIO import StringIO
import gzip
import os
import logging

def read_params(args):
    parser = argparse.ArgumentParser(description='')
    arg = parser.add_argument
    arg( '--refseq_virus_gbff', type=str,
         default = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz",
         help="The gbff NCBI file for viruses (default is ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz)")
    arg( '--taxonomy', metavar='taxonomy', required = True, type=str,
         help="The taxonomy file")
    arg( '--out_dir', required = True, default = None, type = str, 
         help="The output folder")
    arg( '--out_summary', required = True, default = None, type = str, 
         help="The output summary file")
    return vars(parser.parse_args())


if __name__ == '__main__':
    par = read_params(sys.argv)
    logging.basicConfig(level=logging.INFO, stream=sys.stdout, 
                        format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logger = logging.getLogger(sys.argv[0])
    logger.info('Reading the taxonomy from '+par['taxonomy']+'... ')
    with open(par['taxonomy']) as inpf:
        taxids2taxonomy = dict([l.strip().split('\t')[1:] for l in inpf])
    logger.info('Done.')
    
    
    logger.info('Downloading and reading '+par['refseq_virus_gbff']+'... ')
    table = urllib2.urlopen( par['refseq_virus_gbff'] )
    f = gzip.GzipFile(fileobj= StringIO( table.read())  )
    logger.info('Done.')



    if not os.path.exists( par['out_dir'] ):
        os.mkdir( par['out_dir'] )
    outdir = par['out_dir']+ "/"

    summary = {}

    nrec = 0
    for seq_record in SeqIO.parse( f, "genbank"):
        accession = seq_record.annotations['accessions'][0]
        logger.info('Processing accession '+accession)

        fna = seq_record
        ffn, faa = [], []
        taxid = None
   
        for feat in seq_record.features:
            if feat.type == 'source':
                if 'db_xref' in feat.qualifiers:
                    for tn in feat.qualifiers['db_xref']:
                        if 'taxon' in tn:
                            taxid = tn.split(":")[-1]
    
            if feat.type == 'gene':
                if 'db_xref' in feat.qualifiers:
                    gene_id = feat.qualifiers['db_xref'][0]
                    gene = feat.location.extract(fna)
                    gene.id = gene_id
                    ffn.append( gene )
            if feat.type == 'CDS':
                if 'translation' in feat.qualifiers:
                    prot = SeqRecord( Seq(feat.qualifiers['translation'][0])  )
                    if 'protein_id' in feat.qualifiers:
                        prot.id = feat.qualifiers['protein_id'][0]
                    if 'product' in feat.qualifiers:
                        prot.description = feat.qualifiers['product'][0]
                    faa.append( prot )
       
        if taxid is None:
            continue
        
        if taxid not in taxids2taxonomy:
            logger.info('TaxId not found in taxonomy:  '+taxid)
            continue

        nrec += 1
        
        SeqIO.write( fna, outdir+accession+".fna", "fasta" )
        SeqIO.write( ffn, outdir+accession+".ffn", "fasta" )
        SeqIO.write( faa, outdir+accession+".faa", "fasta" )
  
        summary[accession] = {  'taxonomy':taxids2taxonomy[taxid],
                                'lfna': outdir+accession+".fna", 
                                'lffn': outdir+accession+".ffn", 
                                'lfaa': outdir+accession+".faa" } 
        #print summary

        logger.info('Files exported for '+accession+' [record # '+str(nrec)+'] '+",".join([outdir+accession+".fna",outdir+accession+".ffn",outdir+accession+".faa"]))

    logger.info('Writing summary file to: '+par['out_summary'])
    with open( par['out_summary'], "w" ) as outf:
        for k,v in sorted(summary.items(),key=lambda y:y[0]):
            outf.write( "\t".join( [k,v['taxonomy']] ) +"\n"  )
    logger.info('Done. Exiting.')
    

