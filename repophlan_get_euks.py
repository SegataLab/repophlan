#!/usr/bin/env python
import urllib2
import argparse
import sys
from ftplib import FTP
from time import time, sleep
import pickle
import subprocess as sb
import multiprocessing as mp
from Bio import SeqIO
import math
import logging
import os
import tarfile
import StringIO
import collections

FTP_prot = "ftp://"
NCBI_ftp = "ftp.ncbi.nlm.nih.gov"
NCBI_assembly_folder = '/genomes/Fungi/'
NCBI_prokaryotes_file = '/genomes/GENOME_REPORTS/eukaryotes.txt'
NCBI_euks_ids = '/genomes/IDS/Eukaryota.ids'

logging.basicConfig(level=logging.INFO, stream=sys.stdout, 
                    format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
logger = logging.getLogger(sys.argv[0])

def read_params(args):
    parser = argparse.ArgumentParser(description='')
    arg = parser.add_argument
    arg( '--nproc', default = 15, type = int, 
         help="The number of parallel download processes")
    arg( '--out_dir', required = True, type = str, 
         help="The output folder")
    arg( '--taxonomy', metavar='taxonomy', required = True, type=str,
         help="The taxonomy file")
    arg( '--out_summary', required = True, type = str, 
         help="The output summary file")
    return vars(parser.parse_args())

def add_protocol( fn ):
    if fn.startswith("ftp"):
        return "ftp://"+fn
    return "http://"+fn


def retry(tries, delay=3, backoff=2):
    if backoff <= 1:
        raise ValueError("backoff must be greater than 1")
    
    tries = math.floor(tries)
    if tries < 0:
        raise ValueError("tries must be 0 or greater")
    
    if delay <= 0:
        raise ValueError("delay must be greater than 0")
    
    def deco_retry(f):
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay # make mutable
            #rv = f(*args, **kwargs) # first attempt
            while mtries > 1:
                try:
                    return f(*args, **kwargs) 
                except Exception, e:
                    if "No such file or directory" in e.reason and "550" in e.reason:
                        msg = "No remote file found (some ffn and faa are know to be missing remotely). Aborting the download of %s. %s" % (str(args[0]), str(e) )
                        logger.warning(msg)
                        raise e 
                    else:
                        msg = "%s: %s, Retrying in %d seconds..." % (str(args[0]), str(e), mdelay)
                        logger.warning(msg)
                sleep( mdelay )
                mtries -= 1  # consume an attempt
                mdelay *= backoff  # make future wait longer
            return f(*args, **kwargs) # Ran out of tries :-(
        return f_retry # true decorator -> decorated function
    return deco_retry  # @retry(arg[, ...]) -> true decorator

@retry(tries=8, delay=20, backoff=2)
def get_remote_file_with_retry( url, info_str = "", compr_seq = False ):
    loc = urllib2.urlopen( url, timeout = 200 )

    logger.info('Parsing of '+info_str+' '+url)

    if url.endswith( ".tgz" ):
        compressedFile = StringIO.StringIO(loc.read()) 
        tarf = tarfile.open(fileobj=compressedFile)
        seqs = []
        for m in tarf.getmembers():
            #seqs += list(StringIO.StringIO(tarf.extractfile( m )))
            seqs += list(SeqIO.parse(tarf.extractfile( m ),"fasta"))
    else:
        seqs = StringIO.StringIO(loc.read())
    return seqs

def euk_genome_dwl( info ):
    seqtypes = ['fna','ffn','frn','faa']
    try:
        res = dict([(s,[]) for s in seqtypes]) 
        #dwl_folder = NCBI_ftp + NCBI_assembly_folder + info['dwl_folder']

        for chromosome in info['dwl_contigs']:
            ch_name = chromosome.split(".")[0]

            for file_type in seqtypes:
                try:
                    ret = get_remote_file_with_retry( add_protocol( NCBI_ftp + NCBI_assembly_folder + chromosome+"."+file_type ), "FINAL"  )
                except Exception, e:
                    logger.error('Error in downloading of '+NCBI_ftp + NCBI_assembly_folder + chromosome+"."+file_type+' '+str(e)) 
                    break
                res[file_type] += list(SeqIO.parse(ret,'fasta'))
        
        for file_type in seqtypes:
            if len(res[file_type]):
                out_fna = info['outdir']+"/"+file_type+"/"+info['converted_assn']+"."+file_type
                if not os.path.exists( info['outdir']+"/"+file_type+"/"):
                    os.mkdir( info['outdir']+"/"+file_type+"/" )
                SeqIO.write( res[file_type], out_fna, "fasta")
            else:
                print "empty "+file_type
    except Exception, e2: 
        logger.error("Something wrong in retrieving information for "+info['converted_assn']+": "+str(e2) )

def final_genome_dwl( info ):
    seqtypes = ['fna','ffn','frn','faa']
    try:
        res = dict([(s,[]) for s in seqtypes]) 
        dwl_folder = NCBI_ftp + NCBI_assembly_folder + info['dwl_folder']

        for chromosome in info['Chromosomes/RefSeq'].split(","):
            ch_name = chromosome.split(".")[0]

            for file_type in seqtypes:
                try:
                    ret = get_remote_file_with_retry( add_protocol( dwl_folder+"/"+ch_name+"."+file_type ), "FINAL"  )
                except Exception, e:
                    logger.error('Error in downloading of '+dwl_folder+"/"+ch_name+"."+file_type+' '+str(e)) 
                    break
                res[file_type] += list(SeqIO.parse(ret,'fasta'))
        
        for file_type in seqtypes:
            if len(res[file_type]):
                out_fna = info['outdir']+"/"+file_type+"/"+info['converted_assn']+"."+file_type
                if not os.path.exists( info['outdir']+"/"+file_type+"/"):
                    os.mkdir( info['outdir']+"/"+file_type+"/" )
                SeqIO.write( res[file_type], out_fna, "fasta")
            else:
                print "empty "+file_type
    except Exception, e2: 
        logger.error("Something wrong in retrieving information for "+info['converted_assn']+": "+str(e2) )


def draft_genome_dwl( info ):
    seqtypes = ['fna','ffn','frn','faa']
    try:
        res = dict([(s,[]) for s in seqtypes]) 
        dwl_folder = NCBI_ftp + NCBI_assembly_folder + info['dwl_folder']

        for file_type in seqtypes:
            try:
                rname = dwl_folder+"/NZ_"+info['WGS'][:4]+'00000000.scaffold.'+file_type+'.tgz'
                ret = get_remote_file_with_retry( add_protocol( rname ), "DRAFT", compr_seq = True   )
            except Exception, e:
                logger.error('Error in downloading of '+rname+' '+str(e))
                break
            res[file_type] += ret

        for file_type in seqtypes:
            if len(res[file_type]):
                out_fna = info['outdir']+"/"+file_type+"/"+info['converted_assn']+"."+file_type
                if not os.path.exists( info['outdir']+"/"+file_type+"/"):
                    os.mkdir( info['outdir']+"/"+file_type+"/" )
                SeqIO.write( res[file_type], out_fna, "fasta")
            else:
                print "empty "+file_type

    except Exception, e2: 
        logger.error("Something wrong in retrieving information for "+info['converted_assn']+": "+str(e2) )


def get_table_by_assn( remote_file, key, sep = '\t' ):
    table = urllib2.urlopen( remote_file )
    table = [sline.strip().split(sep) for sline in table] 
    table_header, table_content = [v.replace('#','') for v in table[0]],table[1:] 
    ret_table = dict([(d[key],d) for d in [dict(zip(table_header,l)) for l in table_content]
                if key in d])
    return ret_table

def get_ids( remote_file, sep = '\t' ):
    table = urllib2.urlopen( remote_file )
    table = [sline.strip().split(sep) for sline in table]

    retd = collections.defaultdict(list)
    for t in table:
        if "chromosome" in t[6]:
            retd[t[0]].append( t[1] )

    return retd


def get_contigs( remote_dir ):
    ftp = FTP(NCBI_ftp) 
    ftp.login()
    ftp.cwd( remote_dir )
    dirs = ftp.nlst("*")

    retd = {} 
    for d in dirs:
        contig = ".".join( d.split("/")[-1].split(".")[:-1] )
        retd[contig] = d.split("/")[0] 
    return retd 

def get_assn2folder( remote_dir ):
    ftp = FTP(NCBI_ftp)
    ftp.login()
    ftp.cwd( remote_dir )
    dirs = ftp.nlst()

    outd = {}
    for d in dirs:
        if "uid" in d:
            acc = d.split("uid")[-1]
            outd[acc] = d 
    return outd 
    
def conv_ass_format( s ):
    s = s.replace('GCA','GCF')
    return s.split(".")[0]

if __name__ == '__main__':
    par = read_params(sys.argv)
   
    logger.info('Downloading and reading the NCBI eukaryotic assembly list')
    ass_list = get_table_by_assn( add_protocol(NCBI_ftp + NCBI_prokaryotes_file), key = 'Assembly Accession' )
    logger.info('NCBI assembly list read: '+str(len(ass_list))+' assemblies')
    
    logger.info('Downloading and reading the NCBI eukaryotic assembly folders')
    ass2folder = get_assn2folder( NCBI_assembly_folder )
    logger.info('NCBI assembly folder scanned: '+str(len(ass2folder))+' folders')

    
    logger.info('Downloading and reading the NCBI eukaryotic assembly list')
    ass2entries = get_ids( add_protocol(NCBI_ftp + NCBI_euks_ids))
    logger.info('NCBI assembly list read: '+str(len(ass_list))+' assemblies')


    logger.info('Reading the taxonomy from '+par['taxonomy']+'... ')
    with open(par['taxonomy']) as inpf:
        taxids2taxonomy = dict([l.strip().split('\t')[1:] for l in inpf])
    logger.info('Done.')

    """
    for k,v in ass2entries.items():
        if k in taxids2taxonomy:
            print k,taxids2taxonomy[k],v
    """

    contig2folders = get_contigs( NCBI_assembly_folder )
    
    
    if not os.path.exists( par['out_dir'] ):
        os.mkdir( par['out_dir'] )
        logger.warning(par['out_dir']+" does not exist. Creating it.")
    
    nproc = par['nproc'] 
    logger.info('Initializating the pool of '+str(nproc)+" downloaders")
    pool = mp.Pool( nproc )
    
    for assn, info in ass_list.items():
        if not 'TaxID' in info:
            continue

        #if not info['TaxID'] in contig2folders:
        #    continue

        if not info['TaxID'] in taxids2taxonomy:
            continue

        if not info['TaxID'] in ass2entries:
            continue

        ass_list[assn]['dwl_contigs'] = [contig2folders[c]+"/"+c for c in ass2entries[info['TaxID']] if c in contig2folders]
        if len(ass_list[assn]['dwl_contigs']) < 1:
            continue
        
        convass = conv_ass_format(assn)
        ass_list[assn]['outdir'] = par['out_dir']
        ass_list[assn]['converted_assn'] = convass
        ass_list[assn]['taxonomy'] = taxids2taxonomy[info['TaxID']]
        pool.map_async( euk_genome_dwl, [ass_list[assn]] )

    pool.close()
    pool.join()
    
    seqtypes = ['fna','ffn','frn','faa']
    out_sum = []
    for assn, info in ass_list.items():
        if 'converted_assn' not in info:
            continue
        all_files = [os.path.exists( par['out_dir']+"/"+st+"/"+info['converted_assn']+"."+st) for st in seqtypes]
        if all( all_files ):
            out_sum.append( [info['converted_assn'],info['taxonomy']] )
    with open(par['out_summary'],"w") as outf:
        for s in out_sum:
            outf.write( "\t".join(s) +"\n" )

