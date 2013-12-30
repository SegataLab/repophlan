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

FTP_prot = "ftp://"
NCBI_ftp = "ftp.ncbi.nlm.nih.gov"
NCBI_assembly_folder = '/genomes/ASSEMBLY_BACTERIA/'
NCBI_prokaryotes_file = '/genomes/GENOME_REPORTS/prokaryotes.txt'
GENBANK_draft_bacteria = '/genbank/genomes/Bacteria_DRAFT/'

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

def draft_genome_genbank_dwl( info ):
    seqtypes = ['fna','ffn','frn','faa']
    try:
        res = dict([(s,[]) for s in seqtypes]) 
        dwl_folder = NCBI_ftp + GENBANK_draft_bacteria + info['genbank_folder']

        for file_type in seqtypes:
            try:
                rname = dwl_folder+"/"+info['WGS'][:4]+'00000000.contig.'+file_type+'.tgz'
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
                logger.error( 'Error in downloading  '+info['WGS']+": empty "+file_type )

    except Exception, e2: 
        logger.error("Something wrong in retrieving information for "+info['converted_assn']+": "+str(e2) )


def get_table_by_assn( remote_file, key, sep = '\t' ):
    table = urllib2.urlopen( remote_file )
    table = [sline.strip().split(sep) for sline in table] 
    table_header, table_content = [v.replace('#','') for v in table[0]],table[1:] 
    ret_table = dict([(d[key],d) for d in [dict(zip(table_header,l)) for l in table_content]
                if key in d])
    return ret_table

def get_assn2folder( remote_dir ):
    ftp = FTP(NCBI_ftp)
    ftp.login()
    ftp.cwd( remote_dir )
    dirs = ftp.nlst('*')
    return dict( (d.split("/")[::-1] for d in dirs if "/" in d) )
    
def conv_ass_format( s ):
    s = s.replace('GCA','GCF')
    return s.split(".")[0]

def get_assn2folder_genbank( remote_dir ):
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

if __name__ == '__main__':
    par = read_params(sys.argv)
   
    logger.info('Downloading and reading the NCBI microbial assembly list')
    ass_list = get_table_by_assn( add_protocol(NCBI_ftp + NCBI_prokaryotes_file), key = 'Assembly Accession' )
    logger.info('NCBI assembly list read: '+str(len(ass_list))+' assemblies')
    
    logger.info('Downloading and reading the GENBANK draft bacterial genomes list')
    uid2genbankfolder = get_assn2folder_genbank( GENBANK_draft_bacteria )
    logger.info('NCBI assembly list read: '+str(len(ass_list))+' assemblies')
    
    logger.info('Downloading and reading the NCBI microbial assembly folders')
    ass2folder = get_assn2folder( NCBI_assembly_folder )
    logger.info('NCBI assembly folder scanned: '+str(len(ass2folder))+' folders')
    
    logger.info('Reading the taxonomy from '+par['taxonomy']+'... ')
    with open(par['taxonomy']) as inpf:
        taxids2taxonomy = dict([l.strip().split('\t')[1:] for l in inpf])
    logger.info('Done.')

    for assn, info in ass_list.items():
        convass = conv_ass_format(assn)
        ass_list[assn]['converted_assn'] = convass
        if convass in ass2folder:
            ass_list[assn]['dwl_folder'] = ass2folder[conv_ass_format(assn)]+"/"+convass
        if 'TaxID' in info:
            if info['TaxID'] in taxids2taxonomy:
                ass_list[assn]['taxonomy'] = taxids2taxonomy[info['TaxID']]
        if 'taxonomy' not in ass_list[assn]:
            logger.warning(info['TaxID']+" is without taxonomic info and it is thus removed")
            del ass_list[assn]

    nproc = par['nproc'] 
   
    logger.info('Initializating the pool of '+str(nproc)+" downloaders")
    pool = mp.Pool( nproc )

    if not os.path.exists( par['out_dir'] ):
        os.mkdir( par['out_dir'] )
        logger.warning(par['out_dir']+" does not exist. Creating it.")

    for assn, info in ass_list.items():
        if 'dwl_folder' not in info:
            if info['BioProject ID'] in uid2genbankfolder:
                logger.warning(assn+" ("+info['Organism/Name']+" ) does not have an assembly folder, will try to download it from Genbank")
                info['genbank_folder'] = uid2genbankfolder[info['BioProject ID']]
                info['outdir'] = par['out_dir'] 
                info['assn'] = assn
                pool.map_async( draft_genome_genbank_dwl, [info] )
            else:
                logger.warning(assn+" ("+info['Organism/Name']+" ) does not have an assembly folder, and it will not be retrieved")
            continue
        if info['WGS'] == '-':
            logger.info("FINAL")
            info['outdir'] = par['out_dir']
            info['assn'] = assn
            pool.map_async( final_genome_dwl, [info] )
        else:
            logger.info("DRAFT")
            info['outdir'] = par['out_dir'] 
            info['assn'] = assn
            pool.map_async( draft_genome_dwl, [info] ) 

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






