#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pickle
from Bio.Seq import Seq
import urllib2
import sys
import argparse
from ftplib import FTP
from StringIO import StringIO
from time import time, sleep
import gzip
import math
import os
import logging
import collections

NCBI_ftp = "ftp.ncbi.nlm.nih.gov"
#NCBI_assembly_folder = '/genomes/Fungi/'
#GENBANK_assembly_folder = 'genbank/genomes/Fungi/'
NCBI_eukaryotes_file = '/genomes/GENOME_REPORTS/eukaryotes.txt'
NCBI_fungi_assemblies = '/genomes/ASSEMBLY_REPORTS/Eukaryotes/fungi/'
NCBI_protozoa_assemblies = '/genomes/ASSEMBLY_REPORTS/Eukaryotes/protozoa/'

logging.basicConfig(level=logging.INFO, stream=sys.stdout, 
                    format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
logger = logging.getLogger(sys.argv[0])

def read_params(args):
    parser = argparse.ArgumentParser(description='')
    arg = parser.add_argument
    arg( '--refseq_fungi_folder', type=str,
         default = "/refseq/release/fungi/",
         help="The gbff RefSeq file for fungi (default is /refseq/release/fungi/)")
    arg( '--refseq_protozoa_folder', type=str,
         default = "/refseq/release/protozoa/",
         help="The gbff RefSeq file for protozoa (default is /refseq/release/protozoa/)")
    arg( '--taxonomy', metavar='taxonomy', required = True, type=str,
         help="The taxonomy file")
    arg( '--out_dir_fungi', required = True, default = None, type = str, 
         help="The output folder")
    arg( '--out_dir_protozoa', required = True, default = None, type = str, 
         help="The output folder")
    arg( '--out_summary_fungi', required = True, default = None, type = str, 
         help="The output summary file")
    arg( '--out_summary_protozoa', required = True, default = None, type = str, 
         help="The output summary file")
    return vars(parser.parse_args())

def get_table_by_assn( remote_file, key, sep = '\t' ):
    table = urllib2.urlopen( remote_file )
    table = [sline.strip().split(sep) for sline in table] 
    table_header, table_content = [v.replace('#','') for v in table[0]],table[1:] 
    ret_table = dict([(d[key].split(".")[0].replace( "GCF","GCA" ),d) for d in [dict(zip(table_header,l)) for l in table_content]
                if key in d])
    return ret_table

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
                    print e
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

def add_protocol( fn ):
    if fn.startswith("ftp"):
        return "ftp://"+fn
    return "http://"+fn

@retry(tries=8, delay=20, backoff=2)
def get_remote_gbff_with_retry( url ):
    loc = urllib2.urlopen( url, timeout = 200 )

    logger.info('Parsing of '+url)
    f = gzip.GzipFile(fileobj= StringIO( loc.read())) 
    return SeqIO.parse( f, "genbank")

@retry(tries=8, delay=20, backoff=2)
def get_assembly_info_with_retry( url ):
    loc = urllib2.urlopen( add_protocol( url ), timeout = 200 )

    logger.info('Parsing of '+url)

    refseq_contigs, genbank_contigs = [], []
    for l in loc:
        if l[0] == '#':
            continue
        line = l.strip().split('\t')
        if line[-2] and line[-2] != 'na':
            refseq_contigs.append( line[-2] )
        if line[-4] and line[-4] != 'na':
            genbank_contigs.append( line[-4] )
    return refseq_contigs, genbank_contigs

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
    if url.endswith( ".gz" ):
        f = gzip.GzipFile(fileobj= StringIO( loc.read()))
        return SeqIO.parse( f, "fasta")
    else:
        seqs = StringIO.StringIO(loc.read())
    return seqs

def read_gbff( info ):
    contigs2faa = collections.defaultdict( list )
    contigs2ffn = collections.defaultdict( list )

    try:
        ret = get_remote_gbff_with_retry( add_protocol( info['refseq_gbff'] ) )
    except Exception, e:
        logger.error('Error in downloading of '+info['refseq_gbff']+' '+str(e))
    
    for seq_record in ret:
        accession = seq_record.annotations['accessions'][0]
        if accession in info['accessions']:
            logger.info('Processing accession '+accession)

            for feat in seq_record.features:
                gGI, gTI, gPO = "", "", ""
                if feat.type == 'mRNA':
                    if 'transcript_id' in feat.qualifiers:
                        gTI = feat.qualifiers['transcript_id'][0]
                    if 'db_xref' in feat.qualifiers:
                        for q in feat.qualifiers['db_xref']:
                            if "GI:" in q:
                                gGI = q.split("GI:")[-1]
                    contigs2ffn[accession].append( "|".join(["gi",gGI,"ref",gTI,""]) )
                if feat.type == 'CDS':
                    if 'protein_id' in feat.qualifiers:
                        gPI = feat.qualifiers['protein_id'][0]
                    if 'db_xref' in feat.qualifiers:
                        for q in feat.qualifiers['db_xref']:
                            if "GI:" in q:
                                gGI = q.split("GI:")[-1]
                    contigs2faa[accession].append("|".join(["gi",gGI,"ref",gPI,""]) )
    return contigs2ffn, contigs2faa

def get_euks_assemblies( remote_dir ):
    ftp = FTP( NCBI_ftp )
    ftp.login()
    ftp.cwd( remote_dir )
    dirs = ftp.nlst('*/*')
    ret = {}
    for assembly in dirs:
        if 'representatives' not in assembly:
            continue
        if 'assembly' not in assembly:
            continue
        ass = ".".join(assembly.split("/")[-1].split(".")[:1])
        refseq, genbank = get_assembly_info_with_retry( NCBI_ftp + remote_dir + assembly )
        ret[ass] = {'refseq_contigs': refseq, 'genbank_contigs': genbank}

    return ret


def get_contigs_from_folders( remote_dir ):
    ftp = FTP( NCBI_ftp )
    ftp.login()
    ftp.cwd( remote_dir )
    dirs = ftp.nlst('*')
    contigs2files = {}
    for f in dirs:
        if not f.endswith(".fna"):
            continue
        contig = f.split("/")[-1].split(".")[0]
        contigs2files[contig] = f.split(".fna")[0]
    return contigs2files

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

def get_ffnfaa2contigs( folder  ):
    ftp = FTP( NCBI_ftp )
    ftp.login()
    ftp.cwd( folder )
    rem_files = ftp.nlst('*')
    
    gbffs, contigs = [], []

    for f in rem_files:
        if f.endswith(".genomic.gbff.gz"):
            try:
                ret = get_remote_gbff_with_retry( add_protocol( NCBI_ftp + folder + f ) )
            except Exception, e:
                logger.error('Error in downloading of '+info['refseq_gbff']+' '+str(e))
            gbffs.append( ret )

    for gbff in gbffs:
        for seq_record in gbff:
            accession = seq_record.annotations['accessions'][0] 
            if accession in contigs2assembly:
                logger.info('Processing accession '+accession+" ("+contigs2assembly[accession]+")")
                for feat in seq_record.features:
                    gGI, gTI, gPO = "", "", ""
                    if feat.type == 'mRNA':
                        if 'transcript_id' in feat.qualifiers:
                            gTI = feat.qualifiers['transcript_id'][0]
                        if 'db_xref' in feat.qualifiers:
                            for q in feat.qualifiers['db_xref']:
                                if "GI:" in q:
                                    gGI = q.split("GI:")[-1]
                        contigs.append( accession )
                        ffn2contigs[ "|".join(["gi",gGI,"ref",gTI,""]) ] = accession
                    if feat.type == 'CDS':
                        if 'protein_id' in feat.qualifiers:
                            gPI = feat.qualifiers['protein_id'][0]
                        if 'db_xref' in feat.qualifiers:
                            for q in feat.qualifiers['db_xref']:
                                if "GI:" in q:
                                    gGI = q.split("GI:")[-1]
                        faa2contigs[ "|".join(["gi",gGI,"ref",gPI,""]) ] = accession
    faas, ffns, fnas = collections.defaultdict(list), collections.defaultdict(list), collections.defaultdict(list)
    
    for f in rem_files:
        if f.endswith(".1.genomic.fna.gz"):
            try:
                ret = get_remote_file_with_retry( add_protocol( NCBI_ftp + folder + f ) )
            except Exception, e:
                logger.error('Error in downloading of '+ NCBI_ftp + folder + f  ) 
                break
            for s in ret:
                cn = s.id.split("|")[3].split('.')[0] 
                if cn in contigs:
                    fnas[contigs2assembly[cn]].append( s )
   
    for f in rem_files:
        if f.endswith(".rna.fna.gz"):
            try:
                ret = get_remote_file_with_retry( add_protocol( NCBI_ftp + folder + f ) )
            except Exception, e:
                logger.error('Error in downloading of '+ NCBI_ftp + folder + f  ) 
                break
            for s in ret:
                if s.id in ffn2contigs:
                    ffns[contigs2assembly[ffn2contigs[s.id]]].append( s )
        if f.endswith(".protein.faa.gz"):
            try:
                ret = get_remote_file_with_retry( add_protocol( NCBI_ftp + folder + f ) )
            except Exception, e:
                logger.error('Error in downloading of '+ NCBI_ftp + folder + f  ) 
                break
            for s in ret:
                if s.id in faa2contigs:
                    faas[contigs2assembly[faa2contigs[s.id]]].append( s )
    return fnas, ffns, faas

"""
def euk_genome_dwl( info ):
    seqtypes = ['fna','ffn','faa']
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
"""

if __name__ == '__main__':
    par = read_params(sys.argv)
    
    logging.basicConfig(level=logging.INFO, stream=sys.stdout, 
                        format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logger = logging.getLogger(sys.argv[0])
    
    logger.info('Downloading and reading the NCBI eukaryotes assembly list')
    ass_list = get_table_by_assn( add_protocol(NCBI_ftp + NCBI_eukaryotes_file), key = 'Assembly Accession' )
    logger.info('NCBI assembly list read: '+str(len(ass_list))+' assemblies')
    
    logger.info('Reading the taxonomy from '+par['taxonomy']+'... ')
    with open(par['taxonomy']) as inpf:
        taxids2taxonomy = dict([l.strip().split('\t')[1:] for l in inpf])
    logger.info('Done.')
    
    #logger.info('Downloading and reading the NCBI fungal contigs')
    #ncbi_contigs = get_contigs_from_folders( NCBI_assembly_folder )
    #logger.info('NCBI assembly folder scanned: '+str(len(ncbi_contigs))+' contigs')
    
    #logger.info('Downloading and reading the GenBank fungal contigs')
    #genbank_contigs = get_contigs_from_folders( GENBANK_assembly_folder )
    #logger.info('NCBI assembly folder scanned: '+str(len(genbank_contigs))+' contigs')
   
    ass2contigs_protozoa = get_euks_assemblies( NCBI_protozoa_assemblies ) 
    ass2contigs_fungi = get_euks_assemblies( NCBI_fungi_assemblies ) 
    ass2contigs = dict(ass2contigs_fungi.items() + ass2contigs_protozoa.items())

    #output = open('data.pkl', 'wb')
    #pickle.dump(ass2contigs, output)
    #ass2contigs = pickle.load(open('data.pkl'))

    to_get_from_refseq = {}
    contigs2assembly, contigs2ffn, contigs2faa = {}, collections.defaultdict(list), collections.defaultdict(list) 
    ffn2contigs, faa2contigs = {}, {}
    for k,v in ass2contigs.items():
        k = k.replace("GCF","GCA")
        #print "PREPROC "+k, k in ass_list
        if k in ass_list and ass_list[k]['TaxID'] in taxids2taxonomy:
            if v['refseq_contigs']:
                to_get_from_refseq[k] = v['refseq_contigs']
                for c in v['refseq_contigs']:
                    contigs2assembly[c.split('.')[0]] = k
                logger.info( "Will try to retrieve "+k+" ("+taxids2taxonomy[ass_list[k]['TaxID']]+") from RefSeq" ) 
            """
            else:
                print "PROC "+k
                print "genbank_contigs in ncbi_contigs", v['genbank_contigs'], ncbi_contigs
                available_contigs = sum([int(c.split(".")[0] in ncbi_contigs) for c in v['genbank_contigs']])
                print k, "ncbi", available_contigs, len(v['genbank_contigs'])
                print "genbank_contigs in genbank_contigs", v['genbank_contigs'], genbank_contigs
                available_contigs = sum([int(c.split(".")[0] in genbank_contigs) for c in v['genbank_contigs']])
                print k, "genbank", available_contigs, len(v['genbank_contigs'])
            """
   

    for fold in ['refseq_protozoa_folder','refseq_fungi_folder']:
        fnas, ffns, faas = get_ffnfaa2contigs( par[fold] )
       
        outdir = par['out_dir_fungi'] if fold == 'refseq_fungi_folder' else par['out_dir_protozoa']
        if not os.path.exists( outdir ):
            os.mkdir( outdir )
            logger.warning( outdir+" does not exist. Creating it.")
    
        outsum = collections.defaultdict( int )

        logger.info( "Saving fnas, ffn, faas, into "+outdir )
        for k,v in fnas.items():
            outsum[k] += 1
            SeqIO.write( v, outdir+"/"+k+".fna", "fasta" )
        for k,v in ffns.items():
            outsum[k] += 1
            SeqIO.write( v, outdir+"/"+k+".ffn", "fasta" )
        for k,v in faas.items():
            outsum[k] += 1
            SeqIO.write( v, outdir+"/"+k+".faa", "fasta" )
       
        outfn = par['out_summary_fungi'] if fold == 'refseq_fungi_folder' else par['out_summary_protozoa']
        logger.info( "Writing summary into "+outfn )
        with open(outfn,"w") as outf:
            for k,v in outsum.items():
                outf.write( "\t".join([str(k),str(taxids2taxonomy[ass_list[k]['TaxID']])]) +"\n" )




            #print taxids2taxonomy[ass_list[k]['TaxID']] 
    #print read_gbff( {'refseq_gbff': 'ftp.ncbi.nlm.nih.gov/refseq/release/fungi/fungi.15.genomic.gbff.gz', 'accessions':['NC_007194']}  )
