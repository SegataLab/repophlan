#!/usr/bin/env python
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
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
import itertools

FTP_prot = "ftp://"
NCBI_ftp = "ftp.ncbi.nlm.nih.gov"
NCBI_assembly_folder = '/genomes/ASSEMBLY_BACTERIA/'
NCBI_assrep_folder = '/genomes/ASSEMBLY_REPORTS/Bacteria'
NCBI_prokaryotes_file = '/genomes/GENOME_REPORTS/prokaryotes.txt'
NCBI_ASREFSEQ_file = '/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt'
NCBI_ASGENBANK_file = '/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt'
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

def is_plasmid( ppt ):
    if 'plasmid' in ppt[0]:
        return True
    return False

def is_contig( ppt ):
    if 'plasmid' in ppt[0]:
        return False
    if 'chromosome' in ppt[0]:
        return True
    if 'complete genome' in ppt[0]:
        return True
    return False

def assrep_entrez_dwl( info ):
    seqtypes = ['fna','ffn','frn','faa']
    try:
        res = dict([(s,[]) for s in seqtypes]) 

        ret = list(get_remote_file_with_retry( add_protocol( info['assrep'] ), "FINAL"  ) )

        todown = []
        for line in (r.strip().split('\t') for r in ret if r.strip()[0] != '#'):
            if line[3] == 'Chromosome':
                todown.append( line[4]  )

        for ent in todown:

            logger.info( "Downloading from ENTREZ: eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id="+ent+"&rettype=gb" )
            
            retgb = get_remote_file_with_retry( add_protocol( "eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id="+ent+"&rettype=gb" ), "FINAL"  )
            gb = SeqIO.parse(retgb, "genbank")



            #print "++++++++++++++",ent
            for sr in gb:
                #print "===========",ent,sr.id
                gi = sr.annotations['gi']
                ref = sr.annotations['accessions'][0]
                source = sr.annotations['source']
            
                sr.id = "|".join(["gi",gi,"ref",ref,""])
                sr.description = source
                res['fna'].append( sr )
                for feat in sr.features:
                    if feat.type == 'CDS':
                        prot = SeqRecord( Seq(feat.qualifiers['translation'][0])  )
                        protid = feat.qualifiers['protein_id'][0]
                        protgi = feat.qualifiers['db_xref'] [0]
                        prot.id = "|".join(["gi",protgi,"ref",protid,""])
                        prot.description = feat.qualifiers['product'][0]
                        res['faa'].append( prot )
                        fr,to = str(feat.location).split("]")[0].split("[")[1].split(":")
                        sign = str(feat.location).split(")")[0].split("(")[1]
                        fr = fr.replace("<","").replace(">","")
                        to = to.replace("<","").replace(">","")
                        locs = fr+"-"+to if sign == '+' else "c"+to+"-"+fr
                        gene = feat.location.extract(sr)
                        gene.id = "|".join(["gi",gi,"ref",ref,":"+locs])
                        gene.description = source
                        res['ffn'].append( gene )
                    if feat.type == 'tRNA' or feat.type == 'rRNA':
                        frn_l = feat.location.extract(sr)
                        fr,to = str(feat.location).split("]")[0].split("[")[1].split(":")
                        sign = str(feat.location).split(")")[0].split("(")[1]
                        fr = fr.replace("<","").replace(">","")
                        to = to.replace("<","").replace(">","")
                        locs = fr+"-"+to if sign == '+' else fr+"-"+to
                        frn_l.id = "|".join(["gi",gi,":"+locs,feat.qualifiers['product'][0],""])
                        frn_l.description = "["+feat.qualifiers['locus_tag'][0]+"]"
                        res['frn'].append( frn_l )

        for file_type in seqtypes:
            if len(res[file_type]):
                out_fna = info['outdir']+"/"+file_type+"/"+info['assn_nv']+"."+file_type
                if not os.path.exists( info['outdir']+"/"+file_type+"/"):
                    os.mkdir( info['outdir']+"/"+file_type+"/" )
                SeqIO.write( res[file_type], out_fna, "fasta")
            else:
                logger.error("Empty "+file_type+" for "+info['assn_nv'] )


    except Exception, e2: 
        logger.error("Something wrong in retrieving information for "+info['converted_assn']+" using ENTREZ: "+str(e2) )
def full_genome_dwl( info ):
    seqtypes = ['fna','ffn','frn','faa']
    seqtypes_plasmids = ['fna','ffn','faa']
    try:
        res = dict([(s,[]) for s in seqtypes]) 
        res_plasmids = dict([(s,[]) for s in seqtypes_plasmids]) 
        #dwl_folder = NCBI_ftp + NCBI_assembly_folder + info['dwl_folder']

        for chid,chfiles in info['chfin'].items():
            ret = list(get_remote_file_with_retry( add_protocol( chfiles['ptt'] ), "FINAL"  ) )
            if is_contig( ret ):
                for file_type in seqtypes:
                    try:
                        ret = get_remote_file_with_retry( add_protocol( chfiles[file_type] ), "FINAL"  ) 
                    except Exception, e:
                        logger.error('Error in downloading of '+chfiles[file_type]+' '+str(e))
                        break
                    res[file_type] += list(SeqIO.parse(ret,'fasta'))
            elif is_plasmid( ret ):
                for file_type in seqtypes_plasmids:
                    try:
                        ret = get_remote_file_with_retry( add_protocol( chfiles[file_type] ), "FINAL"  ) 
                    except Exception, e:
                        logger.error('Error in downloading of '+chfiles[file_type]+' '+str(e))
                        break
                    res_plasmids[file_type] += list(SeqIO.parse(ret,'fasta'))

        for file_type in seqtypes:
            if len(res[file_type]):
                out_fna = info['outdir']+"/"+file_type+"/"+info['assn_nv']+"."+file_type
                if not os.path.exists( info['outdir']+"/"+file_type+"/"):
                    os.mkdir( info['outdir']+"/"+file_type+"/" )
                SeqIO.write( res[file_type], out_fna, "fasta")
            else:
                logger.error("Empty "+file_type+" for "+info['assn_nv'] )
        
        for file_type in seqtypes_plasmids:
            if len(res_plasmids[file_type]):
                out_fna = info['outdir']+"/"+file_type+"/plasmids/"+info['assn_nv']+"."+file_type
                if not os.path.exists( info['outdir']+"/"+file_type+"/plasmids/"):
                    os.mkdir( info['outdir']+"/"+file_type+"/plasmids/" )
                SeqIO.write( res_plasmids[file_type], out_fna, "fasta")
            else:
                logger.error("Empty plasmid "+file_type+" for "+info['assn_nv'] )

        """
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
        """
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



def wgs_genome_dwl( info ):
    seqtypes = ['fna','ffn','frn','faa']
    try:
        res = dict([(s,[]) for s in seqtypes]) 
        dwl_folder = NCBI_ftp + NCBI_assembly_folder + info['dwl_folder']

        for file_type in seqtypes:
            try:
                rname = dwl_folder+"/NZ_"+info['wgs_master'].split(".")[0]+'.scaffold.'+file_type+'.tgz'
                ret = get_remote_file_with_retry( add_protocol( rname ), "DRAFT", compr_seq = True   )
            except Exception, e:
                logger.error('Error in downloading of '+rname+' '+str(e))
                break
            res[file_type] += ret

        for file_type in seqtypes:
            #print res
            if len(res[file_type]):
                out_fna = info['outdir']+"/"+file_type+"/"+info['assn_nv']+"."+file_type
                if not os.path.exists( info['outdir']+"/"+file_type+"/"):
                    os.mkdir( info['outdir']+"/"+file_type+"/" )
                SeqIO.write( res[file_type], out_fna, "fasta")
            else:
                rname = dwl_folder+"/NZ_"+info['wgs_master'].split(".")[0]+'.scaffold.'+file_type+'.tgz'
                logger.error("Empty "+file_type+" for "+info['assn_nv']+" ["+rname+"]")

    except Exception, e2: 
        logger.error("Something wrong in retrieving information for "+info['assn_nv']+": "+str(e2) )
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
                out_fna = info['outdir']+"/"+file_type+"/"+info['assn_nv']+"."+file_type
                if not os.path.exists( info['outdir']+"/"+file_type+"/"):
                    os.mkdir( info['outdir']+"/"+file_type+"/" )
                SeqIO.write( res[file_type], out_fna, "fasta")
            else:
                print "empty "+file_type

    except Exception, e2: 
        logger.error("Something wrong in retrieving information for "+info['converted_assn']+": "+str(e2) )

def draft_genbank_dwl( info ):
    seqtypes = ['fna','ffn','frn','faa']
    try:
        res = dict([(s,[]) for s in seqtypes]) 
        dwl_folder = NCBI_ftp + GENBANK_draft_bacteria 

        for file_type in seqtypes:
            try:
                rname = dwl_folder+"/"+info['gb_draft'][file_type]
                ret = get_remote_file_with_retry( add_protocol( rname ), "DRAFT", compr_seq = True   )
            except Exception, e:
                logger.error('Error in downloading of '+rname+' '+str(e))
                break
            res[file_type] += ret

        for file_type in seqtypes:
            if len(res[file_type]):
                out_fna = info['outdir']+"/"+file_type+"/"+info['assn_nv']+"."+file_type
                if not os.path.exists( info['outdir']+"/"+file_type+"/"):
                    os.mkdir( info['outdir']+"/"+file_type+"/" )
                #print "writing on "+out_fna
                SeqIO.write( res[file_type], out_fna, "fasta")
            else:
                logger.error( 'Error in downloading  '+info['assn_nv']+": empty "+file_type )

    except Exception, e2: 
        logger.error("Something wrong in retrieving information for "+info['converted_nv']+": "+str(e2) )

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
    table_header, table_content = [v.replace('#','').strip().lower() for v in table[0]],table[1:] 
    ret_table = dict([(d[key],d) for d in [dict(zip(table_header,l)) for l in table_content]
                if key in d])
    return ret_table

def get_assn2folderd( remote_dir ):
    ftp = FTP(NCBI_ftp)
    ftp.login()
    ftp.cwd( remote_dir )
    dirs = ftp.nlst('*/*')
   
    assn2files = {}
    for d in dirs:
        if ".scaffold." not in d:
            continue
        fold,fn = d.split("/") 
        assid = fn.split(".")[0]
        if not assid in assn2files:
            assn2files[assid] = {}
        ft = fn.split(".")[-2]
        assn2files[assid][ft] = d 
    return assn2files


def get_assn2folder( remote_dir ):
    ftp = FTP(NCBI_ftp)
    ftp.login()
    ftp.cwd( remote_dir )
    dirs = ftp.nlst('*/*')
    #print dirs
    #sys.exit()

    assn2folder = {}
    assn2files = {}
    for d in dirs:
        org,ass,fil = d.split("/")
        assn2folder[ass] = org
        fn = fil.split(".")[0]
        if ass not in assn2files:
            assn2files[ass] = {}
        if fn not in assn2files[ass]:
            assn2files[ass][fn] = {}
        assn2files[ass][fn][fil.split(".")[1]] = NCBI_ftp +"/" + remote_dir  + "/"+ d 

    #return dict( (d.split("/")[::-1] for d in dirs if "/" in d) )
    return assn2folder, assn2files

def get_assn2assrep( remote_dir ):
    ftp = FTP(NCBI_ftp)
    ftp.login()
    ftp.cwd( remote_dir )
    dirs = ftp.nlst('*/latest')
    return dict( ((".".join(d.split("/")[2].split(".")[:2]).replace('GCA','GCF'),NCBI_ftp+remote_dir+"/"+d) for d in dirs if d.endswith("assembly.txt")) )
    
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

    # ENSEMBLE ???? ftp://ftp.ensemblgenomes.org/pub/current/species.txt

    par = read_params(sys.argv)
    nproc = par['nproc'] 
  
    logger.info('Retriving the pairs assembly_id - assembly info from '+NCBI_assrep_folder)
    ass2assrep = get_assn2assrep( NCBI_assrep_folder )
    logger.info('Completed reading '+str(len(ass2assrep))+'assembly_id - assembly info from '+NCBI_assrep_folder)

    
    logger.info('Retriving the pairs assembly_id remote draft scaffold files from '+GENBANK_draft_bacteria)
    ass2folderd = get_assn2folderd( GENBANK_draft_bacteria )
    logger.info('Completed reading '+str(len(ass2folderd))+'assembly_id - remote draft scaffold files from '+ GENBANK_draft_bacteria)

    
    logger.info('Downloading and reading the NCBI REFSEQ microbial assembly list')
    ass_list = get_table_by_assn( add_protocol(NCBI_ftp + NCBI_ASREFSEQ_file), key = 'assembly_id' )
    logger.info('NCBI assembly list read: '+str(len(ass_list))+' assemblies')
    
    logger.info('Downloading and reading the NCBI GENBANK microbial assembly list')
    ass_list_gb = get_table_by_assn( add_protocol(NCBI_ftp + NCBI_ASGENBANK_file), key = 'assembly_id' )
    logger.info('NCBI GENBANK assembly list read: '+str(len(ass_list))+' assemblies')

    logger.info('Downloading and reading the NCBI microbial assembly folders')
    ass2folder, assn2files = get_assn2folder( NCBI_assembly_folder )
    logger.info('NCBI assembly folder scanned: '+str(len(ass2folder))+' folders')
    
    logger.info('Reading the taxonomy from '+par['taxonomy']+'... ')
    with open(par['taxonomy']) as inpf:
        taxids2taxonomy = dict([l.strip().split('\t')[1:] for l in inpf])
    logger.info('Done.')

    if not os.path.exists( par['out_dir'] ):
        os.mkdir( par['out_dir'] )
        logger.warning(par['out_dir']+" does not exist. Creating it.")
    
    
    logger.info('Initializating the pool of '+str(nproc)+" downloaders for refseq genomes download")
    pool = mp.Pool( nproc )
    
    for assn, info in ass_list.items():
        
        if info['version_status'] == 'replaced':
            continue
        
        assn_nv = assn.split(".")[0]
        ass_list[assn]['assn_nv'] = assn_nv
        ass_list[assn]['outdir'] = par['out_dir']
        
        if 'taxid' in info:
            if info['taxid'] in taxids2taxonomy:
                ass_list[assn]['taxonomy'] = taxids2taxonomy[info['taxid']]
        if 'taxonomy' not in ass_list[assn]:
            logger.warning(info['taxid']+" is without taxonomic info and it is thus removed")
            del ass_list[assn]
            continue
        #if 'dwl_folder' in info:
        if assn_nv in ass2folder:
            ass_list[assn]['dwl_folder'] = ass2folder[conv_ass_format(assn)]+"/"+assn_nv 
            if 'wgs_master' in info and info['wgs_master'] != 'na':
                pool.map_async( wgs_genome_dwl, [info] )
            else:
                info['chfin'] =  assn2files[info['assn_nv']]
                pool.map_async( full_genome_dwl, [info] )
    
    pool.close()
    pool.join()
    
    logger.info('Initializating the pool of '+str(nproc)+" downloaders for GenBank genomes download")
    pool = mp.Pool( nproc )
    
    for assn, info in ass_list_gb.items():
        if info['version_status'] == 'replaced':
            continue
        if 'taxid' in info:
            if info['taxid'] in taxids2taxonomy:
                ass_list_gb[assn]['taxonomy'] = taxids2taxonomy[info['taxid']]
        if 'taxonomy' not in ass_list_gb[assn]:
            logger.warning(info['taxid']+" is without taxonomic info and it is thus removed")
            del ass_list_gb[assn]
            continue
        assn_nv = assn.split(".")[0]
        ass_list_gb[assn]['assn_nv'] = assn_nv
        ass_list_gb[assn]['outdir'] = par['out_dir']
        if not os.path.exists( par['out_dir'] + "/fna/"+conv_ass_format(assn_nv)+".fna" ):
            if info['wgs_master'] != 'na' and info['wgs_master'].split(".")[0] in ass2folderd:
                ass_list_gb[assn]['gb_draft'] = ass2folderd[info['wgs_master'].split(".")[0]]
                pool.map_async( draft_genbank_dwl, [ass_list_gb[assn]] )
            elif info['wgs_master'] == 'na' and assn.replace('GCA','GCF') in ass2assrep:
                #print "todown " + ass2assrep[assn.replace('GCA','GCF')]
                ass_list_gb[assn]['assrep'] = ass2assrep[assn.replace('GCA','GCF')]
                pool.map_async( assrep_entrez_dwl, [ass_list_gb[assn]] )
        else:
            logger.info("Avoiding the download of "+assn+" because it was already sucessfully retrieved from RefSeq "+par['out_dir'] + "/fna/"+conv_ass_format(assn_nv)+".fna")

    pool.close()
    pool.join()

    
    seqtypes = ['fna','ffn','frn','faa']
    out_sum = []

    for assn, info in itertools.chain(ass_list.items(),ass_list_gb.items()):
        if 'assn_nv' not in info:
            continue
        all_files = [os.path.exists( par['out_dir']+"/"+st+"/"+info['assn_nv']+"."+st) for st in seqtypes]
        if all( all_files ):
            out_sum.append( [info['assn_nv'],info['taxonomy']] )
    with open(par['out_summary'],"w") as outf:
        for s in out_sum:
            outf.write( "\t".join(s) +"\n" )
            
