from Bio import SeqIO
from collections import defaultdict as ddict
import sys
import re
import bz2
from subprocess import call
from tempfile import NamedTemporaryFile
import argparse
import os
import subprocess as sb
import logging
import multiprocessing as mp

aa = ['Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val']
rrnas = {'5S':(100,120), '16S':(1450,1700), '23S':(2900,3500)}
prog = re.compile('nnnnnnnnnn+')

logging.basicConfig(level=logging.INFO, stream=sys.stdout, 
                    format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
logger = logging.getLogger(sys.argv[0])

def read_params(args):
    parser = argparse.ArgumentParser(description='')
    arg = parser.add_argument
    arg( '--nproc', default = 15, type = int, 
         help="The number of parallel download processes")
    arg( '--in_summary', required = True, type = str, 
         help="The input summary file")
    arg( '--out_summary', required = True, type = str, 
         help="The output summary file")
    return vars(parser.parse_args())

def get_102( faa_bz2 ):
    f = NamedTemporaryFile(delete=False)
    f2 = NamedTemporaryFile(delete=False)
    f3 = NamedTemporaryFile(delete=False)
    f.close()
    with open(f.name,"w") as outf:
        call(["bzcat",faa_bz2],stdout=outf)
    call(["hmmscan","--tblout",f2.name,"-o",f3.name,"pfam/102.hmm",f.name])

    p102 = set()
    with open(f2.name) as inp:
        for line in (l.strip().split() for l in inp if not l.startswith('#')):
            if float(line[5]) > 10.0:
                p102.add( line[1].split('.')[0] )
    os.remove( f.name )
    os.remove( f2.name )
    os.remove( f3.name )
    return max(1.0-(102-len(p102))*0.01,0.1)

def get_fna_score( fna_bz2 ):
    ncontigs, ngood, nbad, n10n = 0, 0, 0, 0
    with bz2.BZ2File(fna_bz2,"r") as inp:
        for s in SeqIO.parse(inp,"fasta"):
            ncontigs += 1
            chars = set(s.seq)
            nonnt = 0
            for c in chars:
                if c not in 'atcgATCG':
                    nonnt += s.seq.count(c)
                if c not in 'atcgnATCGN':
                    nbad += s.seq.count(c)

            ngood += (len(s)-nonnt)

            sl = str(s.seq).lower()
            hits = prog.findall(sl)
            if hits:
                n10n += len(hits)

    
    fna_score = float(ngood) / ( ngood + nbad + 10000 * (ncontigs - 1) + 10000 * n10n  )
    return fna_score

def get_rna_scores( frn_bz2 ):
    found_aa = ddict(int)
    rrna_scores = ddict(set)
    for r in rrnas:
        rrna_scores[r].add(0.0)
    with bz2.BZ2File(frn_bz2,"r") as inp:
        for s in SeqIO.parse(inp,"fasta"):
            for a in aa:
                if "tRNA-"+a in s.id:
                    found_aa[a] += 1
            for r in rrnas:
                if '_'+r+'_' in s.id:
                    l = len(s)
                    if rrnas[r][0] < l < rrnas[r][1]:
                        rrna_scores[r].add( 0.3  )
                    elif rrnas[r][0] * 0.5 < l:
                        rrna_scores[r].add( 0.2  )
                    else:
                        rrna_scores[r].add( 0.1  )
    trna_score = max(0.1,1.0-(20.0-len(found_aa))*0.1)
    rrna_score = 0.1 + sum([max(v) for v in rrna_scores.values()])
    return trna_score,rrna_score                

def get_scores( info ):
    logger.info('Processing '+info[0])
    
    lfna, lfaa, lfrn = info[1]['fna_lname'], info[1]['faa_lname'], info[1]['frn_lname']
    info[1]['score_fna'] = get_fna_score( lfna  ) if lfna else 0.0
    info[1]['score_faa'] = get_102( lfaa ) if lfaa else 0.0
    info[1]['score_trna'], info[1]['score_rrna'] = get_rna_scores( lfrn ) if lfrn else (0.0, 0.0)
    for f in ['score_fna','score_faa','score_trna','score_rrna']:
        info[1][f] = str(round(info[1][f],3))
    return info

if __name__ == '__main__':
    par = read_params(sys.argv)
    nproc = par['nproc'] 

    with open(par['in_summary']) as inp:
        lines = [l.strip().strip('#').split('\t') for l in inp] 
        header = lines[0]
        mic = dict([(l[0],dict(zip(header,l))) for l in lines[1:]])

    results = []
    pool = mp.Pool( nproc )
    for k,v in mic.items():
        pool.map_async( get_scores, [(k,v)], callback=results.append ) 
    pool.close()
    pool.join()
    print results

    for r in results:
        k,v = r[0]
        mic[k] = v

    with open(par['out_summary'],"w") as outf:
        for i,(k,v) in enumerate(mic.items()):
            if i == 0:
                outf.write( "\t".join(["#genome"] + list(sorted(v.keys()))) +"\n" )
            outf.write( "\t".join([k] + [v[kk] for kk in sorted(v.keys())]) +"\n" )
    
    #for res in results:
    #    res.get(timeout=1000)    

    #print res.get()
    #fna_score =get_fna_score( sys.argv[2] )
    #p102_score = get_102( sys.argv[3] )
    #trna_score, rrna_score = get_rna_scores( sys.argv[1] )
    #print trna_score, rrna_score, p102_score, fna_score


