#!/usr/bin/env python
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import collections
import argparse 
import os

def read_params(args):
    parser = argparse.ArgumentParser(description='')
    arg = parser.add_argument
    arg( '--gb', metavar='gb', required = True, type=str,
         help="The (compressed) GenBank file")
    arg( '--fna', metavar='fna', required = True, type=str,
         help="The contig multi-fasta file")
    arg( '--taxonomy', metavar='taxonomy', required = True, type=str, 
         help="The taxonomy file")
    arg( '--output', metavar='out', default = None, type = str, 
         help="The output file (and path for the sequence files)")
    return vars(parser.parse_args())


def op( fn ):
    if fn.endswith( ".gz" ):
        import gzip
        return gzip.open(fn)

    return open( fn )

def run( pars ):
    #input_handle is name of genbank flatfile (annotation)
    #fna is a dictionary with fasta accessions pointing to the genomes
    input_handle  = pars['gb']
    fna = SeqIO.index( pars['fna'], "fasta")
    print fna.keys()
    #open taxonomy file. get all fields for each line. tids2accs is dictionary
    #for taxon id -> accession. if four-letter code, it's draft.
    #if not four letter code, it's final or draft. 
    tids2accs = {} 
    for name,status,accs,tid,tax in (l.strip().split('\t') for l in open( pars['taxonomy'] )):
        if len(accs) == 4:
            tids2accs[tid] = {'acc': [accs]}
        else:
            tids2accs[tid] = {'acc': accs.split( "," )}
    
    #create reverse dictionary (accessions to taxon ids)    
    accs2tids = {}
    for k,v in tids2accs.items():
        for vv in v['acc']:
            accs2tids[vv] = k
    
    #
    missed_ids = set()
    for seq_record in SeqIO.parse( op( input_handle ), "genbank") :
        for seq_record_id in [seq_record.id,seq_record.id[3:7]]:
            if seq_record_id in accs2tids:
                if 'genes' not in tids2accs[ accs2tids[seq_record_id] ]:
                    tids2accs[ accs2tids[seq_record_id] ]['genes'] = {} 
                    tids2accs[ accs2tids[seq_record_id] ]['gene_des'] = {} 
                    tids2accs[ accs2tids[seq_record_id] ]['proteins'] = [] 
                    tids2accs[ accs2tids[seq_record_id] ]['fna_genome'] = [] 
                    
                genes = tids2accs[ accs2tids[seq_record_id] ]['genes']
                gene_des = tids2accs[ accs2tids[seq_record_id] ]['gene_des']
                proteins = tids2accs[ accs2tids[seq_record_id] ]['proteins']
                fna_genome = tids2accs[ accs2tids[seq_record_id] ]['fna_genome']
        
                gi = seq_record.annotations['gi']
                n = 0
                id = "|".join(["gi",gi,"ref",seq_record.id,""])
		    try:
			    fna_g = fna[id]
                fna_genome.append( fna_g )
                for seq_feature in seq_record.features :
                    if seq_feature.type == "gene" :
                        gene = seq_feature.extract( fna_g  ) 
                        gene.id = seq_record.id + "_" + str(n)
                        gene.description = str(seq_feature.location)+" "
                        genes[seq_feature.qualifiers['locus_tag'][0]] = gene
                        if 'pseudo' in seq_feature.qualifiers:
                            gene.description += 'pseudo'+ (" "+seq_feature.qualifiers['gene'][0] if 'gene' in seq_feature.qualifiers else "")
                            gene.description += " ["+seq_record.annotations['organism']+"]" 
                        n += 1
                    elif seq_feature.type == "CDS":
                        idp = "|".join(["gi",gi,"ref",seq_feature.qualifiers['protein_id'][0],""])
                        name = seq_feature.qualifiers['product'][0]
                        name = name + (" " + seq_feature.qualifiers['locus_tag'][0] if name == "hypothetical protein" else "")
                        name += " ["+seq_record.annotations['organism']+"]"
                        record = SeqIO.SeqRecord(   Seq(seq_feature.qualifiers['translation'][0]),
                                                    idp,
                                                    '',
                                                    name)
                        proteins.append( record )
                        gene_des[ seq_feature.qualifiers['locus_tag'][0] ] = name
                    elif seq_feature.type == "rRNA":
                        name = seq_feature.qualifiers['product'][0] if 'product' in seq_feature.qualifiers else ""
                        name += " ["+seq_record.annotations['organism']+"]"
                        gene_des[ seq_feature.qualifiers['locus_tag'][0] ] = name
                    elif seq_feature.type == "tRNA":
                        name = seq_feature.qualifiers['product'][0] if 'product' in seq_feature.qualifiers else ""
                        name += " ["+seq_record.annotations['organism']+"]"
                        gene_des[ seq_feature.qualifiers['locus_tag'][0] ] = name
                    elif seq_feature.type == "ncRNA":
                        name = "ncRNA"+(" "+seq_feature.qualifiers['gene'][0] if 'gene' in seq_feature.qualifiers else "")
                        name += " ["+seq_record.annotations['organism']+"]"
                        gene_des[ seq_feature.qualifiers['locus_tag'][0] ] = name
            except:
			    print id
			    missed_ids.add(id)

    outdir = os.path.dirname( pars['output'] )
    if outdir:
        outdir += "/"

    with open( pars['output'], "w" ) as out:
        out.write( "Genomes extracted\n" )
        for k,v in tids2accs.items():
            if 'genes' not in v:
                continue
            all_genes = []
            for kk,vv in v['genes'].items():
                if kk in v['gene_des']:
                    vv.description += v['gene_des'][kk]
                all_genes.append( vv )
    
            stat = 'D'
            if all([('NC_' in g.id) for g in v['fna_genome']]):
                stat = 'F'

            SeqIO.write( all_genes, outdir + "t"+k+stat+".ffn", "fasta")
            SeqIO.write( v['proteins'], outdir +"t"+k+stat+".faa", "fasta")
            SeqIO.write( v['fna_genome'], outdir +"t"+k+stat+".fna", "fasta")
        
            out.write( k+"\n" )


if __name__ == '__main__':
    par = read_params(sys.argv)
    run( par )


