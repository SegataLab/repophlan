import sfleoo
import sys
import collections
import matplotlib
import numpy as np
import os
import glob
import gzip

Import( "*" )
Decider( "MD5-timestamp" )
# EndImport

oo = sfleoo.ooSfle( fileDirInput = "input", fileDirOutput = "output", fileDirTmp = "tmp" )

# ===== GET TAXONOMY AND REFSEQ GENOMES ================
i_ncbi_nodes = oo.fin( "ncbi/nodes.dmp" )
i_ncbi_names = oo.fin( "ncbi/names.dmp" )
i_refseq_cat = oo.fin( "refseq/RefSeq-release59.catalog" )
i_prokaryotes = oo.fin( "ncbi/prokaryotes.txt" )

if_refseq = "refseq/"

o_tax = oo.fout( "full_tax.txt" )
o_red_tax = oo.fout( "full_tax.red.txt" )

t_scaf_all = oo.ftmp( "scaffolds.txt" )

oo.ex( [], t_scaf_all, "src/scafs.py", refseq_dir = "input/refseq/", output = t_scaf_all )
oo.ex( [i_ncbi_nodes, i_ncbi_names, i_refseq_cat, t_scaf_all], [o_tax, o_red_tax], 
        "src/generate_taxonomy.py", 
        nodes = i_ncbi_nodes, names = i_ncbi_names, 
        catalog = i_refseq_cat, prokaryotes = i_prokaryotes,
        scaffs_in_complete = t_scaf_all,
        output = o_tax, output_red = o_red_tax ) 

for fna_in in oo.glob( if_refseq+"complete.*genomic.fna" ):
    gbf = oo.fin( "refseq/"+".".join(os.path.basename( fna_in ).split(".")[:2]) + ".genomic.gbff.gz" )
    if os.path.exists( gbf ):
        out = oo.fout( "genomes/" + oo.rebase( fna_in, "fna", "log" ) )
        oo.ex( [gbf, fna_in, o_red_tax], out, "src/get_genome.py",
               gb = gbf, fna = fna_in, taxonomy = o_red_tax,
               output = out )
        Default( out )

for fna_in in oo.glob( if_refseq+"completeNZ_*genomic.fna" ):
    gbf = oo.fin( "refseq/"+".".join(os.path.basename( fna_in ).split(".")[:1]) + ".genomic.gbff.gz" )
    if os.path.exists( gbf ):
        out = oo.fout( "genomes/" + oo.rebase( fna_in, "fna", "log" ) )
        oo.ex( [gbf, fna_in, o_red_tax], out, "src/get_genome.py",
               gb = gbf, fna = fna_in, taxonomy = o_red_tax,
               output = out )
        Default( out )
Default( o_tax, o_red_tax )
# ===== GET TAXONOMY AND REFSEQ GENOMES ================
















