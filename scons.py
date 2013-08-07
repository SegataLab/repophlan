import sys
import sfleoo
import collections
import matplotlib
import numpy as np
import os
import glob
import gzip

Import( "*" )
Decider( "MD5-timestamp" )
# EndImport
#sfleoo = oo interface for sfle
oo = sfleoo.ooSfle( fileDirInput = "input", fileDirOutput = "output", fileDirTmp = "tmp" )

# ===== GET TAXONOMY AND REFSEQ GENOMES ================
i_ncbi_nodes = oo.fin( "ncbi/nodes.dmp" ) #specify file used initially by the pipeline
i_ncbi_names = oo.fin( "ncbi/names.dmp" )
i_refseq_cat = oo.fin( "refseq/RefSeq-release60.catalog" )
i_prokaryotes = oo.fin( "ncbi/prokaryotes.txt" )
i_eukaryotes = oo.fin( "ncbi/eukaryotes.txt" )
i_viruses = oo.fin( "ncbi/viruses.txt" )

if_refseq = "refseq/"

o_tax = oo.fout( "full_tax.txt" ) #specify output files (not necessarily at end of pipeline)
o_red_tax = oo.fout( "full_tax.red.txt" )
#temporary file (not final output)
t_scaf_all = oo.ftmp( "scaffolds.txt" )

#ex = execute. this take three basic paramaters: list of input files, list of output, 
#script itself. can be a python script or linux script. anything you can calll from linux command 
#line. everything following the script name is just the arguments pertaining to the script.
#notice that second oo.ex uses output of first as its own input
oo.ex( [], t_scaf_all, "src/scafs.py", refseq_dir = "refseq/", output = t_scaf_all )
oo.ex( [i_ncbi_nodes, i_ncbi_names, i_refseq_cat, t_scaf_all], [o_tax, o_red_tax], 
        "src/generate_taxonomy.py", 
        nodes = i_ncbi_nodes, names = i_ncbi_names, 
        catalog = i_refseq_cat, prokaryotes = i_prokaryotes,
        scaffs_in_complete = t_scaf_all,
        output = o_tax, output_red = o_red_tax ) 

#oo.glob is like "ls" command linux + the pipeline automatically recognizes them as input files
for fna_in in oo.glob( if_refseq+"complete.*genomic.fna" ):
    gbf = oo.fin( "refseq/"+".".join(os.path.basename( fna_in ).split(".")[:2]) + ".genomic.gbff.gz" )
    if os.path.exists( gbf ):
        #oo.rebase = substitution
        out = oo.fout( "genomes/" + oo.rebase( fna_in, "fna", "log" ) )
        #get_genome doesn't know the names of the three files a priori, so we create a log file after
        #each run with a predictable name
        oo.ex( [gbf, fna_in, o_red_tax], out, "src/get_genome.py",
               gb = gbf, fna = fna_in, taxonomy = o_red_tax,
               output = out )
        # make sfle aware of what the final target of the pipeline is
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
















