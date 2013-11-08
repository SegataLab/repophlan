This is pretty much a work in progress. The pipeline currently consists of the
following three scripts (reported with default settings):

$ ./generate_taxonomy.py --output taxonomy.txt --output_red taxonomy_reduced.txt --pickle taxonomy.pkl | tee generate_taxonomy.log

[ ~3 hours ]. This script downloads, processes, and saves the NCBI taxonomy in 
a standard format with a fixed number of taxonomic levels.


$ ./repophlan_get_viruses.py --taxonomy taxonomy.txt --out_dir repophlan_viruses --out_summary repophlan_viruses.txt | tee repophlan_viruses.log
[ ~1 hour ]. Downloads and save all sequences available for viruses (from
RefSeq). The files are saved in the 'repophlan_viruses' folder and the taxonomy
of each downloaded set of files is in 'repophlan_viruses.txt'

$ ./repophlan_get_microbes.py --taxonomy taxonomy_reduced.txt --out_dir microbes --nproc 20 --out_summary repophlan_microbes.txt | tee repophlan_microbes.log
[ ~6 hours ]. Downloads and save all available sequences for microbes. It uses
multiple parallel connections to speed up the retrieval: common exceptions 
including issues with the NCBI ftp which temporarily rejects connections for
ftp load problems are handled gracefully with a 'delay and retry' policy. 
IMPORTANT: specifying more than 20 processors to run in parallel causes serious
problems in terms of exceeding the allowed number of connections by NCBI causing
long delays.

As of Nov 8 2013, RepoPhlAn retrieves:
* 4958 viruses (each with fna, ffn, faa files)
* 12277 microbes (each with fna, ffn, frn, faa files)

For microbes, four files are generated for each genome:
* .fna : the genome in multifasta format (one or more contigs)
* .ffn : all the protein coding genes (no rRNA, tRNA)
* .faa : all the proteins
* .frn : all the non-coding genes (rRNA, tRNA)

Several aspects can be optimized and with a bit of refining about 20% more 
microbes could be retrieved. 


Known issues:

* Eukaryotes are missing. Need to find a principled means to get them from NCBI
* For almost all the assembly additional informative file can be automatically
  downloaded:
  - .asn (probably not useful)
  - .gbk (genbank file for the assembly and for each scaffold)
  - .gff (annotations in almost free text format)
  - .ptt (tab-separated file including COG assignments and product names)
  - .rnt (annotation of rRNA and tRNAs)
  - .rpt (few metadata)
  - .val (binary file, not sure what's inside)
* about ~30 strains are missing in the generated taxonomy file. These are
  the representatives of sets of strains with multiple assemblies. It should
  be an easy fix in generate_taxonomy.py
* for some microbes only a subset of the four types of files are present. These
  files are downloaded but the assembly will not appear in the 
  repophlan_microbes.txt for consistency. However, for some downstream analyses
  not all the four files are actually needed.
* several inconsistencies are present in the following file:
  ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt
  entries with inconsistencies are skipped (with a 'try catch' policy), but few
  of these inconsistencies could actually be handled.
* need extensive testing of both the taxonomy and the retrieved files






