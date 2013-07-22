#!/usr/bin/env python
import sys
import collections
import re
import copy 
import argparse

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree as BTree
from Bio.Phylo.BaseTree import Clade as BClade


def read_params(args):
    parser = argparse.ArgumentParser(description='')
    arg = parser.add_argument
    arg( '--nodes', metavar='nodes', required = True, type=str,
         help="The nodes.dmp from the NCBI taxonomy")
    arg( '--names', metavar='names', required = True, type=str,
         help="The names.dmp from the NCBI taxonomy")
    arg( '--prokaryotes', metavar='prokaryotes', required = True, type=str, 
         help="The prokaryotes NCBI file \"prokaryotes.txt\"")
    arg( '--eukaryotes', metavar='eukaryotes', required = True, type=str, 
         help="The eukaryotes NCBI file \"eukaryotes.txt\"")
    arg( '--viruses', metavar='viruses', required = True, type=str, 
         help="The viruses NCBI file \"viruses.txt\"")
    arg( '--catalog', metavar='catalog', required = True, type=str, 
         help="The catalog RefSeq file \"RefSeq-releaseXX.catalog\"")
    arg( '--scaffs_in_complete', metavar='scaffs_in_complete', required = True, type=str, 
         help="The scaffolds entries for draft genomes included in the \"complete\" files")
    arg( '--output', metavar='out', default = None, type = str, 
         help="The output taxonomy")
    arg( '--output_red', metavar='out_red', default = None, type = str, 
         help="The output taxonomy with reduced and fixed number of taxonomic levels")

    return vars(parser.parse_args())




class Nodes:
    #
    # Format of nodes.dmp from RefSeq documentation
    #
    # ---------
    # 
    # This file represents taxonomy nodes. The description for each node includes 
    # the following fields:
    # 
    # 	tax_id					-- node id in GenBank taxonomy database
    #  	parent tax_id				-- parent node id in GenBank taxonomy database
    #  	rank					-- rank of this node (superkingdom, kingdom, ...) 
    #  	embl code				-- locus-name prefix; not unique
    #  	division id				-- see division.dmp file
    #  	inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
    # 	genetic code id				-- see gencode.dmp file
    # 	inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
    # 	mitochondrial genetic code id		-- see gencode.dmp file
    # 	inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
    # 	GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
    # 	hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
    # 	comments				-- free-text comments and citations
    #

    def __init__( self ):
        pass

    def __init__( self, nodes_dmp_file, tax_ids_to_names = None, accessions = None):
        tmp_nodes = {}
        #Read all fields from nodes.dmp file
        #tmp_nodes related tax_id to each clade
        with open( nodes_dmp_file ) as inpf:
            for line in (l.strip().split('\t') for l in inpf):
                ( tax_id, parent_tax_id, rank, embl_code, division_id, inherited_div_flag,
                genetic_code_id, inherited_GC_flag, mitochondrial_genetic_code, inherited_MGC_flag,
                GenBank_hidden_flag, hidden_subtree_root_flag, comments ) = line[::2]

                name = (tax_ids_to_names[int(tax_id)] if tax_ids_to_names else None)

                clade = BClade( clades = [], name = name )
                clade.parent_tax_id = int(parent_tax_id)
                clade.rank = re.sub(r'\W+', '', rank).strip("_")
                clade.tax_id = int(tax_id)
                
                clade.accession = accessions[clade.tax_id] if clade.tax_id in accessions else []
                if clade.tax_id in accessions:
                    clade.sequence_data = True
                    clade.status = clade.accession['status']
            
                tmp_nodes[clade.tax_id] = clade 
                
                # can add any other info in node.dmp

        self.tree = BTree()
        for node in tmp_nodes.values():
            # node = parent is the trick from NCBI to identify the root
            if node.tax_id == node.parent_tax_id:
                self.tree.root = node
                continue
            parent = tmp_nodes[node.parent_tax_id]
            parent.clades.append( node )


    def add_internal_accessions( self, clade = None ):
        if not clade:
            clade = self.tree.root

        clade.all_accessions = [] + ([clade.accession] if clade.accession else [])

        for child in clade.clades:
            clade.all_accessions += self.add_internal_accessions( child )
        return clade.all_accessions
       

    def remove_subtrees_without_accessions( self, clade = None ):
        if not clade:
            clade = self.tree.root
        clade.clades = [c for c in clade.clades if len(c.all_accessions)]
        for c in clade.clades:
            self.remove_subtrees_without_accessions( c )
    
    def remove_plasmids( self, clade = None ):
        if not clade:
            clade = self.tree.root
        clade.clades = [c for c in clade.clades if 'plasmid' not in c.name]
        for c in clade.clades:
            self.remove_plasmids( c )
   
    def print_tree( self, out_file_name, reduced = False ):

        tree = self.reduced_tree if reduced else self.tree 

        to_print = tree.find_clades({"sequence_data": True})
      
        ranks2code = {'superkingdom':'k','phylum':'p','class':'c','order':'o','family':'f','genus':'g','species':'s','taxon':'t'}

        def red_rank( rank ):
            if reduced and rank in ranks2code:
                return ranks2code[rank]
            return rank

        with open(out_file_name,"w") as outf:
            for t in to_print:
                tax = "|".join([red_rank(p.rank)+'__'+p.name for p in tree.get_path( t )])
                outf.write("\t".join( [ t.name,
                                        t.accession['status'],
                                        ",".join(t.accession['gen_seqs']),
                                        str(t.tax_id),
                                        #t.accession['code'],
                                        #",".join(t.accession['accession']),
                                        #str(t.accession['len']),
                                        tax
                                        ]    )+"\n")
        
    
    def get_tree_with_reduced_taxonomy( self, superkingdom = "Bacteria" ):
        reduced_tax_levels = ['superkingdom','phylum','class','order','family','genus','species']
        self.reduced_tree = copy.deepcopy(self.tree)
       
        def remove_noranks( clade ):

            run = True

            while run:
                run = False
                new_clades = []
                for c in clade.clades:
                    #if len(c.clades) and c.rank not in reduced_tax_levels:
                    if not hasattr(c,"sequence_data") and c.rank not in reduced_tax_levels:
                        run = True
                        new_clades += c.clades
                    else:
                        new_clades.append(c)
                    if hasattr(c,"sequence_data") and c.rank not in reduced_tax_levels:
                        c.rank = "norank"
                clade.clades = new_clades
            for c in clade.clades:
                if len(c.clades):
                    remove_noranks( c )

        def add_taxa( clade ):
            if clade.rank == "norank" and hasattr(clade,"sequence_data"):
                clade.rank = "taxon"

            if not len(clade.clades) and clade.accession:
                if clade.rank == 'species':
                    newclade = copy.deepcopy( clade )
                    clade.accession = [] 
                    clade.sequence_data = False
                    newclade.rank = "taxon"
                    clade.clades = [newclade]

            for c in clade.clades:
                add_taxa( c )


        def add_internal_missing_levels( clade, lev = 1 ):
            if clade.rank == "taxon":
                return
            cur_lev = reduced_tax_levels.index(clade.rank) if lev > 0  else 0

            jumps, to_add = [], ""
            for i,c in enumerate(clade.clades):
                if c.rank == 'taxon':
                    continue
                next_lev = reduced_tax_levels.index(c.rank)
                if next_lev == cur_lev + 1 or c.rank == 'superkingdom':
                    add_internal_missing_levels( c )
                    continue

                for i,l in enumerate(reduced_tax_levels[:-1]):
                    if clade.rank == l and c.rank != reduced_tax_levels[i+1]:
                        jumps.append( c )
                        to_add =  reduced_tax_levels[i+1]
            if jumps:
                #print len(jumps), len(clade.clades)
                children_ok = [c for c in clade.clades if c not in jumps]
                newclade = copy.deepcopy( clade )
                newclade.clades = jumps
                clade.clades = [newclade]+children_ok
                newclade.rank = to_add 
                newclade.name = clade.name if "_noname" in clade.name else clade.name+"_noname"
                add_internal_missing_levels( newclade )

        def reduce_double_taxa( clade ):
            if clade.rank == 'species' and len(clade.clades):
                torem = []
                for c in clade.clades:
                    if c.rank == 'taxon':
                        if len(c.clades):
                            clade.clades += c.clades
                            torem.append( c )
                clade.clades = [c for c in clade.clades if c not in torem]
                return
            for c in clade.clades:
                reduce_double_taxa( c )

        remove_noranks( self.reduced_tree.root )
        add_taxa( self.reduced_tree.root )
        add_internal_missing_levels( self.reduced_tree.root, lev = -1 )
        reduce_double_taxa( self.reduced_tree.root )

    #def save( self, out_file_name ):
    #    self.tree = self.tree.as_phyloxml()
    #    Phylo.write( self.tree, out_file_name, "phyloxml")

class Names:
    #
    # Format of names.dmp from RefSeq documentation
    #
    # ---------
    #
    # Taxonomy names file has these fields:
    #
    #   tax_id                  -- the id of node associated with this name
    #   name_txt                -- name itself
    #   unique name             -- the unique variant of this name if name not unique
    #   name class              -- (synonym, common name, ...)
    #

    def __init__( self, names_dmp_file ):
        #Read from file names.dmp, get information in every field
        self.tax_ids_to_names = {}
        with open( names_dmp_file ) as inpf:
            for line in (l.strip().split('\t') for l in inpf):
                tax_id, name_txt, unique, name_class = line[::2]
                
                # extracting scientific names only (at least for now) which are unique!
                if name_class == "scientific name":
                    name = re.sub(r'\W+', '_', name_txt).strip("_")
                    self.tax_ids_to_names[ int(tax_id) ] = name

            
    def get_tax_ids_to_names( self ):
        return self.tax_ids_to_names


class Accessions:
    #
    # Format of RefSeq-releaseXX.catalog
    #
    # Content: Tab-delimited listing of all accessions included in the current 
    # RefSeq release.
    # 
    # Columns:
    #  1. taxonomy ID
    #  2. species name
    #  3. accession.version
    #  4. gi
    #  5. refseq release directory accession is included in
    #       complete + other directories
    #       '|' delimited
    #  6. refseq status
    #       na - not available; status codes are not applied to most genomic records
    #       INFERRED
    #       PREDICTED
    #       PROVISIONAL
    #       VALIDATED
    #       REVIEWED
    #       MODEL
    #       UNKNOWN - status code not provided; however usually is provided for 
    #                 this type of record
    #  7. length  
    #
    # CODES FOR ACCESSIONS:
    # see: http://www.ncbi.nlm.nih.gov/projects/RefSeq/key.html
    #
    # EXTRACTING FOR NOW: NC and NZ

    def __init__( self, prokaryotes, eukaryotes, viruses, scaffs_in_complete, names_dmp_file ):
        self.accessions = {}

        # Create dictionary with scaffold as key and list of contig accessions as value.
        # Simply reads scaffs.txt.
        scafs = dict([(l[0][3:7],l[1].strip().split(',')) for l in (
                        line.split('\t') for line in open(scaffs_in_complete))])
        
        #Populate acc_ok with only chromosome DNA accessions (exclude plastid,plasmid,mitochondria)
        acc_ok = set() 
        with open( names_dmp_file ) as inpf:
            for line in (l.strip().split('\t') for l in inpf):
                code = line[2].split("_")[0]
                if 'plasmid' in line[4] or "mitochondrion" in line[4] or 'plastid' in line[4]:
                    continue
                if code == "NC":
                    acc_ok.add( line[2] )
                if code == "NZ":
                    acc_ok.add( line[2].split("_")[1][:4] ) 
        
        ncbi_files = [prokaryotes,eukaryotes]
       	#The following are the indexes for prokaryotes and eukaryotes NCBI file and what they pertain to depending on organism type.
	#We treat euks and proks separately from viruses because the fields in viruses do not coincide with the fields in euks and proks.
	#Organism/name		1
	#Chromosoms/refseq	8
	#WGS			12
	#Status			19
        for nf in ncbi_files:
	    print nf
            #Parse NCBI prokaryotes and eukaryotes file
            with open( nf ) as inpf:
                for line in (l.strip().split('\t') for l in inpf):
                    #Ignore line if comment, if status of genome = no data, or if both Chromosome/Refseq and WGS fields are empty
                    if line[0][0] == '#':
                        continue
                    if line[19] == "No data":
                        continue
                    if line[8] == '-' and line[12] == '-':
                        continue
                
                    #Get name, taxon id, and status of organism. Status is final 
                    #if WGS field is empty. If it contains four-letter code, 
                    #status is draft.
                    name = line[0]
                    taxid = int(line[1])
                    status = "final" if line[12] == '-' else "draft"
            
                    #If status is final, get NCBI accession number(s). Otherwise, get four letter
                    #code refering to refseq file name in which the draft genome is contained
                    if status == "final":
                        gen_seqs_tmp = line[8].split(",")
                        gen_seqs = list(set([gs for gs in gen_seqs_tmp if gs in acc_ok]))
                    else:
                        if line[12][:4] in scafs:
                            if line[12][:4] in acc_ok:
                                gen_seqs = scafs[line[12][:4]]
                            else:
                                gen_seqs = []
                        else:
                            gen_seqs_tmp = [l[:4] for l in line[12].split(",")]
                            gen_seqs = list(set([gs for gs in gen_seqs_tmp if gs in acc_ok]))


                    if not gen_seqs:
                        #print gen_seqs_tmp
                        continue

                    self.accessions[taxid] = { 'name' : name,
                                           'status' : status,
                                           'gen_seqs' : gen_seqs }
                    
        
       	

        """
        with open( names_dmp_file ) as inpf:
            for line in (l.strip().split('\t') for l in inpf):
                code = line[2].split("_")[0]
                if code == "NC":
                    if 'plasmid' in line[4] or "mitochondrion" in line[4] or 'plastid' in line[4]:
                        continue
                    acc = int(line[0])
                    if acc in self.accessions:
                        self.accessions[int(line[0])]['accession'].append( line[2] )
                        self.accessions[int(line[0])]['gi'].append( line[3] )
                        self.accessions[int(line[0])]['len'].append( line[6] )
                    else:
                        self.accessions[int(line[0])] = { 'name': line[1],
                                                          'accession': [line[2]],
                                                          'gi': [line[3]],
                                                          'refseq_dir': line[4],
                                                          'refseq_status': line[5],
                                                          'len': [line[6]],
                                                          'code': code,
                                                          'status': "final"}
                if code == "NZ":
                    if 'plasmid' in line[4] or "mitochondrion" in line[4]:
                        continue
                    if '00000000' not in line[2]:
                        continue
                    self.accessions[int(line[0])] = { 'name': line[1],
                                                      'accession': [line[2]],
                                                      'gi': [line[3]],
                                                      'refseq_dir': line[4],
                                                      'refseq_status': line[5],
                                                      'len': line[6],
                                                      'code': code,
                                                      'status': "draft"}
    
        """
    def get_accessions( self ):
        return self.accessions
    

if __name__ == '__main__':
    par = read_params(sys.argv)

    accessions = Accessions( par['prokaryotes'], par['eukaryotes'], par['viruses'], par['scaffs_in_complete'], par['catalog'] )
    names = Names( par['names'] )
    tax_tree = Nodes( par['nodes'], names.get_tax_ids_to_names(), accessions.get_accessions() ) 
    tax_tree.add_internal_accessions()
    tax_tree.remove_plasmids()
    tax_tree.remove_subtrees_without_accessions()
    tax_tree.get_tree_with_reduced_taxonomy()

    if par['output']:
        tax_tree.print_tree( par['output'] )
    if par['output_red']:
        tax_tree.print_tree( par['output_red'], reduced = True)



