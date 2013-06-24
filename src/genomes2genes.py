import sys
from Bio import SeqIO
from Bio.Seq import Seq


input_handle  = open("completeNZ_AGIG.genomic.gbff", "r")
fna = SeqIO.index("completeNZ_AGIG.1.genomic.fna", "fasta")

genes, proteins = [],[] 

for seq_record in SeqIO.parse(input_handle, "genbank") :
    gi = seq_record.annotations['gi']
    n = 0
    for seq_feature in seq_record.features :
        if seq_feature.type == "gene" :
            id = "|".join(["gi",gi,"ref",seq_record.id,""])
            seq_record.ref = fna[id]
            gene = seq_feature.extract( fna[id]  ) 
            gene.id = seq_record.id + "_" + str(n)
            gene.description = ""
            genes.append( gene )
            n += 1
        elif seq_feature.type == "CDS":
            id = "|".join(["gi",gi,"ref",seq_feature.qualifiers['protein_id'][0],""])
            name = seq_feature.qualifiers['product'][0]
            name = name + (" " + seq_feature.qualifiers['locus_tag'][0] if name == "hypothetical protein" else "")
            name += " ["+seq_record.annotations['organism']+"]"
            record = SeqIO.SeqRecord(   Seq(seq_feature.qualifiers['translation'][0]),
                                        id,
                                        '',
                                        name)
            proteins.append( record )

SeqIO.write(genes, "my_example.fna", "fasta")
SeqIO.write(proteins, "my_example.faa", "fasta")
            
