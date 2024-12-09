#!/bin/bash

# Input fasta and gff3 were retrieved from NCBI

# Convert the gff3 genes to bed
awk 'BEGIN {OFS="\t"} $3 == "gene" { 
    split($9, attr, ";");
    for (i in attr) {
        if (attr[i] ~ /^Name=/) {
            name = substr(attr[i], 6);
        }
    }
    print $1, $4-1, $5, name, ".", $7;
}' CP001095.1_genes.gff3 > CP001095.1_genes.bed

# Extract the gene sequences from the fasta file
bedtools getfasta -fi CP001095.1_genome.fasta -bed CP001095.1_genes.bed -fo CP001095.1_gene_sequences.fasta -nameOnly

