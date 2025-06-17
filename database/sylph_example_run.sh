#!/bin/bash

# Tell the script where to find the genomes you want to run Sylph on
genomes_path=/pl/active/olm-data2/Projects/2024_Bifidobacteria_Database/bifido_genomes/Database_alpha/fna

# Tell the script where to find the fastq files you want to run Sylph on
fastq_path=/pl/active/olm-data2/Data/Public/15_day/PRJEB73511

# Make a directory to store the genome sketches (optional)
mkdir -p sylph_genome_sketches
cd sylph_genome_sketches


# Run Sylph Sketch on each set of genomes to prepare for rapid Sylph mapping
# https://github.com/bluenote-1577/sylph/wiki/sylph-cookbook#database-sketching-options
sylph sketch $genomes_path/*.fna -t 5 -o my_genomes
cd ..
# This will make a .syldb file

# Sketch (index) fastq files
# https://github.com/bluenote-1577/sylph/wiki/sylph-cookbook#read-sketching-options

# For paired-end reads, use the -1 and -2 options to specify the forward and reverse reads
# -t 5 indicates that 5 cpu threads will be used
sylph sketch -1 $fastq_path/*_1.fastq.gz -2 $fastq_path/*_2.fastq.gz -d fastq_sketches -t 5

# For single-end reads, optionally use -r or just input the read files directly
sylph sketch $fastq_path/*.fastq.gz -d fastq_sketches -t 5

cd ..
# This will make a .sylsp file for each fastq file


# Once you have a genome database (*.syldb) and a fastq sketch database (*.syldb), you can query the fastq sketches against the genome sketches

mkdir -p sylph_genome_queries
cd sylph_genome_queries

# For sylph query, just give it a set of *.sylsp files and a *.syldb file to query against
# Specify an output TSV file with -o
sylph query fastq_sketches/*.sylsp sylph_genome_sketches/my_genomes.syldb -o sylph_genome_queries/genome_query.tsv

# sylph profile behaves very similarly to query, but its output is further processed to be a limited set of genomes per sample, with proportions adding to 100%
sylph profile fastq_sketches/*.sylsp sylph_genome_sketches/my_genomes.syldb -o sylph_genome_queries/genome_profile.tsv

cd ..

