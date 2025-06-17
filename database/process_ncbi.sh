#!/bin/bash

# We already selected only the bifidobacteriaceae genus from NCBI...

ncbi_path=~/2024_Bifidobacteria_Database/bifido_genomes/NCBI_Genomes

cd $ncbi_path

# Run dRep to dereplicate the genomes
dRep dereplicate drep_genomes -g $ncbi_path/ncbi_dataset/data/*/*.fna

# Create a scaffold-to-bin table so that sequences can be mapped back to genomes after concatenation
python parse_stb.py --reverse -f $ncbi_path/drep_genomes/dereplicated_genomes/*.fna -o scaffold_to_bin.tsv



# Sylph pre-mapping of genomes
# sylph sketch /home/bebr1814/olm-data2/Data/Public/15_day/PRJEB73511/ERR13710794_1.fastq.gz $ncbi_path/drep_genomes/dereplicated_genomes/*.fna -o slyph_db
sylph sketch -1 /home/bebr1814/olm-data2/Data/Public/15_day/PRJEB73511/*_1.fastq.gz -2 /home/bebr1814/olm-data2/Data/Public/15_day/PRJEB73511/*_2.fastq.gz -d PRJEB73511_sketches -t 5

sylph sketch $ncbi_path/drep_genomes/dereplicated_genomes/*.fna -t 5 -o drep_genomes

# Get sequence similarities
sylph query PRJEB73511_sketches/*.sylsp drep_genomes.syldb -o PRJEB73511_query.tsv

# Get abundances of genomes in the sample
sylph profile PRJEB73511_sketches/*.sylsp drep_genomes.syldb -o PRJEB73511_profile.tsv

