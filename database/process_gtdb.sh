#!/bin/bash

gtdb_path=$1
cd $gtdb_path

wget https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz > $gtdb_path/gtdb_genomes_reps.tar.gz

tar -xvf gtdb_genomes_reps.tar.gz

grep "Bifido" gtdb_genomes_reps/data/release*/taxonomy/gtdb_taxonomy.tsv > Bifidobacteriaceae_taxonomy.tsv

# Grab the FNA files that match each row in the taxonomy file
# file looks like GCA_000008085.1_genomic.fna.gz
# item looks like GB_GCA_000008085.1
for i in $(cat Bifidobacteriaceae_taxonomy.tsv | cut -f1); do
	match_i=$(echo $i | cut -d'_' -f2-)
	echo $match_i
	find $gtdb_path/data/release*/skani/database/ -name "*$match_i*" -exec cp {} genomes/ \;
	# find $gtdb_path/data/release*/skani/database/ -name "*$match_i*" -exec echo {} \;
done

gunzip genomes/*.gz

