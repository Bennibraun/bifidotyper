#!/bin/bash


ncbi_path=~/2024_Bifidobacteria_Database/bifido_genomes/NCBI_Genomes

gtdb_path=~/2024_Bifidobacteria_Database/bifido_genomes/GTDB_Genomes

genomes_path=~/2024_Bifidobacteria_Database/bifido_genomes/All_Genomes


# # Download GTDB
# export GTDBTK_DATA_PATH="$gtdb_path/data"
# nohup download-db.sh &

# # Grab Bifidobacteriaceae genus
# cd $gtdb_path
# grep "Bifido" $gtdb_path/data/release220/taxonomy/gtdb_taxonomy.tsv > $gtdb_path/Bifidobacteriaceae_taxonomy.tsv

# # Grab the FNA files that match each row in the taxonomy file
# # file looks like GCA_000008085.1_genomic.fna.gz
# # item looks like GB_GCA_000008085.1
# for i in $(cat $gtdb_path/Bifidobacteriaceae_taxonomy.tsv | cut -f1); do
# 	match_i=$(echo $i | cut -d'_' -f2-)
# 	echo $match_i
# 	find $gtdb_path/data/release*/skani/database/ -name "*$match_i*" -exec cp {} genomes/ \;
# 	# find $gtdb_path/data/release*/skani/database/ -name "*$match_i*" -exec echo {} \;
# done
# gunzip $gtdb_path/genomes/*.gz


# # Download NCBI
# # Requires the ncbi-datasets-cli package to be installed:
# # conda install -c conda-forge ncbi-datasets-cli
# datasets download genome taxon bifidobacteriaceae --assembly-level complete --dehydrated
# unzip ncbi_dataset.zip
# datasets rehydrate --directory .



# Run dRep
cd $genomes_path

# Copy all genomes to a shared directory
cp $ncbi_path/ncbi_dataset/data/*/*.fna $genomes_path/fna/
cp $gtdb_path/genomes/*.fna $genomes_path/fna/


# Run dRep on all of these genomes to establish a consensus set

# For experimental purposes, we want to try different thresholds for secondary ANI (95%, 98%)
dRep dereplicate drep_genomes_95 -g $genomes_path/fna/*.fna --S_ani 0.95
dRep dereplicate drep_genomes_98 -g $genomes_path/fna/*.fna --S_ani 0.98

# Run Sylph Sketch on each set of genomes to prepare for rapid Sylph mapping
mkdir -p sylph_genome_sketches
cd sylph_genome_sketches
sylph sketch $genomes_path/fna/*.fna -t 5 -o raw_genomes
sylph sketch ../drep_genomes_95/dereplicated_genomes/*.fna -t 5 -o drep_genomes_95
sylph sketch ../drep_genomes_98/dereplicated_genomes/*.fna -t 5 -o drep_genomes_98
cd ..

# Pre-map fastq files as needed
mkdir -p sylph_fastq_sketches
cd sylph_fastq_sketches
sylph sketch -1 /home/bebr1814/olm-data2/Data/Public/15_day/PRJEB73511/*_1.fastq.gz -2 /home/bebr1814/olm-data2/Data/Public/15_day/PRJEB73511/*_2.fastq.gz -d PRJEB73511_sketches -t 5
cd ..

# Query and Profile the genomes with the fastq sketches
mkdir -p sylph_genome_queries
cd sylph_genome_queries

sylph query sylph_fastq_sketches/PRJEB73511_sketches/*.sylsp sylph_genome_sketches/raw_genomes.syldb -o sylph_genome_queries/PRJEB73511_raw_genomes_query.tsv
sylph profile sylph_fastq_sketches/PRJEB73511_sketches/*.sylsp sylph_genome_sketches/raw_genomes.syldb -o sylph_genome_queries/PRJEB73511_raw_genomes_profile.tsv

sylph query sylph_fastq_sketches/PRJEB73511_sketches/*.sylsp sylph_genome_sketches/drep_genomes_95.syldb -o sylph_genome_queries/PRJEB73511_drep_genomes_95_query.tsv
sylph profile sylph_fastq_sketches/PRJEB73511_sketches/*.sylsp sylph_genome_sketches/drep_genomes_95.syldb -o sylph_genome_queries/PRJEB73511_drep_genomes_95_profile.tsv

sylph query sylph_fastq_sketches/PRJEB73511_sketches/*.sylsp sylph_genome_sketches/drep_genomes_98.syldb -o sylph_genome_queries/PRJEB73511_drep_genomes_98_query.tsv
sylph profile sylph_fastq_sketches/PRJEB73511_sketches/*.sylsp sylph_genome_sketches/drep_genomes_98.syldb -o sylph_genome_queries/PRJEB73511_drep_genomes_98_profile.tsv

cd ..

