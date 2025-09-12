#!/bin/bash


# This script downloads data from GTDB and NCBI
# and runs dRep to establish a set of representative genomes
# For the genus Bifidobacteriaceae

# To run this script, you will need the following tools installed:
# - GTDB-Tk (https://github.com/Ecogenomics/GTDBTk)
# - NCBI Datasets CLI (https://github.com/ncbi/datasets)
# - dRep (https://github.com/MrOlm/drep)
# - sylph (https://github.com/bluenote-1577/sylph)
# You can install them via conda:
# 	conda install -c bioconda gtdbtk ncbi-datasets-cli drep sylph

# If running as sbatch, update the name of your conda env below
# eval "$(conda shell.bash hook)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bif_test

# Usage example:
# sbatch -N 1 -c 1 -p short -t 0-11:00 --mem=10G src/util/recreate_sylph_genome_db.sh /Users/bebr1814/scratch/bifidotyper/testing/genome_db /Users/bebr1814/scratch/bifidotyper/testing/genome_db/genomes.csv

# You can specify a data directory for the downloads
# as the first argument or just use the current dir
if [ -n "$1" ]; then
	data_dir="$1"
else
	data_dir="."
fi

# You can also specify the path to an existing genomes.csv, which comes with Bifidotyper and adds some convenient annotation to each genome
# Or you can just generate a new one by leaving this blank
if [ -n "$2" ]; then
	genomes_csv="$2"
else
	genomes_csv="genomes.csv"
fi

ncbi_path="$data_dir/NCBI_Genomes"
gtdb_path="$data_dir/GTDB_Genomes"
genomes_path="$data_dir/All_Genomes"

# Create directories if they don't exist
mkdir -p "$ncbi_path" "$gtdb_path" "$genomes_path"


# Download NCBI with the Datasets CLI
datasets download genome taxon bifidobacteriaceae --assembly-level complete --dehydrated
mv ncbi_dataset.zip "$ncbi_path/"
unzip "$ncbi_path/ncbi_dataset.zip" -d "$ncbi_path/"
datasets rehydrate --directory "$ncbi_path"

# Download GTDB with https://github.com/Ecogenomics/GTDBTk
export GTDBTK_DATA_PATH="$gtdb_path/data"
# echo "Downloading GTDB data..."
# wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz -O "$gtdb_path/gtdbtk_data.tar.gz" -q
# mkdir -p "$gtdb_path/data"
# tar xvzf "$gtdb_path/gtdbtk_data.tar.gz" -C "$gtdb_path/data" --strip 1
# # Remove the tar file after extraction
# rm "$gtdb_path/gtdbtk_data.tar.gz"

# Grab Bifidobacteriaceae genus from GTDB
grep "Bifido" "$gtdb_path/data/taxonomy/gtdb_taxonomy.tsv" > "$gtdb_path/Bifidobacteriaceae_taxonomy.tsv"
# Grab the FNA files that match each row in the taxonomy file
# file looks like GCA_000008085.1_genomic.fna.gz
# item looks like GB_GCA_000008085.1
mkdir -p "$gtdb_path/genomes"
for i in $(cat $gtdb_path/Bifidobacteriaceae_taxonomy.tsv | cut -f1); do
	match_i=$(echo $i | cut -d'_' -f2-)
	# echo $match_i
	find $gtdb_path/data/skani/database/ -name "*$match_i*" -exec cp {} $gtdb_path/genomes/ \;
	# find $gtdb_path/data/release*/skani/database/ -name "*$match_i*" -exec echo {} \;
done
gunzip $gtdb_path/genomes/*.gz


# Copy all genomes to a shared directory
mkdir -p $genomes_path/fna
cp $ncbi_path/ncbi_dataset/data/*/*.fna $genomes_path/fna/
cp $gtdb_path/genomes/*.fna $genomes_path/fna/

# Run dRep on all of these genomes to establish a consensus set
dRep dereplicate $genomes_path/drep_genomes -g $genomes_path/fna/*.fna --S_ani 0.95

# Finally, run Sylph Sketch on the dereplicated genome set
sylph sketch $genomes_path/drep_genomes/dereplicated_genomes/*.fna -t 5 -o bifidobacteria_sketches

# Bifidotyper needs a genomes.csv to go along with the sketches. We can update the existing one with this script
# This script will create a new genomes.csv file in the current directory or modify one that's already there
# So if you have one, make sure to provide it as $2
python update_bifidotyper_genomes_csv.py $genomes_path/drep_genomes/dereplicated_genomes/ $genomes_csv

echo "Done! Now you can drop the new bifidobacteria_sketches.syldb and genomes.csv into Bifidotyper under data/reference."
echo "

