import pandas as pd
import os
import glob

# This script takes as arguments a path to a set of .fna files and a path to a CSV file that may or may not exist
# And updates the CSV file with the IDs, names, sizes, and colors for each genome that's not already present

# parse CLI arguments in order
# first one is path, second is csv
import argparse
parser = argparse.ArgumentParser(description='Update bifidotyper genomes CSV file with genome information.')
parser.add_argument('genomes_path', type=str, help='Path to the directory containing .fna files.')
parser.add_argument('csv_path', type=str, help='Path to the CSV file to update.')
args = parser.parse_args()

if os.path.exists(csv_path):
	df = pd.read_csv(csv_path)
else:
	df = pd.DataFrame(columns=['Genome_file','Label','Genome_size','Color'])

# Get all .fna files in the specified directory
genome_files = glob.glob(os.path.join(args.genomes_path, '*.fna'))

if not genome_files:
	print(f"No .fna files found in the directory: {args.genomes_path}")

# Iterate through each genome file and extract information
for genome_file in genome_files:
	# Check if it's in the DataFrame already
	name = os.path.basename(genome_file)
	if name in df['Genome_file'].values:
		print(f"Genome {name} already exists in the CSV. Skipping.")
		continue
	
	# get label and size from fna
	with open(genome_file, 'r') as f:
		first_line = f.readline().strip()
		# First line is ">label"
		label = first_line[1:]
		# Rest is the sequence
		# Go through and count the length of any line that doesn't start with ">"
		sequence_length = 0
		for line in f:
			if not line.startswith('>'):
				sequence_length += len(line.strip())
	
	# Pick a random color
	import random
	color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
	# Append the new genome information to the DataFrame
	new_row = {
		'Genome_file': name,
		'Label': label,
		'Genome_size': sequence_length,
		'Color': color
	}
	df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)

# Save the updated DataFrame back to the CSV file
df.to_csv(args.csv_path, index=False)
