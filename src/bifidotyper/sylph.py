import os
import glob
import subprocess
import typing
from pathlib import Path
from .logger import logger
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

class SylphUtils:
    def __init__(self,args,sylph_executable):
        self.args = args
        self.sylph_executable = sylph_executable
        self.genome_sketch_dir = 'sylph_genome_sketches'
        self.fastq_sketch_dir = 'sylph_fastq_sketches'
        self.genome_query_dir = 'sylph_genome_queries'

        # Ensure directories exist
        for dir_path in [self.genome_sketch_dir, 
                         self.fastq_sketch_dir, 
                         self.genome_query_dir]:
            os.makedirs(dir_path, exist_ok=True)
    
    def _run_command(self, command: typing.List[str]) -> subprocess.CompletedProcess:
        try:
            if self.args.verbose:
                logger.info(f"Running command: {' '.join(command)}")
            result = subprocess.run(command, check=True, text=True, capture_output=True, shell=False)
            return result
        except subprocess.CalledProcessError as e:
            logger.error(f"Command '{' '.join(command)}' failed with error: {e.stderr}")
            raise Error(f"Command '{' '.join(command)}' failed with error: {e.stderr}")
    
    def sketch_reads(self,
                     fastq_se: list = None,
                     fastq_r1: list = None,
                     fastq_r2: list = None,
                     threads: int = 1):
        missing_files = []
        command = []

        def find_existing_sylsp(fastq_files):
            """Find existing .sylsp files for the given FASTQ files."""
            existing_files = []
            for fastq in fastq_files:
                base_name = os.path.basename(fastq).replace('.fastq.gz', '').replace('.fastq', '').replace('.fq.gz', '').replace('.fq', '')
                sylsp_pattern = f"{base_name}*.sylsp"  # Match any .sylsp file with the base name
                sylsp_files = glob.glob(os.path.join(self.fastq_sketch_dir, sylsp_pattern)) + glob.glob(sylsp_pattern)
                existing_files.extend(sylsp_files)
            return existing_files

        if fastq_se:
            existing_files = find_existing_sylsp(fastq_se)
            missing_files = [f for f in fastq_se if not any(os.path.basename(f).replace('.fastq.gz', '').replace('.fastq', '').replace('.fq.gz', '').replace('.fq', '') in ef for ef in existing_files)]
            if len(existing_files) > 0:
                logger.info(f"      Found {len(existing_files)} existing .sylsp files. Sketching the remaining {len(missing_files)} samples.")
            else:
                logger.info(f"      No existing .sylsp files found. Sketching all {len(missing_files)} samples.")
            if missing_files:
                command = [self.sylph_executable, 'sketch', *missing_files, '-t', str(threads)]
        elif fastq_r1 and fastq_r2:
            existing_files_r1 = find_existing_sylsp(fastq_r1)
            missing_files_r1 = [f for f in fastq_r1 if not any(os.path.basename(f) in ef for ef in existing_files_r1)]
            missing_files_r2 = [f.replace(self.args.r1_suffix,self.args.r2_suffix) for f in missing_files_r1]
            if len(existing_files_r1) > 0:
                logger.info(f"      Found {len(existing_files_r1)} existing .sylsp files. Sketching the remaining {len(missing_files_r1)} samples.")
            else:
                logger.info(f"      No existing .sylsp files found. Sketching all {len(missing_files_r1)} samples.")
            if missing_files_r1:
                command = [self.sylph_executable, 'sketch', '-1', *missing_files_r1, '-2', *missing_files_r2, '-t', str(threads)]
        else:
            raise Error("Either fastq_se or fastq_r1 and fastq_r2 must be provided")

        if command:
            logger.info(f"Missing .sylsp files: {missing_files}. Running Sylph sketch command.")
            self._run_command(command)
        else:
            logger.info("       All .sylsp files are present. Skipping Sylph sketch command.")

        # Move any .sylsp files in the base dir to the fastq_sketch_dir
        for sylsp in glob.glob('*.sylsp'):
            os.rename(sylsp, os.path.join(self.fastq_sketch_dir, sylsp))

        # Return paths to all existing .sylsp files
        return glob.glob(os.path.join(self.fastq_sketch_dir, '*.sylsp'))
    
    def query_genomes(self, 
                      sylsp_files: typing.List[str], 
                      syldb_file: str, 
                      output_name: str = 'genome_query.tsv') -> str:
        
        # Construct sylph query command
        command = [self.sylph_executable, 'query'] + sylsp_files + [syldb_file, '-o', output_name]
        self._run_command(command)

        # Move the output file into the genome query directory
        os.rename(output_name, os.path.join(self.genome_query_dir, output_name))

        return os.path.join(self.genome_query_dir, output_name)
    
    def profile_genomes(self, 
                        sylsp_files: typing.List[str], 
                        syldb_file: str, 
                        output_name: str = 'genome_profile.tsv') -> str:
        
        # Construct sylph profile command
        command = [self.sylph_executable, 'profile'] + sylsp_files + [syldb_file, '-o', output_name]
        self._run_command(command)

        # Move the output file into the genome query directory
        os.rename(output_name, os.path.join(self.genome_query_dir, output_name))
        
        return os.path.join(self.genome_query_dir, output_name)





# Example usage
def main():
    try:
        # Initialize Sylph utility
        sylph = SylphUtils()
        
        # Sketch genomes
        genome_db = sylph.sketch_genomes('/path/to/genomes/*.fna')
        
        # Sketch reads (paired-end example)
        read_sketches = sylph.sketch_reads('/path/to/reads/*_1.fastq.gz', is_paired_end=True)
        
        # Query genomes
        query_result = sylph.query_genomes(read_sketches, genome_db)
        
        # Profile genomes
        profile_result = sylph.profile_genomes(read_sketches, genome_db)
        
        logger.info(f"Query result: {query_result}")
        logger.info(f"Profile result: {profile_result}")
    
    except Error as e:
        logger.info(f"Sylph processing error: {e}")
        raise Error(f"Sylph processing error: {e}")

if __name__ == '__main__':
    main()

