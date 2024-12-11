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
            logger.info(f"Running command: {' '.join(command)}")
            result = subprocess.run(command, check=True, text=True, capture_output=True, shell=False)
            return result
        except subprocess.CalledProcessError as e:
            raise SylphError(f"Command '{' '.join(command)}' failed with error: {e.stderr}")
    
    def sketch_genomes(self, 
                       genomes: list, 
                       output_name: str = 'my_genomes', 
                       threads: int = 1) -> str:
                
        # Construct sylph sketch command
        command = [self.sylph_executable, 'sketch',*genomes,'-t', str(threads), '-o', output_name]
        self._run_command(command)

        # Move all output files into the genome sketch directory
        syldb = os.path.join(self.genome_sketch_dir, f"{output_name}.syldb")
        os.rename(f"{output_name}.syldb", syldb)

        return syldb
    
    def sketch_reads(self,
                     fastq_se: list = None,
                     fastq_r1: list = None,
                     fastq_r2: list = None,
                     threads: int = 1):
        
        # Construct sylph sketch command
        if fastq_se:
            command = [self.sylph_executable, 'sketch', *fastq_se, '-t', str(threads)]
        elif fastq_r1 and fastq_r2:
            command = [self.sylph_executable,'sketch','-1',*fastq_r1,'-2',*fastq_r2,'-t',str(threads)]
        else:
            raise SylphError("Either fastq_se or fastq_r1 and fastq_r2 must be provided")
        
        self._run_command(command)
        
        # Move all output files into the fastq sketch directory
        for sylsp in glob.glob('*.sylsp'):
            os.rename(sylsp, os.path.join(self.fastq_sketch_dir, sylsp))
        
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
    
    except SylphError as e:
        logger.info(f"Sylph processing error: {e}")
        raise SylphError(f"Sylph processing error: {e}")

if __name__ == '__main__':
    main()

