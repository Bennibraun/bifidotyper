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

class SylphError(Exception):
    """Custom exception for Sylph-related errors."""
    pass

class SylphUtils:
    def __init__(self,
                 args,
                 genome_sketch_dir: str = 'sylph_genome_sketches', 
                 fastq_sketch_dir: str = 'sylph_fastq_sketches', 
                 genome_query_dir: str = 'sylph_genome_queries',
                 ):
        """
        Initialize Sylph utility with configurable directory paths.
        
        Args:
            genome_sketch_dir (str): Directory to store genome sketches
            fastq_sketch_dir (str): Directory to store fastq sketches
            genome_query_dir (str): Directory to store genome query results
        """
        self.args = args
        self.genome_sketch_dir = genome_sketch_dir
        self.fastq_sketch_dir = fastq_sketch_dir
        self.genome_query_dir = genome_query_dir
        
        # Ensure directories exist
        for dir_path in [self.genome_sketch_dir, 
                         self.fastq_sketch_dir, 
                         self.genome_query_dir]:
            os.makedirs(dir_path, exist_ok=True)
    
    def _run_command(self, command: typing.List[str]) -> subprocess.CompletedProcess:
        """
        Execute a Sylph command with error handling.
        
        Args:
            command (List[str]): Command to execute
        
        Returns:
            subprocess.CompletedProcess: Completed process object
        
        Raises:
            SylphError: If the command fails
        """
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
        """
        Sketch genome files using Sylph.
        
        Args:
            genomes_path (str): Path to genome files (e.g., '*.fna')
            output_name (str): Base name for output .syldb file
            threads (int): Number of CPU threads to use
        
        Returns:
            str: Path to the generated .syldb file
        """
                
        # Construct sylph sketch command
        command = ['sylph', 'sketch',*genomes,'-t', str(threads), '-o', output_name]
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
        """
        Sketch read files using Sylph.
        
        Args:
            fastq_path (str): Path to FASTQ files (e.g., '*.fastq.gz')
            is_paired_end (bool): Whether the reads are paired-end
            threads (int): Number of CPU threads to use
        
        Returns:
            List[str]: Paths to the generated .sylsp files
        """
        
        # Construct sylph sketch command
        if fastq_se:
            command = ['sylph', 'sketch', *fastq_se, '-t', str(threads)]
        elif fastq_r1 and fastq_r2:
            command = ['sylph','sketch','-1',*fastq_r1,'-2',*fastq_r2,'-t',str(threads)]
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
        """
        Query genomes using Sylph.
        
        Args:
            sylsp_files (List[str]): List of .sylsp files to query
            syldb_file (str): Path to the .syldb file to query against
            output_name (str): Name of the output TSV file
        
        Returns:
            str: Path to the generated query TSV file
        """
        
        # Construct sylph query command
        command = ['sylph', 'query'] + sylsp_files + [syldb_file, '-o', output_name]
        self._run_command(command)

        # Move the output file into the genome query directory
        os.rename(output_name, os.path.join(self.genome_query_dir, output_name))

        return os.path.join(self.genome_query_dir, output_name)
    
    def profile_genomes(self, 
                        sylsp_files: typing.List[str], 
                        syldb_file: str, 
                        output_name: str = 'genome_profile.tsv') -> str:
        """
        Profile genomes using Sylph.
        
        Args:
            sylsp_files (List[str]): List of .sylsp files to profile
            syldb_file (str): Path to the .syldb file to profile against
            output_name (str): Name of the output TSV file
        
        Returns:
            str: Path to the generated profile TSV file
        """
        
        # Construct sylph profile command
        command = ['sylph', 'profile'] + sylsp_files + [syldb_file, '-o', output_name]
        self._run_command(command)

        # Move the output file into the genome query directory
        os.rename(output_name, os.path.join(self.genome_query_dir, output_name))
        
        return os.path.join(self.genome_query_dir, output_name)

    def plot_sylph_results(self, profile_tsv: str, query_tsv: str):
        """
        Plot Sylph query and profile results using matplotlib.
        
        Args:
            profile_tsv (str): Path to the profile TSV file
            query_tsv (str): Path to the query TSV file
        """
        
        # Read TSV file into a pandas DataFrame
        pdf = pd.read_csv(profile_tsv,sep='\t')
        qdf = pd.read_csv(query_tsv,sep='\t')

        # Make a strain label that cuts off the contig name at 30 chars
        qdf['Strain'] = qdf['Contig_name'].apply(lambda x: x[:30])
        pdf['Strain'] = pdf['Contig_name'].apply(lambda x: x[:30])

        qdf['Sample'] = qdf['Sample_file'].apply(lambda x: os.path.basename(x).replace('.fastq.gz','').replace(self.args['r1_suffix'],'').replace(self.args['r2_suffix'],''))
        pdf['Sample'] = pdf['Sample_file'].apply(lambda x: os.path.basename(x).replace('.fastq.gz','').replace(self.args['r1_suffix'],'').replace(self.args['r2_suffix'],''))

        # Set up a color scheme to uniquely identify each strain
        strains = pdf['Strain'].unique()




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