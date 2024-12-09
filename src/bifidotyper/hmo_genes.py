import os
import glob
import subprocess
import typing
from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from .logger import logger

class HMOError(Exception):
    """Custom exception for HMO-related errors."""
    pass

class HMOUtils:
    def __init__(self,
                 args,
                 salmon_executable: str,
                 sample_name: str,
                 genes_fasta: str, 
                 hmo_annotations: str,
                 fastq_se: str = None,
                 fastq_pe1: str = None,
                 fastq_pe2: str = None,
                 output_dir: str = 'hmo_quantification',
                 threads: int = 1):
        """
        Initialize HMO utility with input file paths.

        Args:
            genome_fasta (str): Path to genome FASTA file
            gene_annotations_gff3 (str): Path to gene annotations GFF3 file
            hmo_annotations (str): Path to HMO annotations file
        """
        self.args = args
        self.salmon_executable = salmon_executable
        self.sample_name = sample_name
        self.genes_fasta = genes_fasta
        self.hmo_annotations = hmo_annotations
        self.fastq_se = fastq_se
        self.fastq_pe1 = fastq_pe1
        self.fastq_pe2 = fastq_pe2
        self.output_dir = output_dir
        self.threads = threads

        # ensure that either se or pe1/2 are provided, but not both
        if self.fastq_se and (self.fastq_pe1 or self.fastq_pe2):
            raise HMOError("Please provide either single-end or paired-end fastq files, not both")
        elif not self.fastq_se and not (self.fastq_pe1 and self.fastq_pe2):
            raise HMOError("Please provide either single-end or paired-end fastq files")

        
        # Ensure all the files exist
        valid_files = [self.genes_fasta, self.hmo_annotations]
        if self.fastq_se:
            valid_files.append(self.fastq_se)
        else:
            valid_files.extend([self.fastq_pe1, self.fastq_pe2])

        for file_path in valid_files:
            if not os.path.exists(file_path):
                raise HMOError(f"File '{file_path}' does not exist")
        
        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

        # Run Salmon
        self.run_salmon()

        # Process gene counts
        self.process_gene_counts()

    
    def _run_command(self, command: typing.List[str], logfile:str = None, stdout=None) -> subprocess.CompletedProcess:
        """
        Execute an HMO command with error handling.
        
        Args:
            command (List[str]): Command to execute
        
        Returns:
            subprocess.CompletedProcess: Completed process object
        
        Raises:
            HMOError: If the command fails
        """

        if logfile:
            logfile = open(logfile, 'w')
        else:
            logfile = subprocess.DEVNULL
        
        if stdout:
            stdout = open(stdout, 'w')
        else:
            stdout = logfile
        
        try:
            logger.info(f"Running command: {' '.join(command)}")
            subprocess.run(command, check=True, text=True, capture_output=False, shell=False, stdout=stdout, stderr=logfile)
        except subprocess.CalledProcessError as e:
            raise HMOError(f"Command '{' '.join(command)}' failed with error: {e.stderr}")
    
    def run_salmon(self):
        
        if not os.path.exists(os.path.join(self.output_dir, 'B_longum_salmon_index')):
            index_cmd = [self.salmon_executable, 'index', '-t', self.genes_fasta, '-i', os.path.join(self.output_dir, 'B_longum_salmon_index')]
            self._run_command(index_cmd, logfile=os.path.join(self.output_dir, 'salmon_index.log'))

        if self.fastq_se:
            quant_cmd = [self.salmon_executable, 'quant', '-i', os.path.join(self.output_dir, 'B_longum_salmon_index'), '-l', 'A', '-r', self.fastq_se, '-p', str(self.threads), '--validateMappings', '-o', os.path.join(self.output_dir, self.sample_name+'_salmon')]
        else:
            quant_cmd = [self.salmon_executable, 'quant', '-i', os.path.join(self.output_dir, 'B_longum_salmon_index'), '-l', 'A', '-1', self.fastq_pe1, '-2', self.fastq_pe2, '-p', str(self.threads), '--validateMappings', '-o', os.path.join(self.output_dir, self.sample_name+'_salmon')]
        self._run_command(quant_cmd, logfile=os.path.join(self.output_dir, self.sample_name+'_salmon.log'))


    def process_gene_counts(self):
        
        salmon_counts = pd.read_csv(os.path.join(self.output_dir, self.sample_name+'_salmon','quant.sf'), sep='\t')

        hmo = pd.read_csv(self.hmo_annotations,sep=';')[['Blon','Cluster']]
        hmo = hmo.rename(columns={'Blon':'Name'})
        # Some cells have multiple "Blon" ids, we need to split them into separate rows
        hmo = hmo.assign(Name=hmo.Name.str.split(' ')).explode('Name')
        # Now remove any row that doesn't match the format "Blon_XXXX"
        hmo = hmo[hmo['Name'].str.match(r'Blon_\d+')]

        salmon_counts = pd.merge(salmon_counts,hmo,left_on='Name',right_on='Name',how='left')
        salmon_counts.dropna(inplace=True)
        salmon_counts['Present'] = salmon_counts['NumReads'] > 0
        salmon_counts.to_csv(os.path.join(self.output_dir,self.sample_name+'.salmon_counts_annotated.tsv'), index=False, sep='\t')
        logger.info('Saved annotated salmon counts to {}'.format(os.path.join(self.output_dir,self.sample_name+'.salmon_counts_annotated.tsv')))

        clusters = salmon_counts.groupby('Cluster')['Present'].all().reset_index()
        clusters.to_csv(os.path.join(self.output_dir,self.sample_name+'.cluster_presence.tsv'), index=False, sep='\t')
        logger.info('Saved cluster presence table to {}'.format(os.path.join(self.output_dir,self.sample_name+'.cluster_presence.tsv')))


# Example usage
def main():
    pass

if __name__ == '__main__':
    main()