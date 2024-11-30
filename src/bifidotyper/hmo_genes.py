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
                 genome_fasta: str, 
                 gene_annotations_gff3: str, 
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
        self.genome_fasta = genome_fasta
        self.gene_annotations_gff3 = gene_annotations_gff3
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
        valid_files = [self.genome_fasta, self.gene_annotations_gff3, self.hmo_annotations]
        if self.fastq_se:
            valid_files.append(self.fastq_se)
        else:
            valid_files.extend([self.fastq_pe1, self.fastq_pe2])

        for file_path in valid_files:                            
            if not os.path.exists(file_path):
                raise HMOError(f"File '{file_path}' does not exist")
        
        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

        # Run bowtie2
        sam_file = self.run_bowtie2()

        # Run samtools postprocessing
        bam_file = self.run_samtools_postprocessing(sam_file)

        # Run gene quantification
        gene_counts = self.run_gene_quantification(bam_file)

        # Process gene counts
        self.process_gene_counts(gene_counts)

    
    def _run_command(self, command: typing.List[str], logfile:str = None) -> subprocess.CompletedProcess:
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
        
        try:
            logger.info(f"Running command: {' '.join(command)}")
            subprocess.run(command, check=True, text=True, capture_output=False, shell=False, stdout=logfile, stderr=logfile)
        except subprocess.CalledProcessError as e:
            raise HMOError(f"Command '{' '.join(command)}' failed with error: {e.stderr}")
    
    def run_bowtie2(self) -> str:
        """
        Run Bowtie2 to align reads to the genome.

        Returns:
            str: Path to the generated SAM file

        Raises:
            HMOError: If the command fails
        """

        genome_index = os.path.join(self.output_dir, 'genome_index')

        # Construct bowtie2-build command
        command = ['bowtie2-build', self.genome_fasta, genome_index]
        self._run_command(command, logfile=os.path.join(self.output_dir, 'bowtie2_build.log'))
        
        # Construct bowtie2 command
        if self.fastq_se:
            sam_file = os.path.join(self.output_dir, os.path.basename(self.fastq_se)+'.B.longum.sam')
            command = ['bowtie2','--no-unal', '-x', genome_index, '-U', self.fastq_se, '-S', sam_file, '-p', str(self.threads)]
        else:
            sam_file = os.path.join(self.output_dir, os.path.basename(self.fastq_pe1)+'.B.longum.sam')
            command = ['bowtie2','--no-unal', '-x', genome_index, '-1', self.fastq_pe1, '-2', self.fastq_pe2, '-S', sam_file, '-p', str(self.threads)]
        
        self._run_command(command, logfile=os.path.join(self.output_dir, 'bowtie2.log'))

        return sam_file

    def run_samtools_postprocessing(self, sam_file: str) -> str:

        # Construct samtools view command
        bam_file = sam_file.replace('.sam', '.bam')
        command = ['samtools', 'view', '-bS', sam_file, '-o', bam_file, '-@', str(self.threads)]
        self._run_command(command)

        os.remove(sam_file)

        # Construct samtools sort command
        sorted_bam_file = bam_file.replace('.bam', '.sorted.bam')
        command = ['samtools', 'sort', bam_file, '-o', sorted_bam_file, '-@', str(self.threads)]
        self._run_command(command)

        os.remove(bam_file)

        # Construct samtools index command
        command = ['samtools', 'index', sorted_bam_file, '-@', str(self.threads)]
        self._run_command(command)

        return sorted_bam_file

    def run_gene_quantification(self, bam_file: str) -> str:

        # Construct htseq-count command
        gene_counts_file = os.path.join(self.output_dir, 'gene_counts.txt')
        gene_counts_log = os.path.join(self.output_dir, 'htseq_count.log')
        command = ['htseq-count', '-f', 'bam', '-r', 'pos', '-s', 'no', '-t', 'CDS', '-i', 'ID', bam_file, self.gene_annotations_gff3]
        self._run_command(command,logfile=gene_counts_log)

        os.remove(bam_file)
        os.remove(bam_file+'.bai')

        return gene_counts_file

    def process_gene_counts(self, gene_counts_file: str):
        
        df = pd.read_csv(gene_counts_file, header=None, names=['gene','count'], sep='\t')
        # drop last 5 rows which are metadata
        df = df[:-5]

        df['Present'] = np.where(df['count'] > 0, True,False)
        df['gene_id'] = df['gene'].str.replace('cds-', '')

        HMOs = pd.read_csv(self.hmo_annotations,sep=';')
        HMOs['gene_id'] = HMOs['HMOgenes'].str.split('_cds_').str[1].str.split('_').str[0]

        df = df.merge(HMOs, on='gene_id', how='left').dropna()

        clusters = df.groupby("Cluster")["Present"].all().reset_index()
        clusters.to_csv(os.path.join(self.output_dir,'cluster_presence.csv'), index=False)

        logger.info('Saved cluster presence table to {}'.format(os.path.join(self.output_dir,'cluster_presence.csv')))

        # Plot "Present" genes per cluster as bar plot
        fig,ax = plt.subplots(figsize=(3,2), dpi=300)
        clust_ct = df.groupby('Cluster').sum()['Present'].reset_index()
        # Get total genes in each cluster
        clust_ct['Total'] = df.groupby('Cluster').count()['gene'].values
        clust_ct['Percent'] = clust_ct['Present'] / clust_ct['Total'] * 100
        sns.barplot(x='Cluster', y='Percent', data=clust_ct, ax=ax)
        plt.xticks(rotation=90)
        plt.title('Percent of HMO genes detected\nin each cluster')
        plt.ylabel('Percent')
        for i in range(0,101,25):
            plt.axhline(i, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=0)
        sns.despine()
        plt.savefig(os.path.join(self.output_dir,'HMO_percent_gene_detection.pdf'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir,'HMO_percent_gene_detection.png'), dpi=300, bbox_inches='tight')

        logger.info('Saved HMO gene detection plot to {}'.format(os.path.join(self.output_dir,'HMO_percent_gene_detection.pdf')))

        return


# Example usage
def main():
    pass

if __name__ == '__main__':
    main()