import os
import glob
import subprocess
import typing
from pathlib import Path

class SylphError(Exception):
    """Custom exception for Sylph-related errors."""
    pass

class SylphUtils:
    def __init__(self, 
                 genome_sketch_dir: str = 'sylph_genome_sketches', 
                 fastq_sketch_dir: str = 'sylph_fastq_sketches', 
                 genome_query_dir: str = 'sylph_genome_queries'):
        """
        Initialize Sylph utility with configurable directory paths.
        
        Args:
            genome_sketch_dir (str): Directory to store genome sketches
            fastq_sketch_dir (str): Directory to store fastq sketches
            genome_query_dir (str): Directory to store genome query results
        """
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
            # print(f"Running command: {' '.join(command)}")
            result = subprocess.run(command, check=True, text=True, capture_output=True, shell=False)
            # print(result)
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
        # Ensure we're in the genome sketch directory
        os.chdir(self.genome_sketch_dir)
                
        # Construct sylph sketch command
        command = ['sylph', 'sketch',*genomes,'-t', str(threads), '-o', output_name]
        self._run_command(command)

        os.chdir('..')
        
        return os.path.join(self.genome_sketch_dir, f"{output_name}.syldb")
    
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
        # Ensure we're in the fastq sketch directory
        os.chdir(self.fastq_sketch_dir)
        
        # Construct sylph sketch command
        if fastq_se:
            command = ['sylph', 'sketch', *fastq_se, '-t', str(threads)]
        elif fastq_r1 and fastq_r2:
            command = ['sylph','sketch','-1',*fastq_r1,'-2',*fastq_r2,'-t',str(threads)]
        else:
            raise SylphError("Either fastq_se or fastq_r1 and fastq_r2 must be provided")
        
        self._run_command(command)

        os.chdir('..')
        
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
        # Ensure we're in the genome query directory
        os.chdir(self.genome_query_dir)
        
        # Construct sylph query command
        command = ['sylph', 'query'] + sylsp_files + [syldb_file, '-o', output_name]
        self._run_command(command)

        os.chdir('..')
        
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
        # Ensure we're in the genome query directory
        os.chdir(self.genome_query_dir)
        
        # Construct sylph profile command
        command = ['sylph', 'profile'] + sylsp_files + [syldb_file, '-o', output_name]
        self._run_command(command)
        
        os.chdir('..')
        
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
        
        print(f"Query result: {query_result}")
        print(f"Profile result: {profile_result}")
    
    except SylphError as e:
        print(f"Sylph processing error: {e}")

if __name__ == '__main__':
    main()