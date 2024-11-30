import os
import sys
import pkg_resources
import pathlib
from .logger import logger

class ReferenceManager:
    """Manages access to reference files for the bifidotyper package"""
    
    def __init__(self):
        # Get the package's root directory
        self.package_root = pathlib.Path(__file__).parent.parent.absolute()
        self.reference_dir = self.package_root / 'data' / 'reference'
        
        # Ensure reference directory exists
        if not self.reference_dir.exists():
            raise FileNotFoundError(f"Reference directory not found at {self.reference_dir}")
        
        # Dictionary to store reference file paths
        self._reference_files = {
            'humann2_hmo': self.reference_dir / 'humann2_HMO_annotation.csv',
            'bl_genome': self.reference_dir / 'CP001095.1_genome.fasta',
            'bl_coding_seqs': self.reference_dir / 'CP001095.1_coding_seqs.gb',
            'bl_genes': self.reference_dir / 'CP001095.1_genes.gff3',
            'genomes_dir': self.reference_dir / 'genomes_ncbi_gtdb_drep95',
        }
        
        # Validate all reference files exist
        self._validate_references()

        logger.info("References validated.")
    
    def _validate_references(self):
        """Ensure all required reference files are present"""
        for name, path in self._reference_files.items():
            if not path.exists():
                raise FileNotFoundError(f"Required reference file '{name}' not found at {path}")
    
    def get_reference_path(self, reference_name):
        """Get the path to a specific reference file"""
        if reference_name not in self._reference_files:
            raise ValueError(f"Unknown reference '{reference_name}'. Available references: {list(self._reference_files.keys())}")
        return str(self._reference_files[reference_name])
    
    @property
    def available_references(self):
        """List all available reference files"""
        return list(self._reference_files.keys())

# Usage in your main program:
def main():
    try:
        ref_manager = ReferenceManager()
        # Use reference files in your program
        ref1_path = ref_manager.get_reference_path('ref1')
        logger.info(f"Using reference file: {ref1_path}")
        
    except FileNotFoundError as e:
        logger.info(f"Error: {e}")
        logger.info("Please ensure the package is properly installed with reference files.")
        sys.exit(1)