import os
import sys
import pkg_resources
import pathlib
from .logger import logger

class ReferenceManager:
    
    def __init__(self):
        # Get the package's root directory
        self.package_root = pathlib.Path(__file__).parent.parent.absolute()
        self.reference_dir = self.package_root / 'bifidotyper' / 'data' / 'reference'
        
        # Ensure reference directory exists
        if not self.reference_dir.exists():
            raise FileNotFoundError(f"Reference directory not found at {self.reference_dir}")
        
        # Dictionary to store reference file paths
        self._reference_files = {
            'humann2_hmo': self.reference_dir / 'humann2_HMO_annotation.csv',
            'bl_genome': self.reference_dir / 'CP001095.1_genome.fasta',
            'bl_genes': self.reference_dir / 'CP001095.1_gene_sequences.fasta',
            'genomes_df': self.reference_dir / 'genomes.csv',
            'bifidobacteria_sketches': self.reference_dir / 'bifidobacteria_sketches.syldb',
        }

        # Flag to track validation status
        self._references_validated = False

        # Validate all reference files exist
        self._validate_references()
        logger.info("References validated.")  # Log only once after validation

    def _validate_references(self):
        if self._references_validated:  # Skip validation if already done
            return
        for name, path in self._reference_files.items():
            if not path.exists():
                raise FileNotFoundError(f"Required reference file '{name}' not found at {path}")
        self._references_validated = True  # Mark as validated
    
    def get_reference_path(self, reference_name):
        if reference_name not in self._reference_files:
            raise ValueError(f"Unknown reference '{reference_name}'. Available references: {list(self._reference_files.keys())}")
        return str(self._reference_files[reference_name])

    def get_reference_dir(self):
        return str(self.reference_dir)

    @property
    def available_references(self):
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