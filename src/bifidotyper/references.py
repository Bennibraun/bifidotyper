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
        self.reference_dir = self.package_root / 'bifidotyper' / 'data' / 'reference'
        self.bin_dir = self.package_root / 'bifidotyper' / 'data' / 'bin'
        
        # Ensure reference directory exists
        if not self.reference_dir.exists():
            raise FileNotFoundError(f"Reference directory not found at {self.reference_dir}")
        
        # Dictionary to store reference file paths
        self._reference_files = {
            'humann2_hmo': self.reference_dir / 'humann2_HMO_annotation.csv',
            'bl_genome': self.reference_dir / 'CP001095.1_genome.fasta',
            'bl_genes': self.reference_dir / 'CP001095.1_gene_sequences.fasta',
            'bifidobacteria_sketches': self.reference_dir / 'bifidobacteria_sketches.syldb',
            'genomes_df': self.reference_dir / 'genomes.csv',
        }

        self._bin_files = {
            'sylph_aarch64-any-linux': self.bin_dir / 'sylph_aarch64-any-linux',
            'sylph_arm64-any-darwin': self.bin_dir / 'sylph_arm64-any-darwin',
            'sylph_x86_64-any-linux': self.bin_dir / 'sylph_x86_64-any-linux',
            'sylph_x86_64-any-darwin': self.bin_dir / 'sylph_x86_64-any-darwin',
        }
        
        # Validate all reference files exist
        self._validate_references()

        logger.info("References validated.")
    
    def _validate_references(self):
        """Ensure all required reference files are present"""
        for name, path in self._reference_files.items():
            if not path.exists():
                raise FileNotFoundError(f"Required reference file '{name}' not found at {path}")
    
    def _validate_bins(self):
        """Ensure all required binary files are present"""
        for name, path in self._bin_files.items():
            if not path.exists():
                raise FileNotFoundError(f"Required binary file '{name}' not found at {path}")

    def get_reference_path(self, reference_name):
        """Get the path to a specific reference file"""
        if reference_name not in self._reference_files:
            raise ValueError(f"Unknown reference '{reference_name}'. Available references: {list(self._reference_files.keys())}")
        return str(self._reference_files[reference_name])

    def get_bin_path(self, bin_name):
        """Get the path to a specific binary file"""
        if bin_name not in self._bin_files:
            raise ValueError(f"Unknown binary '{bin_name}'. Available binaries: {list(self._bin_files.keys())}")
        return str(self._bin_files[bin_name])

    def get_reference_dir(self):
        """Get the path to the reference directory"""
        return str(self.reference_dir)

    def get_bin_dir(self):
        """Get the path to the binary directory"""
        return str(self.bin_dir)
    
    @property
    def available_references(self):
        """List all available reference files"""
        return list(self._reference_files.keys())
    
    @property
    def available_bins(self):
        """List all available binary files"""
        return list(self._bin_files.keys())

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