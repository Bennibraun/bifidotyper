"""Module for processing and validating input files and samples.

This module provides utilities for handling FASTQ files and organizing them into
sample groups for both single-end and paired-end sequencing data.
"""

import os
import re
from typing import List, Dict, Optional, Union
from pathlib import Path
from collections import defaultdict
from .logger import logger

VALID_EXTENSIONS = ('.fastq', '.fastq.gz', '.fq', '.fq.gz')

def validate_files(files: List[str]) -> None:
    """Validate that all files exist and have correct extensions.
    
    Args:
        files: List of file paths to validate.
        
    Raises:
        FileNotFoundError: If any file is missing.
        ValueError: If any file has an invalid extension.
    """
    for file in files:
        if not os.path.isfile(file):
            raise FileNotFoundError(f"File not found: {file}")
            
        if not any(file.lower().endswith(ext) for ext in VALID_EXTENSIONS):
            raise ValueError(
                f"Invalid file extension for {file}. Must be one of: {', '.join(VALID_EXTENSIONS)}"
            )

def get_base_name(filename: str, r1_suffix: str, r2_suffix: str) -> str:
    """Extract the base sample name from a FASTQ filename.
    
    Args:
        filename: The FASTQ file name/path.
        r1_suffix: Suffix identifying R1 files in paired-end mode.
        r2_suffix: Suffix identifying R2 files in paired-end mode.
        
    Returns:
        The base sample name with all suffixes removed.
    """
    basename = os.path.basename(filename)
    # Remove common suffixes first
    basename = re.sub(r'\.f(ast)?q(\.gz)?$', '', basename, flags=re.IGNORECASE)
    # Remove R1/R2 patterns using provided suffixes
    basename = basename.replace(r1_suffix, '').replace(r2_suffix, '')
    return basename

def make_absolute_path(file: Union[str, Path], base_dir: Union[str, Path]) -> str:
    """Convert a file path to an absolute path if it isn't already.
    
    Args:
        file: File path to convert.
        base_dir: Base directory to use for relative paths.
        
    Returns:
        An absolute path to the file.
    """
    file_path = Path(file)
    if not file_path.is_absolute():
        file_path = Path(base_dir) / file_path
    return str(file_path.resolve())

def build_sample_dict(
    single_end: Optional[List[str]] = None, 
    paired_end: Optional[List[str]] = None,
    r1_suffix: str = "_R1",
    r2_suffix: str = "_R2"
) -> Dict[str, Dict[str, Union[str, Dict[str, str]]]]:
    """Build a dictionary mapping sample names to their input files.
    
    Args:
        single_end: List of single-end FASTQ files.
        paired_end: List of paired-end FASTQ files (both R1 and R2).
        r1_suffix: Suffix identifying R1 files in paired-end mode.
        r2_suffix: Suffix identifying R2 files in paired-end mode.
        
    Returns:
        A dictionary mapping sample names to their file information.
        Format: {
            'sample_name': {
                'type': 'single-end' | 'paired-end',
                'files': {
                    'R1': '/path/to/R1.fastq',
                    'R2': '/path/to/R2.fastq'  # Only in paired-end mode
                }
            }
        }
        
    Raises:
        Exception: With descriptive message if inputs are invalid.
    """
    base_dir = os.getcwd()
    sample_dict: Dict[str, Dict[str, Union[str, Dict[str, str]]]] = {}
    
    # Validate input parameters
    if not single_end and not paired_end:
        raise ValueError("Must provide either single_end or paired_end files")
    if single_end and paired_end:
        raise ValueError("Cannot provide both single_end and paired_end files")
    if not isinstance(r1_suffix, str) or not isinstance(r2_suffix, str):
        raise TypeError("r1_suffix and r2_suffix must be strings")
    if r1_suffix == r2_suffix:
        raise ValueError("r1_suffix and r2_suffix must be different")

    try:
        if single_end:
            # Process single-end files
            logger.info(f"Processing {len(single_end)} single-end files")
            single_end = [make_absolute_path(file, base_dir) for file in single_end]
            validate_files(single_end)
            
            # Check for duplicate sample names
            sample_names = [get_base_name(f, r1_suffix, r2_suffix) for f in single_end]
            duplicate_names = {name for name in sample_names if sample_names.count(name) > 1}
            if duplicate_names:
                raise ValueError(f"Found duplicate sample names: {', '.join(duplicate_names)}")
            
            # Build sample dictionary
            for file in single_end:
                sample_name = get_base_name(file, r1_suffix, r2_suffix)
                sample_dict[sample_name] = {
                    'type': 'single-end',
                    'files': {'R1': file}
                }
                logger.debug(f"Added single-end sample {sample_name}: {file}")
                
        elif paired_end:
            # Process paired-end files
            logger.info(f"Processing {len(paired_end)} paired-end files")
            paired_end = [make_absolute_path(file, base_dir) for file in paired_end]
            validate_files(paired_end)
            
            # Group paired-end files by their base name
            r1_files: Dict[str, str] = {}
            r2_files: Dict[str, str] = {}
            unmatched_files: List[str] = []
            
            # First pass: sort files into R1 and R2 groups
            for file in paired_end:
                basename = get_base_name(file, r1_suffix, r2_suffix)
                if r1_suffix in file and r2_suffix in file:
                    raise ValueError(
                        f"File contains both R1 and R2 suffixes: {file}. "
                        f"Check suffix settings: R1='{r1_suffix}', R2='{r2_suffix}'"
                    )
                elif r1_suffix in file:
                    if basename in r1_files:
                        raise ValueError(f"Duplicate R1 file found for sample {basename}")
                    r1_files[basename] = file
                elif r2_suffix in file:
                    if basename in r2_files:
                        raise ValueError(f"Duplicate R2 file found for sample {basename}")
                    r2_files[basename] = file
                else:
                    unmatched_files.append(file)
            
            if unmatched_files:
                raise ValueError(
                    f"The following files do not match R1 ({r1_suffix}) or R2 ({r2_suffix}) "
                    f"patterns: {', '.join(unmatched_files)}"
                )
            
            # Second pass: match pairs and validate
            all_samples = set(r1_files.keys()) | set(r2_files.keys())
            if not all_samples:
                raise ValueError(
                    f"No valid paired-end files found. Ensure files contain '{r1_suffix}' "
                    f"or '{r2_suffix}' in their names."
                )
                
            # Check for missing pairs
            unpaired_r1 = set(r1_files.keys()) - set(r2_files.keys())
            unpaired_r2 = set(r2_files.keys()) - set(r1_files.keys())
            if unpaired_r1:
                raise ValueError(f"Missing R2 files for samples: {', '.join(unpaired_r1)}")
            if unpaired_r2:
                raise ValueError(f"Missing R1 files for samples: {', '.join(unpaired_r2)}")
            
            # Build final sample dictionary
            for sample in sorted(all_samples):
                sample_dict[sample] = {
                    'type': 'paired-end',
                    'files': {
                        'R1': r1_files[sample],
                        'R2': r2_files[sample]
                    }
                }
                logger.debug(
                    f"Added paired-end sample {sample}: "
                    f"R1={r1_files[sample]}, R2={r2_files[sample]}"
                )
                
        logger.info(f"Successfully processed {len(sample_dict)} samples")
        return sample_dict
        
    except Exception:
        # Log full traceback to aid debugging in batch environments
        logger.exception("Error processing input files")
        raise
