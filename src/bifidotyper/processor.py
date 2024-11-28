import os
from collections import defaultdict

def validate_files(files):
    for file in files:
        assert os.path.isfile(file), f"File not found: {file}"

def get_base_name(filename, r1_suffix, r2_suffix):
    """Extract base name from filename, handling various R1/R2 patterns"""
    basename = os.path.basename(filename)
    # Remove common suffixes first
    basename = basename.replace('.fastq.gz', '').replace('.fastq', '').replace('.fq.gz', '').replace('.fq', '')
    # Remove R1/R2 patterns using provided suffixes
    basename = basename.replace(r1_suffix, '').replace(r2_suffix, '')
    return basename

def make_absolute_path(file, base_dir):
    """Convert a file path to an absolute path if it is not already"""
    if not file.startswith(("~", "/")):
        return os.path.join(base_dir, file)
    return file

def build_sample_dict(single_end=None, paired_end=None, r1_suffix="_R1", r2_suffix="_R2"):
    base_dir = os.getcwd()
    sample_dict = {}
    
    if single_end:
        single_end = [make_absolute_path(file, base_dir) for file in single_end]
        validate_files(single_end)
        for file in single_end:
            sample_name = get_base_name(file, r1_suffix, r2_suffix)
            sample_dict[sample_name] = {
                'type': 'single-end',
                'files': {'R1': file}
            }
    elif paired_end:
        paired_end = [make_absolute_path(file, base_dir) for file in paired_end]
        validate_files(paired_end)
        # Group paired-end files by their base name
        r1_files = {}
        r2_files = {}
        
        # First pass: sort files into R1 and R2 groups
        for file in paired_end:
            basename = get_base_name(file, r1_suffix, r2_suffix)
            if r1_suffix in file:
                r1_files[basename] = file
            elif r2_suffix in file:
                r2_files[basename] = file
            else:
                raise ValueError(f"Cannot determine if file is R1 or R2 using suffixes {r1_suffix} and {r2_suffix}: {file}")
        
        # Second pass: match pairs
        all_samples = set(r1_files.keys()) | set(r2_files.keys())
        for sample in all_samples:
            if sample not in r1_files:
                raise ValueError(f"Missing R1 file for sample: {sample}")
            if sample not in r2_files:
                raise ValueError(f"Missing R2 file for sample: {sample}")
            
            sample_dict[sample] = {
                'type': 'paired-end',
                'files': {
                    'R1': r1_files[sample],
                    'R2': r2_files[sample]
                }
            }
            
        if not sample_dict:
            raise ValueError(f"No valid paired-end files found. Ensure files contain '{r1_suffix}'/'{r2_suffix}' in their names.")
            
    return sample_dict
