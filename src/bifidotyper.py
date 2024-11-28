import os
import argparse
from collections import defaultdict
def parse_args():
    parser = argparse.ArgumentParser(description="Process FASTQ files for bifidotyper.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-se', '--single-end', nargs='+', help="Single-end FASTQ files.")
    group.add_argument('-pe', '--paired-end', nargs='+', help="Paired-end FASTQ files (R1 and R2 files, supports wildcards).")
    
    # Add suffix arguments as a mutually inclusive group
    suffix_group = parser.add_argument_group('paired-end options')
    suffix_group.add_argument('--r1-suffix', help="Suffix for R1 files (only for paired-end mode)")
    suffix_group.add_argument('--r2-suffix', help="Suffix for R2 files (only for paired-end mode)")
    
    args = parser.parse_args()
    
    # Validate suffix arguments are only used with paired-end mode
    if (args.r1_suffix or args.r2_suffix) and not args.paired_end:
        parser.error("--r1-suffix and --r2-suffix can only be used with paired-end mode (-pe)")
    
    # Validate suffix arguments are used together
    if bool(args.r1_suffix) != bool(args.r2_suffix):
        parser.error("--r1-suffix and --r2-suffix must be used together")
    
    # Set default suffixes if none provided in paired-end mode
    if args.paired_end and not args.r1_suffix:
        args.r1_suffix = "_R1"
        args.r2_suffix = "_R2"
    
    return args

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

def build_sample_dict(single_end=None, paired_end=None, r1_suffix="_R1", r2_suffix="_R2"):
    sample_dict = {}
    if single_end:
        validate_files(single_end)
        for file in single_end:
            sample_name = get_base_name(file, r1_suffix, r2_suffix)
            sample_dict[sample_name] = {
                'type': 'single-end',
                'files': {'R1': file}
            }
    elif paired_end:
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

def main():
    args = parse_args()
    
    # Build the sample dictionary
    if args.single_end:
        sample_dict = build_sample_dict(single_end=args.single_end, 
                                      r1_suffix=args.r1_suffix, 
                                      r2_suffix=args.r2_suffix)
    elif args.paired_end:
        sample_dict = build_sample_dict(paired_end=args.paired_end,
                                      r1_suffix=args.r1_suffix,
                                      r2_suffix=args.r2_suffix)
    
    # Print the resulting dictionary
    print(sample_dict)

    # Get the reference files
    refs = 


if __name__ == "__main__":
    main()