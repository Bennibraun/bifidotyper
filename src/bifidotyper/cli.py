import argparse
import sys
from .processor import build_sample_dict
from .references import ReferenceManager
import glob
import os
from .sylph import SylphUtils, SylphError

def parse_args():
    parser = argparse.ArgumentParser(description="Process FASTQ files for bifidotyper.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-se', '--single-end', nargs='+', help="Single-end FASTQ files.")
    group.add_argument('-pe', '--paired-end', nargs='+', help="Paired-end FASTQ files (R1 and R2 files, supports wildcards).")
    
    # Add suffix arguments as a mutually inclusive group
    suffix_group = parser.add_argument_group('paired-end options')
    suffix_group.add_argument('--r1-suffix', help="Suffix for R1 files (only for paired-end mode)")
    suffix_group.add_argument('--r2-suffix', help="Suffix for R2 files (only for paired-end mode)")
    
    parser.add_argument('-t', '--threads', type=int, default=1, help="Number of threads to use for parallel processing.")
    
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

def get_reference_files():
    try:
        ref_manager = ReferenceManager()
        references = ref_manager.available_references
        reference_files = {ref: ref_manager.get_reference_path(ref) for ref in references}
        return reference_files
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please ensure the package is properly installed with reference files.")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

def main():
    args = parse_args()
    
    if args.single_end:
        sample_dict = build_sample_dict(single_end=args.single_end, 
                                      r1_suffix=args.r1_suffix, 
                                      r2_suffix=args.r2_suffix)
    elif args.paired_end:
        sample_dict = build_sample_dict(paired_end=args.paired_end,
                                      r1_suffix=args.r1_suffix,
                                      r2_suffix=args.r2_suffix)
    
    # Get the reference files
    refs = get_reference_files()
    
    # Initialize Sylph utility
    sylph = SylphUtils()
    
    try:
        # Sketch genomes
        genome_db = sylph.sketch_genomes(genomes=glob.glob(os.path.join(refs['genomes_dir'], '*.fna.gz')), output_name='genome_sketches', threads=args.threads)
        
        if args.single_end:
            fastq_files = []
            for sample in sample_dict.values():
                fastq_files.append(sample['files'].values())
            read_sketches = sylph.sketch_reads(fastq_se=fastq_files, threads=args.threads)
        else:
            # grab all R1 and R2 files separately, in order
            fastq_files_r1, fastq_files_r2 = [], []
            for sample in sample_dict.values():
                fastq_files_r1.append(sample['files']['R1'])
                fastq_files_r2.append(sample['files']['R2'])
            read_sketches = sylph.sketch_reads(fastq_r1=fastq_files_r1, fastq_r2=fastq_files_r2, threads=args.threads)
        
        print(read_sketches)
        # Query genomes
        query_result = sylph.query_genomes(read_sketches, genome_db)
        
        # Profile genomes
        profile_result = sylph.profile_genomes(read_sketches, genome_db)
        
        print(f"Query result: {query_result}")
        print(f"Profile result: {profile_result}")
    
    except SylphError as e:
        print(f"Sylph processing error: {e}")

if __name__ == "__main__":
    main()