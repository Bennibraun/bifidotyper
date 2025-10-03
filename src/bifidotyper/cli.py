import argparse
import sys
from .processor import build_sample_dict
from .references import ReferenceManager
import glob
import platform
import os
import subprocess
import shutil
from .sylph import SylphUtils
from .hmo_genes import HMOUtils
from .plotting import PlotUtils
from .logger import logger
from .phylogenetic import PhylogeneticUtils
import tqdm

import warnings
warnings.filterwarnings("ignore")


disable_tqdm = not sys.stdout.isatty()  # Disable if output is redirected

# Ensure uncaught exceptions are logged with full tracebacks
def _log_uncaught_exceptions(exctype, value, tb):
    try:
        # Log full traceback
        logger.exception("Uncaught exception", exc_info=(exctype, value, tb))
    finally:
        # Also print the traceback to stderr so batch systems like SLURM capture it
        import traceback
        traceback.print_exception(exctype, value, tb, file=sys.stderr)

sys.excepthook = _log_uncaught_exceptions

def parse_args():
    parser = argparse.ArgumentParser(description="Process FASTQ files for bifidotyper.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-se', '--single-end', nargs='+', help="Single-end FASTQ files.")
    group.add_argument('-pe', '--paired-end', nargs='+', help="Paired-end FASTQ files (R1 and R2 files, supports wildcards).")
    
    # Add suffix arguments as a mutually inclusive group
    suffix_group = parser.add_argument_group('paired-end options')
    suffix_group.add_argument('--r1-suffix', help="Suffix for R1 files (only for paired-end mode)")
    suffix_group.add_argument('--r2-suffix', help="Suffix for R2 files (only for paired-end mode)")

    # parser.add_argument('-l', '--read-length', type=int, default=None, help="Read length for accurate plotting (only affects some plots).")
    
    parser.add_argument('-r', '--rpm-threshold', type=float, default=10, help="Minimum RPM threshold for HMO genes to be considered present (default: 10).")

    parser.add_argument('-t', '--threads', type=int, default=1, help="Number of threads to use for parallel processing.")

    parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose output.", default=False)
    
    args = parser.parse_args()
    
    # Validate suffix arguments are used together if provided
    if (args.r1_suffix and not args.r2_suffix) or (args.r2_suffix and not args.r1_suffix):
        parser.error("--r1-suffix and --r2-suffix must be used together")

    # Always have a concrete suffix value to pass into processors
    if not args.r1_suffix:
        args.r1_suffix = "_R1"
    if not args.r2_suffix:
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

    print(r'''
┏━━┓ ━━ ┏━┓ ━  ┏┓ ━━━ ┏┓ ━━━━  ━━━━━  ━━━
┃┏┓┃ ━  ┃┏┛  ━ ┃┃ ━  ┏┛┗┓  ━━━━━━━   ━━━━
┃┗┛┗┓┏┓┏┛┗┓┏┓┏━┛┃┏━━┓┗┓┏┛┏┓ ┏┓┏━━┓┏━━┓┏━┓
┃┏━┓┃┣┫┗┓┏┛┣┫┃┏┓┃┃┏┓┃ ┃┃ ┃┃ ┃┃┃┏┓┃┃┏┓┃┃┏┛
┃┗━┛┃┃┃ ┃┃ ┃┃┃┗┛┃┃┗┛┃ ┃┗┓┃┗━┛┃┃┗┛┃┃┃━┫┃┃
┗━━━┛┗┛ ┗┛ ┗┛┗━━┛┗━━┛ ┗━┛┗━┓┏┛┃┏━┛┗━━┛┗┛
━━  ━━━  ━━━━  ━━━━━  ━━ ┏━┛┃ ┃┃ ━━━━━ ━━
━━━  ━━━━━  ━━━   ━━━━━  ┗━━┛ ┗┛  ━━━━━━━
    ''')

    args = parse_args()

    print('Loading software and reference data...')

    # Expand any wildcards manually (Windows shells won't expand globs)
    se_files = None
    pe_files = None
    if args.single_end:
        se_files = []
        for pattern in args.single_end:
            expanded = glob.glob(pattern)
            se_files.extend(expanded if expanded else [pattern])
        sample_dict = build_sample_dict(single_end=se_files,
                                        r1_suffix=args.r1_suffix,
                                        r2_suffix=args.r2_suffix)
    elif args.paired_end:
        pe_files = []
        for pattern in args.paired_end:
            expanded = glob.glob(pattern)
            pe_files.extend(expanded if expanded else [pattern])
        sample_dict = build_sample_dict(paired_end=pe_files,
                                        r1_suffix=args.r1_suffix,
                                        r2_suffix=args.r2_suffix)
    
    if args.single_end:
        fastq_files = [file for sample in sample_dict.values() for file in sample['files'].values()]
    else:
        fastq_files_r1 = [sample['files']['R1'] for sample in sample_dict.values()]
        fastq_files_r2 = [sample['files']['R2'] for sample in sample_dict.values()]
    
    # Get the reference files
    refs = get_reference_files()

    logger.info("Processing FASTQ files with Sylph...")
    print('Processing FASTQ files with Sylph...')

    # Resolve external executables
    sylph_exec = shutil.which('sylph')
    salmon_exec = shutil.which('salmon')
    if not sylph_exec:
        logger.error("Could not find 'sylph' in PATH. Please install Sylph and ensure it is available.")
        print("Error: 'sylph' not found in PATH. Install via conda: conda install -c bioconda sylph")
        sys.exit(1)
    if not salmon_exec:
        logger.error("Could not find 'salmon' in PATH. Please install Salmon and ensure it is available.")
        print("Error: 'salmon' not found in PATH. Install via conda: conda install -c bioconda salmon")
        sys.exit(1)

    # Initialize Sylph utility
    sylph_u = SylphUtils(args=args, sylph_executable=sylph_exec)

    try:
        genome_db = refs['bifidobacteria_sketches']

        if args.single_end:
            read_sketches = sylph_u.sketch_reads(fastq_se=fastq_files, threads=args.threads)
        else:
            read_sketches = sylph_u.sketch_reads(fastq_r1=fastq_files_r1, fastq_r2=fastq_files_r2, threads=args.threads)

        print('Querying samples against genomes...')
        query_result = sylph_u.query_genomes(read_sketches, genome_db)
        profile_result = sylph_u.profile_genomes(read_sketches, genome_db)

        logger.info(f"Query result: {query_result}")
        logger.info(f"Profile result: {profile_result}")

    except Exception as e:
        logger.info(f"Sylph processing error: {e}")
        raise Exception(f"Sylph processing error: {e}")

    # Now run HMO quantification
    print('Detecting HMO genes...')

    def get_sample_name(fastq):
        name = os.path.basename(fastq)
        for ext in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
            if name.endswith(ext):
                name = name[: -len(ext)]
                break
        # Remove suffixes if present
        name = name.replace(args.r1_suffix, '').replace(args.r2_suffix, '')
        return name

    if args.single_end:
        missing_samples = []
        for fastq_se in fastq_files:
            sample_name = get_sample_name(fastq_se)
            if not os.path.exists(f'hmo_quantification/{sample_name}.salmon_counts_annotated.tsv') or not os.path.exists(f'hmo_quantification/{sample_name}.cluster_presence.tsv'):
                missing_samples.append(fastq_se)
    
        if len(missing_samples) > 0:
            for fastq_se in tqdm.tqdm(missing_samples, desc="Quantifying HMO genes", unit="samples", total=len(missing_samples), disable=disable_tqdm):
                sample_name = get_sample_name(fastq_se)
                HMOUtils(args=args,
                         salmon_executable=salmon_exec,
                         sample_name=sample_name,
                         genes_fasta=refs['bl_genes'],
                         hmo_annotations=refs['humann2_hmo'],
                         fastq_se=fastq_se,
                         output_dir='hmo_quantification',
                         threads=args.threads)
    else:
        missing_samples = []
        for fastq_r1, fastq_r2 in zip(fastq_files_r1, fastq_files_r2):
            sample_name = get_sample_name(fastq_r1)
            if not os.path.exists(f'hmo_quantification/{sample_name}.salmon_counts_annotated.tsv') or not os.path.exists(f'hmo_quantification/{sample_name}.cluster_presence.tsv'):
                missing_samples.append((fastq_r1, fastq_r2))
        
        if len(missing_samples) > 0:
            for fastq_r1, fastq_r2 in tqdm.tqdm(missing_samples, desc="Quantifying HMO genes", unit="samples", total=len(missing_samples), disable=disable_tqdm):
                sample_name = get_sample_name(fastq_r1)
                HMOUtils(args=args,
                         salmon_executable=salmon_exec,
                         sample_name=sample_name,
                         genes_fasta=refs['bl_genes'],
                         hmo_annotations=refs['humann2_hmo'],
                         fastq_pe1=fastq_r1,
                         fastq_pe2=fastq_r2,
                         output_dir='hmo_quantification',
                         threads=args.threads)
    
    # Run plotting
    print('Plotting results...')
    logger.info("Plotting results...")

    plot_u = PlotUtils(args=args,
                        sylph_profile='sylph_genome_queries/genome_profile.tsv',
                        sylph_query='sylph_genome_queries/genome_query.tsv',
                        hmo_genes='hmo_quantification/*.salmon_counts_annotated.tsv',
                        genomes_df=refs['genomes_df'],
                        output_dir='plots')
    
    plot_u.plot_hmo_genes()

    plot_u.plot_sylph_profile()

    plot_u.plot_sylph_query()

    # Generate and plot phylogenetic tree
    # print('Generating phylogenetic tree...')
    # logger.info("Generating phylogenetic tree...")
    try:
        phylo_utils = PhylogeneticUtils(
            genomes_df=refs['genomes_df'],
            sylph_profile='sylph_genome_queries/genome_profile.tsv',
            output_dir='plots'
        )
        tree = phylo_utils.generate_phylogenetic_tree()
        phylo_utils.plot_cladogram(tree)
    except Exception as e:
        logger.warning(f"Skipping phylogenetic tree due to error: {e}")

    print('Done!')

if __name__ == "__main__":
    main()