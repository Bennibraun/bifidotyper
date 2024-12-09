import argparse
import sys
from .processor import build_sample_dict
from .references import ReferenceManager
import glob
import platform
import os
import subprocess
from .sylph import SylphUtils, SylphError
from .hmo_genes import HMOUtils, HMOError
from .plotting import PlotUtils
from .logger import logger
import tqdm

import warnings
warnings.filterwarnings("ignore")


disable_tqdm = not sys.stdout.isatty()  # Disable if output is redirected

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
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-g', '--genome-dir', type=str, default=None, help="(Optional) Directory containing genomes to use in place of the provided Bifidobacterial genomes. Please provide genomes as .fna or .fna.gz files.")
    group.add_argument('-s', '--genome-sketch', type=str, default=None, help="(Optional) Path to a pre-sketched genome database (.syldb) to use in place of the provided Bifidobacterial genomes. You can use `sylph sketch` to generate this.")

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

def get_bin_files():
    try:
        ref_manager = ReferenceManager()
        bin_files = {bin: ref_manager.get_bin_path(bin) for bin in ref_manager.available_bins}
        return bin_files
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
    bins = get_bin_files()

    # Check for Sylph
    try:
        subprocess.run(['sylph', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sylph = 'sylph'
    except FileNotFoundError:

        # We need to scan the cpu architecture and use the right binary
        system = platform.system()
        arch = platform.machine()

        if system == "Darwin":  # macOS
            if arch == "x86_64":
                sylph = bins['sylph_x86_64-any-darwin']
            elif arch == "arm64":
                sylph = bins['sylph_arm64-any-darwin']
        elif system == "Linux":
            if arch == "x86_64":
                sylph = bins['sylph_x86_64-any-linux']
            elif arch == "aarch64":
                sylph = bins['sylph_aarch64-any-linux']
        else:
            logger.error(f"Unfortunately, we do not have a precompiled Sylph binary for your architecture ({arch}). Please install Sylph manually and ensure it is in your PATH, then try again.")
            raise RuntimeError(f"Unfortunately, we do not have a precompiled Sylph binary for your architecture ({arch}). Please install Sylph manually and ensure it is in your PATH, then try again.")
        
        try:
            subprocess.run([sylph, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except:
            logger.error(f"Sylph failed to run. Please install it manually and ensure it si in your PATH, then try again.")
            raise RuntimeError(f"Sylph failed to run. Please install it manually and ensure it si in your PATH, then try again.")

    # Check for Salmon
    try:
        subprocess.run(['salmon', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        salmon = 'salmon'
    except:
        # Need to download Salmon because it's huge
        # wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz
        # tar -xvf salmon-1.10.0_linux_x86_64.tar.gz
        # cp salmon-1.10.0_linux_x86_64/bin/salmon 
        ref_manager = ReferenceManager()
        bin_dir = ref_manager.get_bin_dir()
        salmon = os.path.join(bin_dir,'salmon')
        if not os.path.exists(salmon):
            subprocess.run(['wget', 'https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz'], cwd=bin_dir, capture_output=True, text=True)
            logger.info(subprocess.run(['tar', '-xvf', 'salmon-1.10.0_linux_x86_64.tar.gz'], cwd=bin_dir, capture_output=True, text=True).stdout)
            logger.info(subprocess.run(['mv', 'salmon-latest_linux_x86_64/bin/salmon', bin_dir], cwd=bin_dir, capture_output=True, text=True).stdout)
            logger.info(subprocess.run(['mv', 'salmon-latest_linux_x86_64/lib', os.path.dirname(bin_dir)], cwd=bin_dir, capture_output=True, text=True).stdout)
            logger.info(subprocess.run(['rm', '-rf', 'salmon-latest_linux_x86_64', 'salmon-1.10.0_linux_x86_64.tar.gz'], cwd=bin_dir, capture_output=True, text=True).stdout)
            assert os.path.exists(salmon), "Salmon binary not found after download."

        try:
            subprocess.run([salmon, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except:
            logger.error(f"Salmon failed to run. Please install it manually and ensure it is in your PATH, then try again.")
            raise RuntimeError(f"Salmon failed to run. Please install it manually and ensure it is in your PATH, then try again.")
        
    if args.single_end:
        sample_dict = build_sample_dict(single_end=args.single_end, 
                                      r1_suffix=args.r1_suffix, 
                                      r2_suffix=args.r2_suffix)
    elif args.paired_end:
        sample_dict = build_sample_dict(paired_end=args.paired_end,
                                      r1_suffix=args.r1_suffix,
                                      r2_suffix=args.r2_suffix)
    
    if args.single_end:
        fastq_files = []
        for sample in tqdm.tqdm(sample_dict.values(), desc="Sketching reads", unit="samples", total=len(sample_dict), disable=disable_tqdm):
            fastq_files.append(sample['files'].values())
    else:
        fastq_files_r1, fastq_files_r2 = [], []
        for sample in tqdm.tqdm(sample_dict.values(), desc="Sketching reads", unit="samples", total=len(sample_dict), disable=disable_tqdm):
            fastq_files_r1.append(sample['files']['R1'])
            fastq_files_r2.append(sample['files']['R2'])
    
    # Get the reference files
    refs = get_reference_files()

    logger.info("Processing FASTQ files with Sylph...")
    print('Processing FASTQ files with Sylph...')

    run_sylph = True
    output_files = ['sylph_genome_queries/genome_profile.tsv', 'sylph_genome_queries/genome_query.tsv']
    # If all of these files exist, skip Sylph processing
    if all([os.path.exists(f) for f in output_files]):
        print("Sylph output files already exist. Skipping Sylph processing.")
        print("To force re-run, delete the output files and re-run the script.")
        run_sylph = False
    
    if run_sylph:
        # Initialize Sylph utility
        sylph = SylphUtils(args=args,sylph_executable=sylph)

        try:
            # Sketch genomes
            if args.genome_dir:
                print('Sketching genomes...')
                genomes = glob.glob(os.path.join(args.genome_dir, '*.fna.gz'))
                if len(genomes) == 0:
                    genomes = glob.glob(os.path.join(args.genome_dir, '*.fna'))
                    if len(genomes) == 0:
                        raise FileNotFoundError(f"No genome files found in {args.genome_dir}")
                genome_db = sylph.sketch_genomes(genomes=genomes, output_name='genome_sketches', threads=args.threads)
            elif args.genome_sketch:
                genome_db = args.genome_sketch
                assert os.path.exists(genome_db), f"Genome sketch file {genome_db} not found."
            else:
                genome_db = refs['bifidobacteria_sketches']
            
            if args.single_end:
                read_sketches = sylph.sketch_reads(fastq_se=fastq_files, threads=args.threads)
            else:
                read_sketches = sylph.sketch_reads(fastq_r1=fastq_files_r1, fastq_r2=fastq_files_r2, threads=args.threads)
            
            print('Querying the genome database...')

            # Query genomes
            query_result = sylph.query_genomes(read_sketches, genome_db)
            # Profile genomes
            profile_result = sylph.profile_genomes(read_sketches, genome_db)
            
            logger.info(f"Query result: {query_result}")
            logger.info(f"Profile result: {profile_result}")
        
        except SylphError as e:
            logger.info(f"Sylph processing error: {e}")
            raise SylphError(f"Sylph processing error: {e}")

    # Now run HMO quantification
    print('Detecting HMO genes...')

    def get_sample_name(fastq):
        return os.path.basename(fastq).replace('.fastq.gz','').replace(args.r1_suffix,'').replace(args.r2_suffix,'')

    if args.single_end:
        for fastq_se in tqdm.tqdm(fastq_files, desc="Quantifying HMO genes", unit="samples", total=len(fastq_files), disable=disable_tqdm):
            sample_name = get_sample_name(fastq_se)
            if all(os.path.exists(f) for f in [f'hmo_quantification/{sample_name}.salmon_counts_annotated.tsv',f'hmo_quantification/{sample_name}.cluster_presence.tsv']):
                print(f"Skipping {sample_name} as output files already exist.")
                continue
            HMOUtils(args=args,
                        salmon_executable=salmon,
                        sample_name=sample_name,
                        genes_fasta=refs['bl_genes'],
                        hmo_annotations=refs['humann2_hmo'],
                        fastq_se=fastq_se,
                        output_dir='hmo_quantification',
                        threads=args.threads)
    else:
        for fastq_r1, fastq_r2 in tqdm.tqdm(zip(fastq_files_r1, fastq_files_r2), desc="Quantifying HMO genes", unit="samples", total=len(fastq_files_r1), disable=disable_tqdm):
            sample_name = get_sample_name(fastq_r1)
            if all(os.path.exists(f) for f in [f'hmo_quantification/{sample_name}.salmon_counts_annotated.tsv',f'hmo_quantification/{sample_name}.cluster_presence.tsv']):
                print(f"Skipping {sample_name} as output files already exist.")
                continue
            HMOUtils(args=args,
                        salmon_executable=salmon,
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



    print('Done!')

if __name__ == "__main__":
    main()