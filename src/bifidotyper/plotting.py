import os
import glob
from pathlib import Path
from .logger import logger
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import gzip
import sys

class PlotUtils:
    def __init__(self,args,sylph_profile:str,sylph_query:str,hmo_genes:str,output_dir:str='plots'):
        """
        Initialize Plotting utility with configurable directory paths.
        
        Args:
            args (argparse.Namespace): Command line arguments
            sylph_profile (str): Sylph profile result (tsv)
            sylph_query (str): Sylph query result (tsv)
            hmo_genes (str): HMO gene quantification table
        """
        self.args = args
        self.sylph_profile = sylph_profile
        self.sylph_query = sylph_query
        self.hmo_genes = glob.glob(hmo_genes)
        self.read_counts = self.retrieve_bowtie2_read_counts()
        self.output_dir = output_dir

        os.makedirs(self.output_dir,exist_ok=True)
        
        # Ensure files exist
        for file_path in [self.sylph_profile,
                            self.sylph_query]:
                if not os.path.exists(file_path):
                    raise FileNotFoundError(f"File not found: {file_path}") 
    
    def retrieve_bowtie2_read_counts(self):
        # Results should be at hmo_quantification/*.bowtie2.log
        bowtie2_logs = glob.glob('hmo_quantification/*.bowtie2.log')
        if len(bowtie2_logs) == 0:
            raise FileNotFoundError("Bowtie2 log files not found")
        read_counts = pd.DataFrame(columns=['Sample','Total_Reads','Aligned_Reads'])
        for bowtie2_log in bowtie2_logs:
            sample_name = os.path.basename(bowtie2_log).replace('.bowtie2.log','')
            with open(bowtie2_log) as f:
                for line in f:
                    if 'reads; of these' in line:
                        total_reads = int(line.split(' ')[0])
                    elif 'overall alignment rate' in line:
                        aligned_reads = int(float(line.split(' ')[0].replace('%',''))/100 * total_reads)
            if self.args.paired_end:
                # Bowtie2 reports for pairs, not reads, so double
                total_reads *= 2
                aligned_reads *= 2
            read_counts.loc[len(read_counts)] = [sample_name,total_reads,aligned_reads]
        
        return read_counts

    def plot_sylph_profile(self):
        """
        Plot Sylph profile results using matplotlib.
        """
        
        pdf = pd.read_csv(self.sylph_profile,sep='\t')

        pdf['Strain'] = pdf['Contig_name'].apply(lambda x: x[:40]+'...' if len(x) > 40 else x)
        pdf['Sample'] = pdf['Sample_file'].apply(lambda x: os.path.basename(x).replace('.fastq.gz','').replace(self.args.r1_suffix,'').replace(self.args.r2_suffix,''))
        strains = pdf['Strain'].unique()
        strain_colors_dict = {strain: color for strain, color in zip(strains, plt.cm.hsv([i / len(strains) for i in range(len(strains))]))}

        ### Taxonomic Abundance Heatmap ###
        heatmap_data = pdf.pivot_table(
            index='Sample', 
            columns='Strain', 
            values='Taxonomic_abundance', 
            aggfunc='first'
        )
        
        # Determine optimal figure size based on data dimensions
        num_rows, num_cols = heatmap_data.shape
        base_size = 0.7  # Size per cell
        fig_width = max(8, num_cols * base_size)
        fig_height = max(6, num_rows * base_size)
        
        plt.figure(figsize=(fig_width, fig_height))
        
        # Create heatmap with masked values
        mask = heatmap_data.isna()
        
        # Custom color palette with ascending purple
        cmap = sns.color_palette('Purples', as_cmap=True)
        
        # Generate heatmap
        sns.heatmap(
            heatmap_data, 
            annot=True,  # Show values
            fmt='.2f',   # Two decimal places
            cmap=cmap,   # Color palette
            mask=mask,   # Mask missing values
            cbar_kws={
                'label': 'Taxonomic Abundance (%)',
                'shrink': 0.8
            },
            linewidths=0.5,  # Add grid lines
            square=True,     # Make cells square
        )
        
        plt.title('Taxonomic Abundance Heatmap', fontsize=14, pad=20)
        plt.xlabel('Bifidobacterial Strains', fontsize=10, labelpad=10)
        plt.ylabel('Samples', fontsize=10, labelpad=10)
        
        # Rotate x and y axis labels for readability
        plt.xticks(rotation=45, ha='right', fontsize=8)
        plt.yticks(rotation=0, fontsize=8)
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir,f'taxonomic_abundance_heatmap.pdf'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir,f'taxonomic_abundance_heatmap.png'), dpi=300, bbox_inches='tight')

        ### Taxonomic Abundance ###

        fig,ax = plt.subplots(figsize=(10+pdf['Sample'].nunique()*0.2,6),dpi=300)
        strains = pdf['Strain'].unique()
        strain_colors = dict(zip(strains,sns.color_palette('tab20',len(strains)).as_hex()))

        # Plot with the correct color palette
        pdf.pivot(index='Sample', columns='Strain', values='Taxonomic_abundance').plot(
            kind='bar',
            stacked=True,
            ax=ax,
            color=[strain_colors[strain] for strain in pdf['Strain'].unique()],
        )

        ax.legend(loc='center left', bbox_to_anchor=(1, 0), frameon=False)
        plt.ylabel('Relative Abundance (%)')
        plt.xticks(rotation=45, ha='right')
        plt.title('Relative abundance of strains in each sample')
        plt.tight_layout()
        sns.despine()
        plt.savefig(os.path.join(self.output_dir,f'taxonomic_abundance_profile_barplot.pdf'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir,f'taxonomic_abundance_profile_barplot.png'), dpi=300, bbox_inches='tight')

        # Standardize colors - set similar genomes to similar colors, make lookup table with nice names
        
        ### Taxonomic Abundance Relative to Seq Depth ###

        fig,ax = plt.subplots(figsize=(10,4),dpi=300)
        strains = pdf['Strain'].unique()
        strain_colors = dict(zip(strains,sns.color_palette('tab20',len(strains)).as_hex()))
        rel_pdf = pdf.merge(self.read_counts,on='Sample')
        
        # Refactor normalization
        # Use effective coverage from sylph and genome size for each genome to calculate reads mapped to genome (eff_cov * genome_length)
        # Divide by read size (default 150? user input) to get number of reads mapped to each genome
        # Then normalize taxonomic abundance to that

        def get_read_lengths(fastq_gz):
            # open a gzipped fastq file and read the first 200 lines, checking the length of each read.
            lengths = []
            with gzip.open(fastq_gz, 'rt') as f:
                for i in range(1,200*4):
                    line = f.readline()
                    if (i-1) % 4 == 1:
                        lengths.append(len(line.strip()))
            # return the most common read length
            return max(set(lengths), key=lengths.count)
        
        if self.args.read_length is None:
            rel_pdf['Read_length'] = rel_pdf['Sample_file'].apply(get_read_lengths)
        else:
            rel_pdf['Read_length'] = self.args.read_length
        
        def get_genome_size(fna_gz):
            total_size = 0
            with gzip.open(fna_gz, 'rt') as f:
                for line in f.readlines():
                    if line.startswith('>'):
                        continue
                    total_size += len(line.strip())
            return total_size
    
        genome_sizes = {}
        for fna_gz in rel_pdf['Genome_file'].unique():
            genome_sizes[fna_gz] = get_genome_size(fna_gz)
        
        rel_pdf['Genome_size'] = rel_pdf['Genome_file'].apply(lambda x: genome_sizes[x])

        # rel_pdf['Taxonomic_abundance_norm'] = rel_pdf['Taxonomic_abundance'] * (rel_pdf['Aligned_Reads'] / rel_pdf['Total_Reads'])
        rel_pdf['n_reads_mapped'] = rel_pdf['Eff_cov'] * rel_pdf['Genome_size'] / rel_pdf['Read_length']
        
        rel_pdf['n_reads_mapped_out_of_all_reads'] = rel_pdf['n_reads_mapped'] / rel_pdf['Total_Reads']

        # 1. # reads mapped to each genome
        for sample in rel_pdf['Sample'].unique():
            plot_df = rel_pdf[rel_pdf['Sample'] == sample]
            sorted_data = plot_df.sort_values('n_reads_mapped', ascending=False)
            fig, ax = plt.subplots(figsize=(4,sorted_data.shape[0]*0.5), dpi=300)
            sns.barplot(data=sorted_data, y='Strain', x='n_reads_mapped', hue='Strain', ax=ax, palette=strain_colors)

            ax.legend(loc='center left', bbox_to_anchor=(1, 0), frameon=False)
            plt.xlabel('# of Reads Mapped')
            # plt.xticks(rotation=45, ha='right')
            plt.title(f'{sample}\nNumber of reads mapped to each strain')
            plt.tight_layout()
            sns.despine()

            plt.savefig(os.path.join(self.output_dir, f'{sample}_number_of_reads_per_genome_barplot.pdf'), dpi=300, bbox_inches='tight')
            plt.savefig(os.path.join(self.output_dir, f'{sample}_number_of_reads_per_genome_barplot.png'), dpi=300, bbox_inches='tight')

        # Combined version
        fig,ax = plt.subplots(figsize=(10,4),dpi=300)
        sns.barplot(data=rel_pdf, x='Sample', y='n_reads_mapped', hue='Strain', ax=ax, palette=strain_colors)

        ax.legend(loc='center left', bbox_to_anchor=(1, 0), frameon=False)
        plt.ylabel('# of Reads Mapped')
        plt.xticks(rotation=45, ha='right')
        plt.title('Number of reads mapped to each strain')
        plt.tight_layout()
        sns.despine()
        plt.savefig(os.path.join(self.output_dir,f'number_of_reads_per_genome_combined_barplot.pdf'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir,f'number_of_reads_per_genome_combined_barplot.png'), dpi=300, bbox_inches='tight')

        # 2. # reads mapped to each genome normalized to total reads in sample
        fig,ax = plt.subplots(figsize=(10,4),dpi=300)
        rel_pdf.pivot(index='Sample', columns='Strain', values='n_reads_mapped_out_of_all_reads').plot(
            kind='bar',
            stacked=True,
            ax=ax,
            color=[strain_colors[strain] for strain in rel_pdf['Strain'].unique()],
        )

        ax.legend(loc='center left', bbox_to_anchor=(1, 0), frameon=False)
        plt.ylabel('Relative Abundance of Reads Mapped to Genome')
        plt.xticks(rotation=45, ha='right')
        plt.title('Relative abundance of reads mapped to each genome')
        plt.tight_layout()
        sns.despine()
        plt.savefig(os.path.join(self.output_dir,f'relative_abundance_of_reads_per_genome_barplot.pdf'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir,f'relative_abundance_of_reads_per_genome_barplot.png'), dpi=300, bbox_inches='tight')


        rel_pdf.to_csv('/home/bebr1814/playground/plots/debug.tsv',sep='\t')


    def plot_sylph_query(self):

        qdf = pd.read_csv(self.sylph_query,sep='\t')
        qdf['Strain'] = qdf['Contig_name'].apply(lambda x: x[:40]+'...' if len(x) > 40 else x)
        qdf['Sample'] = qdf['Sample_file'].apply(lambda x: os.path.basename(x).replace('.fastq.gz','').replace(self.args.r1_suffix,'').replace(self.args.r2_suffix,''))
        strains = qdf['Strain'].unique()
        strain_colors_dict = {strain: color for strain, color in zip(strains, plt.cm.hsv([i / len(strains) for i in range(len(strains))]))}

        ### Containment Indices ###
        # This makes a separate plot for each sample
        for sample in qdf['Sample'].unique():
            fig = self.containment_indices_barplot_horiz(qdf[qdf['Sample'] == sample])
            fig.savefig(os.path.join(self.output_dir,f'{sample}_containment_indices.pdf'), dpi=300, bbox_inches='tight')
            fig.savefig(os.path.join(self.output_dir,f'{sample}_containment_indices.png'), dpi=300, bbox_inches='tight')


    def containment_indices_barplot_horiz(self,df):
        n_genomes = df['Strain'].nunique()
        sample = df['Sample'].unique()[0]
        fig,ax = plt.subplots(figsize=(n_genomes/2,4),dpi=300)
        df['numerator'] = df['Containment_ind'].apply(lambda x: x.split('/')[0]).astype(int)
        df['denominator'] = df['Containment_ind'].apply(lambda x: x.split('/')[1]).astype(int)
        df['containment_index_float'] = df['numerator']/df['denominator']
        df.sort_values('containment_index_float',ascending=False,inplace=True)
        color = '#000'
        sns.barplot(data=df,x='Strain',y='denominator',color=color,label='Genome Size (bp)',ax=ax, fill=False)
        sns.barplot(data=df,x='Strain',y='numerator',color=color,label='Genome Coverage in Sample',ax=ax)
        ax.set_ylabel('Genome Size (bp)')
        ax.set_xlabel('Genomes')
        ax.set_title(f'{sample}\nContainment Indices')
        ax.set_xticklabels(ax.get_xticklabels(),rotation=45,ha='right')
        ax.legend(frameon=False,bbox_to_anchor=(1,1))
        sns.despine()
        return fig

