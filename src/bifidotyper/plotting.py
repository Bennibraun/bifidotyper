import os
import glob
from .logger import logger
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.gridspec import GridSpec
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import gzip
import natsort
import palettable

class PlotUtils:
    def __init__(self,args,sylph_profile:str,sylph_query:str,hmo_genes:str,genomes_df:str,output_dir:str='plots'):
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
        self.genomes_df = pd.read_csv(genomes_df)
        self.output_dir = output_dir


        os.makedirs(self.output_dir,exist_ok=True)
        
        # Ensure files exist
        for file_path in [self.sylph_profile,
                            self.sylph_query]:
                if not os.path.exists(file_path):
                    raise FileNotFoundError(f"File not found: {file_path}") 

    def plot_sylph_profile(self):
        """
        Plot Sylph profile results using matplotlib.
        """
        
        pdf = pd.read_csv(self.sylph_profile,sep='\t')

        pdf['Sample'] = pdf['Sample_file'].apply(lambda x: os.path.basename(x).replace('.fastq.gz','').replace(self.args.r1_suffix,'').replace(self.args.r2_suffix,''))
        pdf = pdf.merge(self.genomes_df,how='left',on='Genome_file').drop_duplicates()
        pdf['Strain'] = pdf['Label']
        strains = pdf['Strain'].unique()

        # pdf.to_csv('/home/bebr1814/playground/plotting_df.csv')

        if pdf['Strain'].nunique() < 20:
            # use palettable.tableau.GreenOrange_12.mpl_colormap as the colormap
            cmap = ListedColormap(palettable.tableau.GreenOrange_12.mpl_colors)
            strain_colors = {strain: cmap(i) for i, strain in enumerate(pdf['Strain'].unique())}
        else:
            # Try using the colors from the pdf
            # 'Color' contains a unique hex code for each strain, grouped by ANI clustering
            strain_colors = {strain: pdf[pdf['Strain'] == strain]['Color'].values[0] for strain in pdf['Strain'].unique()}

        ### Taxonomic Abundance Heatmap ###
        heatmap_data = pdf.pivot_table(
            index='Sample', 
            columns='Strain', 
            values='Taxonomic_abundance', 
            aggfunc='first'
        )

        # Remove columns where the highest value is lower than 5
        heatmap_data = heatmap_data.loc[:, heatmap_data.max() >= 5]

        # Determine optimal figure size based on data dimensions
        num_rows, num_cols = heatmap_data.shape
        base_size = 0.7  # Size per cell
        fig_height = 6 + num_rows * 0.2
        fig_width = 4 + num_cols * 0.5

        fig,ax = plt.subplots(figsize=(fig_width,fig_height),dpi=300)

        # Order the columns from highest to lowest total abundance
        heatmap_data = heatmap_data[heatmap_data.sum().sort_values(ascending=False).index]
        
        # Sort by sample
        heatmap_data = heatmap_data.reindex(natsort.natsorted(heatmap_data.index.tolist()))

        # Create heatmap with masked values
        mask = heatmap_data.isna()

        # Custom color palette with ascending purple
        cmap = sns.color_palette('Purples', as_cmap=True)

        # Generate heatmap
        sns.heatmap(
            heatmap_data, 
            annot=True,  # Show values
            fmt='.1f',   # Rounded values
            cmap=cmap,   # Color palette
            mask=mask,   # Mask missing values
            cbar_kws={
                'label': 'Taxonomic Abundance (%)',
                'shrink': 0.2
            },
            linewidths=0.5,  # Add grid lines
            # square=True,     # Make cells square
            ax=ax,
        )

        plt.title('Taxonomic Abundance Per Sample', fontsize=14, pad=20)
        plt.xlabel('Bifidobacterial Strains', fontsize=10, labelpad=10)
        plt.ylabel('Samples', fontsize=10, labelpad=10)

        # Rotate x and y axis labels for readability
        plt.xticks(rotation=45, ha='right', fontsize=8)
        plt.yticks(rotation=0, fontsize=8)

        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir,f'taxonomic_abundance_heatmap.pdf'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir,f'taxonomic_abundance_heatmap.png'), dpi=300, bbox_inches='tight')


        ### Taxonomic Abundance Full ###

        fig,ax = plt.subplots(figsize=(8,10+pdf['Sample'].nunique()*0.2),dpi=300)
        # strains = pdf['Strain'].unique()

        pivot_pdf = pdf.pivot(index='Sample', columns='Strain', values='Taxonomic_abundance')
        pivot_pdf = pivot_pdf.reindex(natsort.natsorted(pivot_pdf.index.tolist(),reverse=True))

        # add edges
        pivot_pdf.plot(
            kind='barh',
            stacked=True,
            ax=ax,
            edgecolor='grey',
            width=0.8,
	        color=[strain_colors[strain] for strain in pivot_pdf.columns],
        )

        ax.legend(loc='center left', bbox_to_anchor=(1, 0.9), frameon=False)
        plt.xlabel('Taxonomic Abundance (%)')
        plt.xticks(rotation=45, ha='right')
        plt.title('Taxonomic Abundance Per Sample')
        plt.tight_layout()
        sns.despine()
        plt.savefig(os.path.join(self.output_dir,f'taxonomic_abundance_profile_barplot.pdf'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir,f'taxonomic_abundance_profile_barplot.png'), dpi=300, bbox_inches='tight')


        ### Taxonomic Abundance per Sample ###
        for sample in pdf['Sample'].unique():
            sample_df = pdf[pdf['Sample'] == sample]
            sample_df = sample_df.sort_values('Taxonomic_abundance', ascending=False)
            plt.figure(figsize=(8, 4), dpi=300)
            sns.barplot(y='Strain', x='Taxonomic_abundance', hue='Strain', data=sample_df, palette=strain_colors)
            plt.title(f'Taxonomic Abundance\n{sample}')
            plt.ylabel('')
            plt.xlabel('Taxonomic Abundance (%)')
            sns.despine()
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir,f'{sample}_taxonomic_abundance.pdf'), dpi=300, bbox_inches='tight')
            plt.savefig(os.path.join(self.output_dir,f'{sample}_taxonomic_abundance.png'), dpi=300, bbox_inches='tight')

    def plot_hmo_genes(self):
        
        dfs = []

        for file in self.hmo_genes:
            df = pd.read_csv(file, sep='\t')
            label = os.path.basename(file).replace('.salmon_counts_annotated.tsv','')
            
            ### Cluster Completion Plots ###
            fig = self.cluster_completion_plot(df, label)
            fig.savefig(os.path.join(self.output_dir,f'{label}_HMO_cluster_completion.pdf'), dpi=300, bbox_inches='tight')
            fig.savefig(os.path.join(self.output_dir,f'{label}_HMO_cluster_completion.png'), dpi=300, bbox_inches='tight')
            
            df = df[['Name', 'Cluster', 'TPM']]
            df = df.rename(columns={'TPM': label})
            dfs.append(df)

        # Merge all DataFrames on 'Name' and 'Cluster'
        salmon_df = dfs[0]
        for df in dfs[1:]:
            salmon_df = pd.merge(salmon_df, df, on=['Name', 'Cluster'], how='outer')
        
        ### "RPM" Heatmap ###

        if salmon_df.shape[1] < 10:
            # Just sort by name, can't cluster with too few samples
            sorted_cols = natsort.natsorted(salmon_df.columns.tolist()[2:])
            salmon_df = salmon_df[['Name', 'Cluster'] + sorted_cols]

        else:
            # K-means clustering
            def find_optimal_k(data, max_k=12):
                silhouette_scores = []
                k_range = range(2, max_k + 1)
                
                for k in k_range:
                    kmeans = KMeans(n_clusters=k, random_state=42)
                    kmeans.fit(data)
                    silhouette_scores.append(silhouette_score(data, kmeans.labels_))
                
                optimal_k = k_range[np.argmax(silhouette_scores)]
                return optimal_k

            transposed = salmon_df.iloc[:,2:].T
            optimal_k = find_optimal_k(transposed)

            # Perform K-means clustering with the optimal K
            kmeans = KMeans(n_clusters=optimal_k, random_state=42)
            transposed['Cluster'] = kmeans.fit_predict(transposed)

            # Reorder the DataFrame based on clusters
            ordered_samples = transposed.sort_values('Cluster').index
            salmon_df = salmon_df[['Name','Cluster',*ordered_samples]]


        # Clusters and layout
        clusters = ['H1', 'H2', 'H3', 'H4', 'H5', 'Urease']
        # fig = plt.figure(figsize=(salmon_df.shape[1] * 0.5, salmon_df.shape[1] * 0.5), dpi=300)
        fig = plt.figure(figsize=(max(5,salmon_df.shape[1]*0.25),max(8,salmon_df.shape[1]*0.19)), dpi=300)

        # Distribute vertical space
        cluster_sizes = salmon_df.groupby('Cluster').size()
        cluster_sizes = cluster_sizes / cluster_sizes.sum()
        gs = GridSpec(len(clusters), 2, width_ratios=[0.98, 0.02], height_ratios=cluster_sizes)

        # Define colors
        color_palette = ListedColormap(palettable.cartocolors.qualitative.Antique_6.mpl_colors)
        cluster_colors = {cluster: color_palette(i) for i, cluster in enumerate(clusters)}

        for i, cluster in enumerate(clusters):
            # Heatmap
            ax = fig.add_subplot(gs[i, 0])
            cluster_df = salmon_df[salmon_df['Cluster'] == cluster].set_index('Name').drop('Cluster', axis=1)
            cmap = sns.light_palette(cluster_colors[cluster], as_cmap=True)
            
            # Compute normalization
            vmin, vmax = cluster_df.min().min(), cluster_df.max().max()
            norm = Normalize(vmin=vmin, vmax=vmax)
            
            # Draw heatmap
            sns.heatmap(
                cluster_df,
                ax=ax,
                cmap=cmap,
                cbar=False,
                norm=norm,
                linewidths=0.2,
                linecolor='white',
            )

            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=6)
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=6)
            ax.set_ylabel(cluster, rotation=0, labelpad=20, verticalalignment='center',horizontalalignment='right',color=cluster_colors[cluster],fontsize=10)
            if i != len(clusters) - 1:
                ax.set_xticks([])

            # Colorbar
            cbar_ax = fig.add_subplot(gs[i, 1])
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, cax=cbar_ax, orientation='vertical', fraction=0.6, pad=0.4)
            cbar.ax.tick_params(labelsize=6, width=0.5)  # Increase font size and tick width
            cbar.set_label('RPM', fontsize=6, labelpad=10)  # Adjust label size and spacing
            cbar.outline.set_linewidth(0.5)  # Set the top border line width
            cbar.outline.set_edgecolor('black')  # Set the top border line color

        plt.suptitle('Presence of HMO Genes by Cluster (RPM)', y=0.99, fontsize=10)
        plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout to account for title and spacing
        plt.savefig(os.path.join(self.output_dir,'hmo_gene_cluster_RPM_heatmap.pdf'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir,'hmo_gene_cluster_RPM_heatmap.png'), dpi=300, bbox_inches='tight')


        ### HMO Gene Cluster Heatmap (Binary Version) ###

        clusters = ['H1', 'H2', 'H3', 'H4', 'H5', 'Urease']
        # fig = plt.figure(figsize=(salmon_df.shape[1] * 0.5, salmon_df.shape[1] * 0.5), dpi=300)
        fig = plt.figure(figsize=(max(5,salmon_df.shape[1]*0.25),max(8,salmon_df.shape[1]*0.19)), dpi=300)

        # Distribute vertical space
        cluster_sizes = salmon_df.groupby('Cluster').size()
        cluster_sizes = cluster_sizes / cluster_sizes.sum()
        gs = GridSpec(len(clusters), 2, width_ratios=[0.98, 0.02], height_ratios=cluster_sizes)

        # Define colors
        color_palette = ListedColormap(palettable.cartocolors.qualitative.Antique_6.mpl_colors)
        cluster_colors = {cluster: color_palette(i) for i, cluster in enumerate(clusters)}

        for i, cluster in enumerate(clusters):
            # Heatmap
            ax = fig.add_subplot(gs[i, 0])
            cluster_df = salmon_df[salmon_df['Cluster'] == cluster].set_index('Name').drop('Cluster', axis=1)
            
            # Make the matrix binary
            cluster_df = cluster_df.map(lambda x: 1 if x > 0 else 0)
            cmap = ListedColormap(['white', cluster_colors[cluster]])
            
            # Draw heatmap
            sns.heatmap(
                cluster_df,
                ax=ax,
                cmap=cmap,
                cbar=False,
                linewidths=0.2,
                linecolor='white',
            )

            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=6)
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=6)
            ax.set_ylabel(cluster, rotation=0, labelpad=20, verticalalignment='center',horizontalalignment='right',color=cluster_colors[cluster],fontsize=10)
            if i != len(clusters) - 1:
                ax.set_xticks([])

        plt.suptitle('Presence of HMO Genes by Cluster (RPM)', y=0.99, fontsize=10)
        plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout to account for title and spacing
        plt.savefig(os.path.join(self.output_dir,'hmo_gene_cluster_RPM_heatmap_binary.pdf'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir,'hmo_gene_cluster_RPM_heatmap_binary.png'), dpi=300, bbox_inches='tight')


        ### Gene Cassette Plots ###

        for sample in salmon_df.columns.tolist()[2:]:
            fig = self.gene_cassette_plots(salmon_df, rpm_col=sample)
            fig.savefig(os.path.join(self.output_dir,f'{sample}_hmo_gene_cluster_RPM.pdf'), dpi=300, bbox_inches='tight')
            fig.savefig(os.path.join(self.output_dir,f'{sample}_hmo_gene_cluster_RPM.png'), dpi=300, bbox_inches='tight')


    def plot_sylph_query(self):

        qdf = pd.read_csv(self.sylph_query,sep='\t')
        qdf['Sample'] = qdf['Sample_file'].apply(lambda x: os.path.basename(x).replace('.fastq.gz','').replace(self.args.r1_suffix,'').replace(self.args.r2_suffix,''))
        qdf = qdf.merge(self.genomes_df,how='left',on='Genome_file').drop_duplicates()
        qdf['Strain'] = qdf['Label']
        strains = qdf['Strain'].unique()
        # strain_colors_dict = {strain: color for strain, color in zip(strains, plt.cm.hsv([i / len(strains) for i in range(len(strains))]))}


        ### Containment Indices ###
        # This makes a separate plot for each sample
        for sample in qdf['Sample'].unique():
            fig = self.containment_indices_barplot_horiz(qdf[qdf['Sample'] == sample])
            fig.savefig(os.path.join(self.output_dir,f'{sample}_containment_indices.pdf'), dpi=300, bbox_inches='tight')
            fig.savefig(os.path.join(self.output_dir,f'{sample}_containment_indices.png'), dpi=300, bbox_inches='tight')




    def containment_indices_barplot_horiz(self,df,max_bars=10):
        n_genomes = df['Strain'].nunique()
        sample = df['Sample'].unique()[0]
        fig,ax = plt.subplots(figsize=(n_genomes/2,4),dpi=300)
        df['numerator'] = df['Containment_ind'].apply(lambda x: x.split('/')[0]).astype(int)
        df['denominator'] = df['Containment_ind'].apply(lambda x: x.split('/')[1]).astype(int)
        df['containment_index_float'] = df['numerator']/df['denominator']
        df.sort_values('containment_index_float',ascending=False,inplace=True)
        color = '#000'
        sns.barplot(data=df.iloc[:max_bars,:],x='Strain',y='denominator',color=color,label='Genome Size (bp)',ax=ax, fill=False)
        sns.barplot(data=df.iloc[:max_bars,:],x='Strain',y='numerator',color=color,label='Genome Coverage in Sample',ax=ax)
        ax.set_ylabel('Genome Size (bp)')
        ax.set_xlabel('Genomes')
        ax.set_title(f'{sample}\nContainment Indices')
        ax.set_xticklabels(ax.get_xticklabels(),rotation=45,ha='right')
        ax.legend(frameon=False,bbox_to_anchor=(1,1))
        sns.despine()
        return fig


    def cluster_completion_plot(self, salmon_df, label, read_threshold=0):
        # Plot "Present" genes per cluster as bar plot
        salmon_df['Present'] = salmon_df['NumReads'] > read_threshold
        clust_ct = salmon_df.groupby('Cluster').sum()['Present'].reset_index()
        # Get total genes in each cluster
        clust_ct['Total'] = salmon_df.groupby('Cluster').count()['Name'].values
        clust_ct['Percent'] = clust_ct['Present'] / clust_ct['Total'] * 100
        # Define colors
        clusters = salmon_df['Cluster'].unique()
        color_palette = ListedColormap(palettable.cartocolors.qualitative.Antique_6.mpl_colors)
        cluster_colors = {cluster: color_palette(i) for i, cluster in enumerate(clusters)}
        fig,ax = plt.subplots(figsize=(3,2), dpi=300)
        sns.barplot(x='Cluster', y='Percent', hue='Cluster', data=clust_ct, ax=ax, palette=cluster_colors)
        plt.xticks(rotation=90)
        plt.title('Percent of HMO genes detected\nin each cluster')
        plt.ylabel('Percent')
        plt.xlabel(label)
        for i in range(0,101,25):
            plt.axhline(i, color='black', linestyle='--', alpha=0.5, linewidth=0.5, zorder=0)
        sns.despine()
        return fig

    
    def gene_cassette_plots(self, salmon_df, rpm_col):
        # Ensure the dataframe is sorted by Cluster and then by Name
        salmon_df_sorted = salmon_df[['Cluster', 'Name', rpm_col]].sort_values(['Cluster', 'Name'])
        salmon_df_sorted = salmon_df_sorted.rename(columns={rpm_col:'RPM'})

        # Get unique clusters in order
        clusters = salmon_df_sorted['Cluster'].unique()

        # Create a color palette with a unique color for each cluster
        # color_palette = ListedColormap(palettable.colorbrewer.qualitative.Set1_6.mpl_colors)
        color_palette = ListedColormap(palettable.cartocolors.qualitative.Antique_6.mpl_colors)
        cluster_colors = {cluster: color_palette(i) for i, cluster in enumerate(clusters)}

        # Get the maximum number of genes per cluster
        max_genes_per_cluster = salmon_df_sorted.groupby('Cluster').size().max()

        # Create the plot
        fig, ax = plt.subplots(len(clusters), 1, figsize=(15, len(clusters)*0.7), dpi=300, sharex=True)
        fig.subplots_adjust(hspace=0.5)

        # Iterate through clusters
        for i, cluster in enumerate(clusters):
            # Filter genes for this cluster
            cluster_genes = salmon_df_sorted[salmon_df_sorted['Cluster'] == cluster]

            # Plot each gene in the cluster
            for j, (_, gene) in enumerate(cluster_genes.iterrows()):
                # Determine color based on cluster
                color = cluster_colors[cluster]
                # Determine hue based on RPM
                cluster_max = cluster_genes['RPM'].max()
                color = sns.light_palette(color, as_cmap=True)(gene['RPM'] / cluster_max)
                norm = Normalize(vmin=0, vmax=cluster_max)


                # Create the pointed rectangle
                triangle_size = 0.15
                verts = [(j, 0), (j+(1-triangle_size), 0), (j+1, 0.5), (j+(1-triangle_size), 1), (j,1), (j+triangle_size, 0.5), (j, 0)]
                codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
                path = Path(verts, codes)
                patch = PathPatch(path, facecolor=color, edgecolor='black', linewidth=1)
                ax[i].add_patch(patch)

                # Label gene name
                ax[i].text(j + 0.5, -0.2, gene['Name'],#.replace('Blon_',''),
                        horizontalalignment='center',
                        verticalalignment='center',
                        fontsize=6)
                # Label RPM
                ax[i].text(j + 0.5, 0.5, int(gene['RPM']),
                        horizontalalignment='center',
                        verticalalignment='center',
                        fontsize=8)

            # Set axis limits and labels for each cluster
            ax[i].set_xlim(0, max_genes_per_cluster * 1.03)
            ax[i].set_ylim(0, 1)
            ax[i].set_xticks([])
            ax[i].set_yticks([])
            ax[i].set_ylabel(cluster, fontweight='bold', rotation=0, labelpad=10, verticalalignment='center',horizontalalignment='right',color=cluster_colors[cluster],fontsize=12)
            sns.despine(ax=ax[i], left=True, bottom=True, right=True, top=True)

            # # Add independent colorbar for each cluster
            # cbar_ax = fig.add_axes([0.90, 0.775 - i * 0.134, 0.02, 0.1])  # Adjust position as needed
            # cmap = sns.light_palette(cluster_colors[cluster], as_cmap=True)
            # sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            # sm.set_array([])
            # cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical')
            # cbar.ax.tick_params(labelsize=6)
            # cbar.set_label('RPM', fontsize=6)

        plt.suptitle('Presence of HMO Gene Clusters (RPM)')
        return fig