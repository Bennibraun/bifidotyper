import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio.Phylo.BaseTree import Clade
from .logger import logger

class PhylogeneticUtils:
    def __init__(self, genomes_df: str, sylph_profile: str, output_dir: str = 'plots'):
        self.genomes_df = pd.read_csv(genomes_df)
        self.sylph_profile = pd.read_csv(sylph_profile, sep='\t')
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def generate_phylogenetic_tree(self):
        # Extract strains and ANI distances
        strains = self.genomes_df['Label'].unique()
        ani_matrix = self._build_ani_matrix(strains)

        # Construct the phylogenetic tree using UPGMA
        constructor = DistanceTreeConstructor()
        logger.info(f"Generating phylogenetic tree using UPGMA method.")
        tree = constructor.upgma(ani_matrix)

        # Save the tree in Newick format
        tree_path = os.path.join(self.output_dir, 'phylogenetic_tree.newick')
        Phylo.write(tree, tree_path, 'newick')
        logger.info(f"Phylogenetic tree saved to {tree_path}")

        return tree

    def _build_ani_matrix(self, strains):
        # Create a distance matrix based on ANI values
        matrix = []
        labels = []
        for strain in strains:
            labels.append(strain)
            distances = []
            for other_strain in strains:
                if strain == other_strain:
                    distances.append(0)
                else:
                    ani = self.genomes_df.loc[
                        (self.genomes_df['Label'] == strain) & 
                        (self.genomes_df['Label'] == other_strain), 
                        'Adjusted_ANI'
                    ]
                    distances.append(1 - ani.values[0] if not ani.empty else 1)
            matrix.append(distances)
        return DistanceMatrix(labels, matrix)

    def plot_cladogram(self, tree):
        # Plot the phylogenetic tree as a cladogram
        fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
        Phylo.draw(tree, axes=ax, do_show=False)
        plt.title('Phylogenetic Tree of Bifido Strains', fontsize=14)
        plt.tight_layout()
        plot_path = os.path.join(self.output_dir, 'phylogenetic_tree_cladogram.png')
        plt.savefig(plot_path, dpi=300)
        logger.info(f"Cladogram saved to {plot_path}")

# Example usage
def main():
    genomes_df = 'src/bifidotyper/data/reference/genomes.csv'
    sylph_profile = 'sylph_genome_queries/genome_profile.tsv'
    phylo_utils = PhylogeneticUtils(genomes_df, sylph_profile)
    tree = phylo_utils.generate_phylogenetic_tree()
    phylo_utils.plot_cladogram(tree)

if __name__ == '__main__':
    main()