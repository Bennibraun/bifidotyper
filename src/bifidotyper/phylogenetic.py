import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio.Phylo.BaseTree import Clade
from Bio.SeqIO import parse
from .logger import logger

class PhylogeneticUtils:
    def __init__(self, genomes_df: str, sylph_profile: str, output_dir: str = 'plots'):
        self.genomes_df = pd.read_csv(genomes_df)
        self.sylph_profile = pd.read_csv(sylph_profile, sep='\t')
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def generate_phylogenetic_tree(self):
        # Extract strains and genome file paths
        strains = self.genomes_df['Label'].unique()
        genome_files = dict(zip(self.genomes_df['Label'], self.genomes_df['Genome_file']))
        ani_matrix = self._build_ani_matrix(strains, genome_files)

        # Construct the phylogenetic tree using UPGMA
        constructor = DistanceTreeConstructor()
        logger.info(f"Generating phylogenetic tree using UPGMA method.")
        tree = constructor.upgma(ani_matrix)

        # Save the tree in Newick format
        tree_path = os.path.join(self.output_dir, 'phylogenetic_tree.newick')
        Phylo.write(tree, tree_path, 'newick')
        logger.info(f"Phylogenetic tree saved to {tree_path}")

        return tree

    def _build_ani_matrix(self, strains, genome_files):
        # Create a lower-triangular distance matrix based on an ANI-like score
        labels = list(strains)
        n = len(labels)
        if n < 2:
            # DistanceMatrix requires at least 2 taxa; duplicate one with zero distance
            labels = labels + [labels[0]]
            matrix = [[0.0], [0.0, 0.0]]
            return DistanceMatrix(labels, matrix)

        matrix = []
        for i in range(n):
            row = []
            for j in range(i + 1):
                if i == j:
                    row.append(0.0)
                else:
                    ani = self._calculate_ani(genome_files[labels[i]], genome_files[labels[j]])
                    # Convert similarity to distance
                    row.append(max(0.0, 1.0 - ani))
            matrix.append(row)
        return DistanceMatrix(labels, matrix)

    def _calculate_ani(self, genome_file_1, genome_file_2, k: int = 8):
        """Approximate ANI using Jaccard index of k-mer sets.

        This avoids heavy alignments and external tools, providing a stable
        fallback that produces values in [0,1].
        """
        try:
            kmers1 = set()
            for rec in parse(genome_file_1, "fasta"):
                s = str(rec.seq).upper()
                kmers1.update(s[i:i+k] for i in range(0, max(0, len(s) - k + 1)))
            kmers2 = set()
            for rec in parse(genome_file_2, "fasta"):
                s = str(rec.seq).upper()
                kmers2.update(s[i:i+k] for i in range(0, max(0, len(s) - k + 1)))
            if not kmers1 or not kmers2:
                return 0.0
            inter = len(kmers1 & kmers2)
            union = len(kmers1 | kmers2)
            return inter / union if union else 0.0
        except Exception as e:
            logger.error(f"Error calculating ANI between {genome_file_1} and {genome_file_2}: {e}")
            return 0.0

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