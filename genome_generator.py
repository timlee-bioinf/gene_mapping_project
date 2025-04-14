import random
import string
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import itertools
import numpy as np
import networkx as nx
from IPython.display import display

# this function generates a genome with 20 pieces of genes distributed on 3 chromosomes
def genomecreate(chromosomes=3, total_genes=20, seed=407):
    random.seed(seed)
    
    def distribute_genes(total_genes, chromosomes, min_genes=6):
        while True:
            sizes = [min_genes] * chromosomes
            remaining = total_genes - sum(sizes)
            for _ in range(remaining):
                sizes[random.randint(0, chromosomes - 1)] += 1
            if sum(sizes) == total_genes:
                return sizes

    gene_counts = distribute_genes(total_genes, chromosomes)

    all_genes = list(string.ascii_uppercase[:total_genes])
    random.shuffle(all_genes)

    chromosomes_data = []
    gene_idx = 0

    for chr_num, gene_count in enumerate(gene_counts, 1):
        chr_length = random.randint(1000, 2000)
        gene_names = all_genes[gene_idx:gene_idx + gene_count]

        min_pos = int(chr_length * 0.02)
        max_pos = int(chr_length * 0.98)

        quarter = (max_pos - min_pos) // 4
        required_positions = [
            random.randint(min_pos + i * quarter + 1, min_pos + (i + 1) * quarter)
            for i in range(4)
        ]

        remaining = gene_count - 4
        position_pool = list(range(min_pos, max_pos + 1))
        remaining_positions = random.sample(position_pool, remaining)

        unique_remaining = list(set(remaining_positions) - set(required_positions))
        while len(unique_remaining) < remaining:
            new_pos = random.randint(1, chr_length)
            if new_pos not in required_positions and new_pos not in unique_remaining:
                unique_remaining.append(new_pos)

        positions = sorted(required_positions + unique_remaining[:remaining])
        genes = list(zip(gene_names, positions))

        chrA = [(g, pos) for g, pos in genes]
        chrB = [(g, pos) for g, pos in genes]

        df = pd.DataFrame({
            "Chromosome": f"Chr{chr_num}",
            "Gene": gene_names,
            "Position": positions,
            "ChrA_Before": [g.upper() for g, _ in chrA],
            "ChrB_Before": [g.lower() for g, _ in chrB],
        })
        chromosomes_data.append((f"Chr{chr_num}", df, chr_length))
        gene_idx += gene_count

    return chromosomes_data


def print_genome(genome_data):
    for chr_name, df, chr_length in genome_data:
        print(f"\n=== {chr_name} | Length: {chr_length} bp ===")
        print(df[["Gene", "Position", "ChrA_Before", "ChrB_Before"]].to_string(index=False))