from IPython.display import display, HTML
import random
import string
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import itertools
import numpy as np
import networkx as nx
from scipy.optimize import basinhopping
from itertools import combinations


def calculate_recombination_frequency(data, i, j):
    recombinants = 0
    total = 0

    for genotype in data:
        genes = genotype.split("/")
        g1 = ''.join(sorted(genes[i], key=str.lower))
        g2 = ''.join(sorted(genes[j], key=str.lower))

        het1 = g1[0] != g1[1]
        het2 = g2[0] != g2[1]

        if het1 != het2:  
            recombinants += 1

        total += 1

    return recombinants / total if total else 0.0

def compute_recombination_matrix(df_genotypes):
    genotypes = df_genotypes["Genotype"].tolist()
    num_loci = len(genotypes[0].split("/"))
    loci_labels = list(string.ascii_uppercase)[:num_loci]

    matrix = pd.DataFrame(index=loci_labels, columns=loci_labels, dtype=float)

    for i in range(num_loci):
        for j in range(num_loci):
            if i == j:
                matrix.iloc[i, j] = 0.0
            elif pd.isna(matrix.iloc[i, j]):
                r = calculate_recombination_frequency(genotypes, i, j)
                matrix.iloc[i, j] = r
                matrix.iloc[j, i] = r  

    return matrix


def group_loci_by_linkage(recomb_matrix, threshold=46):
    G = nx.Graph()
    loci = list(recomb_matrix.columns)

    G.add_nodes_from(loci)

    for i in range(len(loci)):
        for j in range(i + 1, len(loci)):  
            locus1 = loci[i]
            locus2 = loci[j]
            r_val = recomb_matrix.loc[locus1, locus2]
            if r_val <= threshold:
                G.add_edge(locus1, locus2)
                print(f"Linking {locus1}-{locus2} with r={r_val:.3f}")

    groups = list(nx.connected_components(G))
    return groups
    
def get_chromosome_matrices(recomb_matrix, groups):
    chrom_matrices = {}
    for idx, group in enumerate(groups, 1):
        group = sorted(group)
        submatrix = recomb_matrix.loc[group, group]
        chrom_matrices[f"Chr{idx}"] = submatrix.round(1)
    return chrom_matrices
    
def display_chromosome_matrices(chrom_matrices):
    html_blocks = []
    for chr_name, matrix in chrom_matrices.items():
        html = f"<h4>{chr_name}</h4>" + matrix.to_html()
        wrapped = f"<div style='display: inline-block; margin-right: 30px; vertical-align: top'>{html}</div>"
        html_blocks.append(wrapped)

    display(HTML(''.join(html_blocks)))


def map_genes(chrom_matrix):
    loci = chrom_matrix.columns.tolist()
    n = len(loci)
    
    pairs = list(combinations(range(n), 2))
    

    D = np.array([chrom_matrix.iloc[i, j] for i, j in pairs])
    
    def reconstruction_loss(x):
        return sum(abs((abs(x[i] - x[j]) - D[idx])) for idx, (i, j) in enumerate(pairs))

    bounds = [(2, 48)] * n
    x0 = np.random.uniform(2, 48, size=n)

   
    minimizer_kwargs = {"method": "L-BFGS-B", "bounds": bounds}
    result = basinhopping(reconstruction_loss, x0, minimizer_kwargs=minimizer_kwargs, niter=1000)
    x_recovered = result.x

    print("Mapped positions:")
    for gene, pos in zip(loci, np.round(x_recovered, 2)):
         print(f"{gene}: {pos} cM")
    print("\nTotal absolute error:", round(result.fun, 3))
    
    return x_recovered