# Gene Mapping Project
#### Code pieces:
- please use the `recombination_visualization.ipynb` to generate the plot for gene locations on chromosomes.
- please use `plots.ipynb` to generate the contour plot and convexity map for convexity analysis.
- the `experiment.ipynb` is intended for the author’s exploratory analysis and idea testing. Some analyses found in this notebook are not included in the report due to length constraints.
- `solver.ipynb` creates a genome with three chromosomes and twenty genes and simulates reproduction and birthing of two thousand offpsring. Then given the offspring datasets, predicts the position of each gene, return the result with total error.

#### Functions (updating):
Here is an ongoing list of functions used for this project:
##### in `functions.py`:
- `calculate_recombination_frequency` calculates the frequency of recombination happening in the chromosomes
- `compute_recombination_matrix` computes the recombination distance matrix
- `group_loci_by_linkage` creates a distance submatrix for every grouping of genes
- 
    
