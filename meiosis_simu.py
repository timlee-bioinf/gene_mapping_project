import random
from functions import *
from genome_generator import *

# Given a genome, simulates recombination, selection of chromosomes. Returns alleles found on resulting gamete.
def meiosis(chromosomes_data, seed=None):
    if seed is not None:
        random.seed(seed)

    recombinant_data = []

    for chr_name, df, chr_length in chromosomes_data:
        do_recombine = random.random() < 0.5

        if do_recombine:
            crossover = random.randint(1, chr_length)
            swap_above = random.random() < 0.5

            chr1 = list(zip(df["ChrA_Before"], df["Position"]))
            chr2 = list(zip(df["ChrB_Before"], df["Position"]))

            recombinant1 = []
            recombinant2 = []

            for (g1, pos), (g2, _) in zip(chr1, chr2):
                if (swap_above and pos > crossover) or (not swap_above and pos < crossover):
                    recombinant1.append((g2.lower(), pos))
                    recombinant2.append((g1.upper(), pos))
                else:
                    recombinant1.append((g1.upper(), pos))
                    recombinant2.append((g2.lower(), pos))

            df["ChrA_After"] = [g for g, _ in recombinant1]
            df["ChrB_After"] = [g for g, _ in recombinant2]

        else:
            df["ChrA_After"] = df["ChrA_Before"]
            df["ChrB_After"] = df["ChrB_Before"]

        recombinant_data.append((chr_name, None, None, df, chr_length))

    gamete_alleles = []
    for chr_name, _, _, df, _ in recombinant_data:
        selected = df["ChrA_After"] if random.random() < 0.5 else df["ChrB_After"]
        sorted_alleles = [allele for _, allele in sorted(zip(df["Position"], selected))]
        gamete_alleles.extend(sorted_alleles)

    return ",".join(sorted(gamete_alleles, key=str.lower))


# Simulates the fertilization. Given two gametes, returns the resulting offspring genotype
def fertilization(gamete1, gamete2):
    alleles1 = gamete1.split(",")
    alleles2 = gamete2.split(",")

    allele_dict = {}

    for a in alleles1 + alleles2:
        base = a.lower()
        if base not in allele_dict:
            allele_dict[base] = []
        allele_dict[base].append(a)

    child_genotype = []
    for base in sorted(allele_dict.keys()):
        pair = allele_dict[base]
        if len(pair) == 1:
            pair.append(pair[0])
        sorted_pair = sorted(pair, key=lambda x: (x.islower(), x))
        child_genotype.append("".join(sorted_pair))

    return "/".join(child_genotype)