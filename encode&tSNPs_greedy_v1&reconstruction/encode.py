"""
Created on Wed 5 Dec 2018

@author: Sergei Kotelnikov, Jake Crosser
"""

import glob
import json
import numpy as np

regions = ["SORCS3", "PAH", "17q25"]
SNPs_number = 1000


# Jake's function
def encode(snp_str, ref, alt):
    snp_list = snp_str.split('|')
    if ref < alt:
        if snp_list == ['0', '0']:
            val = 1
        elif snp_list == ['1', '1']:
            val = -1
        else:
            val = 0
    else:
        if snp_list == ['0', '0']:
            val = -1
        elif snp_list == ['1', '1']:
            val = 1
        else:
            val = 0
    return val

with open("populations/indexdata.txt") as f:
    all_populations = [line.split(";")[0] for line in f.read().strip().split("\n")]
    population_files = ["populations/{}.txt".format(population) for population in all_populations]
population_dict = {}
for population_file in population_files:
    with open(population_file) as f:
        name, full_name, individuals = f.read().strip().split("\n")
        for individual in individuals.split():
            population_dict[individual] = name

for region in regions:
    print(region)
    vcf_file = glob.glob("{}/*.vcf".format(region))[0]

    with open(vcf_file) as f:
        file = f.read().strip().split("\n")
        comments = [line for line in file if line[:2] == "##"]
        content = [line for line in file if line[:2] != "##"]
        header = content[0]
        content_lines = content[1:]

    all_individuals = header.split("\t", 9)[-1].split("\t")
    population_of_individuals = [population_dict[individual] for individual in all_individuals]

    minority_freqs = []
    for content_line in content_lines:
        chrom, pos, id, ref, alt, qual, filter, info, format, variations = content_line.strip().split("\t", 9)
        info = [x.split("=") for x in info.split(";") if "=" in x]
        info = {key: value for key, value in info}
        info["AF"] = [float(val) for val in info["AF"].split(",")]
        if len(info["AF"]) == 1 and info["VT"] == "SNP":
            minority_freqs.append(info["AF"][0])
    minority_freqs.sort(reverse=True)
    minority_freq_cutoff = minority_freqs[SNPs_number]
    print("Minority frequency cutoff:", minority_freq_cutoff)


    SNPs = {population: [] for population in all_populations}
    for content_line in content_lines:
        chrom, pos, id, ref, alt, qual, filter, info, format, variations = content_line.strip().split("\t", 9)
        variations = variations.split("\t")
        info = [x.split("=") for x in info.split(";") if "=" in x]
        info = {key: value for key, value in info}
        info["AF"] = [float(val) for val in info["AF"].split(",")]
        if info["VT"] == "SNP" and len(info["AF"]) == 1 and info["AF"][0] > minority_freq_cutoff:  # we can add here more criteria
            variations = [encode(variation, ref, alt) for variation in variations]

            SNP_column = {population: [] for population in all_populations}
            for variation, population in zip(variations, population_of_individuals):
                SNP_column[population].append(variation)

            for population in all_populations:
                SNPs[population].append(SNP_column[population])

    for population, population_SNPs in SNPs.items():
        with open("{}/{}_{}_encode.json".format(region, region, population), "w") as f:
            population_SNPs_transposed = np.transpose(np.array(population_SNPs))
            population_SNPs_transposed = [list(x) for x in population_SNPs_transposed]
            json.dump(population_SNPs_transposed, f)