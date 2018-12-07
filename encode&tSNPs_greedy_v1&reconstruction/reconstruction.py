"""
Created on Wed 5 Dec 2018

@author: Sergei Kotelnikov
"""

import json
import numpy as np
from numpy import linalg as la

# 1 - reference population
# 2 - target population

regions = ["PAH", "17q25", "SORCS3"]
errors = [0.1, 0.01]
with open("populations/indexdata.txt") as f:
    populations = [line.split(";")[0] for line in f.read().strip().split("\n")]
    population_pairs = [(population_1, population_2) for population_1 in populations for population_2 in populations]

for region in regions:
    print(region)
    for error in errors:
        output = []
        for population_1, population_2 in population_pairs:
            with open("{}/{}_{}_encode.json".format(region, region, population_1)) as f:
                SNPs_matrix_1 = np.matrix(json.load(f))
            with open("{}/{}_{}_encode.json".format(region, region, population_2)) as f:
                SNPs_matrix_2 = np.matrix(json.load(f))
            with open("{}/{}_{}_{}_tSNPs.json".format(region, region, population_1, error)) as f:
                tSNP_indexes_1 = np.array(json.load(f))

            tSNPs_matrix_1 = SNPs_matrix_1[:, tSNP_indexes_1]
            tSNPs_matrix_2 = SNPs_matrix_2[:, tSNP_indexes_1]

            SNPs_matrix_2_predicted = np.sign(np.around(tSNPs_matrix_2 @ la.pinv(tSNPs_matrix_1) @ SNPs_matrix_1))

            delta = np.count_nonzero(SNPs_matrix_2_predicted - SNPs_matrix_2) / SNPs_matrix_2.size
            output.append([population_1, population_2, delta])
        with open("{}/{}_{}_reconstruction.json".format(region, region, error), "w") as f:
            json.dump(output, f)