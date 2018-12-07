"""
Created on Wed 5 Dec 2018

@author: Sergei Kotelnikov
"""

import json
import numpy as np
from numpy import linalg as la

def project(matrix, subspace):  # Matrix projection onto subspace
    coefs, residuals, rank, singular_values = la.lstsq(subspace, matrix)
    return subspace @ coefs

errors = [0.1, 0.01]
regions = ["PAH", "17q25", "SORCS3"]
with open("populations/indexdata.txt") as f:
    populations = [line.split(";")[0] for line in f.read().strip().split("\n")]

for region in regions:
    for population in populations:
        print(region, population)
        with open("{}/{}_{}_encode.json".format(region, region, population)) as f:
            SNP_matrix = np.matrix(json.load(f))
        for error in errors:
            S = set()
            E = SNP_matrix
            done = False
            n = SNP_matrix.shape[1]
            while not done:
                f = enumerate([1 - la.norm(project(E, E[:, i]))**2 / la.norm(E)**2 for i in range(n)])
                S.add(min(f, key=lambda x: x[1])[0])
                tmp = project(SNP_matrix, SNP_matrix[:, list(S)])
                E = SNP_matrix - tmp
                SNP_matrix_approx = np.sign(np.around(tmp))
                delta = np.count_nonzero(SNP_matrix - SNP_matrix_approx) / SNP_matrix.size
                if delta <= error:
                    done = True
            with open("{}/{}_{}_{}_tSNPs.json".format(region, region, population, error), "w") as f:
                json.dump(list(S), f)
