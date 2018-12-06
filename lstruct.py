# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 18:55:15 2018

@author: Jacob
"""

import numpy as np
import json

def matrix_rounder(matrix):
    matrix_copy = np.copy(matrix)
    with np.nditer(matrix_copy, op_flags =['readwrite']) as iterator:
        for x in iterator:
            c = float(x)        # Rounds each entry to the nearest integer
            x[...] = int(round(c))
            if x[...] < -1:     # These two conditionals correct any outlier 
                x[...] = -1     # values ( >1 or <-1) to the appropriate value
            elif x[...] > 1:    # (1 or -1, respectively).
                x[...] = 1
    return matrix_copy

def lsq_fit_pre(mat):
    mat_star = np.dot(np.linalg.inv(np.dot(mat.transpose(), mat)), mat.transpose())
    return mat_star
# If we want to approximate a column Ai as a linear combination of columns in Uk
# (i.e. Ai = Uk*vect), we can find a least-squares estimate z of vect with 
# z = Inverse(Uk_transpose * Uk) * Uk_transpose * Ai. The above function 
# precomputes the matrix product mat_star to be used in finding the vector z for
# each column Ai.

def approx_error(ref, est):
    c = ref-est
    n = np.count_nonzero(c)
    delta = n/(c.shape[0]*c.shape[1])
    return delta            
# The above function finds the difference between the estimate matrix and the 
# input (data) matrix, counts the number of non-zero entries, and divides it by
# the size of the matrix. This is the reported error of approximation delta.

def lin_struct(matrix, delt):
    u, sig, vT = np.linalg.svd(matrix, full_matrices=False)
    k = 0
    done = False
    output = []
    while not done:
        k += 1
        Uk = u[:,:(k)]
        #print('uk =', Uk)
        Uk_star = lsq_fit_pre(Uk)   # Finds the lsq fit matrix for a subset of matrix
        #print('Uk_star =', Uk_star) # U with the first k columns.
        kth_approx_list= []
        for i in range(matrix.shape[1]):
            z = np.dot(Uk_star, matrix[:,i]) # Finishes the l. approx. of col i by mult. 
            coli_kth_approx = np.dot(Uk, z) # Uk_star by col i, and taking this output vector
            coli_kth_approx = coli_kth_approx.tolist() # z and operating on it again with Uk. The output vector
            kth_approx_list.append(coli_kth_approx) # is appended to a list of the columns of Ak.
        kth_approx_T = np.array(kth_approx_list)
        kth_approx = kth_approx_T.transpose()   # The list of Ak's columns is converted into a matrix
        kth_approx = matrix_rounder(kth_approx) # and transposed to get the appropriate orientation.
        delta = approx_error(matrix,kth_approx)  # The approximation error is calculated and compared
        if delta <= delt: # to the desired threshold value. If the error
            output.append(k)    # is lower than the threshold, the value of k and the
            output.append(delta) # observed error are saved and the while loop is broken.
            print('Threshold met.')
            done = True  
        elif k == len(u):
            print('Max k hit.')
            done = True
    return output


file = open('pops_all.txt')
pop_file = [[item for item in line.split(',')] for line in file]
pops = pop_file[0]
file2 = open('regions.txt')
region_file = [[item for item in line.split(',')] for line in file2]
regions = region_file[0]
delta_threshold = 0.1

for region in regions:
    for population in pops:
        q = region
        r = population
        print((population+'+'+region))
        with open(q+'_'+str(r)+'_encode.json') as f:
            m = json.load(f)
            A = np.array(m)
        print(A.shape)
        #output = lin_struct(A, delta_threshold)
        #dt = str(int(delta_threshold*100))
        #print('output =', output)
        #f = open((q+'_'+str(r)+'_lstruct'+dt+'.txt'),'w')
        #f.write('[k, delta_obs] \n')
        #f.write(str(output))
        #f.close()


#A = A.transpose()
#A = np.array([[1,0,1,1,1],[-1,1,-1,0,0],[0,1,0,0,1]])
#u, sig, vT = np.linalg.svd(A, full_matrices=False)

