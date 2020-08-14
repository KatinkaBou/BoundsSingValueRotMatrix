# -*- coding: utf-8 -*-
"""
    Created on 06/01/2020
    Last Modification 06/02/2020
    
    @creator: Corentin JEUDY
    
    ###################################
    #   Negacyclic Singular Values    #
    ###################################
    """
#-- Example --#
# ..: python3 negacyclic_matrices.py

#-- Import --#
import numpy as np
from scipy.linalg import toeplitz
#------------#

def negacyclic(vec):
	"""
        Generates the negacyclic matrix whose first column is vec
        Input: vec (n-array)
        Output: mat (nxn-array)
    """ 
	row = np.zeros(len(vec))
	row[0] = vec[0]
	row[1:] = -vec[-1:0:-1]
	return toeplitz(vec, row)

def extreme_sing_values(mat):
	"""
        Returns the extreme singular values of mat
        Input: mat (nxn-array)
        Output: smax, smin (float)
	"""
	u,s,v = np.linalg.svd(mat) # SVD over Q
	return np.max(s), np.min(s)

def random_testing(dim, trials, bound):
	"""
        Verification of the Von Neumann bounds for extreme singular values via random testing
        Input: dim_pow (int, dimension), trials (int, number of random matrices tested), 
		bound (int, bound for matrices entries)
        Output: low, up (int, number of failures)
	"""
	smax_average=0
	smin_average=0
	low, up = 0, 0
	for i in range(trials):
		vec = np.random.randint(bound, size=dim) # Generating uniformly random binary vector
		#print(vec)
		if np.linalg.norm(vec, ord=0)==0:
			vec[np.random.randint(dim)] = 1
		matrix = negacyclic(vec)
		#print(matrix)
		smax, smin = extreme_sing_values(matrix)
		#print(smax,smin)
		smax_average=smax_average+smax
		smin_average=smin_average+smin
        
        
		upper_bound = 10*dim**(1/2) # Von Neumann upper bound for smax
		lower_bound = dim**(-1/2)/10 # Von Neumann lower bound for smin
		tight_lower_bound = dim**(-1/2)
        
		if smin < lower_bound: 
			low += 1
		if smax > upper_bound: 
			up += 1
		smax_average=smax_average/trials
		smin_average=smin_average/trials
	print(
	'''Dimension {0}:
        - Number of failures:
        (depass lower bound) {1}/{2}
        (depass upper bound) {3}/{2}
        - Average smin: {4}
        - Average smax: {5}
        - Lower bound: {6}
        - Tight lower bound: {7}
        {8}
	'''.format(dim,low,trials,up,round(smin_average,5),round(smax_average,5),lower_bound,tight_lower_bound, 60*'*'))
	return low, up

if __name__ == "__main__":
	""" Parameters """
	dim_pow = 11 # Matrix dimension up to 2^(dim_pow - 1)
	trials = 1000 # Number of generated matrices
	bound = 2 # Bound for the matrices' entries
    
	lows = np.zeros(dim_pow)
    
	for i in range(dim_pow):
		low, up = random_testing(2**i, trials, bound)
		lows[i] = low

	print(np.min(lows))
