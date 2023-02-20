# -*- coding: utf-8 -*-
"""
	Brief: Heuristics on the smallest singular value of discrete Gaussian 
    matrices over modules over cyclotomic rings. [Constant behavior]
"""

#-- Import --#
from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
from matplotlib import pyplot as plt
from scipy import linalg
import argparse
#-- End Import --#

#-- Parameter parser --#
parser = argparse.ArgumentParser(description='Heuristics on the smallest \
	singular value of discrete Gaussian matrices over modules (over \
	cyclotomic rings).')
parser.add_argument('--nu', metavar='nu', type=int, default=128, 
	help='Integer representing the nu-th cyclotomic field (default: 128).')
parser.add_argument('--d', metavar='d', type=int, default=4, 
	help='Module rank (default: 4).')
parser.add_argument('--nb', metavar='nb', type=int, default=100, 
	help='Number of sampled matrices (default: 100).')

#-- Module Discrete Gaussian Sampler --#
def sample_complex_matrices(n, d, sampler):
	"""
	Sample dxd complex matrices that correspond to the canonical
	embeddings of a discrete Gaussian matrix over a cyclotomic ring
	- input:
		(int)    n 	      -- Cyclotomic ring degree
		(int)	 d 		  -- Module rank
		(   )    sampler  -- Discrete Gaussian Lattice Sampler
	- output:
		(list)   matrices -- List of dxd complex matrices
	"""
	matrices = [zero_matrix(CDF, d, d) for _ in range(n//2)]
	for i in range(d):
		for j in range(d):
			v = sampler()/sqrt(2)
			for k in range(n//2):
				matrices[k][i,j] = v[2*k] + I*v[2*k+1]
			del v
	return matrices

def evaluate_prob(singular_values, epsilon, d):
	"""
	Evaluate the probability P[s_min(A) > epsilon / sqrt(d)]
	- input:
		(list) singular_values -- List of singular values of (embedded) 
								  discrete Gaussian matrices of size dxd
		(flt)  epsilon 		   
		(int)  d 			   -- Module rank (dimension of matrices)
	- output:
		(flt) 				   -- Heuristical value of 
								  P[s_min(A) > epsilon / sqrt(d)]
	"""
	lower = 0
	bound = epsilon*d**(-1/2)
	for smin in singular_values:
		if smin < bound:
			lower += 1
	return lower/len(singular_values)


if __name__ == '__main__':
	args = parser.parse_args()

	# Computing lattice basis of sigma_H(R)	
	K = CyclotomicField(args.nu)
	n = K.degree()
	"""
	Using K.minkowski_embedding(prec=None) does not generate a matrix in RDF
	and K.minkowski_embedding().change_ring(RDF) does not work for DGS for 
	some reason
	"""
	basis = zero_matrix(RDF, n,n)
	basis.set_block(0,0, K.minkowski_embedding())
	del K 

	# Set up parameters
	gamma_min = RDF(sqrt(args.d * log(n,2)))
	gamma_max = 5 * n * args.d
	nb_points = 50
	gamma_list = [gamma_min + k * (gamma_max - gamma_min)/nb_points for k in range(nb_points)]
	
	nb_trials = ceil(args.nb * 2/n)

	# Heuristics
	C_gamma_list, c_gamma_list = [], []
	for gamma in gamma_list:
		# Initializing Discrete Gaussian Sampler
		DGSampler = DiscreteGaussianDistributionLatticeSampler(basis, gamma)

		# Sampling and computing smallest singular values
		singular_values = []
		for _ in range(nb_trials):
			singular_values += [linalg.svd(A, compute_uv=False)[-1] for A in sample_complex_matrices(n, args.d, DGSampler)]

		epsilons = [k * RDF(sqrt(n)) for k in range(100)]
		probabilities = [evaluate_prob(singular_values, eps, args.d) for eps in epsilons]

		# Evaluate constants
		C_gamma_list.append((probabilities[-1] - probabilities[0])/(epsilons[-1]-epsilons[0]))
		c_gamma_list.append(probabilities[0]**(1/args.d))

	plt.plot(gamma_list, C_gamma_list, label='$C_\gamma$')
	plt.plot(gamma_list, c_gamma_list, label='$c_\gamma$')
	plt.plot(gamma_list, [1/gamma**(3/2) for gamma in gamma_list], label='$1/\gamma^{3/2}$')
	plt.xlabel('$\gamma$')
	plt.legend(loc='best')
	plt.show()