# -*- coding: utf-8 -*-
"""
    Brief: Heuristics on the smallest singular value of discrete Gaussian 
    matrices over modules over cyclotomic rings. [Bound testing]
"""

#-- Import --#
from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
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
def sample_multiplication_matrix(n, sampler):
	"""
	Sample a discrete Gaussian vector from sampler and compute its
	multiplication matrix w.r.t. the Minkowski embedding over
	cyclotomic rings.
	- input:
		(int)    n 	      -- Cyclotomic ring degree
		(   )    sampler  -- Discrete Gaussian Lattice Sampler
	- output:
		(matrix) M_sigmaH -- (n x n) multiplication matrix
	"""
	v = sampler()/sqrt(2)
	M_sigmaH = block_diagonal_matrix(
		[matrix(RDF, 2, [[v[2*i], -v[2*i+1]],[v[2*i+1], v[2*i]]]) for i in range(n//2)]
	)
	del v 
	return M_sigmaH

def sample_structured_gaussian_matrix(n, d, sampler):
	return block_matrix(RDF, d, d, 
		[[sample_multiplication_matrix(n, sampler) for _ in range(d)] for _ in range(d)])


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

	# Setting parameters
	q = 8380417 # Modulus chosen in CRYSTALS Dilithium
	gamma = 2**(1/(2*args.d-1)) * n * q**((args.d-1)/(n*(2*args.d-1)))
	lower_bound = gamma * (n*args.d)**(-1/2) / 10

	# Initializing Discrete Gaussian Sampler
	DGSampler = DiscreteGaussianDistributionLatticeSampler(basis, gamma)
	
	# Heuristics
	smin_sum, lower, lowest = 0, 0, n * args.d * gamma
	for _ in range(args.nb):
		# Sampling matrix
		mat = sample_structured_gaussian_matrix(n, args.d, DGSampler)

		# Computing smallest singular value
		smin = linalg.svd(mat, compute_uv=False)[-1]
		del mat
		smin_sum += smin

		if smin < lower_bound: lower += 1
		if smin < lowest: lowest = smin
	
	# Output
	print('###### HEURISTICS RESULTS (%d-th cyclotomic field) ######'%args.nu)
	print('[*] Tested bound: gamma * (n*d)**(-1/2) / 10 for \
gamma = 2^{1 / (2d-1)} * n * q^{(d-1) / (n(2d-1))}')
	print(60*'_')
	print('''- Degree: %d
- Rank: %d 
- Dimension: %d
- Gaussian parameter: %.5f
- Lowest s_min: %.5f
- Average s_min: %.5f
- Tested Lower Bound: %.5f
- Number of failures (below lower bound): %d/%d''' % (n, args.d, n*args.d, gamma, 
	lowest, smin_sum/args.nb, lower_bound, lower, args.nb))
	print(60*'_')