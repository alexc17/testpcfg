#approximateoracle.py
## We use this when we have a full grammar which 
# satisfies none of the conditions and we want to estimate 
# tree KLD if the result. 

# Here the results may not necessarily be greater than the true parameters
## So we need to renomralise and compute the tree kld.
# ALso how can we be sure that we will recover the true kernels.
# Examine this later.



import math

import utility
import wcfg
from collections import defaultdict
import argparse
import os.path
import logging

N_SAMPLES = 100
MAX_SAMPLES = 1000
EPSILON = 1e-8
RENYI = 5


def divergence(p,q):
	"""
	Estimate the divergence using the samples from 
	p and q 
	"""
	pass

class OracleApproximate:
	"""
	Here we take a grammar that doesn't satisfy the conditions and see what happens asymptotically.
	i.e. how do small perturbations affect the work.

	Also use this to validate the Reny divergence approximation.
	Evaluate how ? 

	Compute the estimate for every one. Take the partition function of the result.
	"""
	def __init__(self, target_pcfg):
		self.target_pcfg = target_pcfg
		self.terminals = set(target_pcfg.terminals)
		self.te = target_pcfg.terminal_expectations()
		self.pe = target_pcfg.production_expectations()
		self.sampler = wcfg.Sampler(self.target_pcfg)
		self.insider = wcfg.InsideComputation(self.target_pcfg)
		self.contexts_map = {}
		## map from nts to 

	def find_best_anchor(self, nt):
		posteriors = [ (a,self.pe.get(prod,0) /self.te[a]) for a in self.terminals ]
		return max(posteriors, key = lambda x : x[1])


	def find_approximate_kernel(self):
		self.kernel_map = {}
		posteriors = []
		for nt in self.target_pcfg.nonterminals:
			a,posterior = self.find_best_anchor(nt)
			self.kernel_map[nt] = a
			posteriors.append(posterior)
		## if min(posteriors) == 1 then it is all good.
		return min(posteriors)


	def sample_contexts(self, nt):
		a = self.kernel_map[nt]
		iwcfg = self.target_pcfg.intersect(a)
		iwcfg.renormalise()
		iwcfg.locally_normalise()
		asampler = wcfg.Sampler(iwcfg)
		contexts = []
	
		for _ in range(maxsamples):
			w = asampler.sample_string()
			positions = [ i for i,b in enumerate(w) if b == a ]
			assert len(positions) > 0
			position = numpy.random.choice(positions)
			left = w[:position]
			right = w[position+1:]
			contexts.append((tuple(left),tuple(right)))
		return contexts

	def test_unary_renyi(self, nt,b, contexts):
		a = self.kernel_map[nt]
		## contexts are sampled.
		for l,r in contexts:
			lar = l + (a,) + r
			lbr = l + (b,) + r
			plar = self.insider.inside_log_probability(lar)
			plbr = self.insider.inside_log_probability(lbr)



		


