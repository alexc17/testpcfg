# Learner that uses exact queries to construct a grammar.

# This is useful to get an upper bound on the performance of an algorihthm 
#  using just samples from the distribution.

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


class OracleApproximate:
	"""
	Here we take a grammar that doesn't satisfy the conditions and see what happens asymptotically.
	i.e. how do small perturbations affect the work.

	Also use this to validate the Reny divergence approximation.
	Evaluate how ? FIXME
	"""
	def __init__(self, target_pcfg):
		self.target_pcfg = target_pcfg
		self.terminals = set(target_pcfg.terminals)
		self.te = target_pcfg.terminal_expectations()
		self.pe = target_pcfg.production_expectations()
		self.sampler = wcfg.Sampler(self.target_pcfg)
		self.renyi = RENYI
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


	def test_unary(self, nt, b):
		a = self.kernel_map[nt]
		## sample contexts of a 
		## compute probs.
		## estimate divergence. 



class OracleLearner:
	"""
	We take a PCFG that is tight and with finite expectations.
	We assume that the grammar is sparse (i.e. nontrivial support)

	First test to see if it is anchored.


	Then test to see if the learned grammar is isomorphic to the target
	Then compute the distance.


	Method -- given the kernel, construct a grammar with with only K as terminals.
	Sample contexts from these until we have N per context.
	Compute the true value of each parameter.

	"""
	def __init__(self, target_pcfg):
		self.target_pcfg = target_pcfg
		self.target_wcfg = target_pcfg.convert_parameters_pi2xi()
		self.insider = wcfg.InsideComputation(target_pcfg)
		self.terminals = set(target_pcfg.terminals)
		self.te = target_pcfg.terminal_expectations()
		self.max_length = 2 *   len(self.target_pcfg.nonterminals)
		self.kernel = []
		self.nonterminal_map = {}
		self.parameters = {}
		self.kernel_grammar = None
		self.sampler = None
		self.asymptotic_grammar = None
		self.prod2context = {}

	def test(self):
		"""
		Test the target PCFG.

		return a Wcfg, or None if the grammar is not anchored. 
		"""
		ntmap = self.test_if_anchored()
		if len(ntmap) == 0:
			return None
		print(ntmap)
		self.nonterminal_map = ntmap
		self.kernel_map = { ntmap[a]: a for a in ntmap}
		self.kernel =  ntmap.values()
		print(self.kernel)
		self.kernel_grammar = self.make_kernel_grammar()
		self.sampler = wcfg.Sampler(self.kernel_grammar)
		self.test_all(N_SAMPLES, MAX_SAMPLES)
		return self.make_grammar()

		# delta_lexical = self.asymptotic_grammar.count_lexical() - self.target_pcfg.count_lexical()
		# delta_binary = self.asymptotic_grammar.count_binary() - self.target_pcfg.count_binary()
		# if delta_binary > 0:
		# 	for prod in self.asymptotic_grammar.productions:
		# 		if not prod in self.target_pcfg.productions:
		# 			print("EXTRA",prod)
		# pf = self.asymptotic_grammar.compute_partition_function_fp()
		# if pf[self.target_pcfg.start] == math.inf:
		# 	return ("DIVERGENT",math.inf, delta_lexical, delta_binary)
		# else:
		# 	return ("CONVERGENT",  pf[self.target_pcfg.start], delta_lexical, delta_binary)



	def test_if_anchored(self):
		"""
		Return a map from nonterminals to anchors.
		picking the most likely one.
		Return  empty map if it does not satisfy the condition.
		"""
		lhs_counter = defaultdict(list)
		answer = {}
		for prod in self.target_pcfg.productions:
			if len(prod) == 2:
				lhs_counter[prod[1]].append(prod[0])
		print(lhs_counter)
		for nt in self.target_pcfg.nonterminals:
			candidates = [ a for a in self.target_pcfg.terminals if lhs_counter[a] == [nt] ]
			if len(candidates) == 0:
				return {}
			else:
				answer[nt] =  max(candidates, key = lambda a : self.target_pcfg.parameters[ (nt,a)])
		return answer

	def make_kernel_grammar(self):
		"""
		Create a grammar with terminals only the unambiguous ones.
		This is then used to sample contexts from. since we only need contexts of these terminals
		and these ones will be the best ones.
		"""
		ua = wcfg.WCFG()
		ua.nonterminals = set(self.target_pcfg.nonterminals)
		for a in self.kernel:
			ua.terminals.add(a)
		ua.start = self.target_pcfg.start
		for prod in self.target_pcfg.productions:
			if len(prod) == 3:
				ua.productions.append(prod)
				ua.parameters[prod] = self.target_pcfg.parameters[prod]
			else:
				lhs = prod[0]
				prod2 = (lhs, self.nonterminal_map[lhs])
				if prod2 in ua.productions:
					ua.parameters[prod2] += self.target_pcfg.parameters[prod]
				else:
					ua.productions.append(prod2)
					ua.parameters[prod2] = self.target_pcfg.parameters[prod]
		ua.set_log_parameters()
		return ua

	def test_all(self, nsamples, maxsamples):
		for a in self.kernel:
			print("Testing",a)
			# xi = self.test_both(a, nsamples, maxsamples)
			# for prod in xi:
			# 	self.parameters[prod] = xi[prod]
			x1,x2 = self.test_both_smart(a, nsamples, maxsamples,verbose=True)
			for prod in x1:
				self.parameters[prod] = x1[prod]
			for prod in x2:
				self.parameters[prod] = x2[prod]

	def sample_contexts(self, nsamples, maxsamples):
		""" 
		Test unary rules with a on the lhs.
		"""	
		contexts = { a: set() for a in self.kernel}
		done = set()
		for _ in range(maxsamples):
			w = self.sampler.sample_string()
			if w in done:
				continue
			done.add(w)
			for i,a in enumerate(w):
				left = w[:i]
				right = w[i+1:]
				contexts[a].add((tuple(left),tuple(right)))
			if len(contexts) >= nsamples:
				break
		contexts = list(contexts)
		contexts.sort(key = lambda x : len(x[0]) + len(x[1]))
		return contexts

	def sample_contexts1(self, a, nsamples, maxsamples):
		""" 
		Test unary rules with a on the lhs.
		"""	
		contexts = set()
		for _ in range(maxsamples):
			w = self.sampler.sample_string()
			if len(w) < self.max_length:
				for i,b in enumerate(w):
					if a == b:
						left = w[:i]
						right = w[i+1:]
						contexts.add((tuple(left),tuple(right)))
						if len(contexts) >= nsamples:
							return contexts
		return contexts

	def sample_contexts_smart(self, a, nsamples, maxsamples):
		#print("Sampling contexts for ", a)
		iwcfg = self.target_pcfg.intersect(a)
		iwcfg.renormalise()
		iwcfg.locally_normalise()
		asampler = wcfg.Sampler(iwcfg)
		contexts = set()
	
		for _ in range(maxsamples):
			w = asampler.sample_string()
			for i,b in enumerate(w):
				if a == b:
					left = w[:i]
					right = w[i+1:]
					contexts.add((tuple(left),tuple(right)))
					if len(contexts) >= nsamples:
						return contexts

		return contexts

	def test_both_smart(self, a, n, maxsamples,verbose=False):
		nta = self.kernel_map[a]
		contexts = self.sample_contexts_smart(a,n,maxsamples)
		ll = 0.0
		for l,r in contexts:
			ll += len(l) + len(r)
		if verbose: print("Contexts:",len(contexts), ll/len(contexts))

		xi1,cc  = self.test_unary_smart(a,contexts)
		#print("halfway")
		xi2 = self.test_binary_smart(a, cc, contexts)
		
		
		return (xi1,xi2)


	def test_binary_smart(self, a, cc, contexts):
		# cc is c ahcarcterising context of a 
		# but that doesnt mean it will work all the time.
		xis = {}
		denoms = {}
		for l,r in contexts:
			denoms[(l,r)] = self.insider.inside_log_probability(l + (a,) + r)
		nta = self.kernel_map[a]
		for b in self.kernel:
			for c in self.kernel:
				ntb = self.kernel_map[b]
				ntc = self.kernel_map[c]
				if (nta,ntb,ntc) in self.target_wcfg.parameters:
					target = self.target_wcfg.parameters[ (nta,ntb,ntc)]
				else:
					target = 0.0
				best = math.inf
				try:
					bestcontext = None
					for l,r in contexts:
						denom = denoms[(l,r)]
						num = self.insider.inside_log_probability(l + (b,c) + r)
						ratio = math.exp(num - denom) * self.te[a]/(self.te[b] * self.te[c])
						if ratio < best:
							best = ratio
							bestcontext = (l,r)
						
						
					xis[(a,b,c)] = best
					self.prod2context[(a,b,c)] = bestcontext
				except utility.ParseFailureException:
					# if we catch an exception then we have a zero so we omit it.
					self.prod2context[(a,b,c)] = (l,r)
					
		return xis

	def test_unary_smart(self, a, contexts):
		"""
		Do it more efficiently
		"""

		denoms = {}
		cc = None
		nta = self.kernel_map[a]
		for l,r in contexts:
			n = 0
			for b in self.kernel:
				nt = self.kernel_map[b]
				try:
					lp =  self.insider.inside_log_probability(l + (b,) + r)
					p = math.exp(lp) / self.target_pcfg.parameters[(nt,b)]
					n += 1
				except utility.ParseFailureException:
					p = 0
				denoms[(l,nt,r)] = p
			if n == 1:
				print("Found characterising context",l,r)
				contexts = set()
				contexts.add((l,r))
				cc = (l,r)

				break
			
		#print("COntexts now of size",len(contexts))
		## Now we have a very simple way to compute the lbr probabilities
		def _compute(l,b,r):
			p = 0
			for nt in self.target_pcfg.nonterminals:
				if (nt,b) in self.target_pcfg.parameters:
					p += denoms[(l,nt,r)] * self.target_pcfg.parameters[(nt,b)]
			return p

		## so now e have a set of contexts perhap of cardianlity 1.
		xis ={}
		for b in self.terminals:
			if (nta,b) in self.target_wcfg.parameters:
				target = self.target_wcfg.parameters[ (nta,b)]
			else:
				target = 0.0
			best = math.inf
			bestcontext = None
			for l,r in contexts:
				denom = denoms[(l,nta,r)] * self.target_pcfg.parameters[(nta,a)]
				num = _compute(l,b,r)
				ratio = (num/denom)  * self.te[a]
				if ratio < best:
					best = ratio
					bestcontext = (l,r)
				if best == 0:
					break
				if target > 0 and abs(target - best) < EPSILON:
					#print("Exiting early")
					break
			if best > 0:
				xis[(a,b)] = best
			self.prod2context[(a,b)] = bestcontext
		return (xis, cc)
		
	def test_both(self, a, n, maxsamples):
		nta = self.kernel_map[a]
		contexts = self.sample_contexts_smart(a,n,maxsamples)
		ll = 0.0
		for l,r in contexts:
			ll += len(l) + len(r)
		print("Smart Contexts:",len(contexts), ll/len(contexts))
		denoms = {}
		for l,r in contexts:
			denoms[(l,r)] = self.insider.inside_log_probability(l + (a,) + r)
		xis = {}
		print("Unary")
		for b in self.terminals:
			if (nta,b) in self.target_wcfg.parameters:
				target = self.target_wcfg.parameters[ (nta,b)]
			else:
				target = 0.0
			best = math.inf
			try:
				for l,r in contexts:
					denom = denoms[(l,r)]
					num = self.insider.inside_log_probability(l + (b,) + r)
					ratio = math.exp(num - denom) * self.te[a]
					best = min(ratio,best)
					if target > 0 and abs(target - best) < EPSILON:
						#print("Exiting early")
						break
				xis[(a,b)] = best
			except utility.ParseFailureException:
				pass
		print("Binary")

		for b in self.kernel:
			for c in self.kernel:
				ntb = self.kernel_map[b]
				ntc = self.kernel_map[c]
				if (nta,ntb,ntc) in self.target_wcfg.parameters:
					target = self.target_wcfg.parameters[ (nta,ntb,ntc)]
				else:
					target = 0.0
				best = math.inf
				try:
					for l,r in contexts:
						denom = denoms[(l,r)]
						num = self.insider.inside_log_probability(l + (b,c) + r)
						ratio = math.exp(num - denom) * self.te[a]/(self.te[b] * self.te[c])
						best = min(ratio,best)
				
					xis[(a,b,c)] = best
				except utility.ParseFailureException:
					# if we catch an exception then we have a zero so we omit it.
					pass
		return xis


	def test_unary(self, a, n, maxsamples):
		contexts = self.sample_contexts1(a,n,maxsamples)
		scale = self.te[a]
		denoms = {}
		for l,r in contexts:
			denoms[(l,r)] = self.insider.inside_log_probability(l + (a,) + r)
		xis = {}
		for b in self.terminals:
			best = math.inf
			try:
				for l,r in contexts:
					denom = denoms[(l,r)]
					num = self.insider.inside_log_probability(l + (b,) + r)
					ratio = math.exp(num - denom) * scale
					best = min(ratio,best)
			
				xis[b] = best
			except utility.ParseFailureException:
				pass
		return xis

	def test_binary(self, a, n, maxsamples):
		contexts = self.sample_contexts1(a,n,maxsamples)
		
		denoms = {}
		for l,r in contexts:
			denoms[(l,r)] = self.insider.inside_log_probability(l + (a,) + r)
		xis = {}
		for b in self.kernel:
			for c in self.kernel:
				best = math.inf
				try:
					for l,r in contexts:
						denom = denoms[(l,r)]
						num = self.insider.inside_log_probability(l + (b,c) + r)
						ratio = math.exp(num - denom) * self.te[a]/(self.te[b] * self.te[c])
						best = min(ratio,best)
				
					xis[(b,c)] = best
				except utility.ParseFailureException:
					pass
		return xis

	def make_grammar(self):
		"""
		Create the grammar
		"""
		self.asymptotic_grammar = wcfg.WCFG()
		g = self.asymptotic_grammar
		g.nonterminals = set(self.target_pcfg.nonterminals)
		g.start = self.target_pcfg.start
		g.terminals = set(self.target_pcfg.terminals)
		for prod in self.parameters:
			if len(prod) == 2:
				a,b = prod
				newprod = (self.kernel_map[a],b)
			else:
				a,b,c = prod
				newprod = (self.kernel_map[a],self.kernel_map[b],self.kernel_map[c])
			g.productions.append(newprod)
			g.parameters[newprod] = self.parameters[prod]
		g.set_log_parameters()
		return g
