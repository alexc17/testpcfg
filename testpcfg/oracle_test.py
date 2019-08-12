#oracle_test.py
## This file takes a pcfg and outputs a wcfg that approximates Algorithm A

import argparse
import wcfg
import os.path
import oracle_learner
import random
import numpy.random
import logging

parser = argparse.ArgumentParser(description='Compute asymptotic grammar using oracle.')
parser.add_argument("--nsamples", type=int, default=100, help="Number of distinct contexts for each nonterminal (default 100)")
parser.add_argument("--maxsamples", type=int, default=1000, help="Maximum samples for each nonterminal (default 1000)")

parser.add_argument("--seed",help="Choose random seed",type=int)
parser.add_argument("--verbose",help="Print useful information",action="store_true")

parser.add_argument('input', type=str,  help='location of the target pcfg.')
parser.add_argument('output', type=str, 	help='location of the output wcfg.')


args = parser.parse_args()

if args.seed:
	random.seed(args.seed)
	numpy.random.seed(args.seed)

oracle_learner.N_SAMPLES = args.nsamples
oracle_learner.MAX_SAMPLES = args.maxsamples

i = args.input
target_pcfg = wcfg.load_wcfg_from_file(i)
logging.info("Loaded")
output_wcfg = args.output
	
ol = oracle_learner.OracleLearner(target_pcfg)
og = ol.test()
if og:
	og.store(output_wcfg)
else:
	## create an empty file
	logging.warning("Error: Unanchored, empty wcfg")
	open(output_wcfg, 'a').close()
