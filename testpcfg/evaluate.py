

import argparse
import wcfg
import math
import os.path
import oracle_learner

parser = argparse.ArgumentParser(description='Evaluate the asympotic wcfg versus the original pcfg.')
parser.add_argument('input', type=str,  help='location of the target pcfg.')
parser.add_argument('output', type=str, 	help='location of the output wcfg.')


args = parser.parse_args()
verbose = False

target_pcfg = wcfg.load_wcfg_from_file(args.input)
asymptotic_wcfg = wcfg.load_wcfg_from_file(args.output)

pf = asymptotic_wcfg.compute_partition_function_fast()
print("Log partition function: %e" % math.log(pf[asymptotic_wcfg.start]))

delta_lexical = asymptotic_wcfg.count_lexical() - target_pcfg.count_lexical()
delta_binary = asymptotic_wcfg.count_binary() - target_pcfg.count_binary()
print("Extra lexical productions: %d" % delta_lexical)
print("Extra binary productions: %d" % delta_binary)

asymptotic_pcfg = asymptotic_wcfg.convert_parameters_pi2xi()
kld = 0.0
pe = target_pcfg.production_expectations()
for prod in pe:
	e = pe[prod]
	alpha = target_pcfg.parameters[prod]
	beta = asymptotic_pcfg.parameters[prod]
	delta =  e * math.log(alpha/beta)
	if verbose:
		print(prod,newprod,delta)
		kld +=delta
print("Labeled tree KLD %e" %kld)
