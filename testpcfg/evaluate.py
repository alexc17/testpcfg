

import argparse
import wcfg
import math
import json
import os.path
import oracle_learner

parser = argparse.ArgumentParser(description='Evaluate the asympotic wcfg versus the original pcfg.')
parser.add_argument('input', type=str,  help='location of the target pcfg.')
parser.add_argument('output', type=str, 	help='location of the output wcfg.')
parser.add_argument('--json', type=str, 	help='location of the output json file if needed.')

args = parser.parse_args()
verbose = False

target_pcfg = wcfg.load_wcfg_from_file(args.input)
target_ambiguity = target_pcfg.estimate_ambiguity()
print("Target grammar ambiguity H(T}W): %e" % target_ambiguity)
asymptotic_wcfg = wcfg.load_wcfg_from_file(args.output)

try:
	pf = asymptotic_wcfg.compute_partition_function_fast()
	print("Log partition function: %e" % math.log(pf[asymptotic_wcfg.start]))

	asymptotic_pcfg = asymptotic_wcfg.convert_parameters_xi2pi()
	kld = 0.0
	pe = target_pcfg.production_expectations()
	for prod in pe:
		e = pe[prod]
		alpha = target_pcfg.parameters[prod]
		beta = asymptotic_pcfg.parameters[prod]
		delta =  e * math.log(alpha/beta)
		kld +=delta
	print("Labeled tree KLD %e" %kld)
except ValueError as e:
	pf = math.inf
	kld = math.inf
	print("Divergent")
delta_lexical = asymptotic_wcfg.count_lexical() - target_pcfg.count_lexical()
delta_binary = asymptotic_wcfg.count_binary() - target_pcfg.count_binary()
print("Extra lexical productions: %d" % delta_lexical)
print("Extra binary productions: %d" % delta_binary)


if args.json:
	result_dict = { "partition function" : pf, "extra lexical" : delta_lexical, "ambiguity" : target_ambiguity, "extra binary" : delta_binary }
	with open(args.json,'w') as ff:
		json.dump(result_dict, ff, sort_keys=True, indent=4)
