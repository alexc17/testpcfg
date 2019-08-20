

import argparse
import wcfg
import math
import json
import os.path
import oracle_learner
import uniformsampler

parser = argparse.ArgumentParser(description='Evaluate the asympotic wcfg versus the original pcfg.')
parser.add_argument('input', type=str,  help='location of the target pcfg.')
parser.add_argument('output', type=str, 	help='location of the asympotic wcfg.')
parser.add_argument('--json', type=str, 	help='location of the output json file if needed.')
parser.add_argument("--seed",help="Choose random seed",type=int)
parser.add_argument('--length', type=int, default=10, 	help='length to measure the string density at.')
parser.add_argument('--samples', type=int, default=1000, 	help='samples to measure the string density.')
args = parser.parse_args()


if args.seed:
	random.seed(args.seed)
	numpy.random.seed(args.seed)

verbose = False
result_dict = {}
target_pcfg = wcfg.load_wcfg_from_file(args.input)
target_ambiguity = target_pcfg.estimate_ambiguity()
result_dict["ambiguity"] =  target_ambiguity
print("Target grammar ambiguity H( tree | word): %e" % target_ambiguity)
asymptotic_wcfg = wcfg.load_wcfg_from_file(args.output)


## Check to see if the kernel is in fact identifiable
ol = oracle_learner.OracleLearner(target_pcfg)
ntmap = ol.test_if_anchored()
oops = 0
for nt,a in ntmap.items():
    for ntb in target_pcfg.nonterminals:
        if ntb != nt:
            if (ntb,a) in asymptotic_wcfg.productions:
                print("oops")
                print(ntb,a)
                oops += 1

result_dict["anchor errors"] = oops

## Now try string denisyt using a sensible approach.
us = uniformsampler.UniformSampler(target_pcfg, args.length)
sd = us.string_density(args.length,args.samples)
print("String density: %e" % sd)
result_dict["string density"] = sd
try:
	pf = asymptotic_wcfg.compute_partition_function_fast()
	print("Log partition function: %e" % math.log(pf[asymptotic_wcfg.start]))
	result_dict["partition function"] = pf
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
	result_dict["kld"] = kld
except ValueError as e:
	pf = math.inf
	kld = math.inf
	print("Divergent")
delta_lexical = asymptotic_wcfg.count_lexical() - target_pcfg.count_lexical()
result_dict["extra lexical"] = delta_lexical
delta_binary = asymptotic_wcfg.count_binary() - target_pcfg.count_binary()
result_dict["extra binary"] = delta_binary
print("Extra lexical productions: %d" % delta_lexical)
print("Extra binary productions: %d" % delta_binary)


if args.json:
	with open(args.json,'w') as ff:
		json.dump(result_dict, ff, sort_keys=True, indent=4)
