

import argparse
import wcfg
import math
import os.path
import oracle_learner
import collections
parser = argparse.ArgumentParser(description='Evaluate the asympotic wcfg versus the original pcfg.')
parser.add_argument('input', type=str, nargs="+", help='location of the target pcfgs')



args = parser.parse_args()
verbose = False

extra_rules = collections.Counter()
histogram = [0, 0, 0, 0, 0]
for i in args.input:
	print(i)
	target_pcfg = wcfg.load_wcfg_from_file(i)
	o = i[:-4] + "wcfg"
	asymptotic_wcfg = wcfg.load_wcfg_from_file(o)
	asymptotic_wcfg.compute_all_scores()
	#pf = asymptotic_wcfg.compute_partition_function_fast()
	#print("Partition function: %e" % (asymptotic_wcfg.partition_function))
	print("Log partition function: %e" % math.log(asymptotic_wcfg.partition_function))

	delta_lexical = asymptotic_wcfg.count_lexical() - target_pcfg.count_lexical()
	delta_binary = asymptotic_wcfg.count_binary() - target_pcfg.count_binary()
	print("Extra lexical productions: %d" % delta_lexical)
	print("Extra binary productions: %d" % delta_binary)
	extra_rules[ (delta_lexical > 0, delta_binary > 0)] += 1

	asymptotic_pcfg = asymptotic_wcfg.convert_parameters_xi2pi()
	asymptotic_pcfg.compute_all_scores()
	# print(asymptotic_pcfg.check_local_normalisation())
	# print(asymptotic_pcfg.partition_function)
	kld = 0.0
	pe = target_pcfg.production_expectations()
	for prod in pe:
		e = pe[prod]
		alpha = target_pcfg.parameters[prod]
		beta = asymptotic_pcfg.parameters[prod]
		delta =  e * math.log(alpha/beta)
		kld += delta
	print("Labeled tree KLD %e" % kld)
	if kld < 0.001:
		histogram[0] += 1
	elif kld < 0.01:
		histogram[1] += 1
	elif kld < 0.1:
		histogram[2] += 1
	elif kld < 1:
		histogram[3] += 1
	else:
		histogram[4] += 1
print(extra_rules)
print(histogram)