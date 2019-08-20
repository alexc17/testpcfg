#plot_diagrams.py


# all hard codes just a script really.


import glob
import matplotlib.pyplot as plt
import math
import json
import numpy as np
rootdir = "../data/"

xs = [20,30,40,50,60,70,80]
N = len(xs)

## Box plot of density


values = []
for i,c in enumerate(xs):
	##
	#print(c)
	cdir = rootdir + "test%d/grammar*.json" % c
	print(cdir)
	ys = []
	for f in glob.glob(cdir):
		#print(f)
		with open(f) as json_file:
			data = json.load(json_file)
			# estimate may be > 1 in which case we could cap it
			# ys.append(min(1.0,data["string density"]))
			# better to leave the excess values in to give idea of the variability in the Monte Carlo estimates.
			ys.append(data["string density"])

	values.append(ys)

print(len(values))
#print(values)
plt.boxplot(values,labels = xs)
# p2 = plt.bar(ind, womenMeans, width,
#              bottom=menMeans, yerr=womenStd)

plt.ylabel('String density at length 10')
plt.xlabel('Binary rules')
plt.savefig("../diagrams/boxplot_density.pdf")
plt.yscale('log')
# plt.title('Scores by group and gender')
# plt.xticks(ind, ('G1', 'G2', 'G3', 'G4', 'G5'))
# plt.yticks(np.arange(0, 81, 10))
# plt.legend((p1[0], p2[0]), ('Men', 'Women'))
# plt.ylabel("Too few nonterminals")
#plt.xticks(xs)
# #plt.title("Correct     Convergent     Divergent")
# plt.title("Errors with anchors", y= 1.15)
plt.savefig("../diagrams/boxplot_log_density.pdf")