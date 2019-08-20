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

## Plot graph of the anchor failures


values = np.zeros( N)
totals = np.zeros(N)
for i,c in enumerate(xs):
	##
	#print(c)
	cdir = rootdir + "test%d/grammar*.json" % c
	print(cdir)
	for f in glob.glob(cdir):
		#print(f)
		with open(f) as json_file:
			data = json.load(json_file)
			if data["anchor errors"] > 0:
				values[i] += 1
			totals[i] += 1


#print(values)
plt.plot(xs, values/totals,'.-')
# p2 = plt.bar(ind, womenMeans, width,
#              bottom=menMeans, yerr=womenStd)

plt.xlabel('Binary rules')
plt.ylim(0,1)
# plt.title('Scores by group and gender')
# plt.xticks(ind, ('G1', 'G2', 'G3', 'G4', 'G5'))
# plt.yticks(np.arange(0, 81, 10))
# plt.legend((p1[0], p2[0]), ('Men', 'Women'))
plt.ylabel("Too few nonterminals")
plt.xticks(xs)
#plt.title("Correct     Convergent     Divergent")
plt.title("Errors with anchors", y= 1.15)
plt.savefig("../diagrams/anchors.pdf")