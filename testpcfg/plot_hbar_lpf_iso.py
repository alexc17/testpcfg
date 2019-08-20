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

## First plot
## stacked bar chart  of partition functions per condition: only with isomorphic grammars
bins = 5

values = np.zeros( (N,bins))

for i,c in enumerate(xs):
	##
	#print(c)
	cdir = rootdir + "test%d/grammar*.json" % c
	print(cdir)
	for f in glob.glob(cdir):
		#print(f)
		with open(f) as json_file:
			data = json.load(json_file)
			if data["anchor errors"] == 0 and data["extra binary"] == 0 and data["extra lexical"] == 0:
				pf = data.get("partition function",math.inf)
				#print(pf)
				if pf == math.inf:
					values[i,4] += 1
				else:
					spf = math.log(pf['S'])
					#print(spf)
					if spf < 0.001:
						values[i,0] += 1
					elif spf < 0.01:
						print(f,spf)
						values[i,1] += 1
					elif spf < 0.1:
						values[i,2] += 1
					else:
						values[i,3] += 1


#print(values)
p1 = plt.barh(xs, values[:,0], height = 5, label="< 0.001")
p1 = plt.barh(xs, values[:,1], height = 5, left = values[:,0],label=" (<0.01)")
p1 = plt.barh(xs, values[:,2], height = 5, left = values[:,0] + values[:,1],label=" (<0.1)")
p1 = plt.barh(xs, values[:,3], height = 5, left = values[:,0] + values[:,1] + values[:,2],label="$\geq 0.1$")
p1 = plt.barh(xs, values[:,4], height = 5, left = values[:,0] + values[:,1] + values[:,2] + values[:,3],label="$\infty$")

# p2 = plt.bar(ind, womenMeans, width,
#              bottom=menMeans, yerr=womenStd)

# plt.ylabel('Scores')
# plt.title('Scores by group and gender')
# plt.xticks(ind, ('G1', 'G2', 'G3', 'G4', 'G5'))
# plt.yticks(np.arange(0, 81, 10))
# plt.legend((p1[0], p2[0]), ('Men', 'Women'))
plt.ylabel("Binary rules")
plt.yticks(xs)
#plt.title("Correct     Convergent     Divergent")
plt.legend(loc='upper center',ncol=4,bbox_to_anchor=(0.5,1.15))
plt.title("Log Partition Function", y= 1.15)
plt.savefig("../diagrams/hbar_lpf_iso.pdf")