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

## 
## stacked bar chart  of violations of local ambiguity.
bins = 4

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
			if data["anchor errors"] == 0 or True:
				el = data["local ambiguity errors lexical"]
				eb = data["local ambiguity errors binary"]
				if el == 0 and eb == 0:
					idx = 0
				if el > 0 and eb == 0:
					idx = 1
				if el == 0 and eb > 0:
					idx = 2
				if el > 0 and eb > 0:
					idx = 3
			else:
				idx = 4
			values[i,idx] += 1
#print(values)
p1 = plt.barh(xs, values[:,0], height = 5, label="Correct")
p1 = plt.barh(xs, values[:,1], height = 5, left = values[:,0],label="L")
p1 = plt.barh(xs, values[:,2], height = 5, left = values[:,0] + values[:,1],label="B")
p1 = plt.barh(xs, values[:,3], height = 5, left = values[:,0] + values[:,1] + values[:,2],label="L+B")
#p1 = plt.barh(xs, values[:,4], height = 5, left = values[:,0] + values[:,1] + values[:,2] + values[:,3],label="AF")


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
plt.title("Incorrect parameters", y= 1.15)
plt.savefig("../diagrams/hbar_local_amb.pdf")