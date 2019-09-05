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

infinities = [0] * N
values = []

for i,c in enumerate(xs):
	##
	#print(c)
	cdir = rootdir + "test%d/grammar*.json" % c
	print(cdir)
	y = []
	for f in glob.glob(cdir):
		#print(f)
		with open(f) as json_file:
			data = json.load(json_file)
			if data["anchor errors"] > 0 or not "kld" in data:
				infinities[i] += 1
			else:
				pf = data["kld"]
				y.append(pf)
	values.append(y)


#print(values)

plt.boxplot(values,labels = xs)
# p2 = plt.bar(ind, womenMeans, width,
#              bottom=menMeans, yerr=womenStd)
infy = plt.gca().get_ylim()[1] + 0.1
for i in range(N):
    # y = values[i]
    # x = np.random.normal(1+i, 0.04, size=len(y))
    # plt.plot(x, y, 'r.', alpha=0.2)

    infs = infinities[i]
    print(i,infs)
    x = np.random.normal(1+i, 0.075, size=infs)
    y = np.random.normal(infy, 0.015, size=infs)
    plt.plot(x,y ,'r^',alpha=0.2 )

plt.ylabel('KL Divergence (nats)')
# plt.title('Scores by group and gender')
# plt.xticks(ind, ('G1', 'G2', 'G3', 'G4', 'G5'))
# plt.yticks(np.arange(0, 81, 10))
# plt.legend((p1[0], p2[0]), ('Men', 'Women'))
# plt.xlim(0,100)
plt.xlabel("Binary rules")
# plt.yticks(xs)
# #plt.title("Correct     Convergent     Divergent")
# plt.legend(loc='upper center',ncol=4,bbox_to_anchor=(0.5,1.15))
#plt.title("Labeled Tree KLD")
plt.savefig("../diagrams/boxplot_kld_iso.pdf")