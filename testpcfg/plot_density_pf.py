#plot_diagrams.py


# all hard codes just a script really.


import glob
import matplotlib.pyplot as plt
import math
import json
import numpy as np
rootdir = "../data/"

xs = [20,30,40,50,60,70,80]

## Third plot
## X axis is ambiguity, Y axis pF. 

x = []
y = []
for c in xs:
	##
	#print(c)
	cdir = rootdir + "test%d/grammar*.json" % c
	print(cdir)
	for f in glob.glob(cdir):
		#print(f)
		with open(f) as json_file:
			data = json.load(json_file)

		pf = data.get("partition function",math.inf)
		#print(pf)
		if data["anchor errors"] == 0:
			if pf == math.inf:
				print("Divergent but anchors are correct")
				continue
			spf = math.log(pf['S'])
			amb = data["string density"]
			#if spf > 0.0001:
			x.append(amb)
			y.append(spf)

#print(values)

plt.plot(x,y,'.')

# plt.ylabel('Scores')
# plt.title('Scores by group and gender')
# plt.xticks(ind, ('G1', 'G2', 'G3', 'G4', 'G5'))
# plt.yticks(np.arange(0, 81, 10))
# plt.legend((p1[0], p2[0]), ('Men', 'Women'))
plt.ylabel("Log partition function")
plt.xlabel("String density")
#plt.title("Correct     Convergent     Divergent")
plt.savefig("../diagrams/density_pf.pdf")