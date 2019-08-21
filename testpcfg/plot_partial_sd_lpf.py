#plot_diagrams.py


# all hard codes just a script really.


import glob
import matplotlib.pyplot as plt
import math
import json
import numpy as np
rootdir = "../data/"

xs = [20,40,80]
N = len(xs)

## Box plot of density

colors = ['r','g','b']

for i,c in enumerate(xs):
	##
	x = []
	y = []
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
			if data["anchor errors"] == 0:
				

				pf = data.get("partition function", math.inf)
				if pf != math.inf:
					x.append(data["string density"])
					y.append(math.log(pf['S']))
	plt.plot(x,y,'.' + colors[i],alpha=0.5,label=c)


plt.xlabel('String density at length 10')
plt.ylabel('Partition function')

plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig("../diagrams/partial_sd_lpf.pdf")