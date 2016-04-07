import os
import glob
import re
import numpy
import datetime
import numpy as np
from tabulate import tabulate
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
ax = plt.gca()

results = []
now = datetime.datetime.now()
atm = "%d.%d.%d.%d" % (now.year, now.month, now.day, now.hour)
#os.chdir("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s" % atm)
os.chdir("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/2016.4.7.11")
l = glob.glob('*.txt')

# Sort the files alpabetically and numerically.
convert = lambda text: int(text) if text.isdigit() else text
alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
list_of_files = sorted(l, key = alphanum_key)

# Open and write the files
for resultfile in list_of_files:
	f = open(resultfile, 'r')
	for line in f:
		output = eval(line)
		results.append(output)

ZC = results[0]
x = results[1]
y = results[2]

table = []
headers = [0]
for i in range(len(x)):
	headers.append(x[i])
for i in range(len(ZC)):
	ZC[i].insert(0,y[i])
	#print ZC[i]
	#print len(ZC[i])

print tabulate(ZC,headers, tablefmt="grid")


plot_it = False
if plot_it == True:
	# Plotting
	X,Y= np.meshgrid(np.log(np.array(x)),np.log(np.array(y)))
	fig = plt.figure("Cl")
	ax = fig.gca(projection='3d')
	norm = plt.matplotlib.colors.Normalize(vmin = np.min(0), vmax = np.max(0.01), clip = False)
	ax.scatter(X,Y,(ZC),'ro')
	surf = ax.plot_surface(X,Y,np.array(ZC),cstride=1,rstride=1, alpha=0.6, cmap=cm.jet, linewidth=0)
	fig.colorbar(surf, shrink=0.5, aspect=5)
	#ax.set_zlim(0, 0.008)
	ax.set_xlabel('Circle Resolution')
	ax.set_ylabel('Edge Resolution')
	ax.set_zlabel('Cl Error')

	plt.show()

