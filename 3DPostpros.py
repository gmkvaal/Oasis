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
os.chdir("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/2016.4.7.18")
l = glob.glob('*.txt')

# Sort the files alpabetically and numerically.
convert = lambda text: int(text) if text.isdigit() else text
alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
list_of_files = sorted(l, key = alphanum_key)

# Open and read the files
print ""
for resultfile in list_of_files:
	f = open(resultfile, 'r')
	for line in f:
		output = eval(line)
		results.append(output)

E_matrix  = results[0]
h_matrix  = results[1]
Cd_matrix = results[2]
Cl_matrix  = results[3]
La_matrix = results[4]
dP_matrix = results[5]
ct_matrix = results[6]
x = results[7]
y = results[8]



Cdref = 5.5800
Clref = 0.0106
Lref  = 0.0847
dPref = 0.1174
Cd_E = list((np.array(Cd_matrix)-Cdref).T)
Cl_E = list((np.array(Cl_matrix)-Clref).T)
L_E  = list(np.array(La_matrix)-Lref)
dP_E = list(np.array(dP_matrix)-dPref)

# Making table
table = []; headers = [0]
[headers.append(y[i]) for i in range(len(y))]
# Setting n/a to values not calculated
Cl_E_tab = Cl_E; Cd_E_tab = Cd_E
for i in range(len(Cl_E)):
	Cl_E_tab[i] = ["n/a" if ind==0.9894 else ind for ind in Cl_E[i]] 
	Cl_E_tab[i].insert(0,x[i])
	#print type(Cd_E_tab)
	#Cd_E_tab[i].insert(0,x[i])	

print tabulate(Cl_E_tab,headers, tablefmt="grid")

plot_it = True
if plot_it == True:
	# Plotting
	Z = Cl_E
	X,Y= np.meshgrid(np.log(np.array(x)),np.log(np.array(y)))
	fig = plt.figure("Cl_E")
	ax = fig.gca(projection='3d')
	norm = plt.matplotlib.colors.Normalize(vmin = np.min(0), vmax = np.max(0.01), clip = False)
	ax.scatter(X,Y,Z,'ro')
	surf = ax.plot_surface(X,Y,np.array(Z),cstride=1,rstride=1, alpha=0.6, cmap=cm.jet, linewidth=0)
	fig.colorbar(surf, shrink=0.5, aspect=5)
	#ax.set_zlim(0, 0.008)
	ax.set_xlabel('Edge Resolution')
	ax.set_ylabel('Circle Resolution')
	ax.set_zlabel('Cl Error')

	plt.show()

