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


# Ref: Circle 0.0001 Edge 0.01
Cdref = 5.57933281841
Clref = 0.010618455861
Lref  = 0.08408408408
dPref = 0.11752008300


results = []
now = datetime.datetime.now()
atm = "%d.%d.%d.%d" % (now.year, now.month, now.day, now.hour)
#os.chdir("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s" % atm)
os.chdir("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/2016.4.13.13")
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


Cd_E = list((np.array(Cd_matrix)-Cdref).T)
Cl_E = list((np.array(Cl_matrix)-Clref).T)
L_E  = list((np.array(La_matrix)-Lref).T)
dP_E = list((np.array(dP_matrix)-dPref).T)


#Cl_tab = list((np.array(Cl_matrix)).T)
#Cd_tab = list((np.array(Cd_matrix)).T)
#La_tab = list((np.array(La_matrix)).T)
#dP_tab = list((np.array(dP_matrix)).T)


Cl_tab = [list(i) for i in zip(*Cl_matrix)]
Cd_tab = [list(i) for i in zip(*Cd_matrix)]
La_tab = [list(i) for i in zip(*La_matrix)]
dP_tab = [list(i) for i in zip(*dP_matrix)]

E_matrix = [list(i) for i in zip(*E_matrix)]
h_matrix = [list(i) for i in zip(*h_matrix)]
ct_matrix = [list(i) for i in zip(*ct_matrix)]


# Making table
table = []; headers = [0]
[headers.append(y[i]) for i in range(len(y))]
# Setting n/a to values not calculated
Cl_E_tab = Cl_E; Cd_E_tab = Cd_E; L_E_tab = L_E; dP_E_tab = dP_E


for i in range(len(Cl_E)):
	Cl_E_tab[i] = [ '%.4e' % elem for elem in Cl_E_tab[i] ]
	Cl_E_tab[i] = ["n/a" if ind==9.8938e-01 else ind for ind in Cl_E[i]] 
	Cl_E_tab[i].insert(0,x[i])
	Cd_E_tab[i] = [ '%.4e' % elem for elem in Cd_E_tab[i] ]
	Cd_E_tab[i] = ["n/a" if ind==-4.57933281841 else ind for ind in Cd_E[i]] 
	Cd_E_tab[i].insert(0,x[i])
	L_E_tab[i] = [ '%.4e' % elem for elem in L_E_tab[i] ]
	L_E_tab[i] = ["n/a" if ind==1-0.08408408408 else ind for ind in L_E[i]] 
	L_E_tab[i].insert(0,x[i])
	dP_E_tab[i] = [ '%.4e' % elem for elem in dP_E_tab[i] ]
	dP_E_tab[i] = ["n/a" if ind==1-0.11752008300 else ind for ind in dP_E[i]] 
	dP_E_tab[i].insert(0,x[i])
	E_matrix[i].insert(0,x[i])
	ct_matrix[i].insert(0,x[i])
	Cl_tab[i].insert(0,x[i])
	Cd_tab[i].insert(0,x[i])
	La_tab[i].insert(0,x[i])
	dP_tab[i].insert(0,x[i])

# tables = what to view
#Cl_E_tab = [ '%.2f' % elem for elem in Cl_E_tab ]
tables = ct_matrix
print tabulate(tables,headers, tablefmt="latex", floatfmt=".3e")


plot_it = False
if plot_it == True:
	# Plotting
	Z = Cl_matrix
	X,Y= np.meshgrid(np.log(np.array(x)),np.log(np.array(y)))
	fig = plt.figure("Cl_E")
	ax = fig.gca(projection='3d')
	norm = plt.matplotlib.colors.Normalize(vmin = np.min(0), vmax = np.max(0.01), clip = False)
	ax.scatter(X,Y,np.array(Z),'ro')
	surf = ax.plot_surface(X,Y,np.array(Z),cstride=1,rstride=1, alpha=0.6, cmap=cm.jet, linewidth=0)
	fig.colorbar(surf, shrink=0.5, aspect=5)
	#ax.set_zlim(0, 0.008)
	ax.set_xlabel('Edge Resolution')
	ax.set_ylabel('Circle Resolution')
	ax.set_zlabel('Cl Error')

	plt.show()

