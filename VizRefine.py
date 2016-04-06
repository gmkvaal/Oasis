from Refinement import *
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#from matplotlib import cm


y = []; ZC = []; ZE = []

iterationlist = np.linspace(1,16,16)

#for j in iterationlist:
for j in [1,2,4]:
	initial_edge = 0.05/j
	initial_circle=0.05
	y.append(initial_edge)
	x,CL_E,Elements = runit(initial_edge=initial_edge,initial_circle=initial_circle,iterationlist=iterationlist)
	ZC.append(CL_E)
	ZE.append(Elements)


# Storing the result matrix
now = datetime.datetime.now()
identity = len(ZE[0])
atm = "%d.%d.%d.%d" % (now.year, now.month, now.day, now.hour)
if os.path.isdir("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s" % atm):
	text_file = open("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s/Cl_E_dimention%s.txt" % (atm,identity),"w")
else:
	os.system("mkdir -p /home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s" % atm)
	text_file = open("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s/Cl_E_dimention%s.txt" % (atm,identity),"w")
text_file.write(str(ZC)+"\n")
text_file.write(str(x)+"\n")
text_file.write(str(y))

X,Y= np.meshgrid(np.log(np.array(x)),np.log(np.array(y)))
fig = plt.figure("Cl")
ax = fig.gca(projection='3d')
norm = plt.matplotlib.colors.Normalize(vmin = np.min(0), vmax = np.max(0.01), clip = False)
ax.scatter(X,Y,np.log(ZC),'ro')
surf = ax.plot_surface(X,Y,np.log(ZC),cstride=1,rstride=1, alpha=0.6, cmap=cm.jet, linewidth=0)
fig.colorbar(surf, shrink=0.5, aspect=5)
#ax.set_zlim(0, 0.008)
ax.set_xlabel('Circle Resolution')
ax.set_ylabel('Edge Resolution')
ax.set_zlabel('Cl Error')

#plt.hold(True)

fig = plt.figure("Elements")
ax = fig.gca(projection='3d')
norm = plt.matplotlib.colors.Normalize(vmin = np.min(0), vmax = np.max(0.01), clip = False)
ax.scatter(X,Y,ZE,'ro')
surf = ax.plot_surface(X,Y,ZE,cstride=1,rstride=1, alpha=0.6, cmap=cm.jet, linewidth=0)
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel('Circle Resolution')
ax.set_ylabel('Edge Resolution')
ax.set_zlabel('Elements')
plt.show()

