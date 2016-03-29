from Refinement import *
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#from matplotlib import cm


y = []; ZC = []; ZE = []
#iterationlist = [1,1.5,2,2.5,3,3.5,4,4.5,5]
#iterationlist = [1,1.5,2,2.5,3]
iterationlist = [1,2,4]
for j in iterationlist:
	initial_edge = 0.05/j
	y.append(initial_edge)
	x,CL_E,Elements = runit(initial_edge=initial_edge,initial_circle=0.005,iterationlist=iterationlist)
	ZC.append(CL_E)
	ZE.append(Elements)

	
X,Y= np.meshgrid(np.log(np.array(x)),np.log(np.array(y)))
fig = plt.figure("Cl")
ax = fig.gca(projection='3d')
norm = plt.matplotlib.colors.Normalize(vmin = np.min(0), vmax = np.max(0.01), clip = False)
ax.scatter(X,Y,ZC,'ro')
surf = ax.plot_surface(X,Y,ZC,cstride=1,rstride=1, alpha=0.6, cmap=cm.jet, linewidth=0)
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_zlim(0, 0.008)
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
ax.set_zlabel('Cl Error')
plt.show()

