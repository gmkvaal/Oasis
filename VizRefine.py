from Refinement import *
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#from matplotlib import cm



E_matrix=[]; h_matrix=[]; Cd_matrix=[]; y=[]
Cl_matrix = []; La_matrix=[]; dP_matrix=[]; ct_matrix=[]
iterationlist = np.linspace(1,24,24)
circle_list = 0.07/iterationlist
x = list(circle_list)
for j in [1,2,4,6,8]:
	edgeres = 0.07/j
	y.append(edgeres)
	E, h, Cd, Cl, La, dP, comptime = runit(edgeres,circle_list)
	E_matrix.append(E)
	h_matrix.append(h)
	Cd_matrix.append(Cd)
	Cl_matrix.append(Cl)
	La_matrix.append(La)
	dP_matrix.append(dP)
	ct_matrix.append(comptime)


# Storing the result matrix
now = datetime.datetime.now()
identity = len(Cl)
atm = "%d.%d.%d.%d" % (now.year, now.month, now.day, now.hour)
if os.path.isdir("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s" % atm):
	print "----Path exists, danger of overwriting!----"
else:
	os.system("mkdir -p /home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s" % atm)
text_file1 = open("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s/1E1_%s.txt" % (atm,identity),"w")
text_file2 = open("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s/2h_%s.txt" % (atm,identity),"w")
text_file3 = open("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s/3Cd_%s.txt" % (atm,identity),"w")
text_file4 = open("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s/4Cl_%s.txt" % (atm,identity),"w")
text_file5 = open("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s/5La_%s.txt" % (atm,identity),"w")
text_file6 = open("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s/6dP_%s.txt" % (atm,identity),"w")
text_file7 = open("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s/7ct_%s.txt" % (atm,identity),"w")
text_file8 = open("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s/8X_%s.txt" % (atm,identity),"w")
text_file9 = open("/home/guttorm/Desktop/Master/RefinementData/Re20/2Dref/%s/9Y_%s.txt" % (atm,identity),"w")
text_file1.write(str(E_matrix))
text_file2.write(str(h_matrix))
text_file3.write(str(Cd_matrix))
text_file4.write(str(Cl_matrix))
text_file5.write(str(La_matrix))
text_file6.write(str(dP_matrix))
text_file7.write(str(ct_matrix))
text_file8.write(str(x))
text_file9.write(str(y))


"""
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
"""
