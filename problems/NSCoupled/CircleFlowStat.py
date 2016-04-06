import os
import datetime
from ..NSCoupled import *
import numpy as np
import matplotlib.pyplot as plt
parameters['allow_extrapolation'] = True
set_log_active(False)
from dolfin import *
import time

time_start = time.clock()

NS_parameters.update(
	d = 0.1,
	U = 0.2,
	nu = 0.001,
	plotit = True,
    omega = 1.0,
    max_iter = 1000,
    plot_interval = 10,
    #velocity_degree = 2,
    CFLwrite = False,
    key=1,
    circle = 0.2/1,
	edge = 0.2/1,
	name = "hello",
	makemesh = False,
	resultswrite = False,
	foldername="TempRes",
	)


def mesh(makemesh, name, circle, edge, key, **params):
	if makemesh == True:
		import subprocess
		os.chdir("/home/guttorm/Desktop/Master/Mesh/Circle/Coarse/AutoMesh")
		os.system("python ControlMakeMesh.py %s %f %f" % (name, circle, edge))
		return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Coarse/AutoMesh/CFM%s.xml" % name)
	else:
		if key == 1: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Coarse/1to16ratio/CM715E.xml")
		if key == 2: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Coarse/1to16ratio/CM2256E.xml")
		if key == 3: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Coarse/1to16ratio/CM8041E.xml")
		if key == 4: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Coarse/1to16ratio/CM30461E.xml")
		if key == 5: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Coarse/1to16ratio/CM121589E.xml")


up = "std::abs(x[1]-0.41) < 1e-8"
down = "std::abs(x[1]) < 1e-8"
left = "std::abs(x[0]) < 1e-8"
right = "std::abs(x[0]-2.2) < 1e-8"

class Up(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and near(x[1],0.41)

class Down(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and near(x[1],0)

class Left(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and near(x[0],0)

class Right(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and near(x[0], 2.2)

class Circle(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and not near(x[0],0) \
						   and not near(x[0],2.2) \
						   and not near(x[1],0) \
						   and not near(x[1],0.41)
								
circle = Circle()
upper = Up()
lower= Down()
left = Left()
right = Right()

u0 = Expression(('1.2*x[1]*(0.41-x[1])/(0.41*0.41)','0'))

def create_bcs(VQ,** NS_namespace):

	bc0 = DirichletBC(VQ.sub(0), u0, left)
	bc1 = DirichletBC(VQ.sub(0), (0, 0), circle)
	bc2 = DirichletBC(VQ.sub(0), (0, 0), upper)
	bc3 = DirichletBC(VQ.sub(0), (0, 0), lower)
	bc4 = DirichletBC(VQ.sub(1), 0.0, right)

	return dict(up = [bc0, bc1, bc2, bc3, bc4])



def theend_hook(mesh, q_, p_, u_,u_components, nu, VQ, V, VV, Q, U, d, \
				sys_comp, up_, key, plotit, CFLwrite, resultswrite, foldername, **NS_namespace):
	 
	comptime = (time.clock() - time_start)
	pressure = p_
	boundary = FacetFunction("size_t", mesh)
	boundary.set_all(0)
	circle.mark(boundary, 1)
	ds = Measure("ds", subdomain_data=boundary)
	n = FacetNormal(mesh)
	R = VectorFunctionSpace(mesh, 'R', 0)
	c = TestFunction(R)
	tau = -pressure*Identity(2)+nu*(grad(u_)+grad(u_).T)
	forces = -assemble(dot(dot(tau, n), c)*ds(1)).array()*2/U**2/d

	x_c = np.linspace(0,2.0,1000)
	u_centercut= zeros(len(x_c))
	for i in range(len(x_c)):
		u_centercut[i] = abs(u_[0](array([x_c[i],0.20])))

	# Following is to locate max and min values
	max_list=[]; min_list=[]
	for j in range(0,(len(u_centercut)-2)):
		A = u_centercut[j+2]
		B = u_centercut[j+1]
		C = u_centercut[j]
		if A < B:
			if B > C:
				max_list.append(list(u_centercut).index(B))
		if A > B:
			if B < C:
				min_list.append(list(u_centercut).index(B))
	
	if plotit == True:
		uu = project(u_,V)
		now = datetime.datetime.now()
		identity = mesh.num_cells()
		atm = "%d.%d.%d.%d" % (now.year, now.month, now.day, now.hour)
		if os.path.isdir("/home/guttorm/Desktop/Master/Oasis/results/data/ResultsRe20%s" % atm):
			f = File("/home/guttorm/Desktop/Master/Oasis/results/data/ResultsRe20%s/u_from_CircleStat_E%d.pvd" % (atm,identity))
		else:
			os.system("mkdir -p /home/guttorm/Desktop/Master/Oasis/results/data/ResultsRe20%s" % atm)
			f = File("/home/guttorm/Desktop/Master/Oasis/results/data/ResultsRe20%s/u_from_CircleStat_E%d.pvd" % (atm,identity))
		f << uu

	if resultswrite == True:
		now = datetime.datetime.now()
		#atm = "%d.%d.%d.%d" % (now.year, now.month, now.day, now.hour)
		atm = foldername
		os.system("mkdir -p /home/guttorm/Desktop/Master/RefinementData/Re20/%s" % atm)
		identity = mesh.num_cells()
		text_file = open("/home/guttorm/Desktop/Master/RefinementData/Re20/%s/OutputCircleStationary%d.txt" % (atm,identity), "w")
		text_file.write("-----Result for %dE-----\n" % identity)
		text_file.write("h min = %.11f \n" % mesh.hmin())
		text_file.write("Cd max: %.11f \n" % forces[0])
		text_file.write("Cl max: %.11f \n" % forces[1])
		text_file.write("Resirculation length: %.11f \n" % (x_c[min_list[-1]] - x_c[min_list[-2]]))
		text_file.write("Delta P (front-back): %.11f \n" % (p_(array([0.15,0.20]))- p_(array([0.25,0.20]))))
		text_file.write("Computational time:: %.5f" % comptime)
		text_file.close()


	print ""
	print "----- Results -----"
	print "Cd = {}, CL = {}".format(*forces)
	print "Resirculation length:%.11f" % (x_c[min_list[-1]] - x_c[min_list[-2]])
	print "Delta P (front-back): %.11f" % (p_(array([0.15,0.20]))- p_(array([0.25,0.20])))
	print "Computational time:: %.5f" % comptime
	print ""