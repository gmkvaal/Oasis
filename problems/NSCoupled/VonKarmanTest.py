import os
import datetime
from ..NSCoupled import *
import numpy as np
import matplotlib.pyplot as plt
parameters['allow_extrapolation'] = True
set_log_active(False)
from dolfin import *

NS_parameters.update(
	d = 0.5,
	U = 2/3,
	nu = 10,
	rho = 1,
	plotit = False,
    omega = 1.0,
    max_iter = 1000,
    plot_interval = 10,
    velocity_degree = 2,
    CFLwrite = False,
    key=1,
    circle = 0.01,
	edge = 0.1,
	name = "AutogeneratedMesh",
	makemesh = False,
	resultswrite = False,
	foldername="TempRes",
	)

mesh = Mesh("/home/guttorm/Desktop/MEK4300/Week13/mesh/vankarmanmesh.xml")



class Walls(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and near(x[1],-1) \
						   or  near(x[1],1)

class Outlet(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and near(x[0],6) 

class Inlet(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and near(x[0],0) 

class Circle(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and \
		pow(x[0]-1,2) + pow(x[1],2) - 0.0625 <= DOLFIN_EPS


	    
print mesh

walls = Walls()
outlet = Outlet()
inlet = Inlet()
circle = Circle()

u_inn = Expression(('1-x[1]*x[1]','0'))
u_circle = Expression(('-4*x[1]','4*(x[0]-1)'))

def create_bcs(VQ,** NS_namespace):

	bc1 = DirichletBC(VQ.sub(0),u_inn,walls)
	bc2 = DirichletBC(VQ.sub(0),u_circle,circle)
	bc3 = DirichletBC(VQ.sub(0),u_inn,inlet)

	return dict(up = [bc1,bc2,bc3])



def theend_hook(mesh, q_, p_, u_,u_components, nu, VQ, V, VV, Q, U, d, \
				sys_comp, up_, key, plotit, CFLwrite, resultswrite, foldername, **NS_namespace):
	 
	pressure = p_
	boundary = FacetFunction("size_t", mesh)
	boundary.set_all(0)
	circle.mark(boundary, 1)
	ds = Measure("ds", subdomain_data=boundary)
	n = FacetNormal(mesh)
	R = VectorFunctionSpace(mesh, 'R', 0)
	c = TestFunction(R)
	tau = -pressure*Identity(2)+nu*(grad(u_)+grad(u_).T)
	forces = -assemble(dot(dot(tau, n), c)*ds(1)).array()
	print forces
	plot(u_,interactive=True)