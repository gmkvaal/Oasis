
from ..NSCoupled import *
import numpy as np
import matplotlib.pyplot as plt
parameters['allow_extrapolation'] = True
set_log_active(False)
#mesh = Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM7169E.xml")
#mesh = Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM24927E.xml")
#mesh = Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM94113E.xml")
#mesh = Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM372949E.xml")


NS_parameters.update(
	d = 0.1,
	U = 0.2,
	nu = 0.001,
    omega = 1.0,
    max_iter = 1000,
    plot_interval = 10,
    velocity_degree = 2,
    CFLwrite = False,
    key=1)


def mesh(key,**params):	
	if key == 1: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM14969E.xml")
	if key == 2: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM39883E.xml")
	if key == 3: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM119197E.xml")
	if key == 4: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM372949E.xml")
	if key == 0: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/Old/CFM7169E.xml")




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



def theend_hook(mesh, q_, p_, u_, nu, VQ, V, Q, U, d, sys_comp,  key, CFLwrite, **NS_namespace):

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

	def writefile(**NS_namespace):
		text_file = open("/home/guttorm/Desktop/Master/RefinementData/Re20/OutputCircleStationary%d.txt" %key, "w")
		text_file.write("Cd max: %.11f \n" % forces[0])
		text_file.write("Cl max: %.11f \n" % forces[1])
		text_file.write("Resirculation length: %.11f \n" % (x_c[min_list[-1]] - x_c[min_list[-2]]))
		text_file.write("Delta P (front-back): %.11f" % (p_(array([0.15,0.20]))- p_(array([0.25,0.20]))))
		text_file.close()
	writefile(**NS_namespace)

	print ""
	print "----- Results -----"
	print "Cd = {}, CL = {}".format(*forces)
	print "Resirculation length:%.11f" % (x_c[min_list[-1]] - x_c[min_list[-2]])
	print "Delta P (front-back): %.11f" % (p_(array([0.15,0.20]))- p_(array([0.25,0.20])))
	print ""