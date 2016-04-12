import os
import matplotlib.pyplot as plt
from ..NSfracStep import *
from numpy import sqrt,array
from datetime import date

set_log_active(False)
parameters['allow_extrapolation'] = True

			
#from IPython import embed; embed()

NS_parameters.update(
	U = 1,
	d = 0.1,
	nu = 0.001,
	T = 8,
	dt = 0.001,
    plot_interval = 1000,
    save_step = 10,
    print_intermediate_info = 100,
    L_list = [],
    D_list = [],
    CFL_list = [],
    p_list=[],
    key = 1,
    CFLwrite = False,
    resultswrite =False,
      	
    )

def mesh(key,**params):	
	if key == 1: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM14969E.xml")
	if key == 2: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM39883E.xml")
	if key == 3: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM119197E.xml")
	if key == 4: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM372949E.xml")
	if key == 0: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/Old/CFM7169E.xml")
	if key == 10: return Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/try.xml")


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
		return on_boundary and near(x[0],2.2)

class Circle(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and not near(x[0],0) \
						   and not near(x[0],2.2) \
						   and not near(x[1],0) \
						   and not near(x[1],0.41)

class Front(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and abs(x[0]-0.15) < 1e-3 \
						   and abs(x[1]-0.20) < 1e-3	

class Back(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and abs(x[0]-0.25) < 1e-3 \
						   and abs(x[1]-0.20) < 1e-3
			
circle = Circle()
up = Up()
down = Down()
left = Left()
right = Right()

u_inlet = Expression('4*1.5*x[1]*(0.41-x[1])/(0.41*0.41)')

#u_inlet = Expression('4*1.5*x[1]*(0.41-x[1])/(0.41*0.41)*sin(pi*t/8)',t=0)

def create_bcs(V,Q,** NS_namespace):
	bc_l0 = DirichletBC(V,u_inlet,left)
	bc_l1 = DirichletBC(V,0.0,left)
	bc_r0 = DirichletBC(V,u_inlet,right)
	bc_r1 = DirichletBC(V,0.0,right)
	bc_u0 = DirichletBC(V,0.0,up)
	bc_u1 = DirichletBC(V,0.0,up)
	bc_d0 = DirichletBC(V,0.0,down)
	bc_d1 = DirichletBC(V,0.0,down)

	bc_c = DirichletBC(V,0.0,circle,method="topological",check_midpoint=False)
	bcq_r = DirichletBC(Q,0.0,right)

	return dict(u0=[bc_l0,bc_u0,bc_d0,bc_c],
				u1=[bc_d1,bc_u1,bc_c],
				p=[bcq_r])


def start_timestep_hook(t, u_, n, ds, **NS_namespace):
    u_inlet.t = t
    flux = assemble(dot(u_, n)*ds(2))
    #print flux


def pre_solve_hook(mesh, **NS_namespace):
	g = 1
	h = CellSize(mesh)
	n = FacetNormal(mesh)
	boundary = FacetFunction("size_t", mesh)
	boundary.set_all(0)
	circle.mark(boundary, 1)
	left.mark(boundary,2)
	right.mark(boundary,3)
	#up.mark(boundary,15)
	#plot(boundary,interactive=True)
	ds = Measure("ds", subdomain_data=boundary)
	R = VectorFunctionSpace(mesh, 'R', 0)
	c = TestFunction(R)

	return dict(h=h,n=n,ds=ds,c=c)



def temporal_hook(t, mesh, q_,h, u_, T, nu, dt, plot_interval, \
			      tstep, sys_comp, L_list, D_list, CFL_list, n, ds,  \
			      c, U, d,p_list, **NS_namespace):
	pressure = q_['p']
	tau = -pressure*Identity(2)+nu*(grad(u_)+grad(u_).T)
	forces = -assemble(dot(dot(tau, n), c)*ds(1)).array()*2/U**2/d
	L_list.append(forces[1])
	D_list.append(forces[0])
	print "Cd = {}, CL = {}".format(*forces)
	p_list.append(pressure(array([0.15,0.20]))-pressure(array([0.25,0.20])))
	print pressure(array([0.15,0.20]))-pressure(array([0.25,0.20]))
	#DG = FunctionSpace(mesh, "DG", 0)
	#CFL = project((dot(u_,u_)**0.5)*dt/h, DG)
	#print "CFL max %.6f" % max(CFL.vector().array())
	

def theend_hook(V, Q, U, d, h, q_, mesh, n, L_list, D_list, T, dt, u_, CFL_list, \
				CFLwrite, resultswrite, key, c, ds, p_list, **NS_namespace):

	N = len(L_list)
	Nn = int(9*N/10.)
	DG = FunctionSpace(mesh, "DG", 0)
	L_list_short = array(L_list[Nn:])
	D_list_short = array(D_list[Nn:])
	CFL_list_short = array(CFL_list[Nn:])
	p_list_short = array(p_list[Nn:])
	max_list = []; min_list = []
	for j in range(0,(len(L_list_short)-2)):
		A = L_list_short[j+2]
		B = L_list_short[j+1]
		C = L_list_short[j]
		if A > 0:
			if A < B:
				if B > C:
					max_list.append(list(L_list_short).index(B))
		if A < 0:
			if A > B:
				if B < C:
					min_list.append(list(L_list_short).index(B))
	
	f = d/(U*dt*(max_list[-1]-max_list[-2]))	# Strouhal number				

	def CFLfinder(u_,**NS_namespace):
		mesh1 = Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM14969E.xml")
		mesh2 = Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM39883E.xml")
		mesh3 = Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM119197E.xml")
		mesh4 = Mesh("/home/guttorm/Desktop/Master/Mesh/Circle/Refined/CFM372949E.xml")
		CFL_list = []
		for dt in [0.001,0.0005,0.0001,0.00005,0.000015]:
			print "dt=%.6e" % dt
			for mesh in [mesh1,mesh2,mesh3,mesh4]:
				h = CellSize(mesh)
				DG = FunctionSpace(mesh, "DG", 0)
				CFL = project((dot(u_,u_)**0.5)*dt/h, DG)
				print "CFL max %.6f" % max(CFL.vector().array())
				CFL_list.append(max(CFL.vector().array()))
			text_file = open("/home/guttorm/Desktop/Master/CFL/CFLres/CFLnumbersNew%8f.txt" %dt, "w")
			text_file.write("dt: %.6e \n" % dt)
			for ii in CFL_list:
				text_file.write("CFL max %.8f \n" % ii)

			CFL_list = []
				
		CFLfinder(u_,**NS_namespace)

	def writefile(**NS_namespace):
		import time
		now = time.strftime("%d/%m/%Y")
		os.system("mkdir -p /home/guttorm/Desktop/Master/RefinementData/Re100/ OutputCF%s" % now)
		#os.system("mkdir OutputCF%s" % now)
		text_file = open("/home/guttorm/Desktop/Master/RefinementData/Re100/OutputCF%s/OutputCircleTest%d.txt" % (now, key), "w")
		text_file.write("Cl max: %.8f \n" % max(L_list_short))
		text_file.write("Cd max: %.8f \n" % max(D_list_short))
		text_file.write("Cl min: %.8f \n" % min(L_list_short))
		text_file.write("Cd min: %.8f \n" % min(D_list_short))
		text_file.write("St = %.8f") % f
		text_file.write("Delta P = %.8f") % p_list_short[min_list[-1]]
		text_file.close()
	if resultswrite==True:	
		writefile(**NS_namespace)

	print "----------Results----------"
	print "Cl max=%.6f" % max(L_list_short)
	print "Cd max=%.6f" % max(D_list_short)
	print "Cl min=%.6f" % min(L_list_short)
	print "Cd min=%.6f" % min(D_list_short)
	print "St=%.6f" % f
	print "Delta p at t=t0 + 0.5/f)=%.6f" % p_list_short[min_list[-1]]
	print "Number of timesteps=%.d" % N



