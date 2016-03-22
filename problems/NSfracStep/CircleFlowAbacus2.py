
import matplotlib.pyplot as plt
from ..NSfracStep import *
from numpy import sqrt,array
set_log_active(False)
parameters['allow_extrapolation'] = True

			
#from IPython import embed; embed()

NS_parameters.update(
	nu = 0.001,
	U=1,
	d=0.1,
	T = 8,
	dt = 0.001,
    plot_interval = 1000,
    save_step = 1000,
    use_krylov_solvers=True,
    print_intermediate_info = 100,
    L_list = [],
    D_list = [],
    CFL_list = [],
    key = 1,
    CFLwrite = False,
    resultswrite =False   	
    )

def mesh(key,**params): 
	if key == 1: return Mesh("/uio/hume/student-u61/gmkvaal/Master/Mesh/Circle/Refined/External/CFM14969E.xml")
	if key == 2: return Mesh("/uio/hume/student-u61/gmkvaal/Master/Mesh/Circle/Refined/External/CFM39883E.xml")
	if key == 3: return Mesh("/uio/hume/student-u61/gmkvaal/Master/Mesh/Circle/Refined/External/CFM119197E.xml")
	if key == 4: return Mesh("/uio/hume/student-u61/gmkvaal/Master/Mesh/Circle/Refined/External/CFM372949Eext.xml")
	if key == 10: return Mesh("/uio/hume/student-u61/gmkvaal/Master/Mesh/Circle/Refined/External/CFME7169Eext.xml")
	if key == 20: return Mesh("/uio/hume/student-u61/gmkvaal/Master/Mesh/Circle/Refined/External/CFM24927Eext.xml")
	if key == 30: return Mesh("/uio/hume/student-u61/gmkvaal/Master/Mesh/Circle/Refined/External/CFM94113Eext.xml")
	if key == 40: return Mesh("/uio/hume/student-u61/gmkvaal/Master/Mesh/Circle/Refined/External/CFM372949Eext.xml")

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
			
			
circle = Circle()
up = Up()
down = Down()
left = Left()
right = Right()

u0 = Expression('4*1.5*x[1]*(0.41-x[1])/(0.41*0.41)')

def create_bcs(V,Q,** NS_namespace):
	bc_l0 = DirichletBC(V,u0,left)
	bc_l1 = DirichletBC(V,0.0,left)
	bc_r0 = DirichletBC(V,u0,right)
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

"""
def start_timestep_hook(t, **NS_namespace):
    u0.t = t
"""

def pre_solve_hook(mesh, **NS_namespace):
	h = CellSize(mesh)
	n = FacetNormal(mesh)
	boundary = FacetFunction("size_t", mesh)
	boundary.set_all(0)
	circle.mark(boundary, 1)
	ds = Measure("ds", subdomain_data=boundary)
	R = VectorFunctionSpace(mesh, 'R', 0)
	c = TestFunction(R)
	return dict(h=h,n=n,ds=ds,c=c)

U = 1
d = 0.1

def temporal_hook(u0, t, mesh, q_,h, u_, T, nu, dt, plot_interval, \
			      tstep, sys_comp, L_list, D_list, CFL_list, n, ds, c, \
			      U, d, **NS_namespace):

        
	pressure = q_['p']
	tau = -pressure*Identity(2)+nu*(grad(u_)+grad(u_).T)
	forces = -assemble(dot(dot(tau, n), c)*ds(1)).array()*2/U**2/d
	

        comm = mpi_comm_world()
        mpiRank = MPI.rank(comm)
        print mpiRank
	print "Cd = {}, CL = {}".format(*forces)
	L_list.append(forces[1])
	D_list.append(forces[0])
	
	
	DG = FunctionSpace(mesh, "DG", 0)
	CFL = project((dot(u_,u_)**0.5)*dt/h, DG)
	print "CFL max %.6f" % max(CFL.vector().array())


def theend_hook(V, mesh, L_list, D_list, T, dt, u_, q_, CFL_list,\
				CFLwrite ,resultswrite ,key ,**NS_namespace):
	N = len(L_list)
	n = int(9*N/10.)
	h = CellSize(mesh)
	DG = FunctionSpace(mesh, "DG", 0)
	CFL = project((dot(u_,u_)**0.5)*dt/h, DG)
	print "CFL max %.6f" % max(CFL.vector().array())
	#plot(CFL,interactive=True,range_min=0.0,range_max=0.1,title="CFL field")

	L_list_short = array(L_list[n:])
	D_list_short = array(D_list[n:])
	CFL_list_short = array(CFL_list[n:])
	max_list = []
	for j in range(0,(len(L_list_short)-2)):
		A = L_list_short[j+2]
		B = L_list_short[j+1]
		C = L_list_short[j]
		if A > 0:
			if A < B:
				if B > C:
					max_list.append(list(L_list_short).index(B))
	"""
	def Strouhal(max_list,**NS_namespace):
		print max_list[-1]-max_list[-2]	
		f = 1/(dt*(max_list[-1]-max_list[-2]))
		print f
		d=0.1;U=1
		print d*f/U
		return d*f/U
	Strouhal(max_list,**NS_namespace)	
	"""

	def writefile(**NS_namespace):
		text_file = open("/uio/hume/student-u61/gmkvaal/Master/RefinementData/Re100/Refine%ddt%.6f.txt" %(key,dt), "w")
		text_file.write("Cl max: %.8f \n" % max(L_list_short))
		text_file.write("Cd max: %.8f \n" % max(D_list_short))
		text_file.write("Cl min: %.8f \n" % min(L_list_short))
		text_file.write("Cd min: %.8f \n" % min(D_list_short))
		text_file.close()
	if resultswrite==True:	
		writefile(**NS_namespace)

 
	print "----------Here are the last 10 percents----------"
	print "Cl:"
	print L_list_short
	print "Cd:"
	print D_list_short
	print "Cl last=%e, Cd last=%e" % (L_list_short[-1],D_list_short[-1])
	print "----------Here are Cd and Cl max and min----------"
	print "Cl max=%.6f" % max(L_list_short)
	print "Cd max=%.6f" % max(D_list_short)
	print "Cl min=%.6f" % min(L_list_short)
	print "Cd min=%.6f" % min(D_list_short)
	print "With index respectively %.d and %.d" % \
			(list(L_list_short).index(max(L_list_short)), \
			list(L_list_short).index(max(L_list_short)))
	print "Number of timesteps=%.d" % N
	print max_list
	#plt.plot(L_list_short)
	#plt.show()
	#return u_