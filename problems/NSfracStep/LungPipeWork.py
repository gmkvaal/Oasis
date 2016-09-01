from ..NSfracStep import *
import numpy as np


mesh = Mesh("/home/guttorm/Dropbox/Guttorm/MeshLab/vtp/MeshFiles/LP_coarse_mesh.xml")
#mesh = Mesh("/home/guttorm/Dropbox/Guttorm/MeshLab/vtp/MeshFiles/lp_095_std.xml")
#mesh = Mesh("/usit/abel/u1/gmkvaal/TestDir/blc0tf085slr075vef1_meter.xml")
#mesh = Mesh("/usit/abel/u1/gmkvaal/TestDir/std_cd08.xml")
#mesh = Mesh("/usit/abel/u1/gmkvaal/TestDir/dist_03_1_1d2_noBL.xml")


NS_parameters.update(
	center = [],
	inlet_normal = [],
	nu = 1.48e-5,
	T = 10,
	dt = 0.01,
    plot_interval = 1,
    save_step = 1,
    velocity_degree=1,
    print_intermediate_info = 10, 
    name = "lungtest",
    counter = 0,
    inlet_ID = 1,    
	outlet_ID = 2,
	check_flux = 10,
	folder = "Results"
    )

#def pre_solve_hook(name,T,dt,save_step,** NS_namespace):
#	return dict(folder = "%s_T%s_dt%s_sstep%s" % (name,T,dt,save_step))

def create_bcs(V, Q, mesh, center, inlet_normal, inlet_ID, outlet_ID,** NS_namespace):

	print mesh
	print "DOF=%d" % int(V.dim() + Q.dim())

	D = 3		 	# Dimentions

	p0 = project(Expression("x[0]"), V)
	p1 = project(Expression("x[1]"), V)
	p2 = project(Expression("x[2]"), V)

	# Compute inlet area
	fd = MeshFunction("size_t", mesh, 2, mesh.domains())
	dsi = ds(inlet_ID, domain=mesh, subdomain_data=fd)
	inlet_area = assemble(Constant(1)*dsi)

	# Compute center of inlet coordinate
	for p in [p0, p1, p2]:
		center.append(assemble(p*dsi) / inlet_area)

	# Compute inlet normal vector
	n = FacetNormal(mesh)
	ni = np.array([assemble(n[i]*dsi) for i in xrange(D)])
	n_len = np.sqrt(sum([ni[i]**2 for i in xrange(D)])) # Should always be 1!?
	i_n = -ni / n_len
	for i in range(len(i_n)):
		inlet_normal.append(i_n[i])


	if MPI.rank(mpi_comm_world()) == 0:
		print "-----Inlet info-----"
		print "ID = %d" % inlet_ID
		print "Area = %.2f" % inlet_area
		print "Center = (%.2f,%.2f,%.2f)" % (center[0],center[1],center[2])
		print "Normal = (%.6f,%.6f,%.6f)" % (inlet_normal[0],inlet_normal[1],inlet_normal[2])

	R = sqrt(inlet_area/pi)
	class IP(Expression):
		# Generate a parabolic inlet profile
		def eval(self, values, x):

			r = np.sqrt((center[0]-x[0])**2+(center[1]-x[1])**2+(center[2]-x[2])**2)
			inlet_prof = 2*1.25*(1 - r**2/R**2)

			values[:] = inlet_prof

		def shape_value(self):
			return (1,)

	"""
	# Does same as class IP		
	R = sqrt(inlet_area/np.pi)
	parabolic = Expression("1 - pow(sqrt(pow((%s-x[0]),2)  + pow((%s-x[1]),2) + pow((%s-x[2]),2) ) / %s, 2)" % \
			(center[0],center[1],center[2],R) )
	"""

	ip = IP()
	noslip = Constant(0)

	bc1x = DirichletBC(V,ip*inlet_normal[0],inlet_ID)
	bc1y = DirichletBC(V,ip*inlet_normal[1],inlet_ID)
	bc1z = DirichletBC(V,ip*inlet_normal[2],inlet_ID)
	bc0 = DirichletBC(V,noslip,0)
	bc2p = DirichletBC(Q,noslip,outlet_ID)

	return dict(u0=[bc1x,bc0],
				u1=[bc1y,bc0],
				u2=[bc1z,bc0],
				p=[bc2p])


def pre_solve_hook(velocity_degree, mesh, **NS_namesepace):

    Vv = VectorFunctionSpace(mesh, 'CG', velocity_degree,
                            constrained_domain=constrained_domain)

    uv = Function(Vv)
    n = FacetNormal(mesh)
    normal = FacetNormal(mesh)
    domains = FacetFunction('size_t', mesh, 0)


    return dict(uv=uv,n=n,domains=domains)


def temporal_hook(mesh, u_, V, uv, name, n, domains, inlet_ID, outlet_ID, tstep, check_flux, **NS_namespace):

		if tstep % check_flux == 0:

			[assign(uv.sub(i), u_[i]) for i in range(mesh.geometry().dim())]
			fd = MeshFunction("size_t", mesh, 2, mesh.domains())
			ds1 = ds(outlet_ID, domain=mesh, subdomain_data=fd)
			ds2 = ds(inlet_ID, domain=mesh, subdomain_data=fd)
			ds0 = ds(0, domain=mesh, subdomain_data=fd)

			inlet_flux = assemble(dot(uv, n)*ds2)
			outlet_flux = assemble(dot(uv, n)*ds1)
			walls_flux = assemble(dot(uv, n)*ds0)

			if MPI.rank(mpi_comm_world()) == 0: 
				print ""
				print "inlet flux = %.3f" % inlet_flux
				print "outlet flux = %.3f" % outlet_flux	
				print "wall flux = %.3e" % walls_flux
				print "delta flux = %.6f" % (abs(inlet_flux) - abs(outlet_flux))
				print "leakage = %.6f" % ((abs(inlet_flux) - abs(outlet_flux))/inlet_flux)
				