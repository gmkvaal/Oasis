from ..NSCoupled import *

mesh = Mesh("/home/guttorm/Desktop/Master/MeshLab/vtp/MeshFiles/test_file_mesh.xml")

print mesh






def create_bcs(VQ, mesh, **NS_namespce):
	
	centers = []
	p0 = project(Expression("x[0]"), VQ.sub(0))
	p1 = project(Expression("x[1]"), VQ.sub(0))
	p2 = project(Expression("x[2]"), VQ.sub(0))


	inlet_ID = 1
	fd = MeshFunction("size_t", mesh, 2, mesh.domains())
	#from IPython import embed; embed()
	dsi = ds(inlet_ID, domain=mesh, subdomain_data=fd)
	inlet_area = assemble(Constant(1)*dsi)

	print inlet_area

	center = []
	for p in [p0, p1, p2]:
		center.append(assemble(p*dsi) / inlet_area)

	print center
	# Compute average normal (assuming boundary is actually flat)

	D = 3
	n = FacetNormal(mesh)
	ni = np.array([assemble(n[i]*dsi) for i in xrange(D)])
	n_len = np.sqrt(sum([ni[i]**2 for i in xrange(D)])) # Should always be 1!?
	inlet_normal = -ni / n_len

	print inlet_normal


	bc0 = DirichletBC(VQ.sub(0), (0,0,0), 1)
	bc1 = DirichletBC(VQ.sub(0), (0,0,0), 0)
	bc2 = DirichletBC(VQ.sub(1), 0, 2)
	return dict(up=[bc0, bc1, bc2])  

def theend_hook(u_,**NS_namespace):

	plot(u_,interactive=True)
