
import numpy as np
from ..NSfracStep import *

def mesh(Nx, Ny, Nz, **params):
    return BoxMesh(Point(-pi, -pi, -pi), Point(pi, pi, pi), Nx, Ny, Nz)

def near(x, y, tol=1e-12):
    return bool(abs(x-y) < tol)

class PeriodicDomain(SubDomain):
    
    def inside(self, x, on_boundary):
        return bool((near(x[0], -pi) or near(x[1], -pi) or near(x[2], -pi)) and 
                (not (near(x[0], pi) or near(x[1], pi) or near(x[2], pi))) and on_boundary)

    def map(self, x, y):
        if near(x[0], pi) and near(x[1], pi) and near(x[2], pi):
            y[0] = x[0] - 2.0*pi
            y[1] = x[1] - 2.0*pi
            y[2] = x[2] - 2.0*pi
        elif near(x[0], pi) and near(x[1], pi):
            y[0] = x[0] - 2.0*pi
            y[1] = x[1] - 2.0*pi
            y[2] = x[2]
        elif near(x[1], pi) and near(x[2], pi):
            y[0] = x[0] 
            y[1] = x[1] - 2.0*pi
            y[2] = x[2] - 2.0*pi
        elif near(x[1], pi):
            y[0] = x[0] 
            y[1] = x[1] - 2.0*pi
            y[2] = x[2]
        elif near(x[0], pi) and near(x[2], pi):
            y[0] = x[0] - 2.0*pi
            y[1] = x[1] 
            y[2] = x[2] - 2.0*pi
        elif near(x[0], pi):
            y[0] = x[0] - 2.0*pi
            y[1] = x[1] 
            y[2] = x[2]            
        else: # near(x[2], pi):
            y[0] = x[0] 
            y[1] = x[1]
            y[2] = x[2] - 2.0*pi

constrained_domain = PeriodicDomain()

# Override some problem specific parameters
recursive_update(NS_parameters, dict(
    nu = 1./1000,
    rho = 1,
    T = 10,
    dt = 0.001,
    Nx = 32,
    Ny = 32, 
    Nz = 32,
    folder = "taylorgreen3D_results",
    max_iter = 1,
    velocity_degree = 2,
    pressure_degree = 1,
    save_step = 10000,
    checkpoint = 10000, 
    plot_interval = 100000,
    print_dkdt_info = 10,
    use_krylov_solvers = True,
    kinlist = [],
    dkdtlist = [],
    krylov_solvers = dict(monitor_convergence=False)
  )
)

initial_fields = dict(
        u0='sin(x[0])*cos(x[1])*cos(x[2])',
        u1='-cos(x[0])*sin(x[1])*cos(x[2])',
        u2='0',
        p='1./16.*(cos(2*x[0])+cos(2*x[1]))*(cos(2*x[2])+2)')
    
def initialize(q_, q_1, q_2, VV, initial_fields, OasisFunction, **NS_namespace):
    for ui in q_:
        vv = OasisFunction(Expression((initial_fields[ui])), VV[ui])
        vv()
        q_[ui].vector()[:] = vv.vector()[:]
        if not ui == 'p':
            q_1[ui].vector()[:] = q_[ui].vector()[:]
            q_2[ui].vector()[:] = q_[ui].vector()[:]

kin = zeros(1)
def temporal_hook(u_, p_, tstep, plot_interval, print_dkdt_info, nu, 
                  dt, t, oasis_memory,kinlist,dkdtlist, **NS_namespace):
    #oasis_memory("tmp", True)
    if (tstep % print_dkdt_info == 0 or
        tstep % print_dkdt_info == 1):
        kinetic = assemble(0.5*dot(u_, u_)*dx) / (2*pi)**3
        if tstep % print_dkdt_info == 0:
            kin[0] = kinetic
            kinlist.append(kinetic)
            info_blue("K = {}".format(kinetic))
            #dissipation = assemble(nu*inner(grad(u_), grad(u_))*dx) / (2*pi)**3
        else:
            info_blue("dk/dt = {} at time = {}".format((kinetic-kin[0])/dt, t))
            dkdtlist.append(abs((kinetic-kin[0])/dt))

            
def theend_hook(dkdtlist, kinlist, Nx, nu, dt, **kw):
    if MPI.rank(mpi_comm_world()) == 0:
        print kinlist
        import datetime
        now = datetime.datetime.now()
        atm = "%d.%d.%d.%d" % (now.year, now.month, now.day, now.hour)

        """
        np.savetxt('/home/guttorm/Desktop/Master/TaylorGreen/ReferenceResults/k_ref%sdt%snu%sN%s.txt' \
            % (atm,dt,nu,Nx), kinlist, delimiter=',')
        np.savetxt('/home/guttorm/Desktop/Master/TaylorGreen/ReferenceResults/dkdt_ref%sdt%snu%sN%s.txt' \
            % (atm,dt,nu,Nx), kinlist, delimiter=',')
        """

        np.savetxt('/home/guttorm/Desktop/Master/TaylorGreen/ReferenceResults/k_ref%sdt%snu%sN%s.txt' \
            % (atm,dt,nu,Nx), kinlist, delimiter=',')
        np.savetxt('/home/guttorm/Desktop/Master/TaylorGreen/ReferenceResults/dkdt_ref%sdt%snu%sN%s.txt' \
            % (atm,dt,nu,Nx), dkdtlist, delimiter=',')






