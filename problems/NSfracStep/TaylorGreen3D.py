__author__ = "Mikael Mortensen <mikaem@math.uio.no>"
__date__ = "2013-06-25"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU Lesser GPL version 3 or any later version"

from ..NSfracStep import *
import numpy as np
import datetime


def mesh(N,**params):
    return BoxMesh(Point(-pi, -pi, -pi), Point(pi, pi, pi), N, N, N)

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
    T = 1,
    dt = 0.001,
    N = 10,
    folder = "taylorgreen3D_results",
    max_iter = 1,
    velocity_degree = 2,
    pressure_degree = 1,
    save_step = 10000,
    checkpoint = 10000, 
    plot_interval = 10,
    print_dkdt_info = 10,
    use_krylov_solvers = True,
    krylov_solvers = dict(monitor_convergence=False),
    kinlist = [],
    dkdtlist = [],
    savefile = True,
    path = "give_me_save_path",
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
                  dt, t, oasis_memory, kinlist, dkdtlist, **NS_namespace):
    #oasis_memory("tmp", True)
    if (tstep % print_dkdt_info == 0 or
        tstep % print_dkdt_info == 1):
        kinetic = assemble(0.5*dot(u_, u_)*dx) / (2*pi)**3
        if tstep % print_dkdt_info == 0:
            kin[0] = kinetic
            dissipation = assemble(nu*inner(grad(u_), grad(u_))*dx) / (2*pi)**3
            info_blue("Kinetic energy = {} at time = {}".format(kinetic, t)) 
            info_blue("Energy dissipation rate = {}".format(dissipation))
            kinlist.append(kinetic)
        else:
            dkdt = (kinetic-kin[0])/dt
            info_blue("dk/dt = {} at time = {}".format(dkdt, t))
            dkdtlist.append(abs(dkdt))
    #if tstep % plot_interval == 0:
    #    plot(p_, title='pressure')
    #    plot(u_[0], title='velocity-x')
    #    plot(u_[1], title='velocity-y')
        

def theend_hook(u_, p_ ,dt, nu, N, kinlist, dkdtlist, savefile, velocity_degree, pressure_degree, path, **kw):
    print "End Hook"
    if MPI.rank(mpi_comm_world()) == 0:
        if savefile == True:
            
            now = datetime.datetime.now()
            atm = "%dm%dd%dh" % (now.month, now.day, now.hour)

            np.savetxt('/%s/K_ref_P%sP%s_date%s_dt%sN%s.txt' % \
                        (path, velocity_degree,pressure_degree,atm,dt,N),  kinlist, delimiter=',')
            np.savetxt('/%s/dkdt_ref_P%sP%s_date%s_dt%sN%s.txt' % \
                        (path, velocity_degree,pressure_degree,atm,dt,N),  dkdtlist, delimiter=',')
            

        

