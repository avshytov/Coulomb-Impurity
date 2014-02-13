from numpy import * 
from pylab import * 
import scipy.integrate 
from scipy.interpolate import splev, splrep

Z = 3.0 # Coulomb potential strength
Ef = -0.06
grid = 500
Nfa = 4.0 * 0.9
E = 0.00147485222961
it = 120
r_init = 10.0

data = np.load('data/denssim-Z=%g-N=%d-Nfa=%g-Ef=%g-it=%d.npz'
               %(Z, grid, Nfa, Ef, it))
rgrid = data['r']
Ugrid = data['U']
grid = len(rgrid)
#gcut = int(0.8 * grid)
Uspline = splrep(rgrid, Ugrid)

if False:
    v = linspace(0.01, 2.0*rgrid[-1], 2*grid)  # We extend U beyond rmax
    cut = (np.abs(v - 40.0)).argmin()          # and attach a tail to 
    Udisc = np.zeros((len(v)))                 # prevent spline error  
    Udisc[:cut] = splev(v[:cut], Uspline)      # and repulsive force
    Udisc[cut:] = Udisc[cut-1] * (v[cut-1] / v[cut:])**3
    Uspline = splrep(v,Udisc)
    dUdisc = (Udisc[1:]-Udisc[:-1])/(v[1]-v[0])
    regions = (E-Udisc)**2 - (0.5 / v)**2

    figure()
    plot(rgrid,Ugrid,label='Udata')
    plot(v, Udisc, label='U')
    plot(v, -splev(v,Uspline, der=1), label='Force')
    plot(v[:-1], -dUdisc, label='discrete')
    plot(v, regions, label='regions')
    legend()

#
#   The equations of motion in a Coulomb potential are solved as 
#   follows. Let us introduce the quantity 
#
#        u = gamma * v          gamma = 1 / sqrt(1 - v^2)
#        v = u / gamma = u / sqrt(1 + u^2)
#
#    Then: 
#
#       dx/dt = vx  = ux / sqrt (1 + u^2)
#       dy/dt = vy  = ....
#
#       dux / dt = Fx = A * x / d^3
#       duy / dt = Fy = .... 
#   
#    This is the system of equations solved by this code via
#    scipy.integrate.odeint
#

#
#  Right-hand side of the equations of motion
#

def rhs_coul(u, t):
    x, y, ux, uy = u[0], u[1], u[2], u[3]
    F = zeros((4,))
    d = math.sqrt(x**2 + y**2)
    u2 = ux**2 + uy**2
    invgamma = 1.0 / math.sqrt(1 + u2)
    d = math.sqrt(x**2 + y**2)
    vx = ux * invgamma     #   u = v / sqrt(1 - v^2)
    vy = uy * invgamma     #   v = u / sqrt(1 + u^2)
    F[0] = vx #  dx / dt = vx
    F[1] = vy #  dy / dt = vy
    F[2] = Z * x / d**3  # dux / dt = Fx = A x / d^3
    F[3] = Z * y / d**3  # duy / dt = Fy = .... 
    return F

def rhs_massless(u, t):
    x, y, px, py = u[0], u[1], u[2], u[3]
    F = zeros((4,))
    d = math.sqrt(x**2 + y**2)
    vx = px / np.sqrt(px**2 + py**2)
    vy = py / np.sqrt(px**2 + py**2)
    Ud = - splev(d,Uspline,der=1)
    F[0] = vx #  dx / dt = vx 
    F[1] = vy #  dy / dt = vy
    F[2] = Ud * x / d
    F[3] = Ud * y / d
    if Ud > 0:
        print "x", x, "y", y, "r", d
        print "Fx", Ud * x / d, "Fy", Ud * y / d
    return F

#
#   One step of RK4 method
#

def rk4step(rhs, tau, u):
    u1 = tau * rhs (u, tau)
    u2 = tau * rhs (u + 0.5 * u1, tau)
    u3 = tau * rhs (u + 0.5 * u2, tau)
    u4 = tau * rhs (u + u3, tau)
    return u + (u1 + 2.0 * u2 + 2.0 * u3 + u4) / 6.0;
#
#  
#
def solve (x0, y0, px0, py0, tmin, tmax, tau):
    print "Solve: ", x0, y0, px0, py0
    u0 = zeros ((4,))
    u0[0] = x0
    u0[1] = y0
    u0[2] = px0
    u0[3] = py0
    xs = [] 
    ys = [] 
    ts = []
    t = tmin
    uold = zeros ((4,))
    unew = zeros ((4,))
    uold[:] = u0
    
#    t = arange (tmin, tmax, tau)
#    u = scipy.integrate.odeint (rhs_coul, u0, t)
#    xs = u[:, 0]
#    ys = u[:, 1]
#    return xs, ys, t
    
    while t < tmax:
          unew = rk4step(rhs_massless, tau, uold)
          xs.append (unew[0])
          ys.append (unew[1])
          t += tau
          uold = unew
          ts.append (t)
    return array(xs), array(ys), array(ts)      

figure()
m = 0.5
x0 = r_init
y0 = 0.0
U0 = splev(r_init, Uspline)
print "E-U(r_init) = ", E-U0
py0 = m / r_init #0.23
px0 = np.sqrt((E - U0)**2 - py0**2)
vx0 = px0 / np.sqrt(px0**2 + py0**2)
vy0 = py0 / np.sqrt(px0**2 + py0**2)
#vy0 = 0.0
tmax = 250.0
tau = 0.005
xs, ys, ts = solve (x0, y0, px0, py0, 0.0, tmax, tau)
xlim(-16, 16)
ylim(-16, 16)
plot (xs, ys, 'k', linewidth=2)
plot(xs[:10000], ys[:10000], 'r', linewidth=2)
plot ([0], [0], 'ko')
if (vy0 < 0.0999):
    ax = axes([0.6, 0.2, 0.27, 0.27])
    ax.plot (xs, ys)
    xlim(-0.05, 0.05)
    ylim(-0.05, 0.05)
axes().set_aspect('equal')
show()

if False:
    x0 = -10.0     # Initial position, imitates infinity
    vx0 = 1.0      # x velocity
    vy0 = 0.0      # y velocity
    tmax = 30.0   # max time
    tau = 0.001    # step
    
    #  List of impact parameters. Modify to obtain 
    #  a different set of trajectories. 
    trajectories = [
        0.80, 
        # 0.60,
        0.40, 
        #  0.30,
        0.20, 
        #  0.15, 
        0.10, 
        #  0.05
    ]
    
    figure()
    for y0 in trajectories: 
        xs, ys, ts = solve (x0, y0, vx0, vy0, 0.0, tmax, tau)
        plot (xs, ys)
    plot ([0], [0], 'ko') 
    text (0.0, -0.1, 'The sun or a nucleus', horizontalalignment='center', 
          bbox={})
    xlim(-3, 3)
    ylim(-1.0, 1.0)
    xlabel ('x')
    ylabel ('y')

    ax = axes([0.6, 0.6, 0.27, 0.27])
    for y0 in trajectories:
        xs, ys, ts = solve (x0, y0, vx0, vy0, 0.0, tmax, tau)
        ax.plot (xs, ys)
    xlim(-0.04, 0.04)
    ylim(-0.04, 0.04)
    xticks([-0.04, -0.02, 0.0, 0.02, 0.04])
    yticks([-0.04, -0.02, 0.0, 0.02, 0.04])

show()
