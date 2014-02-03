from numpy import * 
from pylab import * 
import scipy.integrate 

A = -0.1 # Coulomb potential strength

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
    F[2] = A * x / d**3  # dux / dt = Fx = A x / d^3
    F[3] = A * y / d**3  # duy / dt = Fy = .... 
    return F

#
#   One step of RK4 method: not used now
#
def rk4step(rhs, tau, u):
    u1 = tau * rhs (u)
    u2 = tau * rhs (u + 0.5 * u1)
    u3 = tau * rhs (u + 0.5 * u2)
    u4 = tau * rhs (u + u3)
    return u + (u1 + 2.0 * u2 + 2.0 * u3 + u4) / 6.0;

#
#  
#
def solve (x0, y0, vx0, vy0, tmin, tmax, tau):
    u0 = zeros ((4,))
    u0[0] = x0
    u0[1] = y0
    gamma = 1.0 / math.sqrt (1.0 - vx0**2 - vy0**2)
    u0[2] = vx0 * gamma
    u0[3] = vy0 * gamma 
    xs = [] 
    ys = [] 
    ts = []
    t = tmin
    uold = zeros ((4,))
    unew = zeros ((4,))
    uold[:] = u0
    
    t = arange (tmin, tmax, tau)
    u = scipy.integrate.odeint (rhs_coul, u0, t)
    xs = u[:, 0]
    ys = u[:, 1]
    return xs, ys, t
    
    print "Solve: ", x0, y0, vx0, vy0
    while t < tmax:
          unew = rk4step(rhs_coul, tau, uold)
          xs.append (unew[0])
          ys.append (unew[1])
          t += tau
          uold = unew
          ts.append (t)
    return array(xs), array(ys), array(ts)      
figure() 
x0 = 1.0
y0 = 0.0
vx0 = 0.0
vy0 = 0.23 #0.23
#vy0 = 0.0
tmax = 300.0
tau = 0.001
xs, ys, ts = solve (x0, y0, vx0, vy0, 0.0, tmax, tau)
plot (xs, ys)
if (vy0 < 0.0999):
    ax = axes([0.6, 0.2, 0.27, 0.27])
    ax.plot (xs, ys)
    xlim(-0.05, 0.05)
    ylim(-0.05, 0.05)
show()

    
x0 = -10.0     # Initial position, imitates infinity
vx0 = 0.5      # x velocity
vy0 = 0.0      # y velocity
tmax = 100.0   # max time
tau = 0.001    # step

#  List of impact parameters. Modify to obtain 
#  a different set of trajectories. 
trajectories = [
     0.80, 
     0.60,
     0.40, 
     0.30,
     0.20, 
     0.15, 
     0.10, 
     0.05
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
