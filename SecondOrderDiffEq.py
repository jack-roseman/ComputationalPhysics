from numpy import array, arange, sin, pi, cos, sqrt
from pylab import plot, show, legend, xlabel, ylabel


def f(r, t):
    theta = r[0]
    omega = r[1]
    ftheta = omega              ## dtheta/dt = ftheta
    fomega = -g/l * sin(theta)  ##domega/dt = fomega
    return array([ftheta, fomega], float)

## parameters
m = 1.0     ##mass of pendulum
g = 1.0     ##gravitational acceleration
l = 1.0     ##length of pendulum

theta0 = 0.99*pi
omega0 = 0
r = array([theta0, omega0], float)
t0 = 0
tf = 50.0
N = 1000

h = (tf - t0)/N

tpoints = arange(t0, tf, h)
thetapoints = []
omegapoints = []



for t in tpoints:
    thetapoints.append(r[0])        ##appending theta
    omegapoints.append(r[1])        ##appending omega

    ##second order Runge-Kutta
    # rmid = r + h*f(r,t)
    # r += h*f(rmid, t)

    #4th-order Runge-Kutta
    k1 = h*f(r, t)
    k2 = h*f(r + 0.5*k1, t + 0.5*h)
    k3 = h*f(r + 0.5*k2, t + 0.5*h)
    k4 = h*f(r + k3, t + h)
    r += 1/6 * (k1 + 2*k2 + 2*k3 + k4)


#plot results

plot(tpoints, thetapoints, "b-", label="Nonlinear")
plot(tpoints, theta0*cos(sqrt(g/l)*tpoints), "r-", label = "Linear")
xlabel("$t$")
ylabel("$\theta$")
legend()
show()

