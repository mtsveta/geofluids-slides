# Import array and plotting libraries
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np

# Define time array 
t = np.linspace(0, 15, 1000)
# Define initial conditions of y vector
y0 = np.array([10, 15])

# Define coefficients of Lotka-Volterra
a = 1.
b = 0.1
c = 1.5
d = 0.75

def f(y, t=0):
    """ Return the growth rate of predators and prey populations. """
    return np.array([ a*y[0] -   b*y[0]*y[1] ,
                     -c*y[1] + d*b*y[0]*y[1] ])
def f_y(y, t=0):
    """ Return the Jacobian matrix evaluated in X. """
    return np.array([[a -b*y[1],   -b*y[0]     ],
                     [ b*d*y[1],   -c +b*d*y[0]]])

# Solve linearized Lotka-Volterra system
y, info = odeint(f, y0, t, full_output=True)
print("Number of function evaluations: %d, number of Jacobian evaluations: %d" % (info['nfe'][-1], info['nje'][-1]))

# Plot periodic solutions 
plt.plot(t, y[:, 0], label='Predators')
plt.plot(t, y[:, 1], label='Prey')
plt.legend(loc='best')
plt.ylabel('Population [-]')
plt.xlabel('Time [s]')
plt.title('Solution of linearized Lotka-Volterra equation')
plt.savefig('lotka-volterra.pdf', bbox_inches='tight')
