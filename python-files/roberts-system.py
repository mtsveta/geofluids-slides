import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
import numpy as np
import time
from numpy import linalg as la
import math
import matplotlib as mpl

t0, tT = 0, 4 * 1e3             # considered interval
y0 = np.array([1.0, 0.0, 0.0])  # initial condition

neq = 3

t_clas = 0.0

with_jacobian = True

def Jf(y, t=0):
    """ Jacobian of rhs equations for the Roberts problem"""
    return np.array([[- 0.04,                    1e4 * y[2],   1e4 * y[1]],
                     [  0.04, - 1e4 * y[2] - 6 * 1e7 * y[1], - 1e4 * y[1]],
                     [  0.00,                6 * 1e7 * y[1],         0.0]])

def f(y, t=0):
    """ Rhs equations for the Roberts problem"""
    return np.array([- 0.04 * y[0] + 1e4 * y[1] * y[2],
                0.04 * y[0] - 1e4 * y[1] * y[2] - 3 * 1e7 * y[1] ** 2,
                                                 3 * 1e7 * y[1] ** 2])
num = 10000
t = np.linspace(t0, tT, num)

# -------------------------------------------------------------------------------------------------------------------- #
# Solving the primal system
# -------------------------------------------------------------------------------------------------------------------- #

names = ['$n_A$', '$n_B$', '$n_С$']
names_tilde = ['$\\tilde{n}_A$', '$\\tilde{n}_B$', '$\\tilde{n}_С$']
t_start = time.time()

if with_jacobian:   y, info = odeint(f, y0, t, full_output=True, Dfun=Jf)
else:               y, info = odeint(f, y0, t, full_output=True)
t_clas += time.time() - t_start
print("Solving the system of n unknowns:       func.eval. = %d, jacob.eval. = %d" % (info['nfe'][-1], info['nje'][-1]))

# -------------------------------------------------------------------------------------------------------------------- #
# Plotting solutions
# -------------------------------------------------------------------------------------------------------------------- #
colors = ['steelblue', 'aqua',
         'pink', 'maroon',
         'darkgreen', 'limegreen',
         'red', 'coral',
         'pink', 'maroon',
         'orange', 'darkorange',
         'gold', 'brown']
fig, ax = plt.subplots(1, 1, figsize=(6, 3))
for k in range(neq):
    ax.plot(t, y[:, k], label=names[k], color=colors[2 * k])
ax.legend(loc='best')
ax.set_xlabel('Time [s]')
plt.savefig('approximations.pdf', bbox_inches='tight')