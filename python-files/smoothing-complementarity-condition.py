import numpy
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

def inverse(tau, x):
    return tau / x;

x = numpy.linspace(1e-7, 1e-5, 10000001)

tau_1, tau_2, tau_3, tau_4 = 1, 0.5, 0.1, 0.05
z_1 = inverse(tau_1, x)
z_2 = inverse(tau_2, x)
z_3 = inverse(tau_3, x)
z_4 = inverse(tau_4, x)

fig , ax = plt.subplots()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))

plt.xlabel(r'$x_i$')
plt.ylabel(r'$z_i$')
plt.title(r'Smoothing the condition $x_i \, z_i = 0$')
plt.xticks([1e-6, 5e-6, 1e-5])

plt.plot(x, z_1, label=r'$\tau$ = %.2f' % tau_1)
plt.plot(x, z_2, label=r'$\tau$ = %.2f' % tau_2)
plt.plot(x, z_3, label=r'$\tau$ = %.2f' % tau_3)
plt.plot(x, z_4, label=r'$\tau$ = %.2f' % tau_4)

plt.legend(loc='best')
plt.savefig('smoothing-complementarity-condition.pdf', bbox_inches='tight')