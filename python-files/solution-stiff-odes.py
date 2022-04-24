import numpy
import matplotlib.pyplot as plt

def stiff_solution(eps, t):
    return numpy.exp(- t / eps)

t = numpy.linspace(-1.0, 1.0, 201)

eps_1, eps_2, eps_3, eps_4 = 1, 0.1, 0.01, 0.001
x_1 = stiff_solution(eps_1, t)
x_2 = stiff_solution(eps_2, t)
x_3 = stiff_solution(eps_3, t)
x_4 = stiff_solution(eps_4, t)

plt.xlabel('x')
plt.ylabel(r'x(t)')
plt.title(r'Solution of stiff equation $x(t) = e^{-t/\varepsilon}$')
plt.plot(t, x_1, label=r'$\varepsilon$ = %.3f' % eps_1)
plt.legend(loc='best')
plt.savefig('solution-stiff-odes-eps-1.pdf', bbox_inches='tight')

plt.xlabel('x')
plt.ylabel(r'x(t)')
plt.plot(t, x_2, label=r'$\varepsilon$ = %.3f' % eps_2, color='C2')
plt.legend(loc='best')
plt.savefig('solution-stiff-odes-eps-2.pdf', bbox_inches='tight')


plt.xlabel('x')
plt.ylabel(r'x(t)')
plt.plot(t, x_3, label=r'$\varepsilon$ = %.3f' % eps_2, color='C3')
plt.legend(loc='best')
plt.savefig('solution-stiff-odes-eps-3.pdf', bbox_inches='tight')