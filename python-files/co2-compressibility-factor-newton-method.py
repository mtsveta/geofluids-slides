from numpy import *
import matplotlib.pyplot as plt

# Define a function that solves the nonlinear equation f(x)=0,
# where f is the function, fprime its first derivative, and x0 the initial guess for x.
def newton(f, fprime, x0):
    maxiters = 100 # maximum number of iterations
    tolerance = 1e-4 # the tolerance for the convergence
    counter = 0 # the counter of number of iterations
    x = x0 # start with the solution x being the initial guess
    # Perform one or more Newton iterations
    for counter in xrange(maxiters):
        x = x - f(x) / fprime(x) # calculate the new approximation for x
        if abs(f(x)) < tolerance: # check for convergence
            return x # return x if the calculation converged
    # Raise an error if the calculation did not converge.
    raise RuntimeError('Could not calculate the \
        solution of the nonlinear equation in %d iterations.' % counter)

# Define a function to calculate Z using Newton's method,
# where T is temperature, P is pressure, omega is the acentric factor
# of the substance, Tc and Pc are the critical temperature and pressure
# of the substance.
def compressibility_factor_newton(T, P, omega, Tc, Pc):
    # Create variables for the parameters of Peng-Robinson EOS
    epsilon = 1 - sqrt(2.0)
    sigma = 1 + sqrt(2.0)
    Omega = 0.07780
    Psi = 0.45724
    # The reduced temperature, Tr, and reduced pressure, Pr
    Tr = T / Tc
    Pr = P / Pc
    # The evaluation of the alpha(Tr, omega) function
    alpha = (1 + (0.37464 + 1.54226*omega - 0.26992*omega**2)*(1 - sqrt(Tr)))**2
    # The beta and q contants
    beta = Omega * Pr / Tr
    q = Psi/Omega * alpha/Tr
    # Define the function f that represents the nonlinear equation f(x) = 0
    def f(Z):
        return (1 + beta - q*beta*(Z - beta)/((Z + epsilon*beta)*(Z + sigma*beta))) - Z
    # Define the first order derivative of function f'(x)
    def fprime(Z):
        aux = (Z + epsilon*beta)*(Z + sigma*beta)
        return -q*beta/aux * (1.0 - (Z - beta)*(2*Z + (epsilon + sigma)*beta)/aux) - 1
    # Set the initial guess
    Z0 = 1.0
    return newton(f, fprime, Z0) # use newton function to perform the calculation of Z

# The array with temperature values in K
temperatures = linspace(0.0, 300.0, 11) + 273.15

# The array with pressure values in bar
pressures = linspace(1.0, 400.0, 101)

# The critical temperature and pressure of CO2
TcCO2 = 304.2
PcCO2 = 73.83

# The acentric factor of CO2
omegaCO2 = 0.224

# Create a plot
plt.xlabel('Pressure [bar]')
plt.ylabel(r'$Z = \frac{PV}{RT}$')
plt.title('Compressibility Factor of Carbon Dioxide\nUsing Newton\'s Method')

# Plot one curve for each temperature, T, in array temperatures
for T in temperatures:
    # Create a list with the values of Z at current T and pressure P from the array of pressures
    Z = [compressibility_factor_newton(T, P, omegaCO2, TcCO2, PcCO2) for P in pressures]
    # Plot the values of Z over all pressure points and current temperature T
    plt.plot(pressures, Z, label=r'$T=%.0f\;^{\circ}\mathrm{C}$' % (T - 273.15))

# Position the legend and save the figure
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('../figures/compressibility-factor-co2-newton.pdf', bbox_inches='tight')
plt.show()
