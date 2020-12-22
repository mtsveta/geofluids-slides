from numpy import *
import matplotlib.pyplot as plt

# Function that solves the nonlinear equation r(x) = 0
def newton(r, rprime, x0):
    """
    :param r: the function
    :param rprime: the first derivative of f
    :param x0: the initial guess for the root 
    :return x: the root of the equation f(x) = 0
    """

    # Tolerance and max number of iterations
    maxiters = 100 
    tolerance = 1e-4 

    # Counter for the number of iterations
    counter = 0 

    # Initial guess
    x = x0 

    # Perform one or more Newton iterations
    for counter in range(maxiters):

        # Calculate the new approximation for x
        x = x - r(x) / rprime(x) 

        # Check the new value for convergence
        if abs(r(x)) < tolerance: 
            return x # return x if the calculation converged

    # Raise an error if the calculation did not converge.
    raise RuntimeError('Could not calculate the \
        solution of the nonlinear equation in %d iterations.' % counter)

# Function to calculate Z using Newton's method
def compressibility_factor_newton(T, P, omega, Tc, Pc):
    """
    :param T: temperature
    :param P: pressure
    :param omega: the acentric factor
    :param Tc: critical temperature
    :param Pc: critical pressure
    :return Z: compressibility factor
    """

    # Parameters of Peng-Robinson EOS    
    epsilon = 1 - sqrt(2.0)
    sigma = 1 + sqrt(2.0)
    Omega = 0.07780
    Psi = 0.45724

    # Reduced temperature and pressure
    Tr = T / Tc
    Pr = P / Pc

    # Evaluation of the alpha(Tr, omega) function
    alpha = (1 + (0.37464 + 1.54226*omega - 0.26992*omega**2)*(1 - sqrt(Tr)))**2

    # Contants beta and q 
    beta = Omega * Pr / Tr
    q = Psi/Omega * alpha/Tr

    # Residual function r that represents the nonlinear equation r(x) = 0
    def r(Z):
        return (1 + beta - q*beta*(Z - beta)/((Z + epsilon*beta)*(Z + sigma*beta))) - Z

    # First order derivative of function r'(x)
    def rprime(Z):
        aux = (Z + epsilon*beta)*(Z + sigma*beta)
        return -q*beta/aux * (1.0 - (Z - beta)*(2*Z + (epsilon + sigma)*beta)/aux) - 1

    # Initial guess for Z
    Z0 = 1.0

    return newton(r, rprime, Z0) # use newton function to perform the calculation of Z

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

# Plot one curve for each temperature in array temperatures
for T in reversed(temperatures):
    # Create a list with Z values at current T and P from the array pressures
    Z = [compressibility_factor_newton(T, P, omegaCO2, TcCO2, PcCO2) for P in pressures]
    # Plot the values of Z over all pressure points and current temperature T
    plt.plot(pressures, Z, label=r'$T=%.0f\;^{\circ}\mathrm{C}$' % (T - 273.15))

# Position the legend and save the figure
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('co2-compressibility-factor-newton.pdf', bbox_inches='tight')
