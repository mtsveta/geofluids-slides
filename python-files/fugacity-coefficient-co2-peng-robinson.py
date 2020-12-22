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

# Function to calculate the fugacity coefficient of CO2, phiCO2, using Peng-Robinson EOS
def fugacity_coefficient_co2(T, P):
    """
    :param T: temperature
    :param P: pressure
    :return phiCO2: cfugacity coefficient of CO2
    """

    # The critical temperature and pressure of CO2
    Tc = 304.2
    Pc = 73.83

    # The acentric factor of CO2
    omega = 0.224

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
    def r(Z):
        return (1 + beta - q*beta*(Z - beta)/((Z + epsilon*beta)*(Z + sigma*beta))) - Z
    # Define the first order derivative of function f'(x)
    def rprime(Z):
        aux = (Z + epsilon*beta)*(Z + sigma*beta)
        return -q*beta/aux * (1.0 - (Z - beta)*(2*Z + (epsilon + sigma)*beta)/aux) - 1
    
    # Calculate the compressibility factor of CO2 at (T, P)
    Z0 = 1.0 # the initial guess for the compressibility factor
    Z = newton(r, rprime, Z0) # use newton function to perform the calculation of Z
    
    # Calcute theta
    theta = 1.0/(sigma - epsilon) * log((Z + sigma*beta)/(Z + epsilon*beta))
    
    # Calculate phiCO2
    phiCO2 = exp(Z - 1 - log(Z - beta) - q*theta)

    return phiCO2


# The array with temperature values in K
temperatures = linspace(0.0, 300.0, 11) + 273.15

# The array with pressure values in bar
pressures = linspace(1.0, 400.0, 101)

# Create a plot
plt.xlabel('Pressure [bar]')
plt.ylabel(r'$\varphi_\mathrm{CO_2(g)}$')
plt.title('Fugacity Coefficient of Carbon Dioxide\nUsing Peng-Robinson EOS')

# Plot one curve for each temperature, T, in array temperatures
for T in reversed(temperatures):
    # Create a list with the values of Z at current T and pressure P from the array of pressures
    Z = [fugacity_coefficient_co2(T, P) for P in pressures]
    # Plot the values of Z over all pressure points and current temperature T
    plt.plot(pressures, Z, label=r'$T=%.0f\;^{\circ}\mathrm{C}$' % (T - 273.15))

# Position the legend and save the figure
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('fugacity-coefficient-co2-peng-robinson.pdf', bbox_inches='tight')
