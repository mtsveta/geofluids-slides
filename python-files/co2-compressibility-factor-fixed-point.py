from numpy import *
import matplotlib.pyplot as plt

# Function that calculates the fixed point of equation x = f(x)
def fixed_point(f, x0):
    """
    :param f: the function
    :param x0: the initial guess for the fixed point 
    :return x: the fixed point of the equation x = f(x)
    """

    # Tolerance and max number of iterations
    maxiters = 100 
    tolerance = 1e-4 

    # Counter for the number of iterations
    counter = 0
    
    # Initial guess
    x = x0 

    # Perform one or more fixed point iterations
    for counter in range(maxiters):

        # Calculate the new approximation for x
        x_new = f(x)

        # Check the new value for convergence
        if abs(x_new - x) < tolerance: 
            return x # return x if the calculation converged

        # Update x with x_new
        x = x_new 

    # Raise an error if the calculation did not converge.
    raise RuntimeError('Could not calculate the \
        solution of the nonlinear equation in %d iterations.' % counter)
    
# Function to calculate Z using fixed point method
def compressibility_factor_fixedpoint(T, P, omega, Tc, Pc):
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
    alpha = (1 + (0.37464 + 1.54226 * omega - 0.26992 * omega**2) * (1 - sqrt(Tr)))**2

    # Contants beta and q 
    beta = Omega * Pr/Tr
    q = Psi/Omega * alpha/Tr

    # Initial guess for Z
    Z0 = 1.0

    # Define the function f that represents the nonlinear equation x = f(x)
    def f(Z):
        return (1 + beta - q*beta*(Z - beta)/((Z + epsilon*beta)*(Z + sigma*beta)))
    
    return fixed_point(f, Z0) # use fixed point function to perform the calculation of Z

# The array with temperature values in K
temperatures = linspace(60.0, 300.0, 9) + 273.15

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
plt.ylim(0.05, 1.05)
plt.title('Compressibility Factor of CO2(g)\nUsing Fixed Point Method')

# Plot one curve for each temperature in array temperatures
for T in reversed(temperatures):
    # Create a list with Z values at current T and P from the array pressures
    Z = [compressibility_factor_fixedpoint(T, P, omegaCO2, TcCO2, PcCO2) for P in pressures]
    # Plot the values of Z over all pressure points and current temperature T
    plt.plot(pressures, Z, label=r'$T=%.0f\;^{\circ}\mathrm{C}$' % (T - 273.15))

# Position the legend and save the figure
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('co2-compressibility-factor-fixed-point.pdf', bbox_inches='tight')
