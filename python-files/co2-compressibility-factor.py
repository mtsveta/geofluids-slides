from numpy import *
import matplotlib.pyplot as plt

# Define a function to calculate Z using fixed point method,
# where T is temperature, P is pressure, omega is the acentric factor
# of the substance, Tc and Pc are the critical temperature and pressure
# of the substance.
def compressibility_factor_fixedpoint(T, P, omega, Tc, Pc):
    # Create variables for the parameters of Peng-Robinson EOS
    epsilon = 1 - sqrt(2.0)
    sigma = 1 + sqrt(2.0)
    Omega = 0.07780
    Psi = 0.45724
    # The tolerance and max number of iterations
    tolerance = 1.0e-4
    maxiterations = 100
    # The reduced temperature, Tr, and reduced pressure, Pr
    Tr = T / Tc
    Pr = P / Pc
    # The evaluation of the alpha(Tr, omega) function
    alpha = (1 + (0.37464 + 1.54226 * omega - 0.26992 * omega**2) * (1 - sqrt(Tr)))**2
    # The beta and q contants
    beta = Omega * Pr/Tr
    q = Psi/Omega * alpha/Tr
    # The initial guess for Z
    Z = 1.0
    # A counter variable for the number of iterations
    counter = 0
    # Perform the fixed point iterations 
    for counter in xrange(maxiterations):
        Znew = 1 + beta - q * beta * (Z - beta)/((Z + epsilon*beta)*(Z + sigma*beta))
        if abs(Znew - Z) < tolerance: # Check for convergence
            return Znew # Return the new value of Z in case of convergence
        Z = Znew # Update Z with Znew
    # Raise an error if the fixed point iteration did not converge within max iterations
    raise RuntimeError('Could not calculate the ' + \
        'compressibility factor in %d iterations at %f C and %f bar.' % (counter, T, P))

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
plt.title('Compressibility Factor of Carbon Dioxide\nUsing Fixed Point Method')

# Plot one curve for each temperature, T, in array temperatures
for T in temperatures:
    # Create a list with the values of Z at current T and pressure P from the array of pressures
    Z = [compressibility_factor_fixedpoint(T, P, omegaCO2, TcCO2, PcCO2) for P in pressures]
    # Plot the values of Z over all pressure points and current temperature T
    plt.plot(pressures, Z, label=r'$T=%.0f\;^{\circ}\mathrm{C}$' % (T - 273.15))

# Position the legend and save the figure
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
