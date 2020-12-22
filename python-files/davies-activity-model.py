import numpy
import matplotlib.pyplot as plt

# Function calculating activity coeffitient according to the Davis model
def gamma_davies(Zi, I):
	# Debye-Huckel parameter 
    Agamma = 0.5095
    # Vector of ionic strength in the power 1/2
    sqrtI = numpy.sqrt(I)
    return 10**(-Agamma * Zi**2 * (sqrtI / (1.0 + sqrtI) - 0.3 * I))

# Create array with values of ionic strength in a range from 0.0 to 1.0
I = numpy.linspace(0.0, 1.0, 101)

# Calculate activity coefficients for different charges of the aqueous species
gammaZ1 = gamma_davies(1.0, I)
gammaZ2 = gamma_davies(2.0, I)
gammaZ3 = gamma_davies(3.0, I)

# Plot an activity coeffitient as a function an ionic strength for different charges
plt.xlabel('Ionic Strength [molal]')
plt.ylabel('Activity Coefficient')
plt.plot(I, gammaZ1, label=r'$Z_i=1$')
plt.plot(I, gammaZ2, label=r'$Z_i=2$')
plt.plot(I, gammaZ3, label=r'$Z_i=3$')

# Position the legend and save the figure
plt.legend(loc='center right')
plt.savefig('activity-coefficient-davies.pdf', bbox_inches='tight')