# Import array and plotting libraries
import numpy
import matplotlib.pyplot as plt

# Set matplotlib plotting parameters
plt.style.use('ggplot')
plt.rc('font', size=16)
plt.rc('font', family='serif', serif=['New Century Schoolbook'])
plt.rc('text', usetex=True)
plt.rc('lines', linewidth=2)
plt.rc('figure', autolayout=True)

# Set auxiliary constants
Agamma = 0.5095
ln10 = numpy.log(10.0)

# Function calculating activity coeffitient of water according to the ideal model
def activity_water_ideal(I, xH2O):
    return numpy.exp(-(1 - xH2O)/xH2O)

# Function calculating activity coeffitient of water according to the Davis model
def activity_water_davies(I, xH2O):
    sqrtI = numpy.sqrt(I)
    return numpy.exp(ln10/55.5084*Agamma*(2*(I + 2*sqrtI)/(1 + sqrtI) - 4 * numpy.log(1 + sqrtI) - 0.3 * I**2) - (1 - xH2O)/xH2O)

# Create array with values of ionic strength in a range from 0.0 to 1.0
I = numpy.linspace(0.0, 1.0, 101)

# Assuming H2O+NaCl, with 1 kg H2O and I moles of Na+ and Cl-, so that
# xH2O = 55.5084/(55.5084 + nNa+ + nCl-) = 55.5084/(55.5084 + 2*I)
xH2O = 55.5084/(55.5084 + 2*I)

# Calculate activity coefficients for Davis and ideal model
aH2O_davies = activity_water_davies(I, xH2O)
aH2O_ideal = activity_water_ideal(I, xH2O)

# Plot activity coeffitients of Davis vs ideal model as a function an ionic strength
plt.figure()
plt.ylim(0.96, 1.0)
plt.xlabel('Ionic Strength [molal]')
plt.title(r'Water Activity, $a_\mathrm{H_2O(l)}$')
plt.plot(I, aH2O_ideal, label='Ideal Model')
plt.plot(I, aH2O_davies, label='Davies Model')
plt.legend(loc='upper right')
plt.savefig('activity-water-davies-vs-ideal.pdf')

# Plot the difference in activity coeffitients of Davis and ideal model as a function an ionic strength
plt.figure()
plt.xlabel('Ionic Strength [molal]')
plt.title('Relative Difference Between Ideal and \n Davies Models for Water Activity')
plt.plot(I, numpy.abs(aH2O_davies-aH2O_ideal)/aH2O_davies)
plt.savefig('activity-water-davies-vs-ideal-diff.pdf')
