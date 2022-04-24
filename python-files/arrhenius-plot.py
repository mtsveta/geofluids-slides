# Import array and plotting libraries
import numpy
import math
import matplotlib.pyplot as plt

# Do not show any warnings
import matplotlib as mpl
mpl.rc_file("matplotlibrc")

# Set auxiliary constants for methan CH4
A = 2.21E6
lnA = numpy.log(A)
EA = 3.55E3
n = 2

# Function calculating ln(k) in Arrhenius equation
def arrhenius_eq(invT):
    return lnA - EA*invT

# Function calculating ln(k) in Arrhenius equation
def extended_arrhenius_eq(invT):
    return lnA - n*numpy.log(invT) - EA*invT

# Create an array with values of temperature between 220 K and 320 K
T = numpy.linspace(220.0, 320.0, 101)
invT = 1/T

# Calculate H2O(aq) activity for Davies and ideal model
lnk = arrhenius_eq(invT)

# Plot activity of Davies vs ideal model as a function an ionic strength
plt.figure()
plt.xlabel(r'1/T [1/K]')
plt.ylabel(r'$\ln k$')
plt.plot(invT, lnk, label='Arrhenius plot')
plt.legend(loc='best')
plt.savefig('arrhenius-plot.pdf', bbox_inches='tight')


# Create an array with values of temperature between 220 K and 320 K
T = numpy.linspace(300.0, 2220.0, 201)
invT = 1/T
ln_ext = extended_arrhenius_eq(invT)

# Plot the difference in activity of Davies and ideal model as a function an ionic strength
plt.figure()
plt.xlabel(r'1/T [1/K]')
plt.ylabel(r'$\ln k$')
#plt.title('Relative Difference Between Ideal and \n Davies Models for Water Activity')
plt.plot(invT, ln_ext, label='Extended Arrhenius plot')
plt.legend(loc='best')
plt.savefig('extended-arrhenius-plot.pdf', bbox_inches='tight')
