# Import array and plotting libraries
import numpy
import matplotlib.pyplot as plt

# Function calculating activity of CO2(aq) according to the Drummond model 
def gamma_drummond_co2(T, I):
    c1, c2, c3, c4, c5 = -1.0312, 1.2806e-3, 255.9, 0.4445, -1.606e-3
    return numpy.exp((c1 + c2*T + c3/T)*I - (c4 + c5*T)*I/(1 + I))

# Create an array with values of ionic strength in a range from 0.0 to 1.0
I = numpy.linspace(0.0, 1.0, 101)

# Define the temperatures for activity calculation
T1, T2, T3 = 50, 150, 250 # in celsius
# Calculate CO2(aq) activity for Drummond model for selected temperatures (in Kelvin)
gammaT1 = gamma_drummond_co2(T1 + 273.15, I)
gammaT2 = gamma_drummond_co2(T2 + 273.15, I)
gammaT3 = gamma_drummond_co2(T3 + 273.15, I)

# Plot CO2(aq) activity for Drummond model as a function an ionic strength for different T
plt.xlabel('Ionic Strength [molal]')
plt.ylabel('Activity Coefficient CO2(aq)')
plt.plot(I, gammaT1, label=r'T = %d$^\circ$C' % T1)
plt.plot(I, gammaT2, label=r'T = %d$^\circ$C' % T2)
plt.plot(I, gammaT3, label=r'T = %d$^\circ$C' % T3)
plt.legend(loc='best')
plt.savefig('activity-coefficient-co2-drummond.pdf', bbox_inches='tight')
