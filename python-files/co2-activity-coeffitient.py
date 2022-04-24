import numpy
import matplotlib.pyplot as plt

def gamma_drummond_co2(T, I):
    c1 = -1.0312
    c2 = 1.2806e-3
    c3 = 255.9
    c4 = 0.4445
    c5 = -1.606e-3
    return numpy.exp((c1 + c2*T + c3/T)*I - (c4 + c5*T)*I/(1 + I))

I = numpy.linspace(0.0, 1.0, 101)

T1, T2, T3 = 50, 150, 250
gammaT1 = gamma_drummond_co2(T1 + 273.15, I)
gammaT2 = gamma_drummond_co2(T2 + 273.15, I)
gammaT3 = gamma_drummond_co2(T3 + 273.15, I)

plt.xlabel('Ionic Strength [molal]')
plt.ylabel('Activity Coefficient CO2(aq)')
plt.plot(I, gammaT1, label='T = %d C' % T1)
plt.plot(I, gammaT2, label='T = %d C' % T2)
plt.plot(I, gammaT3, label='T = %d C' % T3)
plt.legend(loc='best')
plt.savefig('activity-coefficient-co2-drummond.pdf')