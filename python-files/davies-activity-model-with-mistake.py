import numpy
def gammadavies(Zi, I):
    Agamma = 0.5095
    sqrtI = numpy.sqrt(I)
    return 10**( - Agamma * Zi * 2 * ( I / 1.0 + sqrtI - 0.3 * I ) )
print("gamma(1, 1.5) = ", gammadavies(1.0, 1.5))