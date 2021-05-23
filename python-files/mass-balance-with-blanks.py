# Import the Python package Numpy so we can perform linear algebra calculations
import numpy

# Create a list of species name
species = ['H2O(aq)', 'H+', 'OH-', 'CO3--', 'HCO3-', 'CO2(aq)', 'CO2(g)', 'H2O(g)']

# Create a list of rows of the formula matrix A
A = [[],
     [],
     [],
     []]

# Create a list with the amounts of species (in moles)
n = []

# Transform Python lists A and n into Numpy arrays
A = numpy.array(A)
n = numpy.array(n)

# Multiply matrix A and vector n to calculate the amounts of elements, b
b = A.dot(n)

# Create a list with the names of the elements
elements = ['H', 'O', 'C', 'Z']

# Create a list with the molar masses of H, O, C, Z 
molar_masses_elements = []

# Loop over all elements, their amounts, and their molar masses
for element, amount, molarmass in zip(elements, b, molar_masses_elements):
    print('Element %s has %f moles and %f grams' % (element, amount, amount*molarmass))

print("") # Just to skip one line in the output

# Alternative way of looping over all species using an index i
for i in range(len(species)):
    # Calculate the molar mass of current species
    molar_mass = A[:, i].dot(molar_masses_elements)
    # Calculate the mass of current species
    mass = 
    print('Species {} has {} moles and {} grams.'.format(species[i], n[i], mass))