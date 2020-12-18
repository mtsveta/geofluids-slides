# Import the Python package Numpy so we can perform linear algebra calculations
import numpy

# Create a list of species name
species = ['H2O(aq)', 'H+', 'OH-', 'CO3--', 'HCO3-', 'CO2(aq)', 'CO2(g)', 'H2O(g)']

# Create a list of rows of the formula matrix A
A = [[2,  1,  1,  0,  1,  0,  0,  2],
     [1,  0,  1,  3,  3,  2,  2,  1],
     [0,  0,  0,  1,  1,  1,  1,  0],
     [0,  1, -1, -2, -1,  0,  0,  0]]

# Create a list with the amounts of species (in moles)
n = [55.4551, 1.23485e-4, 8.39739e-11, 4.93648e-11, 1.23484e-4, 0.032861, 1.96702, 0.0531732]

# Transform Python lists A and n into Numpy arrays
A = numpy.array(A)
n = numpy.array(n)

# Multiply matrix A and vector n to calculate the amounts of elements, b
b = A.dot(n)

# Create a list with the names of the elements
elements = ['H', 'O', 'C', 'Z']

# Create a list with the molar masses of H, O, C, Z
molar_masses_elements = [1.0079, 15.9994, 12.0107, 0.0]

# Loop over all elements, their amounts, and their molar masses
for element, amount, molarmass in zip(elements, b, molar_masses_elements):
    print('Element %s has %f moles and %f grams' % (element, amount, amount*molarmass))

print("") # Just to skip one line in the output

# Alternative way of looping over all species using an index i
for i in range(len(species)):
    # Calculate the molar mass of current species
    molar_mass = A[:, i].dot(molar_masses_elements)
    # Calculate the mass of current species
    mass = n[i] * molar_mass
    print('Species {} has {} moles and {} grams.'.format(species[i], n[i], mass))