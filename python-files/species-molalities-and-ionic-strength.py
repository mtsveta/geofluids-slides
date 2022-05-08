import numpy

# Create a list with names of aqueous species of interest
aqueous_species = ['H2O(aq)', 'H+', 'OH-', 'Na+', 'Cl-', 'CO3--', 'HCO3-', 'CO2(aq)']

# The index of species H2O(aq) used in the calculation of molalities
iH2O = 0 

# Create a list with the electrical charges of each species
Z = [0,  1, -1, 1, -1, -2, -1,  0]

# Create a list with the amounts (in moles) of each species
n = [55.4551, 1.23485e-4, 8.39739e-11, 0.92, 0.92, 4.93648e-11, 1.23484e-4, 0.032861]

# Transform the list n into a Numpy array for numerical calculations
n = numpy.array(n)

# Calculate the molalities of all species
m = 55.508 * n / n[iH2O] #  = n / (18.01528 * 1e-3 * n[iH2O])

# Calculate the ionic strength of the aqueous solution
I = 0.5 * sum(m * Z * Z)

# Print the calculated ionic strength
print('Ionic strength of solution is %f molal.' % I)

# Print the calculated molality of each species
for i in range(len(aqueous_species)):
    print(f'Molality of {aqueous_species[i]:8s} is {m[i]:4.4e} molal.')