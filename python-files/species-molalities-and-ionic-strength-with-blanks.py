import numpy

# Create a list with names of aqueous species of interest
aqueous_species = ['H2O(aq)', 'H+', 'OH-', 'Na+', 'Cl-', 'CO3--', 'HCO3-', 'CO2(aq)']

# The index of species H2O(aq) used in the calculation of molalities
iH2O = 0 

# Create a list with the electrical charges of each species
Z = []

# Create a list with the amounts (in moles) of each species
n = []

# Transform the list n into a Numpy array for numerical calculations
n = numpy.array(n)

# Calculate the molalities of all species
m = 

# Calculate the ionic strength of the aqueous solution
I = 

# Print the calculated ionic strength
print('Ionic strength of solution is %f molal.' % I)

# Print the calculated molality of each species
for i in range(len(aqueous_species)):
    print('Molality of %s is %f molal.' % (aqueous_species[i], m[i]))