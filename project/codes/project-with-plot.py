# Import the Python package Numpy so we can perform linear algebra calculations
from numpy import *
import matplotlib.pyplot as plt

plt.style.use('ggplot')
plt.rc('font', size=16)
plt.rc('font', family='serif', serif=['New Century Schoolbook'])
plt.rc('text', usetex=True)
plt.rc('lines', linewidth=2)
plt.rc('figure', autolayout=True)


# The universal gas constant in units of J/(mol*K)
R = 8.3144598

# The reference pressure P(ref) (in units of Pa)
Pref = 1.0e5

# Create a list of component names
# components = ['H', 'O', 'C', 'Na', 'Cl', 'Ca', 'Z']
components = ['H', 'O', 'C', 'Na', 'Cl', 'Ca']

# Create a list of species names
species = ['H2O(l)',
           'H+(aq)',
           'OH-(aq)',
           'HCO3-(aq)',
           'CO3--(aq)',
           'CO2(aq)',
           'Na+(aq)',
           'Cl-(aq)',
           'Ca++(aq)',
           'CO2(g)',
           'CaCO3(s,calcite)']

# Create a list with the charges of the aqueous species
charges = [0, 1, -1, -1, -2, 0, 1, -1, 2]

# The number of components and species
num_components = len(components)
num_species = len(species)

# The number of aqueous, gaseous, and mineral species
num_species_aqueous = 9
num_species_gaseous = 1
num_species_mineral = 1

slice_aqueous = slice(0, num_species_aqueous)
slice_gaseous = slice(num_species_aqueous, num_species_aqueous + num_species_gaseous)
slice_mineral = slice(num_species_aqueous + num_species_gaseous, -1)

# Store the index of H2O(l)
iH2O = species.index('H2O(l)')

# Store the index of CO2(aq)
iCO2 = species.index('CO2(aq)')

# Create a list with the indices of the charged species
icharged = [idx for idx, charge in enumerate(charges) if charge != 0.0]

# Construct the formula matrix of the chemical system
formula_matrix = [
    [2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],  # H
    [1, 0, 1, 3, 3, 2, 0, 0, 0, 2, 3],  # O
    [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1],  # C
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],  # Na
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],  # Cl
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1],  # Ca
    # [0, 1,-1,-1,-2, 0, 1,-1, 2, 0, 0]   # Z
]

# Define the function that calculates the activities of all aqueous species.
# Parameters:
#   - T is temperature in units of K
#   - P is pressure in units of Pa
#   - nphase is an array with the mole amounts of aqueous species only
# Return:
#   - an array with the ln activities of all aqueous species
def ln_activities_aqueous_species(T, P, nphase):

    # Calculate the molalities of the aqueous species and their natural log
    m = 55.508 * nphase / nphase[iH2O]

    # Calculate the ionic strength of the aqueous phase
    I = 0.5 * sum([mi*Zi**2 for mi, Zi in zip(m, charges)])

    # Calculate the mole fraction of water species, H2O(l)
    xH2O = nphase[iH2O] / sum(nphase)

    # Create an array to store the ln activity coeffs of the aqueous species
    ln_g = zeros(num_species_aqueous)

    ###########################################################################
    # Calculate the ln activity coefficient of CO2(aq)
    # ---- USE DRUMMOND (1981) MODEL ----
    ln_g[iCO2] = 0.0
    ###########################################################################

    ###########################################################################
    # Calculate the ln activity coefficient of each charged species
    # ---- USE DAVIES MODEL FOR EACH CHARGED SPECIES ----
    for i in icharged:
        ln_g[i] = 0.0
    ###########################################################################

    # Calculate the natural log of the molalities of the aqueous species
    ln_m = log(m)

    # Calculate the ln activities of the solute species
    ln_a = ln_g + ln_m

    # Calculate the ln activity of the solvent species H2O(l)
    ln_a[iH2O] = -(1 - xH2O)/xH2O

    # Return the array ln_a with the ln activities of the aqueous species
    return ln_a


# Define the function that calculates the activities of all gaseous species.
# Parameters:
#   - T is temperature in units of K
#   - P is pressure in units of Pa
#   - nphase is an array with the mole amounts of gaseous species only
# Return:
#   - an array with the ln activities of all gaseous species
def ln_activities_gaseous_species(T, P, nphase):
    # Calculate the array of mole fractions of the gaseous species
    x = nphase/sum(nphase)

    # Create an array to store the ln fugacity coeffs of the gaseous species
    ln_phi = zeros(num_species_gaseous)

    ###########################################################################
    # Calculate the ln fugacity coefficient of the only gaseous species, CO2(g)
    # ---- USE PENG-ROBINSON (1976) MODEL ----
    ln_phi[0] = 0.0 # note index 0, since there is only one gaseous species
    ###########################################################################

    # Calculate the array of ln mole fractions of the gaseous species
    ln_x = log(x)

    # Calculate the ln activities of the gaseous species
    ln_a = ln_phi + ln_x + log(P/Pref)

    # Return the array ln_a with the ln activities of the gaseous species
    return ln_a


# Define the function that calculates the activities of the mineral species.
# Parameters:
#   - T is temperature in units of K
#   - P is pressure in units of Pa
#   - nphase is an array with the mole amounts of the mineral species
# Return:
#   - an array with the ln activities of mineral species
def ln_activities_mineral_species(T, P, nphase):
    return zeros(num_species_mineral)


def ln_activities_aqueous_species_ddn(T, P, nphase):
    # The molar amount of H2O(l)
    nH2O = nphase[iH2O]

    # The mole fraction of H2O(l)
    xH2O = nH2O/sum(nphase)

    # Calculate the partial molar derivatives of the solute activities
    ddn = diag(1.0/nphase)
    ddn[:, iH2O] = -1.0/nH2O

    # Calculate the partial molar derivatives of the solvent activity
    ddn[iH2O, :] = -1.0/nH2O * xH2O/(xH2O - 1.0)
    ddn[iH2O, iH2O] = -1.0/nH2O

    return ddn


# Define the function that calculates the activities of all species.
# Parameters:
#   - T is temperature in units of K
#   - P is pressure in units of Pa
#   - n is an array with the mole amounts of all species in the chemical system
# Return:
#   - an array with the ln activities of all species in the chemical system
def ln_activities(T, P, n):
    # Create slices of n corresponding to each phase
    n_aqueous = n[slice_aqueous]
    n_gaseous = n[slice_gaseous]
    n_mineral = n[slice_mineral]

    # Calculate the activities of aqueous, gaseous and mineral species and
    # concatenate them into a single array with the ln activities of all
    # species in the chemical system.
    ln_a = zeros(num_species)
    ln_a[slice_aqueous] = ln_activities_aqueous_species(T, P, n_aqueous)
    ln_a[slice_gaseous] = ln_activities_gaseous_species(T, P, n_gaseous)
    ln_a[slice_mineral] = ln_activities_mineral_species(T, P, n_mineral)

    return ln_a


def ln_activities_ddn(T, P, n):
    # Create an array with the entries in n corresponding to aqueous species
    n_aqueous = n[slice_aqueous]

    # The matrix with partial molar derivatives of the activities
    ln_a_ddn = zeros((num_species, num_species))

    ln_a_ddn[slice_aqueous, slice_aqueous] = \
        ln_activities_aqueous_species_ddn(T, P, n_aqueous)

    return ln_a_ddn


def concentrations(T, P, n):
    c = ones(num_species)
    c[slice_aqueous] = 55.508 * n[slice_aqueous]/n[iH2O]
    c[iH2O] = n[iH2O]/sum(n[slice_aqueous])
    c[slice_gaseous] = n[slice_gaseous]/sum(n[slice_gaseous]) * P/Pref
    c[slice_mineral] = 1
    return c


# Create a dictionary to store the standard chemical potentials of all species
u0_table = {}

# The temperatures where the standard chemical potentials are available (in K)
u0_table['temperatures'] = array([25, 50, 75, 100]) + 273.15

# The pressures where the standard chemical potentials are available (in Pa)
u0_table['pressures'] = array([1, 50, 100]) * 1e5

# The standard chemical potentials of species H2O(l) (in J/mol)
u0_table['H2O(l)'] = array([
    [-237181.72, -237093.28, -237003.23],
    [-239006.71, -238917.46, -238826.59],
    [-240977.57, -240887.12, -240795.03],
    [-999999999, -242992.18, -242898.53],
])

# The standard chemical potentials of species H+(aq) (in J/mol)
u0_table['H+(aq)'] = array([
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
])

# The standard chemical potentials of species OH-(aq) (in J/mol)
u0_table['OH-(aq)'] = array([
    [-157297.48, -157320.08, -157342.16],
    [-156905.20, -156924.41, -156943.13],
    [-156307.81, -156328.46, -156348.60],
    [-999999999, -155558.52, -155583.77],
])

# The standard chemical potentials of species HCO3-(aq) (in J/mol)
u0_table['HCO3-(aq)'] = array([
    [-586939.89, -586820.91, -586698.78],
    [-589371.68, -589247.85, -589120.88],
    [-591758.94, -591634.79, -591507.48],
    [-999999999, -593983.38, -593858.95],
])

# The standard chemical potentials of species CO3--(aq) (in J/mol)
u0_table['CO3--(aq)'] = array([
    [-527983.14, -528011.90, -528039.31],
    [-526461.60, -526494.20, -526525.56],
    [-524473.15, -524514.70, -524555.01],
    [-999999999, -522121.51, -522175.91],
])

# The standard chemical potentials of species CO2(aq) (in J/mol)
u0_table['CO2(aq)'] = array([
    [-385974.00, -385813.47, -385650.13],
    [-389147.04, -388979.61, -388809.39],
    [-392728.17, -392556.66, -392382.38],
    [-999999999, -396484.88, -396307.90],
])

# The standard chemical potentials of species Na+(aq) (in J/mol)
u0_table['Na+(aq)'] = array([
    [-261880.74, -261886.18, -261890.73],
    [-263384.24, -263384.84, -263384.60],
    [-264979.92, -264978.35, -264975.97],
    [-999999999, -266664.71, -266661.71],
])

# The standard chemical potentials of species Cl-(aq) (in J/mol)
u0_table['Cl-(aq)'] = array([
    [-131289.74, -131204.66, -131117.63],
    [-132590.59, -132504.28, -132416.09],
    [-133682.46, -133598.32, -133512.30],
    [-999999999, -134501.87, -134420.82],
])

# The standard chemical potentials of species Ca++(aq) (in J/mol)
u0_table['Ca++(aq)'] = array([
    [-552790.08, -552879.51, -552968.88],
    [-551348.98, -551437.62, -551526.28],
    [-549855.33, -549945.95, -550036.62],
    [-999999999, -548400.71, -548495.68],
])

# The standard chemical potentials of species CO2(g) (in J/mol)
u0_table['CO2(g)'] = array([
    [-394358.74, -394358.74, -394358.74],
    [-399740.71, -399740.71, -399740.71],
    [-405197.73, -405197.73, -405197.73],
    [-410726.87, -410726.87, -410726.87],
])

# The standard chemical potentials of species CaCO3(s,calcite) (in J/mol)
u0_table['CaCO3(s,calcite)'] = array([
    [-1129177.92, -1128996.94, -1128812.26],
    [-1131580.05, -1131399.06, -1131214.39],
    [-1134149.89, -1133968.91, -1133784.23],
    [-1136882.60, -1136701.62, -1136516.94],
])


# Define the interpolation function.
# Parameters:
#   - x is the x-coordinate where the interpolation is performed
#   - y is the y-coordinate where the interpolation is performed
#   - xpoints is an array with the x-coordinates where fvalues are available
#   - ypoints is an array with the y-coordinates where fvalues are available
#   - fvalues is a 2D array with the table of values of a function f
def interpolate(x, y, xpoints, ypoints, fvalues):
    ###########################################################################
    # REPLACE THE LINE BELOW WITH A PROPER INTERPOLATION CALCULATION
    ###########################################################################
    return fvalues[0, 0]  # This only gets the first value in fvalues


# Define the function that calculates the standard chemical potentials of all
# species at given temperature and pressure.
# Parameters:
#   - T is temperature in units of K
#   - P is pressure in units of Pa
# Return:
#   - an array with the standard chemical potentials of all species (in J/mol)
def standard_chemical_potentials(T, P):
    # Create an array to store the standard chemical potentials of all species
    u0 = zeros(num_species)

    # Alias to the temperature coordinates of the interpolation tables
    ts = u0_table['temperatures']

    # Alias to the pressure coordinates of the interpolation tables
    ps = u0_table['pressures']

    # Loop over all indices and species names
    for i, name in enumerate(species):
        # Interpolate the standard chemical potential of species i at (T, P)
        u0[i] = interpolate(T, P, ts, ps, u0_table[name])

    # Return the array with standard chemical potentials of the species
    return u0


# Define the function that calculates the chemical potentials of all species
# at given temperature, pressure, and mole amounts of species.
# Parameters:
#   - T is temperature in units of K
#   - P is pressure in units of Pa
#   - n is an array with the mole amounts of the species
# Return:
#   - an array with the chemical potentials of all species (in J/mol)
def chemical_potentials(T, P, n):
    # Calculate the standard chemical potentials of all species at (T, P)
    u0 = standard_chemical_potentials(T, P)

    # Calculate the activities of all species at (T, P, n)
    ln_a = ln_activities(T, P, n)

    # Return the chemical potentials of all species at (T, P, n)
    return u0 + R*T*ln_a


def chemical_potentials_ddn(T, P, n):
    return R*T*ln_activities_ddn(T, P, n)


def gibbs_energy(T, P, n):
    RT = R*T
    u = chemical_potentials(T, P, n)/RT
    H = chemical_potentials_ddn(T, P, n)/RT
    G = n.dot(u)/RT
    return G, u, H


def gibbs_energy_pure_phases(T, P, n):
    RT = R*T
    u = standard_chemical_potentials(T, P)/RT
    H = zeros((num_species, num_species))
    G = n.dot(u)/RT
    return G, u, H


class ObjectiveResult:
    def __init__(self):
        self.f = None  # the value of the objective function at x
        self.g = None  # the gradient of the objective function at x
        self.H = None  # the Hessian of the objective function at x


class OptimumProblem:
    def __init__(self):
        self.objective = None  # the objective function
        self.A = None  # the matrix A in the equality constraints
        self.b = None  # the vector b in the equality constraints


class OptimumState:
    def __init__(self):
        self.x = None  # the variables x = (x1, ..., xn)
        self.y = None  # the Lagrange multipliers y = (y1, ..., ym)
        self.z = None  # the Lagrange multipliers z = (z1, ..., zn)


def minimize(state, problem, **options):
    """Minimize a function f(x) subject to constraints."""

    imax = options.get('imax', 100)
    mu = options.get('mu', 1.0e-14)
    tau = options.get('tau', 0.99999)
    tol = options.get('tol', 1.0e-6)

    A = problem.A
    b = problem.b

    m = A.shape[0]
    n = A.shape[1]
    p = m + n
    t = m + 2 * n

    if state.x is None:
        state.x = full(n, mu)

    if state.y is None:
        state.y = zeros(m)

    if state.z is None:
        state.z = ones(n)

    x = state.x
    y = state.y
    z = state.z

    F = zeros(t)
    J = zeros((t, t))

    f = 0.0
    g = zeros(n)  # F[:n]
    H = zeros((n, n))  # J[:n, :n]

    for it in xrange(imax):
        f, g, H = problem.objective(x)

        # Assemble the negative of the residual vector -F
        #     [g(x) - tr(A)*y - z]
        # F = [      A*x - b     ]
        #     [    X*Z*e - mu    ]
        F[:n] = g - A.T.dot(y) - z
        F[n:p] = A.dot(x) - b
        F[-n:] = x * z - mu

        # Calculate the optimality, feasibility, complementarity errors
        error_opt = linalg.norm(F[:n], inf)
        error_fea = linalg.norm(F[n:p], inf)
        error_com = linalg.norm(F[-n:], inf)

        # Calculate the current total error
        error = max([error_opt, error_fea, error_com])

        # Check if the calculation has converged
        if error < tol:
            break

        # Assemble the Jacoabian matrix J
        #     [H -tr(A) -I]
        # J = [A    0    0]
        #     [Z    0    X]
        J[:n, :n] = H
        J[:n, n:p] = -A.T
        J[:n, -n:] = -eye(n)
        J[n:p, :n] = A
        J[-n:, :n] = diag(z)
        J[-n:, -n:] = diag(x)

        # Compute the Newton step d = [dx dy dx]
        delta = linalg.solve(J, -F)

        dx = delta[:n]
        dy = delta[n:p]
        dz = delta[-n:]

        # Calculate the new values for x and z
        for i in xrange(n):
            x[i] += dx[i] if x[i] + dx[i] > 0.0 else -tau * x[i]
            z[i] += dz[i] if z[i] + dz[i] > 0.0 else -tau * z[i]

        # Calculate the new values for y
        y += dy

    return (it, it < imax, error)


def component_amounts(kgH2O, molCO2, molNaCl, molCaCO3):
    molH2O = 55.508 * kgH2O
    b = [2*molH2O,                        # H
         molH2O + 2*molCO2 + 3*molCaCO3,  # O
         molCO2 + molCaCO3,               # C
         molNaCl,                         # Na
         molNaCl,                         # Cl
         molCaCO3]                        # Ca
    return b


class EquilibriumState:
    def __init__(self):
        self.T = None  # temperature in K
        self.P = None  # pressure in Pa
        self.n = None  # amounts of species in mol
        self.y = None  # Lagrange multipliers y in J/mol
        self.z = None  # Lagrange multipliers z in J/mol
        self.c = None  # concentrations of the species
        self.a = None  # activities of the species
        self.g = None  # activity coefficients of the species

    def __str__(self):
        bar = '='*125+'\n'
        thinbar = '-'*125+'\n'
        res = bar
        res += '{:25} {:25} {:25} {:25}\n'.format(
            'Temperature [K]', 'Temperature [C]',
            'Pressure [Pa]', 'Pressure [bar]')
        res += thinbar
        res += '{:25} {:25} {:25} {:25}\n'.format(
            str(self.T), str(self.T - 273.15),
            str(self.P), str(self.P * 1e-5))
        res += bar
        res += '{:25} {:25} {:25} {:25} {:25}\n'.format(
            'Species', 'Amount [mol]', 'Concentration',
            'Activity', 'Activity Coefficient')
        res += thinbar
        for i, name in enumerate(species):
            res += '{:25} {:25} {:25} {:25} {:25}\n'.format(
                name, str(state.n[i]), str(state.c[i]),
                str(state.a[i]), str(state.g[i]))
        res += bar
        return res


# Calculate the equilibrium state of the chemical system.
# Parameters:
#   - T is temperature in K
#   - P is pressure in Pa
#   - b is an array with molar amounts of components
# Return:
#   - an state object with members n, y, z, so that the
#     molar amounts of the species is given by state.n,
#     and the Lagrange multipliers by state.y and state.z.
def equilibrate(T, P, b):
    # The state of the variables (x, y, z)
    state = OptimumState()

    # Define the minimization problem
    problem = OptimumProblem()
    problem.A = array(formula_matrix)
    problem.b = array(b)

    # Define the objective function to calculate an initial guess
    problem.objective = lambda x: gibbs_energy_pure_phases(T, P, x)

    # Minimize the Gibbs energy assuming each species is a pure phase
    minimize(state, problem)

    # Define now the objective function with a normalized Gibbs energy function
    problem.objective = lambda x: gibbs_energy(T, P, x)

    # Minimize the Gibbs energy of the chemical system
    minimize(state, problem)

    # Create an alias for the x member that resembles chemical notation
    equilibriumstate = EquilibriumState()
    equilibriumstate.T = T
    equilibriumstate.P = P
    equilibriumstate.n = state.x
    equilibriumstate.y = state.y * R * T
    equilibriumstate.z = state.z * R * T
    equilibriumstate.a = exp(ln_activities(T, P, state.x))
    equilibriumstate.c = concentrations(T, P, state.x)
    equilibriumstate.g = equilibriumstate.a/equilibriumstate.c

    return equilibriumstate


T = 298.15  # temperature input in K
P = 1.0e5   # pressure input in Pa

b = component_amounts(kgH2O=1, molCO2=4, molNaCl=0.2, molCaCO3=0.01)

state = equilibrate(T, P, b)

output = open('result.txt', 'w')

print >>output, state

##############################################################################
# The code below is the example given in the project description document.
##############################################################################
T = 333.15  # temperature in K (=60 C)
P = 100e5  # pressure in Pa (=100 bar)
# Calculate the amounts of elements based on a mixture of
# 1 kg of H2O, 1 mol of CO2, 0.5 moles of NaCl, and 5 moles of CaCO3.
b = component_amounts(kgH2O=1, molCO2=1, molNaCl=0.5, molCaCO3=5)
state = equilibrate(T, P, b)  # perform the equilibrium calculation
iH = species.index('H+(aq)')  # the index of species H+(aq)
pH = -log10(state.a[iH])  # calculate the pH of the aqueous solution
print 'The pH in the calculated equilibrium state is', pH


def calculate_pH(molCO2):
    b = component_amounts(kgH2O=1, molCO2=molCO2, molNaCl=0.3, molCaCO3=1e-6)
    state = equilibrate(T, P, b)  # perform the equilibrium calculation
    iH = species.index('H+(aq)')  # the index of species H+(aq)
    pH = -log10(state.a[iH])  # calculate the pH of the aqueous solution
    return pH

print 'The pH in the calculated equilibrium state is', calculate_pH(1.0)


molCO2vals = linspace(1e-6, 2, 100)
pHvals = [calculate_pH(x) for x in molCO2vals]

plt.xlabel('Amount CO2 [mol]')
plt.ylabel('pH')
plt.plot(molCO2vals, pHvals)
# plt.savefig('../figures/activity-coefficient-co2-drummond.pdf')
plt.show()


