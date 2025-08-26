from solver import solve_schrodinger
from potentials import harmonic_oscillator, particle_in_box, double_well
import visualisation as viz
import numpy as np

# Choose spatial grid
x = np.linspace(-10, 10, 1000)

# Choose potential
V = harmonic_oscillator(x, k=1.0)
# V = particle_in_box(x, L=10)
# V = double_well(x, a=5, b=0.5)

# Solve Schrödinger equation
x, energies, wavefuncs = solve_schrodinger(V, x)

# Plot first 5 eigenstates
viz.plot_eigenstates(x, energies, wavefuncs, n_states=5)
