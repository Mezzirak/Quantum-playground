import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# Import your solvers and potentials
from solver import solve_schrodinger
from solver_2d import solve_schrodinger_2d
from potentials import harmonic_oscillator_1d, particle_in_box_1d, double_well_1d
from potentials import harmonic_oscillator_2d, particle_in_box_2d, double_well_2d

# Streamlit sidebar for user input
st.sidebar.title("Quantum Playground")
dimension = st.sidebar.selectbox("Select dimension", ["1D", "2D"])

# Select potential type
if dimension == "1D":
    potential_type = st.sidebar.selectbox("Select potential", ["Harmonic oscillator", "Particle in a box", "Double well"])
elif dimension == "2D":
    potential_type = st.sidebar.selectbox("Select potential", ["Harmonic oscillator", "Particle in a box", "Double well"])

# Number of eigenstates to compute
num_states = st.sidebar.slider("Number of eigenstates", min_value=1, max_value=10, value=3)

# Grid resolution
grid_points = st.sidebar.slider("Grid points per dimension", min_value=50, max_value=200, value=100)

# Construct grid and potential
if dimension == "1D":
    x = np.linspace(-5, 5, grid_points)
    if potential_type == "Harmonic oscillator":
        potential = harmonic_oscillator_1d(x)
    elif potential_type == "Particle in a box":
        potential = particle_in_box_1d(x)
    elif potential_type == "Double well":
        potential = double_well_1d(x)
elif dimension == "2D":
    x = np.linspace(-5, 5, grid_points)
    y = np.linspace(-5, 5, grid_points)
    X, Y = np.meshgrid(x, y)
    if potential_type == "Harmonic oscillator":
        potential_2d = harmonic_oscillator_2d(X, Y)
    elif potential_type == "Particle in a box":
        potential_2d = particle_in_box_2d(X, Y)
    elif potential_type == "Double well":
        potential_2d = double_well_2d(X, Y)

# Solve and plot
if dimension == "1D":
    energies, wavefuncs = solve_schrodinger(x, potential, num_states=num_states)
    st.subheader("Eigenstates and Probability Densities (1D)")

    for i in range(num_states):
        fig, ax = plt.subplots()
        ax.plot(x, wavefuncs[:, i], label=f'ψ{i+1}(x)')
        ax.plot(x, np.abs(wavefuncs[:, i])**2, label=f'|ψ{i+1}(x)|²')
        ax.set_title(f'Eigenstate {i+1}, Energy = {energies[i]:.3f}')
        ax.legend()
        st.pyplot(fig)

elif dimension == "2D":
    energies, wavefuncs = solve_schrodinger_2d(x, y, potential_2d, num_states=num_states)
    st.subheader("Eigenstates and Probability Densities (2D)")

    for i in range(num_states):
        fig, ax = plt.subplots()
        c = ax.imshow(np.abs(wavefuncs[i])**2,
                      extent=[x.min(), x.max(), y.min(), y.max()],
                      origin='lower',
                      cmap='viridis')
        fig.colorbar(c, ax=ax, label='Probability density |ψ(x,y)|²')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(f'Eigenstate {i+1}, Energy = {energies[i]:.3f}')
        st.pyplot(fig)


