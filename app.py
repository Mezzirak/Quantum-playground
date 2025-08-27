import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

# Import solvers and potentials
from solver import solve_schrodinger
from solver_2d import solve_schrodinger_2d
from potentials import (
    harmonic_oscillator_1d, particle_in_box_1d, double_well_1d,
    harmonic_oscillator_2d, particle_in_box_2d, double_well_2d
)

# Streamlit user interface
st.sidebar.title("Quantum Playground")
dimension = st.sidebar.selectbox("Select dimension", ["1D", "2D"])

# Select potential type
if dimension == "1D":
    potential_type = st.sidebar.selectbox("Select potential", ["Harmonic oscillator", "Particle in a box", "Double well"])
elif dimension == "2D":
    potential_type = st.sidebar.selectbox("Select potential", ["Harmonic oscillator", "Particle in a box", "Double well"])

# Number of eigenstates
num_states = st.sidebar.slider("Number of eigenstates", 1, 10, 3)

# Grid points
grid_points = st.sidebar.slider("Grid points per dimension", 50, 200, 100)

# Construct grid and potential
if dimension == "1D":
    x = np.linspace(-5, 5, grid_points)
    if potential_type == "Harmonic oscillator":
        omega = st.sidebar.slider("Harmonic frequency ω", 0.1, 5.0, 1.0)
        potential = harmonic_oscillator_1d(x, omega=omega)
    elif potential_type == "Particle in a box":
        L = st.sidebar.slider("Box length L", 1.0, 20.0, 10.0)
        potential = particle_in_box_1d(x, L=L)
    elif potential_type == "Double well":
        a = st.sidebar.slider("Quartic coefficient a", 0.1, 5.0, 1.0)
        b = st.sidebar.slider("Quadratic coefficient b", 0.1, 5.0, 1.0)
        potential = double_well_1d(x, a=a, b=b)

elif dimension == "2D":
    x = np.linspace(-5, 5, grid_points)
    y = np.linspace(-5, 5, grid_points)
    X, Y = np.meshgrid(x, y)

    if potential_type == "Harmonic oscillator":
        omega_x = st.sidebar.slider("Harmonic frequency ωx", 0.1, 5.0, 1.0)
        omega_y = st.sidebar.slider("Harmonic frequency ωy", 0.1, 5.0, 1.0)
        potential_2d = harmonic_oscillator_2d(X, Y, omega_x=omega_x, omega_y=omega_y)

    elif potential_type == "Particle in a box":
        Lx = st.sidebar.slider("Box length Lx", 1.0, 20.0, 10.0)
        Ly = st.sidebar.slider("Box length Ly", 1.0, 20.0, 10.0)
        potential_2d = particle_in_box_2d(X, Y, Lx=Lx, Ly=Ly)

    elif potential_type == "Double well":
        a = st.sidebar.slider("Quartic coefficient a", 0.1, 5.0, 1.0)
        b = st.sidebar.slider("Quadratic coefficient b", 0.1, 5.0, 1.0)
        potential_2d = double_well_2d(X, Y, a=a, b=b)

# Solve and plot 1D eigenstates
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

# Solve and plot 2D eigenstates with interactive Plotly
elif dimension == "2D":
    energies, wavefuncs = solve_schrodinger_2d(x, y, potential_2d, num_states=num_states)
    st.subheader("Eigenstates and Probability Densities (2D)")

    for i in range(num_states):
        fig = go.Figure(data=[go.Surface(
            z=np.abs(wavefuncs[i])**2,
            x=X,
            y=Y,
            colorscale='Viridis',
            colorbar=dict(title='|ψ(x,y)|²')
        )])
        fig.update_layout(
            title=f"Eigenstate {i+1}, Energy = {energies[i]:.3f}",
            scene=dict(
                xaxis_title='x',
                yaxis_title='y',
                zaxis_title='Probability density'
            ),
            autosize=False,
            width=700,
            height=600
        )
        st.plotly_chart(fig)



