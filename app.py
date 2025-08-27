import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

from solver import solve_1d
from solver_2d import solve_2d
from potentials import (
    harmonic_oscillator_1d, particle_in_box_1d, double_well_1d,
    harmonic_oscillator_2d, particle_in_box_2d, double_well_2d
)

# Sidebar UI
st.sidebar.title("Quantum Playground")
dimension = st.sidebar.selectbox("Select dimension", ["1D", "2D"])
potential_type = st.sidebar.selectbox(
    "Select potential",
    ["Harmonic oscillator", "Particle in a box", "Double well"]
)
custom_potential = st.sidebar.checkbox("Use a custom potential")
num_states = st.sidebar.slider("Number of eigenstates", 1, 10, 3)
grid_points = st.sidebar.slider("Grid points per dimension", 50, 200, 100)

# Construct grid and potential
if dimension == "1D":
    x = np.linspace(-5, 5, grid_points)
    if custom_potential:
        user_input = st.sidebar.text_area(
            "Enter potential function of x (use 'x' and numpy 'np')",
            "0.5 * x**2"
        )
        try:
            potential = eval(user_input, {"x": x, "np": np})
        except Exception as e:
            st.error(f"Error in custom potential: {e}")
            potential = np.zeros_like(x)
    else:
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
    if custom_potential:
        user_input = st.sidebar.text_area(
            "Enter potential function of X and Y (use 'X', 'Y' and numpy 'np')",
            "0.5 * (X**2 + Y**2)"
        )
        try:
            potential_2d = eval(user_input, {"X": X, "Y": Y, "np": np})
        except Exception as e:
            st.error(f"Error in custom potential: {e}")
            potential_2d = np.zeros_like(X)
    else:
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

# Solve and plot
if dimension == "1D":
    energies, wavefuncs = solve_1d(potential, x, n_eigen=num_states)
    st.subheader("Eigenstates and probability densities (1D)")
    for i in range(num_states):
        fig, ax = plt.subplots()
        ax.plot(x, wavefuncs[i], label=f'ψ{i+1}(x)')
        ax.plot(x, np.abs(wavefuncs[i])**2, label=f'|ψ{i+1}(x)|²')
        ax.set_title(f'Eigenstate {i+1}, Energy = {energies[i]:.3f}')
        ax.legend()
        st.pyplot(fig)

elif dimension == "2D":
    energies, wavefuncs = solve_2d(potential_2d, x, y, n_eigen=num_states)
    st.subheader("Eigenstates and probability densities (2D)")
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
