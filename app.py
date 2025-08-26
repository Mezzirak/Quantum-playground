import streamlit as st
import numpy as np
from solver import solve_schrodinger
from potentials import harmonic_oscillator, particle_in_box, double_well
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import io

st.sidebar.title("Quantum Playground")

# Select potential
potential_type = st.sidebar.selectbox(
    "Choose potential",
    ["Harmonic Oscillator", "Particle in a Box", "Double Well", "Custom Expression", "Custom CSV"]
)

# Spatial grid
x = np.linspace(-10, 10, 1000)
V = np.zeros_like(x)

# --- Potential parameters ---
if potential_type == "Harmonic Oscillator":
    k = st.sidebar.slider("Spring constant k", 0.1, 5.0, 1.0)
    V = harmonic_oscillator(x, k)
elif potential_type == "Particle in a Box":
    L = st.sidebar.slider("Box width L", 1.0, 20.0, 10.0)
    V = particle_in_box(x, L)
elif potential_type == "Double Well":
    a = st.sidebar.slider("Well separation a", 1.0, 10.0, 5.0)
    b = st.sidebar.slider("Barrier height b", 0.1, 2.0, 0.5)
    V = double_well(x, a, b)
elif potential_type == "Custom Expression":
    st.sidebar.write("Enter a Python expression for V(x) using `x` and `np` (numpy).")
    expr = st.sidebar.text_input("V(x) =", "0.5*x**2")
    try:
        V = eval(expr)
    except Exception as e:
        st.error(f"Error in custom potential: {e}")
        V = np.zeros_like(x)
elif potential_type == "Custom CSV":
    uploaded_file = st.sidebar.file_uploader("Upload CSV (x, V(x))", type=["csv"])
    if uploaded_file is not None:
        try:
            data = np.genfromtxt(io.StringIO(uploaded_file.getvalue().decode("utf-8")), delimiter=",", skip_header=0)
            x_csv = data[:,0]
            V_csv = data[:,1]
            f = interp1d(x_csv, V_csv, kind='cubic', fill_value="extrapolate")
            V = f(x)
        except Exception as e:
            st.error(f"Error reading CSV: {e}")
            V = np.zeros_like(x)

# --- Solve Schrödinger equation ---
x, energies, wavefuncs = solve_schrodinger(V, x)

# --- Plot ---
st.title(f"{potential_type} Eigenstates")
n_states = st.slider("Number of eigenstates to display", 1, 10, 5)

fig, ax = plt.subplots(figsize=(8,5))
for i in range(n_states):
    ax.plot(x, wavefuncs[:,i] + energies[i], label=f"n={i}")
ax.set_xlabel("x")
ax.set_ylabel("Energy + ψ(x)")
ax.legend()
st.pyplot(fig)

