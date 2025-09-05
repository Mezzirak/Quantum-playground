import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.fft import fft, fftshift, fftfreq

# --- Import all your solvers and potentials ---
from solver import solve_1d
from solver_tdse import solve_time_1d
from solver_tdse_2d import solve_tdse_2d
from potentials import (
    harmonic_oscillator_1d, particle_in_box_1d, double_well_1d,
    finite_square_well_1d, potential_barrier_1d,
    harmonic_oscillator_2d
)

plt.style.use('seaborn-v0_8-darkgrid')
st.set_page_config(layout="wide")

# --- Caching Functions ---
@st.cache_data
def cached_solve_1d(potential_tuple, x_tuple, n_eigen):
    potential, x = np.array(potential_tuple), np.array(x_tuple)
    return solve_1d(potential, x, n_eigen)

@st.cache_data
def cached_solve_time_1d(potential_tuple, x_tuple, psi0_tuple, dt, n_steps):
    potential, x = np.array(potential_tuple), np.array(x_tuple)
    psi0 = np.array(psi0_tuple, dtype=complex)
    return solve_time_1d(potential, x, psi0, dt, n_steps)

@st.cache_data
def cached_solve_tdse_2d(_potential, _x, _y, _psi0, dt, n_steps):
    return solve_tdse_2d(_potential, _x, _y, _psi0, dt, n_steps)

# --- Sidebar UI ---
st.sidebar.title("Quantum Playground")
mode = st.sidebar.selectbox("Select mode", ["Stationary States", "Time Evolution"])
st.sidebar.markdown("---")

if mode == "Stationary States":
    # --- STATIONARY STATES UI (No changes here) ---
    st.sidebar.header("Stationary State Settings")
    dimension_stationary = st.sidebar.selectbox("Select dimension", ["1D"])
    st.header("1D Stationary Schrödinger Equation")
    potential_type_stationary = st.sidebar.selectbox("Select potential", ["Harmonic oscillator", "Particle in a box", "Double well", "Finite square well"], key="stationary_potential")
    num_states = st.sidebar.slider("Number of eigenstates", 1, 10, 3)
    grid_points_stationary = st.sidebar.slider("Grid points (Accuracy)", 50, 500, 200)
    x = np.linspace(-10, 10, grid_points_stationary)
    if potential_type_stationary == "Harmonic oscillator":
        omega = st.sidebar.slider("Harmonic frequency ω", 0.1, 5.0, 1.0); potential = harmonic_oscillator_1d(x, omega=omega)
    elif potential_type_stationary == "Particle in a box":
        L = st.sidebar.slider("Box length L", 1.0, 20.0, 10.0); potential = particle_in_box_1d(x, L=L)
    elif potential_type_stationary == "Double well":
        a = st.sidebar.slider("Quartic coefficient a", 0.1, 5.0, 1.0); b = st.sidebar.slider("Quadratic coefficient b", 0.1, 5.0, 2.0); potential = double_well_1d(x, a=a, b=b)
    elif potential_type_stationary == "Finite square well":
        V0 = st.sidebar.slider("Potential depth V₀", 1.0, 20.0, 10.0); L = st.sidebar.slider("Well width L", 1.0, 10.0, 5.0); potential = finite_square_well_1d(x, V0=V0, L=L)
    energies, wavefuncs = cached_solve_1d(tuple(potential), tuple(x), n_eigen=num_states)
    colors = cm.viridis(np.linspace(0.1, 0.9, num_states)); fig_main, ax_main = plt.subplots(figsize=(8, 6))
    ax_main.plot(x, potential, color='black', linestyle='--', linewidth=2, label='Potential V(x)')
    for i in range(num_states):
        wavefunc_max = np.max(np.abs(wavefuncs[i]));
        if wavefunc_max > 1e-9:
            scaled_wavefunc = wavefuncs[i] / wavefunc_max; energy_spacing = (energies[1] - energies[0]) if num_states > 1 and len(energies) > 1 else 1
            ax_main.plot(x, scaled_wavefunc * energy_spacing + energies[i], color=colors[i], linewidth=2.5, label=f'$E_{i+1}={energies[i]:.2f}$')
    ax_main.set_title("Eigenstates and Energies", fontsize=16); ax_main.set_xlabel("Position (x)", fontsize=12); ax_main.set_ylabel(r"Energy + $\psi(x)$", fontsize=12); ax_main.legend(); st.pyplot(fig_main)
    st.markdown("---"); st.subheader("Probability Densities"); cols = st.columns(num_states)
    for i, col in enumerate(cols):
        with col:
            fig_prob, ax_prob = plt.subplots(); prob_density = np.abs(wavefuncs[i])**2
            ax_prob.plot(x, prob_density, color=colors[i], linewidth=2); ax_prob.fill_between(x, prob_density, alpha=0.3, color=colors[i])
            ax_prob.set_title(f'State $n={i+1}$'); ax_prob.set_xlabel("x"); ax_prob.set_ylabel(r'$|\psi_{' + str(i+1) + r'}(x)|^2$'); ax_prob.set_yticks([]); st.pyplot(fig_prob)


elif mode == "Time Evolution":
    # --- TIME EVOLUTION UI ---
    st.sidebar.header("Time Evolution Settings")
    dimension_tdse = st.sidebar.selectbox("Select dimension", ["1D", "2D"])
    st.sidebar.warning("Higher 'Grid points' or 'Time steps' will significantly increase calculation time, especially for 2D.")

    if dimension_tdse == "1D":
        potential_type_tdse = st.sidebar.selectbox("Select potential", ["Potential Barrier (Tunneling)", "Free Particle", "Harmonic oscillator"])
        grid_points = st.sidebar.slider("Grid points (Accuracy)", 100, 1000, 250); time_steps = st.sidebar.slider("Time steps (Animation length)", 50, 500, 100)
        st.sidebar.subheader("Initial Wave Packet"); x0 = st.sidebar.slider("Initial position x₀", -5.0, 5.0, -3.0); p0 = st.sidebar.slider("Initial momentum p₀", -10.0, 10.0, 5.0); sigma0 = st.sidebar.slider("Initial width σ₀", 0.1, 2.0, 0.5)
        
        if st.sidebar.button("Calculate and Animate 1D"):
            with st.spinner(f"Calculating {time_steps} 1D steps..."):
                x = np.linspace(-20, 20, grid_points)
                if potential_type_tdse == "Harmonic oscillator": potential = harmonic_oscillator_1d(x, omega=1.0)
                elif potential_type_tdse == "Potential Barrier (Tunneling)": potential = potential_barrier_1d(x, V0=5.0, width=0.5, center=0.0)
                else: potential = np.zeros_like(x)
                psi0 = (1 / (2 * np.pi * sigma0**2))**(1/4) * np.exp(-(x - x0)**2 / (4 * sigma0**2)) * np.exp(1j * p0 * x); psi0 /= np.sqrt(np.sum(np.abs(psi0)**2) * (x[1]-x[0]))
                psi_t = cached_solve_time_1d(tuple(potential), tuple(x), tuple(psi0), dt=0.01, n_steps=time_steps)
            
            st.header("1D Time-Dependent Schrödinger Equation"); fig, (ax_main, ax_momentum) = plt.subplots(2, 1, figsize=(8, 8)); plot_placeholder = st.empty(); colors = cm.plasma(np.linspace(0, 1, time_steps))
            for i in range(time_steps):
                # ... (Animation code remains the same) ...
                ax_main.clear(); ax_momentum.clear(); prob_density = np.abs(psi_t[i])**2; ax_main.plot(x, prob_density, color=colors[i], linewidth=2); ax_main.fill_between(x, prob_density, alpha=0.3, color=colors[i]); potential_max = np.max(np.abs(potential)) if np.max(np.abs(potential)) > 0 else 1.0; prob_max = np.max(np.abs(psi_t[0])**2); ax_main.plot(x, potential / potential_max * prob_max, 'k--', label="V(x)"); ax_main.set_ylim(0, prob_max * 1.2); ax_main.set_title("Position Space"); ax_main.set_ylabel(r'$|\psi(x,t)|^2$'); ax_main.legend(loc="upper right")
                p = fftshift(fftfreq(grid_points, x[1]-x[0]))*2*np.pi; phi_p = fftshift(fft(psi_t[i])); ax_momentum.plot(p, np.abs(phi_p)**2, color=colors[i], linewidth=2); ax_momentum.fill_between(p, np.abs(phi_p)**2, alpha=0.3, color=colors[i]); ax_momentum.set_title("Momentum Space"); ax_momentum.set_ylabel(r'$|\phi(p,t)|^2$'); ax_momentum.set_xlim(-20, 20);
                fig.tight_layout(); plot_placeholder.pyplot(fig)

            # --- NEW: SCATTERING ANALYSIS ---
            if potential_type_tdse == "Potential Barrier (Tunneling)":
                st.subheader("Scattering Analysis")
                psi_final = psi_t[-1]
                prob_density_final = np.abs(psi_final)**2
                dx = x[1] - x[0]
                
                # Define reflected and transmitted regions (barrier is at x=0)
                reflected_mask = x < 0
                transmitted_mask = x >= 0
                
                # Calculate R and T by integrating the probability density
                R = np.sum(prob_density_final[reflected_mask]) * dx
                T = np.sum(prob_density_final[transmitted_mask]) * dx
                
                col1, col2, col3 = st.columns(3)
                col1.metric("Reflection Coefficient (R)", f"{R:.3f}")
                col2.metric("Transmission Coefficient (T)", f"{T:.3f}")
                col3.metric("Total Probability (R+T)", f"{R+T:.3f}")

                if not np.isclose(R+T, 1.0, atol=0.01):
                    st.warning("Note: R+T is not close to 1. This can happen if the simulation time is too short for the wave packet to fully scatter, or if it hits the simulation boundaries.")

        else: st.info("Set 1D parameters and click 'Calculate and Animate' to begin.")

    elif dimension_tdse == "2D":
        # ... (2D code remains the same) ...
        st.sidebar.subheader("2D Settings"); grid_points_2d = st.sidebar.slider("Grid points per side", 20, 100, 50); time_steps_2d = st.sidebar.slider("Time steps", 20, 100, 40)
        st.sidebar.subheader("Initial 2D Wave Packet"); x0_2d = st.sidebar.slider("Initial x₀", -5.0, 5.0, -2.0); y0_2d = st.sidebar.slider("Initial y₀", -5.0, 5.0, 0.0); px0_2d = st.sidebar.slider("Initial momentum pₓ", -5.0, 5.0, 2.0); py0_2d = st.sidebar.slider("Initial momentum pᵧ", -5.0, 5.0, 1.0)
        if st.sidebar.button("Calculate and Animate 2D"):
            with st.spinner(f"Calculating {time_steps_2d} 2D steps (this may take a while)..."):
                span = 8; x = np.linspace(-span, span, grid_points_2d); y = np.linspace(-span, span, grid_points_2d); X, Y = np.meshgrid(x, y)
                V = harmonic_oscillator_2d(X, Y, omega_x=0.1, omega_y=0.1)
                sigma = 0.5; psi0 = (np.exp(-((X-x0_2d)**2 + (Y-y0_2d)**2)/(4*sigma**2)) * np.exp(1j * (px0_2d*X + py0_2d*Y))); psi0 /= np.sqrt(np.sum(np.abs(psi0)**2))
                psi_t = cached_solve_tdse_2d(V, x, y, psi0, dt=0.02, n_steps=time_steps_2d)
            st.header("2D Time-Dependent Schrödinger Equation"); fig, ax = plt.subplots(figsize=(7, 6)); plot_placeholder = st.pyplot(fig); cmap = plt.get_cmap('magma') 
            for i in range(time_steps_2d):
                ax.clear(); prob_density_2d = np.abs(psi_t[i])**2; im = ax.imshow(prob_density_2d.T, cmap=cmap, interpolation='bilinear', extent=[-span, span, -span, span], origin='lower'); ax.set_title(f"2D Probability Density $|\psi(x, y, t)|^2$ at Time Step {i}"); ax.set_xlabel("x"); ax.set_ylabel("y"); plot_placeholder.pyplot(fig)
        else: st.info("Set 2D parameters and click 'Calculate and Animate' to begin.")