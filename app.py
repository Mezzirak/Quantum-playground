import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plotly.graph_objects as go
from scipy.fft import fft, fftshift, fftfreq
from scipy.integrate import solve_ivp
from numpy import trapz
import imageio
import io

# Import all solvers and potentials
from solver import solve_1d
from solver_2d import solve_2d
from solver_tdse import solve_time_1d
from solver_tdse_2d import solve_tdse_2d
from potentials import (
    harmonic_oscillator_1d, particle_in_box_1d, double_well_1d,
    finite_square_well_1d, potential_barrier_1d,
    harmonic_oscillator_2d, particle_in_box_2d, double_well_2d
)

plt.style.use('seaborn-v0_8-darkgrid')
st.set_page_config(layout="wide")

# Caching Functions
@st.cache_data
def cached_solve_1d(potential_tuple, x_tuple, n_eigen):
    potential, x = np.array(potential_tuple), np.array(x_tuple)
    return solve_1d(potential, x, n_eigen)

@st.cache_data
def cached_solve_2d(_potential, _x, _y, n_eigen):
    return solve_2d(_potential, _x, _y, n_eigen)

@st.cache_data
def cached_solve_time_1d(potential_tuple, x_tuple, psi0_tuple, dt, n_steps):
    potential, x = np.array(potential_tuple), np.array(x_tuple)
    psi0 = np.array(psi0_tuple, dtype=complex)
    return solve_time_1d(potential, x, psi0, dt, n_steps)

@st.cache_data
def cached_solve_tdse_2d(_potential, _x, _y, _psi0, dt, n_steps):
    return solve_tdse_2d(_potential, _x, _y, _psi0, dt, n_steps)

# Helper function for expectation values and uncertainties
def calculate_expectation_values(psi, x):
    """Calculates expectation values for position and momentum, and their uncertainties"""
    # Ensure wavefunction is normalised for this calculation
    dx = x[1] - x[0]
    norm = np.sqrt(np.sum(np.abs(psi)**2) * dx)
    psi_norm = psi / norm

    # Position expectation value <x>
    exp_x = np.sum(np.conj(psi_norm) * x * psi_norm).real * dx
    
    # Position squared expectation value <x^2>
    exp_x2 = np.sum(np.conj(psi_norm) * (x**2) * psi_norm).real * dx
    
    # Momentum operator application using np gradient for the derivative
    # p_hat * psi = -i * hbar * d/dx(psi) We use hbar=1
    p_psi = -1j * np.gradient(psi_norm, dx)
    
    # Momentum expectation value <p>
    exp_p = np.sum(np.conj(psi_norm) * p_psi).real * dx

    # Momentum squared expectation value <p^2>
    # p_hat^2 * psi = (-i)^2 * d^2/dx^2(psi) = -1 * d/dx(p_psi / -i)
    p2_psi = -1 * np.gradient(np.gradient(psi_norm, dx), dx)
    exp_p2 = np.sum(np.conj(psi_norm) * p2_psi).real * dx
    
    # Calculate uncertainties (variances first)
    var_x = exp_x2 - exp_x**2
    var_p = exp_p2 - exp_p**2

    # Avoid numerical issues with very small negative numbers
    delta_x = np.sqrt(max(0, var_x))
    delta_p = np.sqrt(max(0, var_p))
    
    return exp_x, exp_p, delta_x, delta_p

#Solver for Ehrenfest Theorem
def solve_classical(potential_func, x0, p0, t_span, t_eval):
    """Solves Newton's equation of motion for a classical particle"""
    def classical_equations(t, y):
        x, p = y
        V_plus = potential_func(x + 1e-6); V_minus = potential_func(x - 1e-6)
        force = -(V_plus - V_minus) / (2e-6)
        dxdt = p; dpdt = force
        return [dxdt, dpdt]
    sol = solve_ivp(classical_equations, t_span, [x0, p0], t_eval=t_eval, dense_output=True)
    return sol.y[0]

# Sidebar UI
st.sidebar.title("Quantum Playground")
mode = st.sidebar.selectbox("Select mode", ["Stationary States", "Time Evolution"])
st.sidebar.markdown("***")

if mode == "Stationary States":
    st.sidebar.header("Stationary State Settings")
    dimension_stationary = st.sidebar.selectbox("Select dimension", ["1D", "2D"])

    if dimension_stationary == "1D":
        st.header("1D Stationary Schrödinger Equation")
        potential_type_stationary = st.sidebar.selectbox("Select potential", ["Harmonic oscillator", "Particle in a box", "Double well", "Finite square well"], key="stationary_potential_1d")
        num_states = st.sidebar.slider("Number of eigenstates", 1, 10, 3); grid_points_stationary = st.sidebar.slider("Grid points (Accuracy)", 50, 500, 200)
        x = np.linspace(-10, 10, grid_points_stationary)
        if potential_type_stationary == "Harmonic oscillator": omega = st.sidebar.slider("Harmonic frequency ω", 0.1, 5.0, 1.0); potential = harmonic_oscillator_1d(x, omega=omega)
        elif potential_type_stationary == "Particle in a box": L = st.sidebar.slider("Box length L", 1.0, 20.0, 10.0); potential = particle_in_box_1d(x, L=L)
        elif potential_type_stationary == "Double well": a = st.sidebar.slider("Quartic coefficient a", 0.1, 5.0, 1.0); b = st.sidebar.slider("Quadratic coefficient b", 0.1, 5.0, 2.0); potential = double_well_1d(x, a=a, b=b)
        elif potential_type_stationary == "Finite square well": V0 = st.sidebar.slider("Potential depth V₀", 1.0, 20.0, 10.0); L = st.sidebar.slider("Well width L", 1.0, 10.0, 5.0); potential = finite_square_well_1d(x, V0=V0, L=L)
        energies, wavefuncs = cached_solve_1d(tuple(potential), tuple(x), n_eigen=num_states)
        colours = cm.viridis(np.linspace(0.1, 0.9, num_states)); fig_main, ax_main = plt.subplots(figsize=(8, 6))
        ax_main.plot(x, potential, color='black', linestyle='--', linewidth=2, label='Potential V(x)'); ax_main.text(x[5], potential[5] + 0.1, 'V(x)', fontsize=12, color='black')
        for i in range(num_states):
            wavefunc_max = np.max(np.abs(wavefuncs[i]));
            if wavefunc_max > 1e-9:
                scaled_wavefunc = wavefuncs[i] / wavefunc_max; energy_spacing = (energies[1] - energies[0]) if num_states > 1 and len(energies) > 1 else 1
                ax_main.plot(x, scaled_wavefunc * energy_spacing + energies[i], color=colours[i], linewidth=2.5, label=f'$E_{i+1}={energies[i]:.2f}$')
        ax_main.set_title("Eigenstates and Energies", fontsize=16); ax_main.set_xlabel("Position (x)", fontsize=12); ax_main.set_ylabel(r"Energy + $\psi(x)$", fontsize=12); ax_main.legend(); st.pyplot(fig_main)
        st.markdown("***"); st.subheader("Probability Densities"); cols = st.columns(num_states)
        for i, col in enumerate(cols):
            with col:
                fig_prob, ax_prob = plt.subplots(); prob_density = np.abs(wavefuncs[i])**2
                ax_prob.plot(x, prob_density, color=colours[i], linewidth=2); ax_prob.fill_between(x, prob_density, alpha=0.3, color=colours[i])
                ax_prob.set_title(f'State $n={i+1}$'); ax_prob.set_xlabel("x"); ax_prob.set_ylabel(r'$|\psi_{' + str(i+1) + r'}(x)|^2$'); ax_prob.set_yticks([]); st.pyplot(fig_prob)
        st.markdown("***"); st.subheader("Check for Eigenstate Orthogonality")
        st.write("A fundamental property of stationary states is that they are orthogonal.")
        state_options = list(range(num_states)); col1, col2 = st.columns(2)
        with col1: n_state = st.selectbox("Select first state (n)", state_options, index=0)
        with col2: m_state = st.selectbox("Select second state (m)", state_options, index=1 if num_states > 1 else 0)
        psi_n = wavefuncs[n_state]; psi_m = wavefuncs[m_state]; psi_n_norm = psi_n / np.sqrt(trapz(np.abs(psi_n)**2, x)); psi_m_norm = psi_m / np.sqrt(trapz(np.abs(psi_m)**2, x))
        integrand = np.conj(psi_n_norm) * psi_m_norm; overlap_integral = trapz(integrand, x); st.latex(r"\int \psi_n^*(x) \psi_m(x) dx = " + f"{np.real(overlap_integral):.4f}")
        if n_state == m_state: st.success("Result is ~1 (normalised).")
        else: st.success("Result is ~0 (orthogonal).")

    elif dimension_stationary == "2D":
        st.header("2D Stationary Schrödinger Equation")
        potential_type_2d = st.sidebar.selectbox("Select potential", ["Harmonic oscillator", "Particle in a box", "Double well"], key="stationary_potential_2d")
        num_states_2d = st.sidebar.slider("Number of eigenstates", 1, 6, 3)
        grid_points_2d = st.sidebar.slider("Grid points per side", 20, 100, 40)
        x = np.linspace(-5, 5, grid_points_2d); y = np.linspace(-5, 5, grid_points_2d); X, Y = np.meshgrid(x, y)
        if potential_type_2d == "Harmonic oscillator":
            omega_x = st.sidebar.slider("Harmonic frequency ωx", 0.1, 5.0, 1.0); omega_y = st.sidebar.slider("Harmonic frequency ωy", 0.1, 5.0, 1.5)
            potential_2d = harmonic_oscillator_2d(X, Y, omega_x=omega_x, omega_y=omega_y)
        elif potential_type_2d == "Particle in a box":
            Lx = st.sidebar.slider("Box length Lx", 1.0, 10.0, 8.0); Ly = st.sidebar.slider("Box length Ly", 1.0, 10.0, 8.0)
            potential_2d = particle_in_box_2d(X, Y, Lx=Lx, Ly=Ly)
        elif potential_type_2d == "Double well":
            a = st.sidebar.slider("Quartic coefficient a", 0.1, 5.0, 1.0); b = st.sidebar.slider("Quadratic coefficient b", 0.1, 5.0, 1.0)
            potential_2d = double_well_2d(X, Y, a=a, b=b)
        
        with st.spinner("Solving 2D Eigenstates..."):
            energies, wavefuncs = cached_solve_2d(potential_2d, x, y, n_eigen=num_states_2d)
        
        st.subheader("2D Probability Densities")
        for i in range(num_states_2d):
            prob_density_2d = np.abs(wavefuncs[i])**2
            fig = go.Figure(data=[go.Surface(z=prob_density_2d, x=X, y=Y, colorscale='Viridis', cmin=0, cmax=np.max(prob_density_2d))])
            fig.update_layout(title=f'Eigenstate n={i+1}, Energy = {energies[i]:.3f}',
                              scene=dict(xaxis_title='x', yaxis_title='y', zaxis_title='|ψ(x,y)|²'),
                              autosize=True, margin=dict(l=50, r=50, b=50, t=50))
            st.plotly_chart(fig, use_container_width=True)

elif mode == "Time Evolution":
    #Time Evolution code
    st.sidebar.header("Time Evolution Settings")
    dimension_tdse = st.sidebar.selectbox("Select dimension", ["1D", "2D"])
    st.sidebar.warning("Higher 'Grid points' or 'Time steps' will significantly increase calculation time, especially for 2D.")
    if dimension_tdse == "1D":
        potential_type_tdse = st.sidebar.selectbox("Select potential", ["Harmonic oscillator", "Potential Barrier (Tunnelling)", "Free Particle"])
        grid_points = st.sidebar.slider("Grid points (Accuracy)", 100, 1000, 250); time_steps = st.sidebar.slider("Time steps (Animation length)", 50, 500, 100)
        st.sidebar.subheader("Initial Wave Packet"); x0 = st.sidebar.slider("Initial position x₀", -5.0, 5.0, -3.0); p0 = st.sidebar.slider("Initial momentum p₀", -10.0, 10.0, 5.0); sigma0 = st.sidebar.slider("Initial width σ₀", 0.1, 2.0, 0.5)
        if st.sidebar.button("Calculate and Animate 1D"):
            dt = 0.02; time_array = np.linspace(0, dt*time_steps, time_steps)
            with st.spinner(f"Calculating {time_steps} 1D steps..."):
                x = np.linspace(-20, 20, grid_points)
                if potential_type_tdse == "Harmonic oscillator": potential = harmonic_oscillator_1d(x, omega=1.0); potential_func = lambda pos: harmonic_oscillator_1d(pos, omega=1.0)
                elif potential_type_tdse == "Potential Barrier (Tunnelling)": potential = potential_barrier_1d(x, V0=5.0, width=0.5, center=0.0); potential_func = lambda pos: potential_barrier_1d(pos, V0=5.0, width=0.5, center=0.0)
                else: potential = np.zeros_like(x); potential_func = lambda pos: 0.0
                psi0 = (1 / (2 * np.pi * sigma0**2))**(1/4) * np.exp(-(x - x0)**2 / (4 * sigma0**2)) * np.exp(1j * p0 * x); psi0 /= np.sqrt(np.sum(np.abs(psi0)**2) * (x[1]-x[0]))
                psi_t = cached_solve_time_1d(tuple(potential), tuple(x), tuple(psi0), dt=dt, n_steps=time_steps)
            
            st.header("1D Time-Dependent Schrödinger Equation"); fig, (ax_main, ax_momentum) = plt.subplots(2, 1, figsize=(8, 8)); plot_placeholder = st.empty(); colours = cm.plasma(np.linspace(0, 1, time_steps))
            frames_1d = []
            for i in range(time_steps):
                ax_main.clear(); ax_momentum.clear(); prob_density = np.abs(psi_t[i])**2; ax_main.plot(x, prob_density, color=colours[i], linewidth=2); ax_main.fill_between(x, prob_density, alpha=0.3, color=colours[i]); potential_max = np.max(np.abs(potential)) if np.max(np.abs(potential)) > 0 else 1.0; prob_max = np.max(np.abs(psi_t[0])**2); ax_main.plot(x, potential / potential_max * prob_max, 'k--', label="V(x)"); ax_main.text(x[-20], potential[-20]/potential_max*prob_max, 'V(x)', fontsize=12, color='black'); ax_main.set_ylim(0, prob_max * 1.2); ax_main.set_title("Position Space"); ax_main.set_ylabel(r'$|\psi(x,t)|^2$'); ax_main.legend(loc="upper right")
                p = fftshift(fftfreq(grid_points, x[1]-x[0]))*2*np.pi; phi_p = fftshift(fft(psi_t[i])); ax_momentum.plot(p, np.abs(phi_p)**2, color=colours[i], linewidth=2); ax_momentum.fill_between(p, np.abs(phi_p)**2, alpha=0.3, color=colours[i]); ax_momentum.set_title("Momentum Space"); ax_momentum.set_ylabel(r'$|\phi(p,t)|^2$'); ax_momentum.set_xlim(-20, 20);
                fig.tight_layout(); buf = io.BytesIO(); fig.savefig(buf, format='png', dpi=100); buf.seek(0); frames_1d.append(imageio.imread(buf)); plot_placeholder.pyplot(fig)
            with st.spinner("Creating GIF..."):
                gif_buf = io.BytesIO(); imageio.mimsave(gif_buf, frames_1d, format='gif', duration=60);
                st.download_button(label="Download Animation as GIF", data=gif_buf.getvalue(), file_name="animation_1d.gif", mime="image/gif")

            # Calculate Expectation Values and Uncertainties over Time
            exp_x_list, exp_p_list, delta_x_list, delta_p_list = [], [], [], []
            for psi in psi_t:
                exp_x, exp_p, delta_x, delta_p = calculate_expectation_values(psi, x)
                exp_x_list.append(exp_x)
                exp_p_list.append(exp_p)
                delta_x_list.append(delta_x)
                delta_p_list.append(delta_p)
            
            # The classical trajectory calculation remains the same
            classical_x = solve_classical(potential_func, x0, p0, (0, dt*time_steps), time_array)

            if potential_type_tdse == "Potential Barrier (Tunnelling)": st.subheader("Scattering Analysis"); psi_final = psi_t[-1]; prob_density_final = np.abs(psi_final)**2; dx = x[1] - x[0]; reflected_mask = x < 0; transmitted_mask = x >= 0; R = np.sum(prob_density_final[reflected_mask]) * dx; T = np.sum(prob_density_final[transmitted_mask]) * dx; col1, col2, col3 = st.columns(3); col1.metric("Reflection Coefficient (R)", f"{R:.3f}"); col2.metric("Transmission Coefficient (T)", f"{T:.3f}"); col3.metric("Total Probability (R+T)", f"{R+T:.3f}");
            
            st.subheader("Ehrenfest's Theorem: Quantum vs Classical Motion")
            fig_ehrenfest, ax_ehrenfest = plt.subplots(figsize=(8,4))
            ax_ehrenfest.plot(time_array, classical_x, 'r--', linewidth=2.5, label='Classical Trajectory')
            ax_ehrenfest.plot(time_array, exp_x_list, 'b-', linewidth=2, label='Quantum Expectation Value <x>')
            ax_ehrenfest.set_xlabel("Time (t)")
            ax_ehrenfest.set_ylabel("Position (x)")
            ax_ehrenfest.set_title("Comparison of Quantum and Classical Paths")
            ax_ehrenfest.legend()
            st.pyplot(fig_ehrenfest)

            st.subheader("Heisenberg's Uncertainty Principle")
            st.markdown(r"A visual confirmation of $\Delta x \Delta p \geq \frac{\hbar}{2}$ Watch how the position uncertainty "
                        r"($\Delta x$) grows for a free particle as the wave packet spreads")
            fig_uncertainty, ax_uncertainty = plt.subplots(figsize=(8, 5))
            
            # Plot uncertainties
            ax_uncertainty.plot(time_array, delta_x_list, color='crimson', linestyle='-', label=r'$\Delta x$ (Position Uncertainty)')
            ax_uncertainty.plot(time_array, delta_p_list, color='mediumblue', linestyle='-', label=r'$\Delta p$ (Momentum Uncertainty)')
            
            # Plot the uncertainty product
            uncertainty_product = np.array(delta_x_list) * np.array(delta_p_list)
            ax_uncertainty.plot(time_array, uncertainty_product, 'k--', linewidth=2.5, label=r'$\Delta x \Delta p$ (Uncertainty Product)')
            
            # Plot the theoretical minimum (hbar/2), hbar is 1 in our units
            ax_uncertainty.axhline(y=0.5, color='gray', linestyle=':', linewidth=2, label=r'Theoretical Minimum $\hbar/2 = 0.5$')
            
            ax_uncertainty.set_xlabel("Time (t)")
            ax_uncertainty.set_ylabel("Value")
            ax_uncertainty.set_title("Time Evolution of Uncertainties")
            ax_uncertainty.legend(loc='best')
            ax_uncertainty.grid(True, linestyle='--', alpha=0.6)
            ax_uncertainty.set_ylim(bottom=0) # Uncertainty cannot be negative
            st.pyplot(fig_uncertainty)

        else: st.info("Set 1D parameters and click 'Calculate and Animate' to begin.")
    elif dimension_tdse == "2D":
        st.sidebar.subheader("2D Settings"); grid_points_2d = st.sidebar.slider("Grid points per side", 20, 100, 50); time_steps_2d = st.sidebar.slider("Time steps", 20, 100, 40)
        st.sidebar.subheader("Initial 2D Wave Packet"); x0_2d = st.sidebar.slider("Initial x₀", -5.0, 5.0, -2.0); y0_2d = st.sidebar.slider("Initial y₀", -5.0, 5.0, 0.0); px0_2d = st.sidebar.slider("Initial momentum pₓ", -5.0, 5.0, 2.0); py0_2d = st.sidebar.slider("Initial momentum pᵧ", -5.0, 5.0, 1.0)
        if st.sidebar.button("Calculate and Animate 2D"):
            with st.spinner(f"Calculating {time_steps_2d} 2D steps (this may take a while)..."):
                span = 8; x = np.linspace(-span, span, grid_points_2d); y = np.linspace(-span, span, grid_points_2d); X, Y = np.meshgrid(x, y)
                V = harmonic_oscillator_2d(X, Y, omega_x=0.1, omega_y=0.1); sigma = 0.5; psi0 = (np.exp(-((X-x0_2d)**2 + (Y-y0_2d)**2)/(4*sigma**2)) * np.exp(1j * (px0_2d*X + py0_2d*Y))); psi0 /= np.sqrt(np.sum(np.abs(psi0)**2))
                psi_t = cached_solve_tdse_2d(V, x, y, psi0, dt=0.02, n_steps=time_steps_2d)
            st.header("2D Time-Dependent Schrödinger Equation"); fig, ax = plt.subplots(figsize=(7, 6)); plot_placeholder = st.pyplot(fig); cmap = plt.get_cmap('magma') 
            frames_2d = []
            for i in range(time_steps_2d):
                ax.clear(); prob_density_2d = np.abs(psi_t[i])**2; im = ax.imshow(prob_density_2d.T, cmap=cmap, interpolation='bilinear', extent=[-span, span, -span, span], origin='lower'); ax.set_title(f"2D Probability Density $|\psi(x, y, t)|^2$ at Time Step {i}"); ax.set_xlabel("x"); ax.set_ylabel("y");
                buf = io.BytesIO(); fig.savefig(buf, format='png', dpi=100); buf.seek(0); frames_2d.append(imageio.imread(buf)); plot_placeholder.pyplot(fig)
            with st.spinner("Creating GIF..."):
                gif_buf = io.BytesIO(); imageio.mimsave(gif_buf, frames_2d, format='gif', duration=100);
                st.download_button(label="Download Animation as GIF", data=gif_buf.getvalue(), file_name="animation_2d.gif", mime="image/gif")
        else: st.info("Set 2D parameters and click 'Calculate and Animate' to begin.")