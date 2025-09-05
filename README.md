## Quantum Playground

Quantum Playground is an interactive Python application for exploring the stationary Schrödinger equation in one and two dimensions. It allows users to study quantum systems under a variety of potentials, including harmonic oscillators, infinite square wells, and double wells, as well as user-defined custom potentials.

This project was developed as a self-study tool based on Griffiths' "Introduction to Quantum Mechanics."

It allows you to **compute eigenvalues, eigenfunctions, and probability densities**, and visualise quantum systems interactively.

## Features

The application is organised into two main modes, each with 1D and 2D capabilities:

### Stationary States
- Solves the time-independent Schrödinger equation for bound states in 1D and 2D.
- **1D**: Calculates and visualises energy eigenvalues, eigenfunctions ($ψ_n(x)$), and probability densities ($|\psi_n(x)|^2$).
- **2D**: Displays interactive 3D surface plots of the probability densities ($|\psi_n(x,y)|^2$) using Plotly.
- **Analysis**: Includes a tool to numerically verify the orthogonality of eigenstates.

### 🏃‍♂️ Time Evolution
- Solves the time-dependent Schrödinger equation to simulate the dynamics of a wave packet in 1D and 2D.
- **Physics Visualised**:
  - **Quantum Tunnelling**: A 1D wave packet scattering off a potential barrier.
  - **Momentum Space**: Simultaneously visualises the wave function's evolution in momentum space ($|\phi(p,t)|^2$).
  - **Ehrenfest's Theorem**: Compares the quantum expectation value `<x>` to the trajectory of a classical particle, demonstrating the correspondence principle.
- **Analysis Tools**:
  - **Scattering Analysis**: Calculates the Reflection (R) and Transmission (T) coefficients from 1D scattering simulations.
  - **GIF Export**: Allows the user to download 1D and 2D animations as GIF files for easy sharing.


## 📚 Background

The project numerically solves the two fundamental forms of the Schrödinger equation.

The **time-independent Schrödinger equation** is treated as a matrix eigenvalue problem to find stationary states:
```math
\left[ -\frac{\hbar^2}{2m}\frac{d^2}{dx^2} + V(x) \right] \psi(x) = E \psi(x)
```
The **time-dependent Schrödinger equation** is solved as an initial value problem using the Crank-Nicolson method to simulate dynamics:

```math
i\hbar \frac{\partial}{\partial t}\Psi(x,t) = \left[ -\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x) \right] \Psi(x,t)
```
**Technology Stack**

- Framework: Streamlit

- Numerical Computation: NumPy & SciPy

- Plotting: Matplotlib & Plotly

- GIF Creation: imageio & Pillow

**Installation**

1. Clone the repository and install dependencies:

</pre>```bash
git clone https://github.com/Mezzirak/quantum-playground.git
cd quantum_playground```</pre>

2. Create a virtual environment and activate it:

</pre>```
python3 -m venv venv
source venv/bin/activate  # Mac/Linux
venv\Scripts\activate     # Windows ```</pre>

3. Install dependincies:
   
</pre>```python3 -m pip install -r requirements.txt```</pre>

If requirements.txt is not present, install directly:

</pre>```pip install streamlit numpy matplotlib plotly scipy```</pre>

**Usage**

Run the application:

</pre>```streamlit run app.py```</pre>

- Select 1D or 2D in the sidebar
- Pick a predefined protential or use a custome potential
- Adjust the number or eigenstates, grid points and potential parameters
- Visualise eigenstates and probability densities interactively

**Project structure**

- app.py                  # Main Streamlit application
- solver.py               # 1D Schrödinger equation solver
- solver_2d.py            # 2D Schrödinger equation solver
- potentials.py           # Definitions of 1D and 2D potentials
- visualisation.py        # 1D plotting helper functions
- venv/                   # Virtual environment
- README.md               # This file
- requriements.txt        # Python dependencies
- .gitignore              # Ignore __pycache__ and .pyc files




