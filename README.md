## 🌌 Quantum Playground

Quantum Playground is an interactive Python application for exploring the stationary Schrödinger equation in one and two dimensions. It allows users to study quantum systems under a variety of potentials, including harmonic oscillators, infinite square wells, and double wells, as well as user-defined custom potentials.

It allows you to **compute eigenvalues, eigenfunctions, and probability densities**, and visualise quantum systems interactively.

## Features

- **1D Potentials**
  - Harmonic oscillator
  - Particle in a box
  - Double well
    
- **2D Potentials**
  - Harmonic oscillator
  - Particle in a box
  - Double well
    
- **Interactive Sliders**
  - Adjust physical parameters like `ω`, `L`, `a`, `b`, `ωx`, `ωy`, `Lx`, `Ly`
    
- **Solver**
  - Computes eigenstates and energies using finite difference method

- **Visualisation**
  - 1D: line plots of wavefunctions and probability densities
  - 2D: interactive 3D surfaces of probability densities with Plotly

## 📚 Background
The time-independent Schrödinger equation:

-(ħ² / 2m) * d²ψ(x)/dx² + V(x) * ψ(x) = E * ψ(x)

is discretized using finite differences and turned into a **matrix eigenvalue problem**:

H ψ = E ψ
where \(H\) is the Hamiltonian. Solving this numerically yields energy eigenvalues and eigenfunctions.

---

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

**Future Improvements**

- Add more predefines potentials (e.g. Morse, finite wells)
- Add eigenstate animations to see time evolution.
- 
- Improve solver performance for larger grids.


