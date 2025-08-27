## 🌌 Quantum Playground

Quantum Playground is an interactive Python application for exploring **quantum wavefunctions** and **eigenstates** in both **1D and 2D potentials**. Users can select the potential type, adjust parameters with sliders, and visualise eigenstates and probability densities. The 2D visualisations are fully interactive using **Plotly**, allowing rotation, zoom, and pan.

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

</pre>```python3 -m pip install streamlit numpy matplotlib plotly```</pre>

**Usage**

Run the application:

</pre>```streamlit run app.py```</pre>

- Use the sidebar to select dimension (1D or 2D) and potential type.
- Adjust sliders to change potential parameters and number of eigenstates.
- Visualisations update in real-time.
- For 2D potentials, you can rotate, zoom, and pan the probability density surface.

**Project structure**

Quantum-playground/
│
├─ app.py                  # Main Streamlit application
├─ solver.py               # 1D Schrödinger equation solver
├─ solver_2d.py            # 2D Schrödinger equation solver
├─ potentials.py           # Definitions of 1D and 2D potentials
├─ visualisation.py        # 1D plotting helper functions
├─ venv/                   # Virtual environment
├─ README.md               # This file
└─ .gitignore              # Ignore __pycache__ and .pyc files

**Future Improvements**

- Allow custom potential functions to be uploaded.
- Add eigenstate animations to see time evolution.
- Add 2D wavefunction phase plots.
- Improve solver performance for larger grids.


