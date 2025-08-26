# 🌌 Quantum Playground

**Quantum Playground** is a computational quantum mechanics project that numerically solves the **1D time-independent Schrödinger equation** for user-defined potentials.  
It allows you to **compute eigenvalues, eigenfunctions, and probability densities**, and visualize quantum systems interactively.

## ✨ Features
- 🔬 Solve the Schrödinger equation with finite-difference methods  
- 📦 Built-in potentials: particle in a box, harmonic oscillator, double well, barriers  
- 🧑‍💻 Define your own custom potential functions  
- 📊 Visualize energy levels, eigenfunctions, and probability densities  

## 📚 Background
The time-independent Schrödinger equation:

-(ħ² / 2m) * d²ψ(x)/dx² + V(x) * ψ(x) = E * ψ(x)

is discretized using finite differences and turned into a **matrix eigenvalue problem**:

H ψ = E ψ
where \(H\) is the Hamiltonian. Solving this numerically yields energy eigenvalues and eigenfunctions.

---

## 🛠 Installation
Clone the repository and install dependencies:

bash
git clone https://github.com/Mezzirak/quantum-playground.git
cd quantum_playground
pip install -r requirements.txt 


🚀 Usage

Run the solver

python main.py

Example: Harmonic Oscillator

from solver import solve_schrodinger
from potentials import harmonic_oscillator
import visualization as viz

x, energies, wavefuncs = solve_schrodinger(harmonic_oscillator)

viz.plot_eigenstates(x, energies, wavefuncs, n_states=5)

This will plot the first few eigenstates of the harmonic oscillator.
