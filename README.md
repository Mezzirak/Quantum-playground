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

\[
-\frac{\hbar^2}{2m} \frac{d^2 \psi(x)}{dx^2} + V(x)\psi(x) = E\psi(x)
\]

is discretized using finite differences and turned into a **matrix eigenvalue problem**:

\[
H \psi = E \psi
\]

where \(H\) is the Hamiltonian. Solving this numerically yields energy eigenvalues and eigenfunctions.

---

## 🛠 Installation
Clone the repository and install dependencies:

```bash
git clone https://github.com/Mezzirak/quantum-playground.git
cd quantum_playground
pip install -r requirements.txt '''



