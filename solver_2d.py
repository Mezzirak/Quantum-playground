import numpy as np
from scipy.sparse import diags, kron, identity
from scipy.sparse.linalg import eigsh

def construct_hamiltonian_2d(x, y, potential_2d):
    """
    Construct the 2D Hamiltonian matrix for a quantum system using
    finite differences and sparse matrices.
    """
    Nx = len(x)
    Ny = len(y)
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    # 1D second derivative operators for x and y
    diag_x = np.full(Nx, -2.0)
    off_x = np.full(Nx - 1, 1.0)
    D2x = diags([off_x, diag_x, off_x], offsets=[-1, 0, 1]) / (dx ** 2)

    diag_y = np.full(Ny, -2.0)
    off_y = np.full(Ny - 1, 1.0)
    D2y = diags([off_y, diag_y, off_y], offsets=[-1, 0, 1]) / (dy ** 2)

    # 2D Laplacian using Kronecker products
    Ix = identity(Nx)
    Iy = identity(Ny)
    Laplacian = kron(D2x, Iy) + kron(Ix, D2y)

    # Flatten potential to match Hamiltonian
    V_flat = potential_2d.flatten()
    Potential = diags(V_flat, 0)

    # Hamiltonian: H = -0.5 * Laplacian + V
    H = -0.5 * Laplacian + Potential

    return H

def solve_schrodinger_2d(x, y, potential_2d, num_states=5):
    """
    Solve the 2D time-independent Schrödinger equation.
    Returns eigenvalues and eigenfunctions (reshaped to 2D grids).
    """
    H = construct_hamiltonian_2d(x, y, potential_2d)

    # Sparse eigenvalue solver: lowest 'num_states' energies
    energies, wavefuncs = eigsh(H, k=num_states, which="SM")

    # Sort by energy
    idx = np.argsort(energies)
    energies = energies[idx]
    wavefuncs = wavefuncs[:, idx]

    # Reshape eigenvectors and normalise
    wavefuncs_2d = []
    for i in range(num_states):
        psi = wavefuncs[:, i].reshape((len(x), len(y)))
        norm = np.trapz(np.trapz(np.abs(psi)**2, y), x)
        psi /= np.sqrt(norm)
        wavefuncs_2d.append(psi)

    return energies, wavefuncs_2d