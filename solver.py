import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

def construct_hamiltonian(x, potential):
    """
    Construct the Hamiltonian matrix for a 1D quantum system using
    finite differences and sparse storage.
    """
    N = len(x)
    dx = x[1] - x[0]

    # Kinetic energy operator (second derivative with central difference)
    diagonal = np.full(N, -2.0)
    off_diagonal = np.full(N - 1, 1.0)
    laplacian = diags([off_diagonal, diagonal, off_diagonal], offsets=[-1, 0, 1]) / (dx ** 2)

    # Potential energy operator (diagonal matrix)
    potential_matrix = diags(potential, 0)

    # Hamiltonian H = -0.5 * Laplacian + V(x)
    H = -0.5 * laplacian + potential_matrix
    return H

def solve_schrodinger(x, potential, num_states=5):
    """
    Solve the time-independent Schrödinger equation in 1D using a sparse solver.
    Returns eigenvalues (energies) and eigenvectors (wavefunctions).
    """
    H = construct_hamiltonian(x, potential)

    # Use eigsh for sparse symmetric matrices (Hermitian)
    energies, wavefuncs = eigsh(H, k=num_states, which="SM")  # "SM" = smallest magnitude eigenvalues

    # Sort results by energy
    idx = np.argsort(energies)
    energies = energies[idx]
    wavefuncs = wavefuncs[:, idx]

    # Normalise wavefunctions
    for i in range(num_states):
        norm = np.trapz(np.abs(wavefuncs[:, i])**2, x)
        wavefuncs[:, i] /= np.sqrt(norm)

    return energies, wavefuncs


