import numpy as np
from scipy.linalg import eigh_tridiagonal

def solve_1d(V, x, n_eigen=3):
    """
    Solve 1D stationary Schrödinger equation for a given potential.

    Parameters
    V : array_like
        Potential array, same size as x
    x : array_like
        Spatial grid
    n_eigen : int, optional
        Number of lowest eigenstates to compute (default is 3)

    Returns
    energies : ndarray
        Eigenvalues (energies)
    wavefunctions : ndarray
        Eigenfunctions, shape (n_eigen, len(x))
    """
    N = len(x)
    dx = x[1] - x[0]

    # Kinetic energy using finite difference (tridiagonal)
    diag = np.full(N, 1.0 / dx**2)
    off_diag = np.full(N-1, -0.5 / dx**2)

    H_diag = diag + V  # Hamiltonian main diagonal
    H_off_diag = off_diag

    # Compute eigenvalues and eigenvectors
    energies, wavefunctions = eigh_tridiagonal(H_diag, H_off_diag)

    # Select the lowest n_eigen states
    return energies[:n_eigen], wavefunctions[:, :n_eigen].T