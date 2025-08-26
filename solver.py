import numpy as np
from scipy.sparse import diags
from scipy.linalg import eigh

def solve_schrodinger(V, x=None):
    """
    Solve 1D time-independent Schrödinger equation using finite differences.

    Parameters:
        V : array-like
            Potential energy array
        x : array-like, optional
            Spatial grid

    Returns:
        x : np.ndarray
            Spatial grid
        energies : np.ndarray
            Eigenvalues
        wavefuncs : np.ndarray
            Eigenfunctions (columns correspond to energies)
    """
    if x is None:
        x = np.linspace(-10, 10, 1000)
    dx = x[1] - x[0]

    # Kinetic energy operator (finite difference)
    diag = np.ones(len(x))
    off_diag = np.ones(len(x)-1)
    laplacian = diags([off_diag, -2*diag, off_diag], [-1,0,1]).toarray() / dx**2

    # Hamiltonian
    H = -(0.5) * laplacian + np.diag(V)

    # Solve eigenvalue problem
    energies, wavefuncs = eigh(H)
    return x, energies, wavefuncs

