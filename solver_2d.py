import numpy as np
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import eigsh

def solve_2d(V, x, y, n_eigen=3):
    """
    Solve 2D stationary Schrödinger equation for a given potential.

    Parameters
    V : 2D array
        Potential grid, shape (len(y), len(x))
    x, y : 1D arrays
        Spatial grids
    n_eigen : int
        Number of lowest eigenstates to compute

    Returns
    energies : ndarray
        Eigenvalues (energies)
    wavefunctions : list of 2D arrays
        Eigenfunctions, each of shape (len(y), len(x))
    """
    Nx = len(x)
    Ny = len(y)
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    # Kinetic energy matrices (sparse)
    Tx = diags([1, -2, 1], [-1, 0, 1], shape=(Nx, Nx)) / dx**2
    Ty = diags([1, -2, 1], [-1, 0, 1], shape=(Ny, Ny)) / dy**2

    # 2D Laplacian using Kronecker products (sparse)
    Lap = kron(eye(Ny), Tx) + kron(Ty, eye(Nx))

    # Flatten potential and make sparse diagonal
    V_flat = V.flatten()
    H = -0.5 * Lap + diags(V_flat)

    # Solve lowest n_eigen states
    energies, wavefunctions_flat = eigsh(H, k=n_eigen, which='SA')

    # Reshape eigenvectors to 2D
    wavefunctions = [wavefunctions_flat[:, i].reshape(Ny, Nx) for i in range(n_eigen)]

    return energies, wavefunctions
