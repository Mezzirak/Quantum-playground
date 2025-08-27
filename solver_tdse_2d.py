import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu

def solve_tdse_2d(V, x, y, psi0=None, dt=0.01, n_steps=100):
    """
    Solve 2D time-dependent Schrödinger equation using Crank-Nicolson.

    Parameters
    V : ndarray
        2D potential array shape (Nx, Ny)
    x, y : ndarray
        1D spatial grids
    psi0 : ndarray
        Initial wavefunction (Nx, Ny). If None, use Gaussian at center
    dt : float
        Time step
    n_steps : int
        Number of time steps

    Returns
    psi_t : ndarray
        Array of wavefunctions at each time step, shape (n_steps, Nx, Ny)
    """
    Nx, Ny = len(x), len(y)
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    # Initialize wavefunction
    if psi0 is None:
        X, Y = np.meshgrid(x, y, indexing='ij')
        psi0 = np.exp(-(X**2 + Y**2))  # Gaussian at center
        psi0 /= np.sqrt(np.sum(np.abs(psi0)**2) * dx * dy)  # Normalise

    psi = psi0.copy()
    psi_t = np.zeros((n_steps, Nx, Ny), dtype=complex)
    psi_t[0] = psi

    # Precompute 1D Crank-Nicolson matrices
    def crank_matrix(N, d):
        main_diag = np.full(N, 1 + 1j * dt / (2 * d**2))
        off_diag = np.full(N-1, -1j * dt / (4 * d**2))
        A = diags([off_diag, main_diag, off_diag], offsets=[-1,0,1], format='csc')
        B = diags([-off_diag, 1 - 1j * dt / (2*d**2), -off_diag], offsets=[-1,0,1], format='csc')
        return A, B

    Ax, Bx = crank_matrix(Nx, dx)
    Ay, By = crank_matrix(Ny, dy)
    LUx = splu(Ax)
    LUy = splu(Ay)

    # Time evolution (alternating direction implicit)
    for t in range(1, n_steps):
        # Half step in x
        for j in range(Ny):
            rhs = Bx.dot(psi[:, j])
            rhs -= 1j * dt/2 * V[:, j] * psi[:, j]
            psi[:, j] = LUx.solve(rhs)
        # Half step in y
        for i in range(Nx):
            rhs = By.dot(psi[i, :])
            rhs -= 1j * dt/2 * V[i, :] * psi[i, :]
            psi[i, :] = LUy.solve(rhs)

        psi_t[t] = psi

    return psi_t
