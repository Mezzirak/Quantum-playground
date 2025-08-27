import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu

def solve_time_1d(V, x, psi0, dt, n_steps):
    """
    Solve the 1D time-dependent Schrödinger equation using Crank-Nicolson.

    Parameters
    V : array_like
        Potential array (len(x))
    x : array_like
        Spatial grid
    psi0 : array_like
        Initial wavefunction (len(x)), should be normalized
    dt : float
        Time step
    n_steps : int
        Number of time steps to evolve

    Returns
    psi_t : ndarray
        Wavefunction at each time step, shape (n_steps, len(x))
    """
    N = len(x)
    dx = x[1] - x[0]
    hbar = 1.0
    m = 1.0

    # Kinetic term coefficient
    r = 1j * hbar * dt / (2 * m * dx**2)

    # Construct the tridiagonal matrices
    main_diag = (1 + 2*r + 1j*dt*V/(2*hbar)) * np.ones(N)
    off_diag = -r * np.ones(N-1)

    A = diags([off_diag, main_diag, off_diag], offsets=[-1,0,1], format='csc')
    B = diags([-off_diag, 2 - main_diag + 1j*dt*V/hbar, -off_diag], offsets=[-1,0,1], format='csc')
    
    # LU decomposition for efficiency
    lu = splu(A)

    psi = psi0.copy()
    psi_t = np.zeros((n_steps, N), dtype=complex)
    psi_t[0] = psi

    for n in range(1, n_steps):
        b = B.dot(psi)
        psi = lu.solve(b)
        # Normalize
        psi /= np.sqrt(np.sum(np.abs(psi)**2) * dx)
        psi_t[n] = psi

    return psi_t

